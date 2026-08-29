#![cfg(feature = "ply")]

//! This module implements reading of PLY (Polygon File Format, sometimes Stanford Triangle Format)
//! files into `MeshData3`.
//!
//! PLY is the only mainstream mesh format which is a generic typed property container rather than a
//! fixed schema, which is why scanning and metrology pipelines standardized on it. A file declares
//! its own elements and, for each element, an arbitrary list of named typed properties. Conventional
//! names exist for the common ones (`x`/`y`/`z`, `nx`/`ny`/`nz`, `red`/`green`/`blue`), but a file is
//! free to carry anything else, and in practice they do: the Stanford scanning repository ships
//! `confidence` and `intensity` per vertex, and MeshLab reads and writes `quality`.
//!
//! This reader preserves all of it. Properties with a recognized meaning are routed to the typed
//! fields on `PointAttrSet3` and `FaceAttrSet3`. Everything else is preserved under its original
//! name in the open attribute maps, so data that the library does not interpret still survives a
//! round trip.
//!
//! # Mapping
//!
//! On the `vertex` element:
//!
//! | Property                | Destination                                    |
//! |-------------------------|------------------------------------------------|
//! | `x`, `y`, `z`           | the point buffer (required)                    |
//! | `nx`, `ny`, `nz`        | `point_normals`, normalized                    |
//! | `red`, `green`, `blue`  | `point_colors`                                 |
//! | `stdev`, `std_dev`      | `point_stdev`                                  |
//! | anything else           | `point_attrs` under its own name               |
//!
//! On the `face` element:
//!
//! | Property                        | Destination                            |
//! |---------------------------------|----------------------------------------|
//! | `vertex_indices`/`vertex_index` | the face buffer (required)             |
//! | `red`, `green`, `blue`          | `face_colors`                          |
//! | `label`, `labels`               | `face_labels`                          |
//! | anything else                   | `face_attrs` under its own name        |
//!
//! Elements other than `vertex` and `face` are parsed but discarded, since their payload has to be
//! consumed to keep the stream in sync.
//!
//! # Known limitations
//!
//! - An `alpha` channel is discarded, because the attribute sets store colors as RGB.
//! - Polygons with more than three indices are fan-triangulated, and any per-face attribute is
//!   replicated across the triangles a polygon expands into.
//! - Per-corner properties (PLY's face `texcoord`, for instance) are not supported because the mesh
//!   types do not yet have a per-corner attribute domain.

use crate::geom3::attributes3::Attr3;
use crate::geom3::attributes3::FaceAttrSet3;
use crate::geom3::attributes3::PointAttrSet3;
use crate::geom3::mesh::MeshView3;
use crate::geom3::mesh::data::MeshData3;
use crate::{Point3, Result, UnitVec3, Vector3};
use ply_rs_bw::parser::{Parser, Reader};
use ply_rs_bw::ply::{
    BeginList, ElementDef, Encoding, Header, PropertyAccess, PropertyAccessResult, PropertyDef,
    PropertyType, ScalarType,
};
use ply_rs_bw::writer::Writer;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Load a triangle mesh from a PLY file, preserving every property the file carries.
///
/// The file must have a `vertex` element with `x`, `y`, and `z` properties. A `face` element is
/// optional; without one the result is a mesh with points and no faces, which is the natural way a
/// PLY point cloud lands in this type.
///
/// # Arguments
///
/// * `path`: the path to the PLY file
///
/// returns: `Result<MeshData3>`
pub fn load_ply_mesh_data(path: &Path) -> Result<MeshData3> {
    let file = File::open(path)?;
    read_ply_mesh_data(BufReader::new(file))
}

/// Read a triangle mesh from any source of PLY data, preserving every property it carries.
///
/// See [`load_ply_mesh_data`] for the property mapping.
///
/// # Arguments
///
/// * `source`: a buffered reader positioned at the start of the PLY data
///
/// returns: `Result<MeshData3>`
pub fn read_ply_mesh_data<R: BufRead>(source: R) -> Result<MeshData3> {
    let mut reader = Reader::new(source);
    let scalar_parser = Parser::<ScalarRow>::new();
    let header = scalar_parser.read_header(&mut reader)?;

    let mut points: Option<Vec<Point3>> = None;
    let mut point_attrs = PointAttrSet3::empty();
    let mut face_attrs = FaceAttrSet3::empty();
    let mut faces: Vec<[u32; 3]> = Vec::new();

    // The payload sections appear in header declaration order, so every element has to be consumed
    // in that order even if its contents are of no interest, or the stream desynchronizes.
    for (name, def) in header.elements.iter() {
        match name.as_str() {
            "vertex" => {
                let rows = scalar_parser.read_payload_for_element(&mut reader, def, &header)?;
                points = Some(read_points(def, &rows, &mut point_attrs)?);
            }
            "face" => {
                let face_parser = Parser::<FaceRow>::new();
                let rows = face_parser.read_payload_for_element(&mut reader, def, &header)?;
                faces = read_faces(def, &rows, &mut face_attrs)?;
            }
            _ => {
                let skip_parser = Parser::<SkipRow>::new();
                skip_parser.read_payload_for_element(&mut reader, def, &header)?;
            }
        }
    }

    let points = points.ok_or("PLY file has no 'vertex' element")?;
    MeshData3::new_with_attrs(points, faces, point_attrs, face_attrs)
}

/// Read the point data of a PLY file, ignoring any connectivity it declares.
///
/// This is the point-domain half of [`read_ply_mesh_data`], and routes the `vertex` element's
/// properties exactly the same way. It is the reader behind `PointCloud3::load_ply`.
///
/// A `face` element is **not** an error here, but it is reported through the returned flag so that
/// a caller loading what it believes to be a point cloud can refuse a file which is really a mesh
/// rather than silently throwing the connectivity away.
///
/// # Arguments
///
/// * `source`: a buffered reader positioned at the start of the PLY data
///
/// returns: `Result<(Vec<Point3>, PointAttrSet3, bool)>` of the points, their attributes, and
/// whether the file declared any faces
pub fn read_ply_points<R: BufRead>(source: R) -> Result<(Vec<Point3>, PointAttrSet3, bool)> {
    let mut reader = Reader::new(source);
    let scalar_parser = Parser::<ScalarRow>::new();
    let header = scalar_parser.read_header(&mut reader)?;

    let mut points: Option<Vec<Point3>> = None;
    let mut attrs = PointAttrSet3::empty();
    let mut has_faces = false;

    for (name, def) in header.elements.iter() {
        match name.as_str() {
            "vertex" => {
                let rows = scalar_parser.read_payload_for_element(&mut reader, def, &header)?;
                points = Some(read_points(def, &rows, &mut attrs)?);
            }
            // The face payload still has to be consumed to keep the stream in sync, even though
            // nothing is done with it.
            "face" => {
                has_faces = def.count > 0;
                let face_parser = Parser::<FaceRow>::new();
                face_parser.read_payload_for_element(&mut reader, def, &header)?;
            }
            _ => {
                let skip_parser = Parser::<SkipRow>::new();
                skip_parser.read_payload_for_element(&mut reader, def, &header)?;
            }
        }
    }

    let points = points.ok_or("PLY file has no 'vertex' element")?;
    attrs.validate(points.len())?;

    Ok((points, attrs, has_faces))
}

/// Load the point data of a PLY file from disk, ignoring any connectivity it declares.
///
/// See [`read_ply_points`] for what the returned flag means.
///
/// # Arguments
///
/// * `path`: the path to the PLY file
///
/// returns: `Result<(Vec<Point3>, PointAttrSet3, bool)>`
pub fn load_ply_points(path: &Path) -> Result<(Vec<Point3>, PointAttrSet3, bool)> {
    let file = File::open(path)?;
    read_ply_points(BufReader::new(file))
}

// ===============================================================================================
// Writing
// ===============================================================================================

/// Options controlling how a mesh is written to a PLY file.
///
/// There is deliberately no attribute-loss flag here, because PLY can represent everything a
/// `MeshData3` carries. Formats which cannot are the ones that need the guard.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct PlyWriteOpts {
    /// Write a binary little-endian payload rather than ascii. Binary is smaller, faster, and
    /// round-trips floating point exactly. This is the default.
    pub binary: bool,

    /// Comment lines to record in the header.
    pub comments: Vec<String>,

    /// The width to write floating point values at. Defaults to [`PlyPrecision::Double`], which is
    /// lossless.
    pub precision: PlyPrecision,
}

impl Default for PlyWriteOpts {
    fn default() -> Self {
        Self {
            binary: true,
            comments: Vec::new(),
            precision: PlyPrecision::Double,
        }
    }
}

/// The width used for every floating point value written: positions, normals, standard deviations,
/// and the `Scalar` and `Vector` open attributes.
///
/// `MeshData3` holds positions as `f64`, so writing `Double` is lossless while writing `Single`
/// narrows and cannot be undone. `Single` is worth choosing when the data was `f32` at the source
/// anyway, since the extra mantissa bits are then all zero and cost roughly a third of the file size
/// for nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PlyPrecision {
    /// 32-bit floats, declared in the header as `float`.
    Single,

    /// 64-bit floats, declared in the header as `double`. Lossless, and the default.
    #[default]
    Double,
}

impl PlyPrecision {
    /// The PLY scalar type this precision declares.
    fn scalar_type(&self) -> ScalarType {
        match self {
            PlyPrecision::Single => ScalarType::Float,
            PlyPrecision::Double => ScalarType::Double,
        }
    }
}

/// Write a mesh to a PLY file, preserving every attribute it carries.
///
/// Positions and floating point attributes are written as `double`, so nothing is lost to
/// narrowing. See [`load_ply_mesh_data`] for the property naming, which this is the exact inverse
/// of.
///
/// # Arguments
///
/// * `path`: the path to write to, which is overwritten if it exists
/// * `mesh`: the mesh to write
/// * `opts`: encoding and header options
///
/// returns: `Result<()>`
pub fn write_ply_mesh<'a>(
    path: &Path,
    mesh: impl Into<MeshView3<'a>>,
    opts: &PlyWriteOpts,
) -> Result<()> {
    let file = File::create(path)?;
    let mut out = BufWriter::new(file);
    write_ply_to(&mut out, mesh, opts)?;
    out.flush()?;
    Ok(())
}

/// Write a mesh as PLY to any sink, preserving every attribute it carries.
///
/// The payload is streamed one element at a time rather than being materialized, so memory use does
/// not scale with the size of the mesh.
///
/// # Arguments
///
/// * `out`: the sink to write to
/// * `mesh`: the mesh to write
/// * `opts`: encoding and header options
///
/// returns: `Result<()>`
pub fn write_ply_to<'a, W: Write>(
    out: &mut W,
    mesh: impl Into<MeshView3<'a>>,
    opts: &PlyWriteOpts,
) -> Result<()> {
    let mesh = mesh.into();
    let point_cols = point_columns(mesh.points(), mesh.point_attrs());
    let face_cols = face_columns(mesh.face_attrs());

    let mut header = new_header(opts);
    header.elements.insert(
        "vertex".to_string(),
        vertex_element(&point_cols, mesh.point_count(), opts),
    );

    let mut face = ElementDef::new("face".to_string());
    face.count = mesh.face_count();
    face.properties.insert(
        "vertex_indices".to_string(),
        PropertyDef::new(
            "vertex_indices".to_string(),
            PropertyType::List(ScalarType::UChar, ScalarType::Int),
        ),
    );
    for (name, col) in face_cols.iter() {
        face.properties.insert(
            name.clone(),
            PropertyDef::new(
                name.clone(),
                PropertyType::Scalar(col.scalar_type(opts.precision)),
            ),
        );
    }
    header.elements.insert("face".to_string(), face);

    let writer = Writer::<Row>::new();
    writer.write_header(out, &header)?;
    write_vertex_rows(&writer, out, &header, point_cols, mesh.point_count())?;

    let face_def = &header.elements["face"];
    let mut row = Row::new();
    row.cols = face_cols;
    for (i, f) in mesh.faces().iter().enumerate() {
        row.index = i;
        for (slot, index) in row.face.iter_mut().zip(f.iter()) {
            *slot = i32::try_from(*index).map_err(|_| {
                format!(
                    "Point index {index} does not fit the `int` type PLY faces are written with"
                )
            })?;
        }
        write_row(&writer, out, &row, face_def, header.encoding)?;
    }

    Ok(())
}

/// Write points and their attributes to a PLY file, as a `vertex` element with no connectivity.
///
/// The vertex section is byte-identical to what [`write_ply_mesh`] produces for the same
/// points and attributes; the file simply has no `face` element. This is the writer behind
/// `PointCloud3::save_ply`.
///
/// # Arguments
///
/// * `path`: the path to write to, which is overwritten if it exists
/// * `points`: the point positions
/// * `attrs`: the per-point attributes, whose arrays must match the point count
/// * `opts`: encoding and header options
///
/// returns: `Result<()>`
pub fn write_ply_points(
    path: &Path,
    points: &[Point3],
    attrs: &PointAttrSet3,
    opts: &PlyWriteOpts,
) -> Result<()> {
    let file = File::create(path)?;
    let mut out = BufWriter::new(file);
    write_ply_points_to(&mut out, points, attrs, opts)?;
    out.flush()?;
    Ok(())
}

/// Write points and their attributes as PLY to any sink, as a `vertex` element with no
/// connectivity.
///
/// # Arguments
///
/// * `out`: the sink to write to
/// * `points`: the point positions
/// * `attrs`: the per-point attributes, whose arrays must match the point count
/// * `opts`: encoding and header options
///
/// returns: `Result<()>`
pub fn write_ply_points_to<W: Write>(
    out: &mut W,
    points: &[Point3],
    attrs: &PointAttrSet3,
    opts: &PlyWriteOpts,
) -> Result<()> {
    attrs.validate(points.len())?;

    let cols = point_columns(points, attrs);

    let mut header = new_header(opts);
    header.elements.insert(
        "vertex".to_string(),
        vertex_element(&cols, points.len(), opts),
    );

    let writer = Writer::<Row>::new();
    writer.write_header(out, &header)?;
    write_vertex_rows(&writer, out, &header, cols, points.len())
}

/// Start a header carrying the encoding and comments the options ask for.
fn new_header(opts: &PlyWriteOpts) -> Header {
    let mut header = Header::new();
    header.encoding = if opts.binary {
        Encoding::BinaryLittleEndian
    } else {
        Encoding::Ascii
    };
    header.comments = opts.comments.clone();
    header
}

/// Declare the `vertex` element for a set of columns.
fn vertex_element(cols: &[(String, Col<'_>)], count: usize, opts: &PlyWriteOpts) -> ElementDef {
    let mut vertex = ElementDef::new("vertex".to_string());
    vertex.count = count;
    for (name, col) in cols.iter() {
        vertex.properties.insert(
            name.clone(),
            PropertyDef::new(
                name.clone(),
                PropertyType::Scalar(col.scalar_type(opts.precision)),
            ),
        );
    }
    vertex
}

/// Stream the `vertex` payload, one row at a time.
fn write_vertex_rows<W: Write>(
    writer: &Writer<Row>,
    out: &mut W,
    header: &Header,
    cols: Vec<(String, Col<'_>)>,
    count: usize,
) -> Result<()> {
    let def = &header.elements["vertex"];
    let mut row = Row::new();
    row.cols = cols;
    for i in 0..count {
        row.index = i;
        write_row(writer, out, &row, def, header.encoding)?;
    }

    Ok(())
}

/// Write one row in whichever encoding the header declares.
fn write_row<W: Write>(
    writer: &Writer<Row>,
    out: &mut W,
    row: &Row,
    def: &ElementDef,
    encoding: Encoding,
) -> Result<()> {
    match encoding {
        Encoding::Ascii => writer.write_ascii_element(out, row, def)?,
        Encoding::BinaryLittleEndian => writer.write_little_endian_element(out, row, def)?,
        Encoding::BinaryBigEndian => writer.write_big_endian_element(out, row, def)?,
    };
    Ok(())
}

/// One column of values to be written, borrowed from the mesh rather than copied.
///
/// Multi-component sources carry the index of the component this column represents, which is how a
/// position or a vector attribute is split across the three PLY properties it needs.
enum Col<'a> {
    Point(&'a [Point3], usize),
    Unit(&'a [UnitVec3], usize),
    Scalar(&'a [f64]),
    Label(&'a [u32]),
    Vector(&'a [Vector3], usize),
    Color(&'a [[u8; 3]], usize),
}

impl Col<'_> {
    /// The PLY scalar type this column is declared as.
    ///
    /// Only the floating point columns are affected by the requested precision; labels and colors
    /// are integers and have a fixed width.
    fn scalar_type(&self, precision: PlyPrecision) -> ScalarType {
        match self {
            Col::Point(..) | Col::Unit(..) | Col::Scalar(..) | Col::Vector(..) => {
                precision.scalar_type()
            }
            Col::Label(..) => ScalarType::UInt,
            Col::Color(..) => ScalarType::UChar,
        }
    }

    /// The value at a row, as a double, for the columns which are declared that way.
    fn double(&self, i: usize) -> Option<f64> {
        match self {
            Col::Point(v, c) => Some(v[i][*c]),
            Col::Unit(v, c) => Some(v[i][*c]),
            Col::Scalar(v) => Some(v[i]),
            Col::Vector(v, c) => Some(v[i][*c]),
            _ => None,
        }
    }

    /// The value at a row, as an unsigned integer, for the columns which are declared that way.
    fn uint(&self, i: usize) -> Option<u32> {
        match self {
            Col::Label(v) => Some(v[i]),
            _ => None,
        }
    }

    /// The value at a row, as a byte, for the columns which are declared that way.
    fn uchar(&self, i: usize) -> Option<u8> {
        match self {
            Col::Color(v, c) => Some(v[i][*c]),
            _ => None,
        }
    }
}

/// Build the columns for the `vertex` element, in the order they will be declared.
///
/// This takes the point buffer and the point-domain attributes rather than a whole container, so
/// that a mesh and a point cloud produce byte-identical vertex sections.
fn point_columns<'a>(points: &'a [Point3], attrs: &'a PointAttrSet3) -> Vec<(String, Col<'a>)> {
    let mut cols = Vec::new();

    for (name, c) in [("x", 0), ("y", 1), ("z", 2)] {
        cols.push((name.to_string(), Col::Point(points, c)));
    }

    if let Some(n) = attrs.normals() {
        for (name, c) in [("nx", 0), ("ny", 1), ("nz", 2)] {
            cols.push((name.to_string(), Col::Unit(n, c)));
        }
    }

    if let Some(colors) = attrs.colors() {
        for (name, c) in [("red", 0), ("green", 1), ("blue", 2)] {
            cols.push((name.to_string(), Col::Color(colors, c)));
        }
    }

    if let Some(stdev) = attrs.stdev() {
        cols.push(("stdev".to_string(), Col::Scalar(stdev)));
    }

    let mut names: Vec<&str> = attrs.attr_names().collect();
    names.sort_unstable();
    for name in names {
        push_open_attr(&mut cols, name, attrs.attr(name).unwrap());
    }

    cols
}

/// Build the columns for the `face` element, excluding the index list which is always present.
fn face_columns(attrs: &FaceAttrSet3) -> Vec<(String, Col<'_>)> {
    let mut cols = Vec::new();

    if let Some(colors) = attrs.colors() {
        for (name, c) in [("red", 0), ("green", 1), ("blue", 2)] {
            cols.push((name.to_string(), Col::Color(colors, c)));
        }
    }

    if let Some(labels) = attrs.labels() {
        cols.push(("label".to_string(), Col::Label(labels)));
    }

    let mut names: Vec<&str> = attrs.attr_names().collect();
    names.sort_unstable();
    for name in names {
        push_open_attr(&mut cols, name, attrs.attr(name).unwrap());
    }

    cols
}

/// Add the columns needed for one open-map attribute.
///
/// Scalars and labels are a single property. Vectors and colors have no single-property form in PLY
/// and are split across three properties sharing a base name, which the reader folds back together.
fn push_open_attr<'a>(cols: &mut Vec<(String, Col<'a>)>, name: &str, attr: &'a Attr3) {
    match attr {
        Attr3::Scalar(v) => cols.push((name.to_string(), Col::Scalar(v))),
        Attr3::Label(v) => cols.push((name.to_string(), Col::Label(v))),
        Attr3::Vector(v) => {
            for (c, suffix) in VECTOR_SUFFIXES.iter().enumerate() {
                cols.push((format!("{name}{suffix}"), Col::Vector(v, c)));
            }
        }
        Attr3::Color(v) => {
            for (c, suffix) in COLOR_SUFFIXES.iter().enumerate() {
                cols.push((format!("{name}{suffix}"), Col::Color(v, c)));
            }
        }
    }
}

/// A cursor over the columns, presented to the writer as one row at a time.
///
/// Nothing is allocated per row: the columns borrow the mesh, and moving to the next row is just an
/// index bump. The writer only ever reads through this, so `new` exists to satisfy the trait and is
/// never used to produce a row that gets written.
#[derive(Default)]
struct Row<'a> {
    cols: Vec<(String, Col<'a>)>,
    index: usize,
    face: [i32; 3],
}

impl Row<'_> {
    /// Find the column a property name refers to.
    fn col(&self, name: &str) -> Option<&Col<'_>> {
        self.cols.iter().find(|(n, _)| n == name).map(|(_, c)| c)
    }
}

impl PropertyAccess for Row<'_> {
    fn new() -> Self {
        Self::default()
    }

    fn get_double(&self, name: &str) -> Option<f64> {
        self.col(name)?.double(self.index)
    }

    fn get_float(&self, name: &str) -> Option<f32> {
        self.col(name)?.double(self.index).map(|v| v as f32)
    }

    fn get_uint(&self, name: &str) -> Option<u32> {
        self.col(name)?.uint(self.index)
    }

    fn get_uchar(&self, name: &str) -> Option<u8> {
        self.col(name)?.uchar(self.index)
    }

    fn get_list_int(&self, _name: &str) -> Option<&[i32]> {
        Some(&self.face)
    }
}

// ===============================================================================================
// Row accumulators
// ===============================================================================================

/// Accumulates the scalar properties of one row.
///
/// The parser visits properties in header declaration order, so values land in `values` in the same
/// order the scalar properties appear in the element definition, and can be matched up positionally
/// afterward. Every PLY scalar type is exactly representable in `f64`, since the format has no
/// 64-bit integer types, so widening here is lossless.
#[derive(Default)]
struct ScalarRow {
    values: Vec<f64>,
}

impl ScalarRow {
    fn push(&mut self, value: f64) -> PropertyAccessResult {
        self.values.push(value);
        PropertyAccessResult::Set
    }
}

impl PropertyAccess for ScalarRow {
    fn new() -> Self {
        Self::default()
    }

    fn set_char(&mut self, _: &str, v: i8) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_uchar(&mut self, _: &str, v: u8) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_short(&mut self, _: &str, v: i16) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_ushort(&mut self, _: &str, v: u16) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_int(&mut self, _: &str, v: i32) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_uint(&mut self, _: &str, v: u32) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_float(&mut self, _: &str, v: f32) -> PropertyAccessResult {
        self.push(v as f64)
    }
    fn set_double(&mut self, _: &str, v: f64) -> PropertyAccessResult {
        self.push(v)
    }
}

/// Accumulates one row of the `face` element: the index list plus any scalar properties.
///
/// The index list is filled in place through `begin_list_int` when the file declares `int` data,
/// which is overwhelmingly the common case, and falls back to the allocating setter path for the
/// other integer types.
#[derive(Default)]
struct FaceRow {
    indices: Vec<i32>,
    values: Vec<f64>,
}

impl FaceRow {
    fn take_list<T: Copy + Into<i32>>(&mut self, values: Vec<T>) -> PropertyAccessResult {
        self.indices = values.into_iter().map(|v| v.into()).collect();
        PropertyAccessResult::Set
    }
}

impl PropertyAccess for FaceRow {
    fn new() -> Self {
        Self::default()
    }

    fn begin_list_int(&mut self, _: &str, _len: usize) -> BeginList<'_, i32> {
        BeginList::Fill(&mut self.indices)
    }

    fn set_list_char(&mut self, _: &str, v: Vec<i8>) -> PropertyAccessResult {
        self.take_list(v)
    }
    fn set_list_uchar(&mut self, _: &str, v: Vec<u8>) -> PropertyAccessResult {
        self.take_list(v)
    }
    fn set_list_short(&mut self, _: &str, v: Vec<i16>) -> PropertyAccessResult {
        self.take_list(v)
    }
    fn set_list_ushort(&mut self, _: &str, v: Vec<u16>) -> PropertyAccessResult {
        self.take_list(v)
    }
    fn set_list_uint(&mut self, _: &str, v: Vec<u32>) -> PropertyAccessResult {
        // A mesh with more than 2^31 points is not something this library can represent anyway,
        // so a value which does not fit is a malformed file rather than a size limitation.
        self.indices = v.into_iter().map(|x| x as i32).collect();
        PropertyAccessResult::Set
    }

    fn set_char(&mut self, _: &str, v: i8) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_uchar(&mut self, _: &str, v: u8) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_short(&mut self, _: &str, v: i16) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_ushort(&mut self, _: &str, v: u16) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_int(&mut self, _: &str, v: i32) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_uint(&mut self, _: &str, v: u32) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_float(&mut self, _: &str, v: f32) -> PropertyAccessResult {
        self.values.push(v as f64);
        PropertyAccessResult::Set
    }
    fn set_double(&mut self, _: &str, v: f64) -> PropertyAccessResult {
        self.values.push(v);
        PropertyAccessResult::Set
    }
}

/// Consumes and discards a row of an element this reader has no interest in.
struct SkipRow;

impl PropertyAccess for SkipRow {
    fn new() -> Self {
        SkipRow
    }
}

// ===============================================================================================
// Vertex element
// ===============================================================================================

/// Pull the point positions out of the vertex rows, routing every other scalar property either to
/// its typed attribute field or into the open per-point attribute map.
fn read_points(
    def: &ElementDef,
    rows: &[ScalarRow],
    attrs: &mut PointAttrSet3,
) -> Result<Vec<Point3>> {
    let scalars = scalar_properties(def);
    let columns = transpose(&scalars, rows.iter().map(|r| &r.values), "vertex")?;
    let n = rows.len();

    let x = take_column(&columns, "x").ok_or("PLY 'vertex' element has no 'x' property")?;
    let y = take_column(&columns, "y").ok_or("PLY 'vertex' element has no 'y' property")?;
    let z = take_column(&columns, "z").ok_or("PLY 'vertex' element has no 'z' property")?;

    let points = (0..n)
        .map(|i| Point3::new(x.values[i], y.values[i], z.values[i]))
        .collect();

    if let (Some(nx), Some(ny), Some(nz)) = (
        take_column(&columns, "nx"),
        take_column(&columns, "ny"),
        take_column(&columns, "nz"),
    ) {
        let mut normals = Vec::with_capacity(n);
        for i in 0..n {
            let v = Vector3::new(nx.values[i], ny.values[i], nz.values[i]);
            let unit = UnitVec3::try_new(v, 1.0e-12).ok_or_else(|| {
                format!("PLY vertex {i} has a normal of zero length, which cannot be a direction")
            })?;
            normals.push(unit);
        }
        attrs.set_normals(Some(normals), n)?;
    }

    if let Some(colors) = take_colors(&columns, n) {
        attrs.set_colors(Some(colors), n)?;
    }

    if let Some(stdev) = take_column(&columns, "stdev").or_else(|| take_column(&columns, "std_dev"))
    {
        attrs.set_stdev(Some(stdev.values.clone()), n)?;
    }

    let (composites, taken) = take_composites(&columns);
    for (name, attr) in composites {
        attrs.insert_attr(&name, attr, n)?;
    }

    for column in remaining(&columns, POINT_CONSUMED, &taken) {
        attrs.insert_attr(&column.name, column.to_attr(), n)?;
    }

    Ok(points)
}

/// Property names on the vertex element which are consumed by a typed field rather than carried
/// into the open attribute map.
const POINT_CONSUMED: &[&str] = &[
    "x", "y", "z", "nx", "ny", "nz", "red", "green", "blue", "alpha", "stdev", "std_dev",
];

// ===============================================================================================
// Face element
// ===============================================================================================

/// Pull the triangles out of the face rows, fan-triangulating any polygon with more than three
/// indices and replicating per-face attributes across the triangles it expands into.
fn read_faces(
    def: &ElementDef,
    rows: &[FaceRow],
    attrs: &mut FaceAttrSet3,
) -> Result<Vec<[u32; 3]>> {
    let mut faces = Vec::with_capacity(rows.len());

    // Records which source row each emitted triangle came from, so that per-face attributes can be
    // expanded to match after a fan triangulation.
    let mut source_rows = Vec::with_capacity(rows.len());

    for (i, row) in rows.iter().enumerate() {
        if row.indices.len() < 3 {
            return Err(format!(
                "PLY face {i} has {} indices, but a face needs at least 3",
                row.indices.len()
            )
            .into());
        }

        let mut fan = Vec::with_capacity(row.indices.len());
        for index in &row.indices {
            let value = u32::try_from(*index)
                .map_err(|_| format!("PLY face {i} has a negative point index {index}"))?;
            fan.push(value);
        }

        for k in 1..fan.len() - 1 {
            faces.push([fan[0], fan[k], fan[k + 1]]);
            source_rows.push(i);
        }
    }

    let scalars = scalar_properties(def);
    let columns = transpose(&scalars, rows.iter().map(|r| &r.values), "face")?;
    let columns: Vec<Column> = columns
        .into_iter()
        .map(|c| c.expanded(&source_rows))
        .collect();
    let n = faces.len();

    if let Some(colors) = take_colors(&columns, n) {
        attrs.set_colors(Some(colors), n)?;
    }

    if let Some(labels) = take_column(&columns, "label").or_else(|| take_column(&columns, "labels"))
    {
        attrs.set_labels(Some(labels.as_labels()?), n)?;
    }

    let (composites, taken) = take_composites(&columns);
    for (name, attr) in composites {
        attrs.insert_attr(&name, attr, n)?;
    }

    for column in remaining(&columns, FACE_CONSUMED, &taken) {
        attrs.insert_attr(&column.name, column.to_attr(), n)?;
    }

    Ok(faces)
}

/// Property names on the face element which are consumed by a typed field or by the face buffer
/// itself rather than being carried into the open attribute map.
const FACE_CONSUMED: &[&str] = &[
    "vertex_indices",
    "vertex_index",
    "red",
    "green",
    "blue",
    "alpha",
    "label",
    "labels",
];

// ===============================================================================================
// Column helpers
// ===============================================================================================

/// One property's worth of values, pulled out of the row-oriented parse into a column, along with
/// the PLY type it was declared as so that it can be given the right `Attr3` variant.
struct Column {
    name: String,
    declared: ScalarType,
    values: Vec<f64>,
}

impl Column {
    /// Choose a `Attr3` variant for this column.
    ///
    /// Floating point properties become scalars. Integer properties become labels, which is what
    /// the identifiers they usually hold (region, scan pass, material index) want to be, unless the
    /// column actually contains a negative value, in which case it is a signed quantity rather than
    /// an identifier and stays a scalar.
    fn to_attr(&self) -> Attr3 {
        match &self.declared {
            ScalarType::Float | ScalarType::Double => Attr3::Scalar(self.values.clone()),
            _ if self.values.iter().any(|v| *v < 0.0) => Attr3::Scalar(self.values.clone()),
            _ => Attr3::Label(self.values.iter().map(|v| *v as u32).collect()),
        }
    }

    /// Interpret this column as labels, failing if it holds a negative value.
    fn as_labels(&self) -> Result<Vec<u32>> {
        if let Some(v) = self.values.iter().find(|v| **v < 0.0) {
            return Err(format!("PLY property '{}' holds a negative label {v}", self.name).into());
        }
        Ok(self.values.iter().map(|v| *v as u32).collect())
    }

    /// Interpret this column as an 8-bit color channel, clamping out-of-range values.
    fn as_channel(&self) -> Vec<u8> {
        self.values
            .iter()
            .map(|v| v.round().clamp(0.0, 255.0) as u8)
            .collect()
    }

    /// Rebuild this column so that it has one value per emitted triangle, given the source row each
    /// triangle came from.
    fn expanded(self, source_rows: &[usize]) -> Self {
        if source_rows.len() == self.values.len() {
            return self;
        }

        let values = source_rows.iter().map(|i| self.values[*i]).collect();
        Self { values, ..self }
    }
}

/// List the scalar (non-list) properties of an element, in declaration order.
fn scalar_properties(def: &ElementDef) -> Vec<(String, ScalarType)> {
    def.properties
        .iter()
        .filter_map(|(name, prop)| match &prop.data_type {
            PropertyType::Scalar(t) => Some((name.clone(), t.clone())),
            PropertyType::List(_, _) => None,
        })
        .collect()
}

/// Turn row-oriented scalar values into one column per property.
fn transpose<'a>(
    scalars: &[(String, ScalarType)],
    rows: impl IntoIterator<Item = &'a Vec<f64>>,
    element: &str,
) -> Result<Vec<Column>> {
    let mut columns: Vec<Column> = scalars
        .iter()
        .map(|(name, declared)| Column {
            name: name.clone(),
            declared: declared.clone(),
            values: Vec::new(),
        })
        .collect();

    for (i, values) in rows.into_iter().enumerate() {
        if values.len() != columns.len() {
            return Err(format!(
                "PLY '{element}' row {i} produced {} scalar values, but the header declares {}",
                values.len(),
                columns.len()
            )
            .into());
        }

        for (column, value) in columns.iter_mut().zip(values.iter()) {
            column.values.push(*value);
        }
    }

    Ok(columns)
}

/// Find a column by name.
fn take_column<'a>(columns: &'a [Column], name: &str) -> Option<&'a Column> {
    columns.iter().find(|c| c.name == name)
}

/// Assemble an RGB color column from the `red`, `green`, and `blue` properties, if all are present.
fn take_colors(columns: &[Column], n: usize) -> Option<Vec<[u8; 3]>> {
    let r = take_column(columns, "red")?.as_channel();
    let g = take_column(columns, "green")?.as_channel();
    let b = take_column(columns, "blue")?.as_channel();

    Some((0..n).map(|i| [r[i], g[i], b[i]]).collect())
}

/// List the columns which were not consumed by a typed field or folded into a composite.
fn remaining<'a>(columns: &'a [Column], consumed: &[&str], taken: &[String]) -> Vec<&'a Column> {
    columns
        .iter()
        .filter(|c| !consumed.contains(&c.name.as_str()) && !taken.contains(&c.name))
        .collect()
}

/// Suffix triples which are folded back into a single multi-component attribute.
///
/// A `Attr3::Vector` or `Attr3::Color` has no single-property representation in PLY, so it is
/// written as three properties sharing a base name. Recognizing them on the way back in is what makes
/// the round trip exact, and it also picks up the same convention when another tool happens to use
/// it, which is common for directional fields.
const VECTOR_SUFFIXES: [&str; 3] = ["_x", "_y", "_z"];
const COLOR_SUFFIXES: [&str; 3] = ["_red", "_green", "_blue"];

/// Fold any complete suffix triple among the columns back into a single multi-component attribute.
///
/// Returns the reassembled attributes along with the names of every column they consumed, so those
/// columns are not also emitted individually.
fn take_composites(columns: &[Column]) -> (Vec<(String, Attr3)>, Vec<String>) {
    let mut composites = Vec::new();
    let mut taken = Vec::new();

    for column in columns.iter() {
        let Some(base) = column.name.strip_suffix(VECTOR_SUFFIXES[0]) else {
            continue;
        };

        let parts: Vec<&Column> = VECTOR_SUFFIXES
            .iter()
            .filter_map(|s| take_column(columns, &format!("{base}{s}")))
            .collect();
        if parts.len() != 3 {
            continue;
        }

        let values = (0..column.values.len())
            .map(|i| Vector3::new(parts[0].values[i], parts[1].values[i], parts[2].values[i]))
            .collect();
        composites.push((base.to_string(), Attr3::Vector(values)));
        taken.extend(parts.iter().map(|p| p.name.clone()));
    }

    for column in columns.iter() {
        let Some(base) = column.name.strip_suffix(COLOR_SUFFIXES[0]) else {
            continue;
        };

        let parts: Vec<&Column> = COLOR_SUFFIXES
            .iter()
            .filter_map(|s| take_column(columns, &format!("{base}{s}")))
            .collect();
        if parts.len() != 3 {
            continue;
        }

        let channels: Vec<Vec<u8>> = parts.iter().map(|p| p.as_channel()).collect();
        let values = (0..column.values.len())
            .map(|i| [channels[0][i], channels[1][i], channels[2][i]])
            .collect();
        composites.push((base.to_string(), Attr3::Color(values)));
        taken.extend(parts.iter().map(|p| p.name.clone()));
    }

    (composites, taken)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::get_test_file_path;
    use approx::assert_relative_eq;
    use std::io::Cursor;

    /// A binary little-endian file from Zeiss Inspect, positions and faces only.
    #[test]
    fn reads_a_binary_file() -> Result<()> {
        let mesh = load_ply_mesh_data(&get_test_file_path("sample-clip.ply"))?;

        assert_eq!(mesh.point_count(), 41706);
        assert_eq!(mesh.face_count(), 82275);

        // Nothing beyond positions and faces was declared, so nothing should have been invented.
        assert!(!mesh.has_attrs());

        Ok(())
    }

    /// An ascii file carrying two per-point scalars which the previous reader discarded.
    #[test]
    fn reads_an_ascii_file_and_keeps_its_extra_properties() -> Result<()> {
        let mesh = load_ply_mesh_data(&get_test_file_path("bun_zipper_res4.ply"))?;

        assert_eq!(mesh.point_count(), 453);
        assert_eq!(mesh.face_count(), 948);

        let confidence = mesh
            .point_attr("confidence")
            .expect("the file declares a 'confidence' property")
            .as_scalar()
            .expect("a float property should become a scalar attribute");
        assert_eq!(confidence.len(), 453);
        assert_relative_eq!(confidence[0], 0.850855, epsilon = 1.0e-6);

        let intensity = mesh
            .point_attr("intensity")
            .expect("the file declares an 'intensity' property")
            .as_scalar()
            .expect("a float property should become a scalar attribute");
        assert_relative_eq!(intensity[0], 0.5, epsilon = 1.0e-6);

        // The positions must not have been polluted by the extra properties.
        assert_relative_eq!(mesh.points()[0].x, -0.0312216, epsilon = 1.0e-6);
        assert_relative_eq!(mesh.points()[0].y, 0.126304, epsilon = 1.0e-6);
        assert_relative_eq!(mesh.points()[0].z, 0.00514924, epsilon = 1.0e-6);

        Ok(())
    }

    /// How far a vertex may sit from the embedded fixture's idea of it.
    ///
    /// The fixture is a quantized file, so it never held the PLY's coordinates to better than its
    /// own tolerance, which is about 1.9 um. It was re-encoded once when the point block format
    /// changed, and re-encoding a quantized file necessarily moves its points by up to that same
    /// tolerance again. This is the sum of the two, rounded up, and it is still four orders of
    /// magnitude below anything a genuine reader fault would produce on a 150 mm part.
    const REFERENCE_EPS: f64 = 4e-6;

    /// The positions and faces must match what the old reader produced, so that swapping the
    /// implementation cannot have moved any geometry.
    ///
    /// The embedded fixture is a tcmesh, which renumbers vertices on write, so the two cannot be
    /// compared index for index. The renumbering is derived from connectivity alone, so recomputing
    /// the same plan from the freshly parsed faces says where each one went. A reader that changed
    /// either the geometry or the vertex order still fails this.
    #[test]
    fn agrees_with_the_embedded_reference_mesh() -> Result<()> {
        let mesh = load_ply_mesh_data(&get_test_file_path("bun_zipper_res4.ply"))?;
        let expected = crate::tests::stanford_bun_4();

        assert_eq!(mesh.point_count(), expected.points().len());
        assert_eq!(mesh.face_count(), expected.faces().len());

        let plan = tol_compress::reorder::optimize(mesh.faces(), mesh.point_count())?;

        for (new, &old) in plan.face_order.iter().enumerate() {
            for corner in 0..3 {
                let a = mesh.points()[mesh.faces()[old as usize][corner] as usize];
                let b = expected.points()[expected.faces()[new][corner] as usize];
                assert_relative_eq!(a, b, epsilon = REFERENCE_EPS);
            }
        }

        for (new, &old) in plan.vertex_order.iter().enumerate() {
            assert_relative_eq!(
                mesh.points()[old as usize],
                expected.points()[new],
                epsilon = REFERENCE_EPS
            );
        }

        Ok(())
    }

    fn ascii_ply(body: &str) -> Result<MeshData3> {
        read_ply_mesh_data(Cursor::new(body.to_string()))
    }

    #[test]
    fn routes_recognized_properties_to_typed_fields() -> Result<()> {
        let mesh = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 3\n\
             property float x\n\
             property float y\n\
             property float z\n\
             property float nx\n\
             property float ny\n\
             property float nz\n\
             property uchar red\n\
             property uchar green\n\
             property uchar blue\n\
             property float stdev\n\
             element face 1\n\
             property list uchar int vertex_indices\n\
             property uchar label\n\
             end_header\n\
             0 0 0 0 0 2 255 0 0 0.01\n\
             1 0 0 0 0 1 0 255 0 0.02\n\
             0 1 0 0 0 1 0 0 255 0.03\n\
             3 0 1 2 7\n",
        )?;

        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);

        // Normals are normalized on the way in, so the length-2 vector becomes a unit vector.
        let normals = mesh.point_normals().expect("normals should be present");
        assert_relative_eq!(normals[0].into_inner(), Vector3::z(), epsilon = 1.0e-12);

        assert_eq!(mesh.point_colors().unwrap()[0], [255, 0, 0]);
        // Declared as `float`, so the value arrives with f32 precision.
        assert_relative_eq!(mesh.point_stdev().unwrap()[0], 0.01, epsilon = 1.0e-8);
        assert_eq!(mesh.face_labels().unwrap(), &[7]);

        // Everything recognized went to a typed field, so the open maps stay empty.
        assert_eq!(mesh.point_attrs().attr_names().count(), 0);
        assert_eq!(mesh.face_attrs().attr_names().count(), 0);

        Ok(())
    }

    #[test]
    fn integer_properties_become_labels_unless_they_go_negative() -> Result<()> {
        let mesh = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 2\n\
             property float x\n\
             property float y\n\
             property float z\n\
             property uint scan_pass\n\
             property int offset\n\
             property float quality\n\
             end_header\n\
             0 0 0 4 -1 0.5\n\
             1 0 0 9 2 0.25\n",
        )?;

        assert_eq!(
            mesh.point_attr("scan_pass").unwrap().as_label().unwrap(),
            &[4, 9]
        );
        assert_eq!(
            mesh.point_attr("offset").unwrap().as_scalar().unwrap(),
            &[-1.0, 2.0]
        );
        assert_eq!(
            mesh.point_attr("quality").unwrap().as_scalar().unwrap(),
            &[0.5, 0.25]
        );

        Ok(())
    }

    #[test]
    fn polygons_are_fan_triangulated_and_face_attributes_follow() -> Result<()> {
        let mesh = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 5\n\
             property float x\n\
             property float y\n\
             property float z\n\
             element face 2\n\
             property list uchar int vertex_indices\n\
             property float quality\n\
             end_header\n\
             0 0 0\n\
             1 0 0\n\
             1 1 0\n\
             0 1 0\n\
             2 2 0\n\
             4 0 1 2 3 0.75\n\
             3 0 1 4 0.25\n",
        )?;

        // The quad becomes two triangles and the triangle stays one, so three faces in total.
        assert_eq!(mesh.face_count(), 3);
        assert_eq!(mesh.faces(), &[[0, 1, 2], [0, 2, 3], [0, 1, 4]]);

        // The quad's attribute value is replicated across both of its triangles.
        assert_eq!(
            mesh.face_attr("quality").unwrap().as_scalar().unwrap(),
            &[0.75, 0.75, 0.25]
        );

        Ok(())
    }

    #[test]
    fn a_file_without_faces_loads_as_points_only() -> Result<()> {
        let mesh = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 2\n\
             property float x\n\
             property float y\n\
             property float z\n\
             end_header\n\
             0 0 0\n\
             1 2 3\n",
        )?;

        assert_eq!(mesh.point_count(), 2);
        assert_eq!(mesh.face_count(), 0);
        assert!(!mesh.is_empty());

        Ok(())
    }

    #[test]
    fn unknown_elements_are_skipped_without_desynchronizing() -> Result<()> {
        let mesh = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 2\n\
             property float x\n\
             property float y\n\
             property float z\n\
             element camera 1\n\
             property float view_px\n\
             property list uchar int extra\n\
             element face 1\n\
             property list uchar int vertex_indices\n\
             end_header\n\
             0 0 0\n\
             1 0 0\n\
             1.5 2 8 9\n\
             3 0 1 1\n",
        )?;

        assert_eq!(mesh.point_count(), 2);
        assert_eq!(mesh.faces(), &[[0, 1, 1]]);

        Ok(())
    }

    #[test]
    fn a_missing_vertex_element_is_an_error() {
        let result = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element face 1\n\
             property list uchar int vertex_indices\n\
             end_header\n\
             3 0 1 2\n",
        );

        assert!(result.is_err());
    }

    #[test]
    fn an_out_of_range_face_index_is_an_error() {
        let result = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 2\n\
             property float x\n\
             property float y\n\
             property float z\n\
             element face 1\n\
             property list uchar int vertex_indices\n\
             end_header\n\
             0 0 0\n\
             1 0 0\n\
             3 0 1 5\n",
        );

        assert!(result.is_err());
    }

    // ===========================================================================================
    // Writing
    // ===========================================================================================

    /// Round trip a mesh through an in-memory PLY and hand back what comes out.
    fn round_trip(mesh: &MeshData3, binary: bool) -> Result<MeshData3> {
        let opts = PlyWriteOpts {
            binary,
            ..Default::default()
        };
        let mut buffer = Vec::new();
        write_ply_to(&mut buffer, mesh, &opts)?;
        read_ply_mesh_data(Cursor::new(buffer))
    }

    /// A mesh carrying every typed field and one open attribute of each `Attr3` variant.
    fn loaded_mesh() -> MeshData3 {
        let mut mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )
        .unwrap();

        mesh.set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 3]))
            .unwrap();
        mesh.set_point_colors(Some(vec![[255, 0, 0], [0, 255, 0], [0, 0, 255]]))
            .unwrap();
        mesh.set_point_stdev(Some(vec![0.001, 0.002, 0.003]))
            .unwrap();
        mesh.set_face_colors(Some(vec![[10, 20, 30]])).unwrap();
        mesh.set_face_labels(Some(vec![42])).unwrap();

        mesh.insert_point_attr("confidence", Attr3::Scalar(vec![0.25, 0.5, 0.75]))
            .unwrap();
        mesh.insert_point_attr("scan_pass", Attr3::Label(vec![1, 2, 3]))
            .unwrap();
        mesh.insert_point_attr(
            "principal_dir",
            Attr3::Vector(vec![Vector3::x(), Vector3::y(), Vector3::z()]),
        )
        .unwrap();
        mesh.insert_point_attr("shade", Attr3::Color(vec![[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
            .unwrap();
        mesh.insert_face_attr("quality", Attr3::Scalar(vec![0.875]))
            .unwrap();

        mesh
    }

    /// Assert that two meshes are identical in geometry and in every attribute.
    fn assert_same(a: &MeshData3, b: &MeshData3) {
        assert_eq!(a.point_count(), b.point_count());
        assert_eq!(a.faces(), b.faces());
        for (p, q) in a.points().iter().zip(b.points()) {
            assert_relative_eq!(p, q, epsilon = 0.0);
        }

        assert_eq!(a.point_colors(), b.point_colors());
        assert_eq!(a.point_stdev(), b.point_stdev());
        assert_eq!(a.face_colors(), b.face_colors());
        assert_eq!(a.face_labels(), b.face_labels());

        match (a.point_normals(), b.point_normals()) {
            (Some(x), Some(y)) => {
                for (p, q) in x.iter().zip(y) {
                    assert_relative_eq!(p.into_inner(), q.into_inner(), epsilon = 1.0e-15);
                }
            }
            (None, None) => {}
            _ => panic!("normals present on only one side"),
        }

        let mut a_names: Vec<&str> = a.point_attrs().attr_names().collect();
        let mut b_names: Vec<&str> = b.point_attrs().attr_names().collect();
        a_names.sort_unstable();
        b_names.sort_unstable();
        assert_eq!(a_names, b_names, "point attribute names differ");

        for name in a_names {
            assert_eq!(
                a.point_attr(name),
                b.point_attr(name),
                "point attr '{name}'"
            );
        }

        let mut a_names: Vec<&str> = a.face_attrs().attr_names().collect();
        let mut b_names: Vec<&str> = b.face_attrs().attr_names().collect();
        a_names.sort_unstable();
        b_names.sort_unstable();
        assert_eq!(a_names, b_names, "face attribute names differ");

        for name in a_names {
            assert_eq!(a.face_attr(name), b.face_attr(name), "face attr '{name}'");
        }
    }

    #[test]
    fn binary_round_trip_preserves_everything() -> Result<()> {
        let mesh = loaded_mesh();
        assert_same(&mesh, &round_trip(&mesh, true)?);
        Ok(())
    }

    #[test]
    fn ascii_round_trip_preserves_everything() -> Result<()> {
        let mesh = loaded_mesh();
        assert_same(&mesh, &round_trip(&mesh, false)?);
        Ok(())
    }

    #[test]
    fn positions_survive_at_full_double_precision() -> Result<()> {
        // A value which cannot be represented in f32, to catch any narrowing on the way through.
        let awkward = 1.0 / 3.0;
        let mesh = MeshData3::new(
            vec![
                Point3::new(awkward, -awkward, 1.0e-9),
                Point3::new(1.0e9 + awkward, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;

        for binary in [true, false] {
            let back = round_trip(&mesh, binary)?;
            for (p, q) in mesh.points().iter().zip(back.points()) {
                // Exactly equal, not merely close.
                assert_eq!(p, q, "binary = {binary}");
            }
        }

        Ok(())
    }

    #[test]
    fn a_bare_mesh_round_trips_without_inventing_attributes() -> Result<()> {
        let mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;

        let back = round_trip(&mesh, true)?;
        assert!(!back.has_attrs());
        assert_same(&mesh, &back);

        Ok(())
    }

    #[test]
    fn a_mesh_with_no_faces_round_trips() -> Result<()> {
        let mut mesh = MeshData3::new(
            vec![Point3::new(1.0, 2.0, 3.0), Point3::new(4.0, 5.0, 6.0)],
            Vec::new(),
        )?;
        mesh.set_point_stdev(Some(vec![0.1, 0.2]))?;

        let back = round_trip(&mesh, true)?;
        assert_eq!(back.point_count(), 2);
        assert_eq!(back.face_count(), 0);
        assert_eq!(back.point_stdev().unwrap(), &[0.1, 0.2]);

        Ok(())
    }

    #[test]
    fn multi_component_attributes_are_written_as_suffixed_triples() -> Result<()> {
        let mesh = loaded_mesh();
        let opts = PlyWriteOpts {
            binary: false,
            ..Default::default()
        };
        let mut buffer = Vec::new();
        write_ply_to(&mut buffer, &mesh, &opts)?;
        let text = String::from_utf8(buffer).unwrap();

        // A Vector has no single-property form in PLY, so it becomes three doubles.
        assert!(text.contains("property double principal_dir_x"), "{text}");
        assert!(text.contains("property double principal_dir_y"), "{text}");
        assert!(text.contains("property double principal_dir_z"), "{text}");

        // A Color becomes three bytes.
        assert!(text.contains("property uchar shade_red"), "{text}");
        assert!(text.contains("property uchar shade_green"), "{text}");
        assert!(text.contains("property uchar shade_blue"), "{text}");

        Ok(())
    }

    #[test]
    fn single_precision_narrows_and_shrinks_the_output() -> Result<()> {
        let mesh = loaded_mesh();

        let mut wide = Vec::new();
        write_ply_to(&mut wide, &mesh, &PlyWriteOpts::default())?;

        let narrow_opts = PlyWriteOpts {
            precision: PlyPrecision::Single,
            ..Default::default()
        };
        let mut narrow = Vec::new();
        write_ply_to(&mut narrow, &mesh, &narrow_opts)?;

        assert!(
            narrow.len() < wide.len(),
            "single precision produced {} bytes against {} for double",
            narrow.len(),
            wide.len()
        );

        // Values come back exactly as the f32 narrowing left them, not as the originals.
        let back = read_ply_mesh_data(Cursor::new(narrow))?;
        for (original, returned) in mesh.points().iter().zip(back.points()) {
            for c in 0..3 {
                assert_eq!(returned[c], original[c] as f32 as f64);
            }
        }

        // Integer-valued attributes are unaffected by the float width.
        assert_eq!(back.point_colors(), mesh.point_colors());
        assert_eq!(back.face_labels(), mesh.face_labels());
        assert_eq!(
            back.point_attr("scan_pass").unwrap(),
            mesh.point_attr("scan_pass").unwrap()
        );

        Ok(())
    }

    #[test]
    fn precision_is_declared_in_the_header() -> Result<()> {
        let mesh = loaded_mesh();

        for (precision, expected) in [
            (PlyPrecision::Double, "property double x"),
            (PlyPrecision::Single, "property float x"),
        ] {
            let opts = PlyWriteOpts {
                binary: false,
                precision,
                ..Default::default()
            };
            let mut buffer = Vec::new();
            write_ply_to(&mut buffer, &mesh, &opts)?;
            let text = String::from_utf8(buffer).unwrap();
            assert!(text.contains(expected), "expected `{expected}` in:\n{text}");
        }

        Ok(())
    }

    #[test]
    fn a_file_that_was_f32_at_source_round_trips_at_single_precision() -> Result<()> {
        // The sample file declares `property float x`, so widening it into f64 gains nothing and
        // writing it back at single precision is lossless while being about a third smaller.
        let mesh = load_ply_mesh_data(&get_test_file_path("sample-clip.ply"))?;

        let opts = PlyWriteOpts {
            precision: PlyPrecision::Single,
            ..Default::default()
        };
        let mut buffer = Vec::new();
        write_ply_to(&mut buffer, &mesh, &opts)?;

        let back = read_ply_mesh_data(Cursor::new(buffer))?;
        assert_same(&mesh, &back);

        Ok(())
    }

    #[test]
    fn the_default_options_are_binary_and_lossless() {
        let opts = PlyWriteOpts::default();
        assert!(opts.binary);
        assert_eq!(opts.precision, PlyPrecision::Double);
    }

    #[test]
    fn comments_are_recorded_in_the_header() -> Result<()> {
        let opts = PlyWriteOpts {
            binary: false,
            comments: vec!["Unit: mm".to_string()],
            ..Default::default()
        };
        let mut buffer = Vec::new();
        write_ply_to(&mut buffer, &loaded_mesh(), &opts)?;

        let text = String::from_utf8(buffer).unwrap();
        assert!(text.contains("comment Unit: mm"), "{text}");

        Ok(())
    }

    #[test]
    fn a_file_read_from_disk_round_trips_through_the_writer() -> Result<()> {
        // The ascii bunny carries `confidence` and `intensity`, so this exercises open attributes
        // that came from a real file rather than ones this test constructed.
        let mesh = load_ply_mesh_data(&get_test_file_path("bun_zipper_res4.ply"))?;
        assert_same(&mesh, &round_trip(&mesh, true)?);
        assert_same(&mesh, &round_trip(&mesh, false)?);
        Ok(())
    }

    #[test]
    fn a_large_binary_file_round_trips() -> Result<()> {
        let mesh = load_ply_mesh_data(&get_test_file_path("sample-clip.ply"))?;
        let back = round_trip(&mesh, true)?;

        assert_eq!(back.point_count(), 41706);
        assert_eq!(back.face_count(), 82275);
        assert_same(&mesh, &back);

        Ok(())
    }

    #[test]
    fn a_zero_length_normal_is_an_error() {
        let result = ascii_ply(
            "ply\n\
             format ascii 1.0\n\
             element vertex 1\n\
             property float x\n\
             property float y\n\
             property float z\n\
             property float nx\n\
             property float ny\n\
             property float nz\n\
             end_header\n\
             0 0 0 0 0 0\n",
        );

        assert!(result.is_err());
    }

    // ===========================================================================================
    // Point-only reading and writing
    // ===========================================================================================

    /// The points and attributes of `loaded_mesh`, without its connectivity.
    fn loaded_points() -> (Vec<Point3>, PointAttrSet3) {
        let mesh = loaded_mesh();
        (mesh.points().to_vec(), mesh.point_attrs().clone())
    }

    /// Round trip points and their attributes through an in-memory PLY.
    fn round_trip_points(
        points: &[Point3],
        attrs: &PointAttrSet3,
        binary: bool,
    ) -> Result<(Vec<Point3>, PointAttrSet3, bool)> {
        let opts = PlyWriteOpts {
            binary,
            ..Default::default()
        };
        let mut buffer = Vec::new();
        write_ply_points_to(&mut buffer, points, attrs, &opts)?;
        read_ply_points(Cursor::new(buffer))
    }

    #[test]
    fn point_round_trip_preserves_every_attribute() -> Result<()> {
        let (points, attrs) = loaded_points();

        for binary in [true, false] {
            let (back_points, back_attrs, has_faces) = round_trip_points(&points, &attrs, binary)?;

            assert!(!has_faces, "a point-only file must declare no faces");
            assert_eq!(back_points.len(), points.len());
            for (p, q) in points.iter().zip(back_points.iter()) {
                assert_relative_eq!(p, q, epsilon = 0.0);
            }

            assert_eq!(back_attrs.colors(), attrs.colors(), "binary = {binary}");
            assert_eq!(back_attrs.stdev(), attrs.stdev(), "binary = {binary}");

            for (a, b) in back_attrs
                .normals()
                .unwrap()
                .iter()
                .zip(attrs.normals().unwrap())
            {
                assert_relative_eq!(a.into_inner(), b.into_inner(), epsilon = 1.0e-15);
            }

            let mut names: Vec<&str> = back_attrs.attr_names().collect();
            let mut expected: Vec<&str> = attrs.attr_names().collect();
            names.sort_unstable();
            expected.sort_unstable();
            assert_eq!(names, expected, "open attribute names differ");

            for name in names {
                assert_eq!(back_attrs.attr(name), attrs.attr(name), "attr '{name}'");
            }
        }

        Ok(())
    }

    /// The vertex section a point cloud writes has to be identical to the one a mesh writes for the
    /// same points, or the two readers would not agree on what a file means.
    #[test]
    fn the_vertex_section_matches_what_a_mesh_writes() -> Result<()> {
        let mesh = loaded_mesh();
        let (points, attrs) = loaded_points();
        let opts = PlyWriteOpts::default();

        let mut mesh_bytes = Vec::new();
        write_ply_to(&mut mesh_bytes, &mesh, &opts)?;

        let mut point_bytes = Vec::new();
        write_ply_points_to(&mut point_bytes, &points, &attrs, &opts)?;

        // The point-only header is a prefix of the mesh header up to where the face element is
        // declared, and the vertex payload which follows it is byte for byte the same.
        let split = point_bytes
            .windows(END_HEADER.len())
            .position(|w| w == END_HEADER)
            .expect("the point-only file must have a header")
            + END_HEADER.len();

        let mesh_split = mesh_bytes
            .windows(END_HEADER.len())
            .position(|w| w == END_HEADER)
            .expect("the mesh file must have a header")
            + END_HEADER.len();

        let point_payload = &point_bytes[split..];
        let mesh_payload = &mesh_bytes[mesh_split..mesh_split + point_payload.len()];
        assert_eq!(point_payload, mesh_payload);

        Ok(())
    }

    const END_HEADER: &[u8] = b"end_header\n";

    /// A file which really is a mesh has to be distinguishable from a point cloud, so that a caller
    /// asking for points can refuse rather than silently discard the connectivity.
    #[test]
    fn reading_points_reports_whether_the_file_had_faces() -> Result<()> {
        let mesh = loaded_mesh();
        let mut buffer = Vec::new();
        write_ply_to(&mut buffer, &mesh, &PlyWriteOpts::default())?;

        let (points, attrs, has_faces) = read_ply_points(Cursor::new(buffer))?;

        assert!(has_faces);
        assert_eq!(points.len(), mesh.point_count());
        assert_eq!(attrs.stdev(), mesh.point_stdev());

        Ok(())
    }

    /// A mesh file whose `face` element is declared but empty is a point cloud in practice.
    #[test]
    fn an_empty_face_element_does_not_count_as_faces() -> Result<()> {
        let mesh = MeshData3::new(
            vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
            Vec::new(),
        )?;

        let mut buffer = Vec::new();
        write_ply_to(&mut buffer, &mesh, &PlyWriteOpts::default())?;

        let (points, _, has_faces) = read_ply_points(Cursor::new(buffer))?;

        assert!(!has_faces);
        assert_eq!(points.len(), 2);

        Ok(())
    }

    #[test]
    fn writing_points_rejects_a_mismatched_attribute() {
        let (points, _) = loaded_points();

        let mut attrs = PointAttrSet3::empty();
        attrs.set_stdev(Some(vec![0.1]), 1).unwrap();

        let mut buffer = Vec::new();
        assert!(
            write_ply_points_to(&mut buffer, &points, &attrs, &PlyWriteOpts::default()).is_err()
        );
    }

    /// Because the writer uses a borrowed view, the accelerated mesh must produce the same bytes as
    /// the plain container, including attributes, without first copying it.
    #[test]
    fn either_container_writes_identical_bytes() -> Result<()> {
        let data = loaded_mesh();
        let mesh = crate::Mesh3::from_data(data.clone(), false)?;

        for binary in [true, false] {
            let opts = PlyWriteOpts {
                binary,
                ..Default::default()
            };

            let mut from_data = Vec::new();
            write_ply_to(&mut from_data, &data, &opts)?;

            let mut from_mesh = Vec::new();
            write_ply_to(&mut from_mesh, &mesh, &opts)?;

            assert_eq!(from_data, from_mesh, "binary = {binary}");
        }

        Ok(())
    }
}

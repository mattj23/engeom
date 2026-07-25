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
//! fields on `MeshAttrSet3`, and everything else is carried through into the open attribute maps
//! under its own name, so that data which the library has no opinion about still survives a round
//! trip.
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
//! - An `alpha` channel is discarded, because `MeshAttrSet3` stores colors as RGB.
//! - Polygons with more than three indices are fan-triangulated, and any per-face attribute is
//!   replicated across the triangles a polygon expands into.
//! - Per-corner properties (PLY's face `texcoord`, for instance) are not supported, since
//!   `MeshAttrSet3` has no per-corner domain yet.

use crate::geom3::mesh::data::{MeshAttr3, MeshAttrSet3, MeshData3};
use crate::{Mesh, Point3, Result, UnitVec3, Vector3};
use ply_rs_bw::parser::{Parser, Reader};
use ply_rs_bw::ply::{
    BeginList, ElementDef, PropertyAccess, PropertyAccessResult, PropertyType, ScalarType,
};
use std::fs::File;
use std::io::{BufRead, BufReader};
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
    let mut attrs = MeshAttrSet3::empty();
    let mut faces: Vec<[u32; 3]> = Vec::new();

    // The payload sections appear in header declaration order, so every element has to be consumed
    // in that order even if its contents are of no interest, or the stream desynchronizes.
    for (name, def) in header.elements.iter() {
        match name.as_str() {
            "vertex" => {
                let rows = scalar_parser.read_payload_for_element(&mut reader, def, &header)?;
                points = Some(read_points(def, &rows, &mut attrs)?);
            }
            "face" => {
                let face_parser = Parser::<FaceRow>::new();
                let rows = face_parser.read_payload_for_element(&mut reader, def, &header)?;
                faces = read_faces(def, &rows, &mut attrs)?;
            }
            _ => {
                let skip_parser = Parser::<SkipRow>::new();
                skip_parser.read_payload_for_element(&mut reader, def, &header)?;
            }
        }
    }

    let points = points.ok_or("PLY file has no 'vertex' element")?;
    MeshData3::new_with_attrs(points, faces, attrs)
}

/// Load a triangle mesh from a PLY file into the accelerated `Mesh` type.
///
/// This is a temporary bridge which exists so that callers written against the old reader keep
/// working. It will be replaced by the conversion from `MeshData3`, at which point every attribute
/// the file carried will survive into the mesh instead of being discarded here.
///
/// # Arguments
///
/// * `path`: the path to the PLY file
///
/// returns: `Result<Mesh>`
pub fn load_ply_mesh(path: &Path) -> Result<Mesh> {
    let data = load_ply_mesh_data(path)?;
    Ok(Mesh::new(
        data.points().to_vec(),
        data.faces().to_vec(),
        false,
    ))
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
    attrs: &mut MeshAttrSet3,
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
        attrs.set_point_normals(Some(normals), n)?;
    }

    if let Some(colors) = take_colors(&columns, n) {
        attrs.set_point_colors(Some(colors), n)?;
    }

    if let Some(stdev) = take_column(&columns, "stdev").or_else(|| take_column(&columns, "std_dev"))
    {
        attrs.set_point_stdev(Some(stdev.values.clone()), n)?;
    }

    for column in remaining(&columns, POINT_CONSUMED) {
        attrs.insert_point_attr(&column.name, column.to_attr(), n)?;
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
    attrs: &mut MeshAttrSet3,
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
        attrs.set_face_colors(Some(colors), n)?;
    }

    if let Some(labels) = take_column(&columns, "label").or_else(|| take_column(&columns, "labels"))
    {
        attrs.set_face_labels(Some(labels.as_labels()?), n)?;
    }

    for column in remaining(&columns, FACE_CONSUMED) {
        attrs.insert_face_attr(&column.name, column.to_attr(), n)?;
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
/// the PLY type it was declared as so that it can be given the right `MeshAttr3` variant.
struct Column {
    name: String,
    declared: ScalarType,
    values: Vec<f64>,
}

impl Column {
    /// Choose a `MeshAttr3` variant for this column.
    ///
    /// Floating point properties become scalars. Integer properties become labels, which is what
    /// the identifiers they usually hold (region, scan pass, material index) want to be, unless the
    /// column actually contains a negative value, in which case it is a signed quantity rather than
    /// an identifier and stays a scalar.
    fn to_attr(&self) -> MeshAttr3 {
        match &self.declared {
            ScalarType::Float | ScalarType::Double => MeshAttr3::Scalar(self.values.clone()),
            _ if self.values.iter().any(|v| *v < 0.0) => MeshAttr3::Scalar(self.values.clone()),
            _ => MeshAttr3::Label(self.values.iter().map(|v| *v as u32).collect()),
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

/// List the columns which were not consumed by a typed field.
fn remaining<'a>(columns: &'a [Column], consumed: &[&str]) -> Vec<&'a Column> {
    columns
        .iter()
        .filter(|c| !consumed.contains(&c.name.as_str()))
        .collect()
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
        assert!(mesh.attrs().is_empty());

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

    /// The positions and faces must match what the old reader produced, so that swapping the
    /// implementation cannot have moved any geometry.
    #[test]
    fn agrees_with_the_embedded_reference_mesh() -> Result<()> {
        let mesh = load_ply_mesh_data(&get_test_file_path("bun_zipper_res4.ply"))?;
        let expected = crate::tests::stanford_bun_4();

        assert_eq!(mesh.point_count(), expected.vertices().len());
        assert_eq!(mesh.face_count(), expected.faces().len());

        for (a, b) in mesh.faces().iter().zip(expected.faces().iter()) {
            assert_eq!(a, b);
        }
        for (a, b) in mesh.points().iter().zip(expected.vertices().iter()) {
            assert_relative_eq!(a, b, epsilon = 0.000002);
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
        assert_eq!(mesh.attrs().point_attr_names().count(), 0);
        assert_eq!(mesh.attrs().face_attr_names().count(), 0);

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
}

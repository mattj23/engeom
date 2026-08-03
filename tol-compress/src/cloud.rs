//! Point cloud containers, conventionally `.tccloud2` and `.tccloud3`.
//!
//! The simplest container there is: a points block and a preamble. No connectivity, no ordering
//! guarantee beyond the one the caller supplied, nothing else.
//!
//! Point order *is* preserved, and once increment 5 introduces spatial reordering that will stop
//! being true unless a caller asks for it.
//!
//! Like every container in this crate, a cloud file is an ordered collection, so several scans or
//! several passes can live in one file. See [`crate::find_by_name`].

use crate::container::{self, Kind, Named, item};
use crate::error::{Error, Result};
use crate::metadata::Metadata;
use crate::points::{read_points, write_points};
use crate::raw::MAX_PREALLOC;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// An unordered set of points in `N` dimensions, with an optional name.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Cloud<const N: usize> {
    /// Optional identifier, preserved through a round trip. Order is meaningful regardless.
    pub name: Option<String>,
    /// The points.
    pub points: Vec<[f64; N]>,
    /// Caller metadata, empty when the cloud carries none. Stored, never interpreted.
    pub metadata: Metadata,
}

/// A point cloud in the x, y plane.
pub type Cloud2 = Cloud<2>;

/// A point cloud in x, y, z space.
pub type Cloud3 = Cloud<3>;

impl<const N: usize> Cloud<N> {
    /// A cloud with no name and no metadata.
    pub fn new(points: Vec<[f64; N]>) -> Self {
        Self {
            name: None,
            points,
            metadata: Metadata::new(),
        }
    }

    /// The same cloud carrying a name.
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The same cloud carrying one metadata entry.
    pub fn with_meta(mut self, key: impl Into<String>, value: impl Into<crate::Value>) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }
}

impl<const N: usize> Named for Cloud<N> {
    fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
}

/// The container kind for a cloud of this dimension.
///
/// # Panics
///
/// Panics outside 2 and 3 dimensions, which have no assigned kind byte.
fn kind_for<const N: usize>() -> Kind {
    match N {
        2 => Kind::Cloud2,
        3 => Kind::Cloud3,
        _ => panic!("cloud containers exist only in 2 and 3 dimensions, not {N}"),
    }
}

/// Write a collection of clouds, every item at the same storage tolerance.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis is too wide to meet `tol`, and
/// [`Error::Malformed`] for a non-finite coordinate.
pub fn write_to<W: Write, const N: usize>(
    writer: &mut W,
    clouds: &[Cloud<N>],
    tol: f64,
) -> Result<()> {
    write_to_with_meta(writer, clouds, tol, &Metadata::new())
}

/// Write a collection of clouds with file-level metadata attached. See [`write_to`].
pub fn write_to_with_meta<W: Write, const N: usize>(
    writer: &mut W,
    clouds: &[Cloud<N>],
    tol: f64,
    file_metadata: &Metadata,
) -> Result<()> {
    let count = u32::try_from(clouds.len())
        .map_err(|_| Error::Malformed("container holds more items than a u32 can count"))?;
    container::write_header(writer, kind_for::<N>(), count, file_metadata)?;

    for cloud in clouds {
        write_item(writer, cloud, tol)?;
    }

    Ok(())
}

/// Write a single cloud as a one-item collection.
pub fn write_one_to<W: Write, const N: usize>(
    writer: &mut W,
    cloud: &Cloud<N>,
    tol: f64,
) -> Result<()> {
    write_to(writer, std::slice::from_ref(cloud), tol)
}

/// Read a collection of clouds.
///
/// # Errors
///
/// [`Error::Malformed`] if the file holds a different dimension than `N`.
pub fn read_from<R: Read, const N: usize>(reader: &mut R) -> Result<Vec<Cloud<N>>> {
    let header = container::read_header(reader, kind_for::<N>())?;

    let mut out = Vec::with_capacity((header.count as usize).min(MAX_PREALLOC));
    for _ in 0..header.count {
        out.push(read_item(reader)?);
    }

    Ok(out)
}

/// Read a container that holds exactly one cloud.
///
/// # Errors
///
/// [`Error::NotASingleItem`] if the container holds any other number.
pub fn read_one_from<R: Read, const N: usize>(reader: &mut R) -> Result<Cloud<N>> {
    let header = container::read_header(reader, kind_for::<N>())?;
    if header.count != 1 {
        return Err(Error::NotASingleItem {
            found: header.count,
        });
    }
    read_item(reader)
}

/// Write a collection of clouds to a file. See [`write_to`].
pub fn write_file<const N: usize>(path: &Path, clouds: &[Cloud<N>], tol: f64) -> Result<()> {
    // Buffering is required, not tidiness: the bit reader and writer move a byte at a time, so an
    // unbuffered File would mean a syscall per byte.
    let mut writer = BufWriter::new(File::create(path)?);
    write_to(&mut writer, clouds, tol)?;
    writer.flush()?;
    Ok(())
}

/// Write a single cloud to a file as a one-item collection.
pub fn write_one_file<const N: usize>(path: &Path, cloud: &Cloud<N>, tol: f64) -> Result<()> {
    write_file(path, std::slice::from_ref(cloud), tol)
}

/// Read a collection of clouds from a file.
pub fn read_file<const N: usize>(path: &Path) -> Result<Vec<Cloud<N>>> {
    let mut reader = BufReader::new(File::open(path)?);
    read_from(&mut reader)
}

/// Read a file holding exactly one cloud. See [`read_one_from`].
pub fn read_one_file<const N: usize>(path: &Path) -> Result<Cloud<N>> {
    let mut reader = BufReader::new(File::open(path)?);
    read_one_from(&mut reader)
}

fn write_item<W: Write, const N: usize>(writer: &mut W, cloud: &Cloud<N>, tol: f64) -> Result<()> {
    item::write_preamble(writer, cloud.name.as_deref(), &cloud.metadata, false)?;
    write_points(writer, &cloud.points, tol)?;
    Ok(())
}

fn read_item<R: Read, const N: usize>(reader: &mut R) -> Result<Cloud<N>> {
    let preamble = item::read_preamble(reader, false)?;
    let points: Vec<[f64; N]> = read_points(reader)?;

    Ok(Cloud {
        name: preamble.name,
        points,
        metadata: preamble.metadata,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corpus;
    use crate::find_by_name;
    use crate::testgen::Rng;
    use std::io::Cursor;

    fn distance<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
        (0..N).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    }

    fn assert_matches<const N: usize>(
        original: &Cloud<N>,
        recovered: &Cloud<N>,
        tol: f64,
        what: &str,
    ) {
        assert_eq!(recovered.name, original.name, "{what}: name");
        assert_eq!(recovered.metadata, original.metadata, "{what}: metadata");
        assert_eq!(
            recovered.points.len(),
            original.points.len(),
            "{what}: point count"
        );

        for (i, (o, r)) in original
            .points
            .iter()
            .zip(recovered.points.iter())
            .enumerate()
        {
            let d = distance(o, r);
            assert!(d <= tol, "{what}: point {i} recovered {d} away, tol {tol}");
        }
    }

    fn round_trip<const N: usize>(cloud: &Cloud<N>, tol: f64) -> Cloud<N> {
        let mut buf = Vec::new();
        write_one_to(&mut buf, cloud, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let back = read_one_from(&mut cursor).unwrap();
        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder left bytes unread"
        );
        back
    }

    #[test]
    fn every_corpus_point_set_round_trips() {
        for case in corpus::all() {
            let cloud = Cloud3::new(case.points.clone());
            assert_matches(&cloud, &round_trip(&cloud, case.tol), case.tol, case.name);
        }
    }

    #[test]
    fn two_dimensional_clouds_round_trip() {
        let mut rng = Rng::new(77);
        let cloud = Cloud2::new(rng.points(5_000, -20.0, 20.0)).named("plane");
        assert_matches(&cloud, &round_trip(&cloud, 1e-4), 1e-4, "2d");
    }

    /// A cloud is the mesh container minus its index block, so it must actually be smaller for the
    /// same points. If it is not, the kind is not earning its place.
    #[test]
    fn a_cloud_is_smaller_than_the_same_points_as_a_mesh() {
        let case = corpus::distant_clusters();

        let mut as_cloud = Vec::new();
        write_one_to(&mut as_cloud, &Cloud3::new(case.points.clone()), case.tol).unwrap();

        let mut as_mesh = Vec::new();
        crate::mesh::write_one_to(
            &mut as_mesh,
            &crate::Mesh3::new(case.points.clone(), Vec::new()),
            case.tol,
        )
        .unwrap();

        assert!(as_cloud.len() < as_mesh.len());
    }

    #[test]
    fn a_collection_preserves_order_and_names() {
        let mut rng = Rng::new(78);
        let items: Vec<Cloud3> = (0..5)
            .map(|i| {
                let c = Cloud3::new(rng.points(200, -5.0, 5.0));
                if i % 2 == 0 {
                    c.named(format!("pass-{i}"))
                } else {
                    c
                }
            })
            .collect();

        let mut buf = Vec::new();
        write_to(&mut buf, &items, 1e-4).unwrap();
        let back: Vec<Cloud3> = read_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(back.len(), 5);
        for (i, (o, r)) in items.iter().zip(back.iter()).enumerate() {
            assert_matches(o, r, 1e-4, &format!("item {i}"));
        }
        assert_eq!(find_by_name(&back, "pass-2").unwrap().points.len(), 200);
        assert!(find_by_name(&back, "pass-1").is_none());
    }

    #[test]
    fn metadata_round_trips() {
        let cloud = Cloud3::new(vec![[1.0, 2.0, 3.0]])
            .with_meta("units", "mm")
            .with_meta("acme.scanner", 7i64);

        let back = round_trip(&cloud, 1e-6);
        assert_eq!(back.metadata["units"].as_text(), Some("mm"));
        assert_eq!(back.metadata["acme.scanner"].as_i64(), Some(7));
    }

    #[test]
    fn degenerate_sizes_round_trip() {
        let tol = 1e-6;

        let empty = Cloud3::new(Vec::new());
        assert_matches(&empty, &round_trip(&empty, tol), tol, "empty");

        let single = Cloud2::new(vec![[1.5, -2.5]]);
        assert_matches(&single, &round_trip(&single, tol), tol, "single");
    }

    #[test]
    fn reading_at_the_wrong_dimension_is_rejected() {
        let mut buf = Vec::new();
        write_one_to(&mut buf, &Cloud2::new(vec![[0.0, 1.0]]), 1e-4).unwrap();

        let wrong: Result<Vec<Cloud3>> = read_from(&mut Cursor::new(&buf));
        assert!(matches!(wrong, Err(Error::Malformed(_))));
    }

    /// A cloud has no notion of closure, so the flag is refused rather than ignored.
    #[test]
    fn the_closed_flag_is_rejected() {
        let mut buf = Vec::new();
        write_one_to(&mut buf, &Cloud3::new(vec![[0.0; 3]]), 1e-4).unwrap();
        buf[12] |= item::CLOSED;

        let back: Result<Cloud3> = read_one_from(&mut Cursor::new(&buf));
        assert!(matches!(back, Err(Error::Malformed(_))));
    }

    #[test]
    fn read_one_refuses_a_multi_item_file() {
        let c = Cloud3::new(vec![[0.0; 3]]);
        let mut buf = Vec::new();
        write_to(&mut buf, &[c.clone(), c.clone()], 1e-4).unwrap();

        let back: Result<Cloud3> = read_one_from(&mut Cursor::new(&buf));
        assert!(matches!(back, Err(Error::NotASingleItem { found: 2 })));
    }

    #[test]
    fn round_trips_through_a_real_file() {
        let mut rng = Rng::new(79);
        let cloud = Cloud3::new(rng.points(2_000, 0.0, 50.0)).named("on disk");

        let path = std::env::temp_dir().join(format!(
            "tol-compress-cloud-{}-{:x}.tccloud3",
            std::process::id(),
            rng.next_u64()
        ));

        write_one_file(&path, &cloud, 1e-3).unwrap();
        let back: Cloud3 = read_one_file(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_matches(&cloud, &back, 1e-3, "file");
    }

    #[test]
    fn truncated_input_is_an_error() {
        let mut rng = Rng::new(80);
        let cloud = Cloud3::new(rng.points(500, 0.0, 10.0)).named("truncate me");

        let mut buf = Vec::new();
        write_one_to(&mut buf, &cloud, 1e-3).unwrap();

        for cut in [0, 5, 12, 13, 30, buf.len() / 2, buf.len() - 1] {
            let back: Result<Cloud3> = read_one_from(&mut Cursor::new(&buf[..cut]));
            assert!(back.is_err(), "truncating to {cut} bytes should fail");
        }
    }

    #[test]
    #[should_panic(expected = "only in 2 and 3 dimensions")]
    fn an_unsupported_dimension_panics() {
        let mut buf = Vec::new();
        let _ = write_one_to(&mut buf, &Cloud::<4>::new(Vec::new()), 1e-3);
    }
}

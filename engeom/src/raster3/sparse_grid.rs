//! A sparse, block-allocated 3D voxel grid, and the signed-distance payload that surface
//! reconstruction stores in one.
//!
//! # Why blocks rather than voxels
//!
//! Reconstruction only cares about voxels near a surface. For a scan of ten million points at a
//! voxel size close to the point spacing, the occupied shell is a few percent of the bounding
//! volume. A dense array would therefore waste most of its storage, while the overhead of a hash
//! entry for each voxel would exceed the eight-byte payload several times over.
//!
//! Blocks of `BLOCK_EDGE`³ voxels balance these costs. One hash lookup reaches 512 voxels, and
//! within a block the neighbors a stencil needs are adjacent in memory instead of scattered across
//! the map. The cost is over-allocation at the edges of the occupied region, where a block holds
//! some voxels that are never written. In a narrow band around a surface, approximately one-third
//! to one-half of the allocated voxels are typically filled. This over-allocation avoids the cost
//! of hashing each voxel.
//!
//! # Why the payload carries a weight
//!
//! [`SdfVoxel`] stores both a value and a weight, and a weight of zero means "unknown." The weight
//! costs four bytes per voxel, although a single-scan reconstruction never reads it back. Fusing
//! several scans into one field requires a weighted running mean over the same grid, so the voxel
//! must contain the accumulator from the beginning. This design permits a future multi-scan
//! populator without changes to the grid or the surface extraction that reads it.
//!
//! # Determinism
//!
//! Hash iteration order is not stable between runs, so anything that numbers voxels or vertices
//! must work from [`SparseGrid3::block_keys`], which sorts. `common::voxel_grid` makes the same
//! requirement for its own output. A mesh whose vertex order changes between runs cannot be
//! compared directly with a stored result.

use crate::common::PCoords;
use crate::{Aabb3, Point3, Result};
use faer::prelude::default;
use parry3d_f64::utils::hashmap::HashMap;

/// The number of voxels along one edge of a block.
pub const BLOCK_EDGE: usize = 8;

/// The number of voxels in a block, which is [`BLOCK_EDGE`] cubed.
pub const BLOCK_VOXELS: usize = BLOCK_EDGE * BLOCK_EDGE * BLOCK_EDGE;

const BLOCK_EDGE_I: i32 = BLOCK_EDGE as i32;

/// The integer coordinates of a single voxel in a [`SparseGrid3`].
///
/// A voxel key identifies a *sample position*, not a cell: the sample for key `k` sits at
/// `origin + k * voxel_size`. Marching cubes then treats the eight keys `k + {0,1}³` as the corners
/// of the cell whose lowest corner is `k`.
pub type VoxelKey3 = [i32; 3];

/// The integer coordinates of a block of voxels in a [`SparseGrid3`].
///
/// This is the voxel key divided by [`BLOCK_EDGE`], rounded toward negative infinity so that the
/// blocks tile the whole of space rather than mirroring around the origin.
pub type BlockKey3 = [i32; 3];

/// A value that can be stored in a [`SparseGrid3`] and knows whether it has been written.
///
/// The grid allocates a whole block at a time and fills it with `Default::default()`, so much of its
/// storage can be untouched. The default value must be unknown, which means that `is_known` must
/// return `false` for it. Only the payload type can distinguish written voxels from untouched ones
/// when the grid counts or trims voxels or extracts a surface.
pub trait VoxelValue: Copy + Default + Send + Sync {
    /// Whether this voxel has been written with a meaningful value.
    fn is_known(&self) -> bool;
}

/// A signed distance sample with the accumulated weight that produced it.
///
/// The value is the signed distance to the surface in the units of the cloud, positive on the side
/// the surface normals point toward. The weight is the total confidence behind it; zero means the
/// voxel was never written and the value is meaningless.
///
/// `f32` is enough for the value because a voxel is only ever written inside a narrow band, so the
/// magnitude stored is a few voxel edges at most and the relative precision of `f32` puts the
/// quantization orders of magnitude below any real measurement noise. Coordinates are never stored
/// this way; they are recomputed in `f64` from the integer key and the grid origin.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct SdfVoxel {
    /// The signed distance to the surface, meaningful only when `weight` is positive.
    pub value: f32,

    /// The total weight accumulated into `value`. Zero means unknown.
    pub weight: f32,
}

impl SdfVoxel {
    /// A voxel holding a single observation of `value` with confidence `weight`.
    ///
    /// A non-finite value or a non-positive or non-finite weight produces an unknown voxel and
    /// discards the observation. This prevents a degenerate observation from affecting a later
    /// [`SdfVoxel::fuse`].
    pub fn new(value: f64, weight: f64) -> Self {
        if !weight.is_finite() || weight <= 0.0 || !value.is_finite() {
            return Self::default();
        }

        Self {
            value: value as f32,
            weight: weight as f32,
        }
    }

    /// The signed distance, or `None` if this voxel is unknown.
    pub fn value(&self) -> Option<f64> {
        if self.is_known() {
            Some(self.value as f64)
        } else {
            None
        }
    }

    /// Fold another observation into this voxel as a running weighted mean.
    ///
    /// This is the update a multi-scan fusion performs: after any sequence of `fuse` calls the
    /// value is the weighted average of everything folded in, subject to `f32` rounding after each
    /// call. An observation with a non-finite value or a non-positive or non-finite weight is
    /// ignored, so a caller can pass through a rejected sample without branching.
    ///
    /// The accumulation runs in `f64` and is stored back to `f32`, so the rounding happens once per
    /// call rather than compounding through the intermediate.
    pub fn fuse(&mut self, value: f64, weight: f64) {
        if !weight.is_finite() || weight <= 0.0 || !value.is_finite() {
            return;
        }

        let w0 = self.weight as f64;
        let w1 = w0 + weight;
        let v = ((self.value as f64) * w0 + value * weight) / w1;

        self.value = v as f32;
        self.weight = w1 as f32;
    }
}

impl VoxelValue for SdfVoxel {
    fn is_known(&self) -> bool {
        self.weight > 0.0
    }
}

/// A dense cube of [`BLOCK_EDGE`]³ voxels, the unit a [`SparseGrid3`] allocates.
///
/// The voxels are stored in `x`-major order, so that varying `z` walks contiguous memory. A stencil
/// which sweeps `z` innermost therefore reads the block in order.
#[derive(Clone, Debug)]
pub struct Block3<V> {
    data: Box<[V; BLOCK_VOXELS]>,
}

impl<V: VoxelValue> Default for Block3<V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<V: VoxelValue> Block3<V> {
    /// A block with every voxel set to its default value, which is unknown as required by
    /// [`VoxelValue`].
    pub fn new() -> Self {
        Self {
            data: Box::new([V::default(); BLOCK_VOXELS]),
        }
    }

    /// The voxel at a local index, which must be less than [`BLOCK_VOXELS`].
    pub fn at(&self, local: usize) -> V {
        self.data[local]
    }

    /// A mutable reference to the voxel at a local index.
    pub fn at_mut(&mut self, local: usize) -> &mut V {
        &mut self.data[local]
    }

    /// Every voxel in the block, in local index order.
    pub fn as_slice(&self) -> &[V] {
        self.data.as_slice()
    }

    /// The number of voxels in this block which have been written.
    ///
    /// This counts on every call instead of tracking a running total. A cached count would become
    /// stale when a caller writes through [`SparseGrid3::get_mut_or_insert`].
    pub fn known_count(&self) -> usize {
        self.data.iter().filter(|v| v.is_known()).count()
    }

    /// Whether no voxel in this block has been written.
    pub fn has_no_known(&self) -> bool {
        !self.data.iter().any(|v| v.is_known())
    }

    /// The local index of a voxel from its position within the block, each coordinate in
    /// `0..BLOCK_EDGE`.
    ///
    /// This method defines the storage order. All code that walks a block by coordinate uses this
    /// method, which keeps the order consistent.
    pub fn local_index_of(x: usize, y: usize, z: usize) -> usize {
        (x * BLOCK_EDGE + y) * BLOCK_EDGE + z
    }
}

/// The local index of a voxel key within its own block.
///
/// `rem_euclid` makes a negative key land in `0..BLOCK_EDGE`. Using `%` could produce a negative
/// remainder and address a location outside the block or the wrong voxel.
fn local_index<V: VoxelValue>(key: &VoxelKey3) -> usize {
    Block3::<V>::local_index_of(
        key[0].rem_euclid(BLOCK_EDGE_I) as usize,
        key[1].rem_euclid(BLOCK_EDGE_I) as usize,
        key[2].rem_euclid(BLOCK_EDGE_I) as usize,
    )
}

/// The offset of a local index from the lowest corner of its block, inverting [`local_index`].
fn local_offset(local: usize) -> [i32; 3] {
    let z = local % BLOCK_EDGE;
    let y = (local / BLOCK_EDGE) % BLOCK_EDGE;
    let x = local / (BLOCK_EDGE * BLOCK_EDGE);
    [x as i32, y as i32, z as i32]
}

/// A sparse 3D grid of voxels, allocated in blocks of [`BLOCK_EDGE`]³ and addressed by integer key.
///
/// See the module documentation for why the storage is blocked and why the payload carries a
/// weight. The grid itself knows nothing about surfaces: it maps integer keys to values and
/// converts between keys and world positions. The populator determines what the grid stores, and
/// the extractor determines how to interpret those values.
#[derive(Clone, Debug)]
pub struct SparseGrid3<V> {
    voxel_size: f64,
    origin: Point3,
    blocks: HashMap<BlockKey3, Block3<V>>,
}

impl<V: VoxelValue> SparseGrid3<V> {
    /// An empty grid with the given voxel spacing and sample origin.
    ///
    /// The sample for key `[0, 0, 0]` sits at `origin`, and the sample for key `k` at
    /// `origin + k * voxel_size`. The origin usually remains at the world origin. Shift it to align
    /// the grid with another object; it does not center the grid on the data.
    ///
    /// # Arguments
    ///
    /// * `voxel_size`: The spacing between samples along each axis. Must be finite and positive.
    /// * `origin`: The world position of the sample with key `[0, 0, 0]`.
    ///
    /// returns: `Result<SparseGrid3<V>>`
    pub fn new(voxel_size: f64, origin: Point3) -> Result<Self> {
        if !voxel_size.is_finite() || voxel_size <= 0.0 {
            return Err(format!("Voxel size must be finite and positive, got {voxel_size}").into());
        }

        Ok(Self {
            voxel_size,
            origin,
            blocks: HashMap::with_hasher(default()),
        })
    }

    /// The spacing between samples along each axis.
    pub fn voxel_size(&self) -> f64 {
        self.voxel_size
    }

    /// The world position of the sample with key `[0, 0, 0]`.
    pub fn origin(&self) -> Point3 {
        self.origin
    }

    /// The key of the voxel sample at or below a world position on every axis.
    ///
    /// # A note on very large coordinates
    ///
    /// Keys are `i32`, so a coordinate beyond roughly `2e9 * voxel_size` from the origin saturates
    /// the cast and collides with everything else out that far. This matches the caveat on
    /// `common::voxel_grid` and is documented rather than checked, for the same reason: the check
    /// would cost a branch per coordinate on a path that runs once per point.
    pub fn key_of(&self, point: &impl PCoords<3>) -> VoxelKey3 {
        let c = point.coords();
        let mut key = [0i32; 3];
        for d in 0..3 {
            key[d] = ((c[d] - self.origin.coords[d]) / self.voxel_size).floor() as i32;
        }
        key
    }

    /// The world position of the sample with the given key.
    pub fn corner_position(&self, key: &VoxelKey3) -> Point3 {
        Point3::new(
            self.origin.x + key[0] as f64 * self.voxel_size,
            self.origin.y + key[1] as f64 * self.voxel_size,
            self.origin.z + key[2] as f64 * self.voxel_size,
        )
    }

    /// The key of the block that owns a voxel key.
    pub fn block_key(key: &VoxelKey3) -> BlockKey3 {
        [
            key[0].div_euclid(BLOCK_EDGE_I),
            key[1].div_euclid(BLOCK_EDGE_I),
            key[2].div_euclid(BLOCK_EDGE_I),
        ]
    }

    /// The voxel key of a local index within a block.
    pub fn voxel_key(block_key: &BlockKey3, local: usize) -> VoxelKey3 {
        let offset = local_offset(local);
        [
            block_key[0] * BLOCK_EDGE_I + offset[0],
            block_key[1] * BLOCK_EDGE_I + offset[1],
            block_key[2] * BLOCK_EDGE_I + offset[2],
        ]
    }

    /// The value at a voxel key, or `None` if its block was never allocated.
    ///
    /// An allocated but unwritten voxel comes back as `Some` holding a default value, which
    /// [`VoxelValue::is_known`] reports as unknown. Callers which only care whether there is a
    /// usable value should check that rather than the `Option`.
    pub fn get(&self, key: &VoxelKey3) -> Option<&V> {
        self.blocks
            .get(&Self::block_key(key))
            .map(|b| &b.data[local_index::<V>(key)])
    }

    /// The value at a voxel key, treating an unallocated or unwritten voxel alike as absent.
    pub fn get_known(&self, key: &VoxelKey3) -> Option<V> {
        self.get(key).copied().filter(|v| v.is_known())
    }

    /// A mutable reference to the value at a voxel key, allocating its block if needed.
    pub fn get_mut_or_insert(&mut self, key: &VoxelKey3) -> &mut V {
        let block = self
            .blocks
            .entry(Self::block_key(key))
            .or_insert_with(Block3::new);
        &mut block.data[local_index::<V>(key)]
    }

    /// The block with the given key, if it has been allocated.
    pub fn block(&self, key: &BlockKey3) -> Option<&Block3<V>> {
        self.blocks.get(key)
    }

    /// Put a block into the grid at the given key, replacing whatever was there.
    pub fn insert_block(&mut self, key: BlockKey3, block: Block3<V>) {
        self.blocks.insert(key, block);
    }

    /// Every allocated block key, in ascending order.
    ///
    /// The sort is what makes anything built from this grid reproducible, as the module
    /// documentation explains. It allocates and sorts on every call, so a caller that needs the
    /// order more than once should retain the result.
    pub fn block_keys(&self) -> Vec<BlockKey3> {
        let mut keys: Vec<BlockKey3> = self.blocks.keys().copied().collect();
        keys.sort_unstable();
        keys
    }

    /// The number of allocated blocks.
    pub fn block_count(&self) -> usize {
        self.blocks.len()
    }

    /// The number of allocated voxels, including written and unwritten voxels.
    pub fn allocated_count(&self) -> usize {
        self.blocks.len() * BLOCK_VOXELS
    }

    /// The number of voxels which have been written.
    ///
    /// This walks every allocated voxel, for the reason given on [`Block3::known_count`].
    pub fn known_count(&self) -> usize {
        self.blocks.values().map(|b| b.known_count()).sum()
    }

    /// Whether no voxel in the grid has been written.
    pub fn has_no_known(&self) -> bool {
        self.blocks.values().all(|b| b.has_no_known())
    }

    /// Drop every block in which nothing was written.
    ///
    /// Population activates blocks before it knows whether they contain samples inside the band.
    /// This method removes the unused blocks and returns the number removed.
    pub fn remove_empty_blocks(&mut self) -> usize {
        let before = self.blocks.len();
        self.blocks.retain(|_, b| !b.has_no_known());
        before - self.blocks.len()
    }

    /// The bounding box of the sample positions of every written voxel, or `None` if there are
    /// none.
    ///
    /// This bounds the *samples*, not the cells they are corners of, so it is smaller than the
    /// region the grid describes by up to one voxel on each face.
    pub fn compute_aabb(&self) -> Option<Aabb3> {
        let mut bounds: Option<(Point3, Point3)> = None;

        for (block_key, block) in self.blocks.iter() {
            for (local, value) in block.data.iter().enumerate() {
                if !value.is_known() {
                    continue;
                }

                let p = self.corner_position(&Self::voxel_key(block_key, local));
                bounds = Some(match bounds {
                    None => (p, p),
                    Some((lo, hi)) => (
                        Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
                        Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
                    ),
                });
            }
        }

        bounds.map(|(lo, hi)| Aabb3::new(lo, hi))
    }

    /// The keys of every block whose voxel samples could fall within `radius` of any of `points`.
    ///
    /// This is the candidate set for a narrow band. Computing it per block is inexpensive, but it
    /// can include a block whose voxels all prove too far from the data. This over-inclusion costs
    /// a field evaluation. An evaluation that finds no support leaves its voxel unknown, and
    /// [`SparseGrid3::remove_empty_blocks`] reclaims the unused allocation.
    ///
    /// The test is a sphere against the block's sample box, so a block that a point's sphere only
    /// passes near the corner of is correctly left out. The returned keys are sorted and free of
    /// duplicates. This method does not modify the grid. The populator allocates the blocks
    /// separately after it determines what to write.
    ///
    /// # Samples, not the region a block tiles
    ///
    /// What this measures against is the box spanned by a block's *samples*, which stops one voxel
    /// short of the region the block tiles, because the samples sit at the low corner of each
    /// voxel. A point in that last slab is therefore inside the block by
    /// [`SparseGrid3::key_of`] and still farther than zero from every sample in it. That is the
    /// correct result for a band defined over samples. Consequently, a radius of zero activates
    /// almost nothing.
    ///
    /// It also means a cell straddling two blocks can have some corners in an activated block and
    /// others outside them. A surface extractor that requires all eight corners will drop such a
    /// cell. A radius at least as large as the cell diagonal prevents this omission. The same
    /// condition makes the band wide enough to hold a crossing.
    ///
    /// # Arguments
    ///
    /// * `points`: The points whose neighborhoods define the band.
    /// * `radius`: The distance from a point within which a sample is worth evaluating. Must be
    ///   finite and non-negative. A radius much larger than a block edge makes the per-point
    ///   candidate scan proportionally more expensive.
    ///
    /// returns: `Result<Vec<BlockKey3>>`
    pub fn activate_blocks_near(
        &self,
        points: &[impl PCoords<3>],
        radius: f64,
    ) -> Result<Vec<BlockKey3>> {
        if !radius.is_finite() || radius < 0.0 {
            return Err(format!("Radius must be finite and non-negative, got {radius}").into());
        }

        let mut found: HashMap<BlockKey3, ()> = HashMap::with_hasher(default());

        for point in points.iter() {
            let c = point.coords();

            // The block range whose sample boxes can reach the sphere around this point. The upper
            // bound uses the same flooring as the lower so that both are block keys.
            let mut lo = [0i32; 3];
            let mut hi = [0i32; 3];
            for d in 0..3 {
                let rel = c[d] - self.origin.coords[d];
                lo[d] = ((rel - radius) / self.voxel_size).floor() as i32;
                hi[d] = ((rel + radius) / self.voxel_size).floor() as i32;
                lo[d] = lo[d].div_euclid(BLOCK_EDGE_I);
                hi[d] = hi[d].div_euclid(BLOCK_EDGE_I);
            }

            for bx in lo[0]..=hi[0] {
                for by in lo[1]..=hi[1] {
                    for bz in lo[2]..=hi[2] {
                        let key = [bx, by, bz];
                        if found.contains_key(&key) {
                            continue;
                        }
                        if self.block_sample_box(&key).distance_to_point(&c) <= radius {
                            found.insert(key, ());
                        }
                    }
                }
            }
        }

        let mut keys: Vec<BlockKey3> = found.into_keys().collect();
        keys.sort_unstable();
        Ok(keys)
    }

    /// The box spanned by the sample positions of a block, from its lowest corner sample to its
    /// highest. This is one voxel short of the region the block tiles, because the samples sit at
    /// the corners.
    fn block_sample_box(&self, key: &BlockKey3) -> SampleBox {
        let lo = self.corner_position(&[
            key[0] * BLOCK_EDGE_I,
            key[1] * BLOCK_EDGE_I,
            key[2] * BLOCK_EDGE_I,
        ]);
        let span = (BLOCK_EDGE_I - 1) as f64 * self.voxel_size;
        SampleBox {
            lo,
            hi: Point3::new(lo.x + span, lo.y + span, lo.z + span),
        }
    }
}

/// An axis-aligned box used only for the distance test in
/// [`SparseGrid3::activate_blocks_near`].
struct SampleBox {
    lo: Point3,
    hi: Point3,
}

impl SampleBox {
    /// The distance from the box to a point, or zero when the point is inside it.
    fn distance_to_point(&self, c: &crate::na::SVector<f64, 3>) -> f64 {
        let mut total = 0.0;
        for d in 0..3 {
            let v = c[d];
            let excess = if v < self.lo.coords[d] {
                self.lo.coords[d] - v
            } else if v > self.hi.coords[d] {
                v - self.hi.coords[d]
            } else {
                0.0
            };
            total += excess * excess;
        }
        total.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn grid(voxel_size: f64) -> SparseGrid3<SdfVoxel> {
        SparseGrid3::new(voxel_size, Point3::origin()).expect("grid creation failed")
    }

    // ===============================================================================================
    // Keys, positions, and the block/voxel mapping
    // ===============================================================================================

    #[test]
    fn grid_rejects_a_nonsense_voxel_size() {
        assert!(SparseGrid3::<SdfVoxel>::new(0.0, Point3::origin()).is_err());
        assert!(SparseGrid3::<SdfVoxel>::new(-1.0, Point3::origin()).is_err());
        assert!(SparseGrid3::<SdfVoxel>::new(f64::NAN, Point3::origin()).is_err());
        assert!(SparseGrid3::<SdfVoxel>::new(f64::INFINITY, Point3::origin()).is_err());
    }

    /// The sample for a key must map back to that key because all other grid operations depend on
    /// this relationship.
    #[test]
    fn key_of_a_corner_position_is_that_key() {
        let g = grid(0.25);

        for x in -20..20 {
            for y in -20..20 {
                for z in -20..20 {
                    let key = [x, y, z];
                    let p = g.corner_position(&key);
                    assert_eq!(g.key_of(&p), key, "round trip failed for {key:?}");
                }
            }
        }
    }

    /// A shifted origin must move the samples with it so that an aligned grid samples the intended
    /// positions.
    #[test]
    fn a_shifted_origin_moves_the_samples() {
        let origin = Point3::new(-1.3, 0.7, 4.25);
        let g = SparseGrid3::<SdfVoxel>::new(0.5, origin).expect("grid creation failed");

        assert_relative_eq!(g.corner_position(&[0, 0, 0]), origin, epsilon = 1e-12);
        assert_relative_eq!(
            g.corner_position(&[2, -1, 3]),
            Point3::new(origin.x + 1.0, origin.y - 0.5, origin.z + 1.5),
            epsilon = 1e-12
        );

        for key in [[0, 0, 0], [3, -4, 5], [-7, 8, -9]] {
            assert_eq!(g.key_of(&g.corner_position(&key)), key);
        }
    }

    /// `key_of` floors, so a position inside a voxel belongs to the sample at or below it, and a
    /// negative coordinate must floor away from zero rather than truncating toward it.
    #[test]
    fn key_of_floors_including_below_the_origin() {
        let g = grid(1.0);

        assert_eq!(g.key_of(&Point3::new(0.4, 0.9, 0.1)), [0, 0, 0]);
        assert_eq!(g.key_of(&Point3::new(1.0, 1.0, 1.0)), [1, 1, 1]);
        assert_eq!(g.key_of(&Point3::new(-0.1, -0.5, -0.9)), [-1, -1, -1]);
        assert_eq!(g.key_of(&Point3::new(-1.0, -1.0, -1.0)), [-1, -1, -1]);
        assert_eq!(g.key_of(&Point3::new(-1.1, -2.5, -0.001)), [-2, -3, -1]);
    }

    /// The block boundaries this pins are the ones a `%` instead of a `rem_euclid` would get
    /// wrong: the step from 7 to 8 and the step from 0 down to -1.
    #[test]
    fn block_keys_tile_space_across_zero() {
        assert_eq!(SparseGrid3::<SdfVoxel>::block_key(&[0, 0, 0]), [0, 0, 0]);
        assert_eq!(SparseGrid3::<SdfVoxel>::block_key(&[7, 7, 7]), [0, 0, 0]);
        assert_eq!(SparseGrid3::<SdfVoxel>::block_key(&[8, 8, 8]), [1, 1, 1]);
        assert_eq!(
            SparseGrid3::<SdfVoxel>::block_key(&[-1, -1, -1]),
            [-1, -1, -1]
        );
        assert_eq!(
            SparseGrid3::<SdfVoxel>::block_key(&[-8, -8, -8]),
            [-1, -1, -1]
        );
        assert_eq!(
            SparseGrid3::<SdfVoxel>::block_key(&[-9, -9, -9]),
            [-2, -2, -2]
        );
    }

    /// Every voxel key must map to a distinct local index in its block, and the inverse must
    /// recover the key. A collision here would have two voxels silently share storage.
    #[test]
    fn voxel_keys_and_local_indices_are_a_bijection_within_a_block() {
        for block_key in [[0, 0, 0], [1, -1, 2], [-3, -3, -3]] {
            let mut seen = vec![false; BLOCK_VOXELS];

            for (local, hit) in seen.iter_mut().enumerate() {
                let key = SparseGrid3::<SdfVoxel>::voxel_key(&block_key, local);

                assert_eq!(
                    SparseGrid3::<SdfVoxel>::block_key(&key),
                    block_key,
                    "local {local} of {block_key:?} escaped its own block"
                );
                assert_eq!(
                    local_index::<SdfVoxel>(&key),
                    local,
                    "local index round trip failed"
                );

                assert!(!*hit, "local index {local} was produced twice");
                *hit = true;
            }

            assert!(
                seen.iter().all(|&s| s),
                "some local index was never reached"
            );
        }
    }

    // ===============================================================================================
    // Storage
    // ===============================================================================================

    #[test]
    fn an_absent_voxel_reads_as_none() {
        let g = grid(1.0);
        assert!(g.get(&[0, 0, 0]).is_none());
        assert!(g.get(&[-5, 12, 3]).is_none());
        assert!(g.get_known(&[0, 0, 0]).is_none());
        assert_eq!(g.block_count(), 0);
        assert_eq!(g.known_count(), 0);
        assert!(g.has_no_known());
    }

    #[test]
    fn get_mut_or_insert_allocates_a_block_and_writes_through_it() {
        let mut g = grid(1.0);

        *g.get_mut_or_insert(&[3, -4, 5]) = SdfVoxel::new(0.5, 2.0);

        assert_eq!(g.block_count(), 1);
        assert_eq!(g.known_count(), 1);
        assert!(!g.has_no_known());
        assert_eq!(g.get_known(&[3, -4, 5]).and_then(|v| v.value()), Some(0.5));

        // A neighbor in the same block is allocated but unwritten, which is not the same as absent.
        let neighbor = g
            .get(&[3, -4, 6])
            .copied()
            .expect("block should be present");
        assert!(!neighbor.is_known());
        assert!(g.get_known(&[3, -4, 6]).is_none());
        assert_eq!(g.block_count(), 1);

        // A key in a different block allocates a second one.
        *g.get_mut_or_insert(&[100, 100, 100]) = SdfVoxel::new(-1.0, 1.0);
        assert_eq!(g.block_count(), 2);
        assert_eq!(g.known_count(), 2);
    }

    /// Writes to keys spread across block boundaries must each land in their own voxel and read
    /// back the value they were given.
    #[test]
    fn writes_across_many_blocks_read_back_independently() {
        let mut g = grid(0.5);

        let keys: Vec<VoxelKey3> = (-12i32..12)
            .map(|i| [i, (i * 7).rem_euclid(19) - 9, -i])
            .collect();

        for (n, key) in keys.iter().enumerate() {
            *g.get_mut_or_insert(key) = SdfVoxel::new(n as f64, 1.0);
        }

        assert_eq!(g.known_count(), keys.len());

        for (n, key) in keys.iter().enumerate() {
            assert_eq!(
                g.get_known(key).and_then(|v| v.value()),
                Some(n as f64),
                "wrong value at {key:?}"
            );
        }
    }

    #[test]
    fn insert_block_replaces_what_was_there() {
        let mut g = grid(1.0);

        *g.get_mut_or_insert(&[0, 0, 0]) = SdfVoxel::new(1.0, 1.0);
        assert_eq!(g.known_count(), 1);

        let mut block = Block3::<SdfVoxel>::new();
        *block.at_mut(0) = SdfVoxel::new(-3.0, 4.0);
        *block.at_mut(1) = SdfVoxel::new(-2.0, 4.0);
        g.insert_block([0, 0, 0], block);

        assert_eq!(g.block_count(), 1);
        assert_eq!(g.known_count(), 2);
        assert_eq!(g.get_known(&[0, 0, 0]).and_then(|v| v.value()), Some(-3.0));
    }

    #[test]
    fn remove_empty_blocks_keeps_only_what_was_written() {
        let mut g = grid(1.0);

        // Three blocks allocated, one of them written to.
        g.get_mut_or_insert(&[0, 0, 0]);
        g.get_mut_or_insert(&[64, 0, 0]);
        *g.get_mut_or_insert(&[-64, 0, 0]) = SdfVoxel::new(0.25, 1.0);

        assert_eq!(g.block_count(), 3);
        assert_eq!(g.allocated_count(), 3 * BLOCK_VOXELS);

        assert_eq!(g.remove_empty_blocks(), 2);
        assert_eq!(g.block_count(), 1);
        assert_eq!(g.known_count(), 1);
        assert_eq!(
            g.get_known(&[-64, 0, 0]).and_then(|v| v.value()),
            Some(0.25)
        );

        // Running it again has nothing left to do.
        assert_eq!(g.remove_empty_blocks(), 0);
        assert_eq!(g.block_count(), 1);
    }

    /// Block order out of the hash map is not stable, so `block_keys` sorts. Anything which
    /// numbers vertices from this order depends on it.
    #[test]
    fn block_keys_come_back_sorted_and_repeatable() {
        let build = || {
            let mut g = grid(1.0);
            for i in 0..40 {
                let k = [(i * 13) % 41 - 20, (i * 29) % 37 - 18, (i * 7) % 43 - 21];
                *g.get_mut_or_insert(&[k[0] * 8, k[1] * 8, k[2] * 8]) = SdfVoxel::new(1.0, 1.0);
            }
            g
        };

        let first = build().block_keys();
        let second = build().block_keys();

        let mut sorted = first.clone();
        sorted.sort_unstable();
        assert_eq!(first, sorted, "block keys were not in ascending order");
        assert_eq!(first, second, "block key order changed between builds");

        let mut deduped = first.clone();
        deduped.dedup();
        assert_eq!(deduped.len(), first.len(), "a block key appeared twice");
    }

    #[test]
    fn compute_aabb_bounds_the_written_samples() {
        let mut g = grid(0.5);
        assert!(g.compute_aabb().is_none());

        // Allocated but unwritten voxels must not count toward the bounds.
        g.get_mut_or_insert(&[100, 100, 100]);
        assert!(g.compute_aabb().is_none());

        *g.get_mut_or_insert(&[-2, 0, 4]) = SdfVoxel::new(0.0, 1.0);
        *g.get_mut_or_insert(&[6, -3, 1]) = SdfVoxel::new(0.0, 1.0);

        let aabb = g.compute_aabb().expect("bounds should exist");
        assert_relative_eq!(aabb.mins, Point3::new(-1.0, -1.5, 0.5), epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs, Point3::new(3.0, 0.0, 2.0), epsilon = 1e-12);
    }

    // ===============================================================================================
    // SdfVoxel
    // ===============================================================================================

    #[test]
    fn a_default_voxel_is_unknown() {
        let v = SdfVoxel::default();
        assert!(!v.is_known());
        assert!(v.value().is_none());
    }

    #[test]
    fn a_degenerate_observation_does_not_become_known() {
        for (value, weight) in [
            (1.0, 0.0),
            (1.0, -1.0),
            (1.0, f64::NAN),
            (f64::NAN, 1.0),
            (f64::INFINITY, 1.0),
        ] {
            let v = SdfVoxel::new(value, weight);
            assert!(!v.is_known(), "({value}, {weight}) should not be known");

            let mut fused = SdfVoxel::default();
            fused.fuse(value, weight);
            assert!(!fused.is_known(), "({value}, {weight}) should not fuse in");
        }
    }

    #[test]
    fn fusing_into_an_unknown_voxel_makes_it_known() {
        let mut v = SdfVoxel::default();
        v.fuse(2.5, 3.0);

        assert!(v.is_known());
        assert_relative_eq!(v.value().unwrap(), 2.5, epsilon = 1e-6);
        assert_relative_eq!(v.weight as f64, 3.0, epsilon = 1e-6);
    }

    /// Repeated fusing must equal the closed-form weighted mean within `f32` rounding, which is the
    /// property that multi-scan accumulation requires.
    #[test]
    fn fusing_equals_the_closed_form_weighted_mean() {
        let samples = [
            (1.0, 2.0),
            (-0.5, 1.0),
            (3.25, 4.0),
            (0.125, 0.5),
            (-2.0, 3.5),
        ];

        let mut v = SdfVoxel::default();
        for (value, weight) in samples {
            v.fuse(value, weight);
        }

        let total_weight: f64 = samples.iter().map(|(_, w)| w).sum();
        let expected: f64 = samples.iter().map(|(x, w)| x * w).sum::<f64>() / total_weight;

        assert_relative_eq!(v.value().unwrap(), expected, epsilon = 1e-6);
        assert_relative_eq!(v.weight as f64, total_weight, epsilon = 1e-6);
    }

    /// Accumulating the samples in reverse order must produce the same weighted mean within `f32`
    /// rounding.
    #[test]
    fn fusing_is_order_independent() {
        let samples = [(1.0, 2.0), (-0.5, 1.0), (3.25, 4.0), (0.125, 0.5)];

        let mut forward = SdfVoxel::default();
        for (value, weight) in samples {
            forward.fuse(value, weight);
        }

        let mut backward = SdfVoxel::default();
        for (value, weight) in samples.iter().rev() {
            backward.fuse(*value, *weight);
        }

        assert_relative_eq!(
            forward.value().unwrap(),
            backward.value().unwrap(),
            epsilon = 1e-6
        );
        assert_relative_eq!(
            forward.weight as f64,
            backward.weight as f64,
            epsilon = 1e-6
        );
    }

    // ===============================================================================================
    // Block activation
    // ===============================================================================================

    #[test]
    fn activation_rejects_a_nonsense_radius() {
        let g = grid(1.0);
        let points = [Point3::origin()];
        assert!(g.activate_blocks_near(&points, -1.0).is_err());
        assert!(g.activate_blocks_near(&points, f64::NAN).is_err());
        assert!(g.activate_blocks_near(&points, f64::INFINITY).is_err());
    }

    /// The candidate set must be the blocks whose sample box comes within the radius of some
    /// point, checked against a brute-force sweep of every block in a generous surrounding range.
    #[test]
    fn activation_matches_brute_force() {
        let g = grid(0.5);

        let points: Vec<Point3> = (0..25)
            .map(|i| {
                let t = i as f64 * 0.83;
                Point3::new(t.sin() * 6.0, t.cos() * 6.0, (t * 0.4).sin() * 3.0 - 1.0)
            })
            .collect();

        for radius in [0.0, 0.3, 1.0, 2.5, 6.0] {
            let actual = g
                .activate_blocks_near(&points, radius)
                .expect("activation failed");

            // Brute force over a range which comfortably contains anything reachable, testing the
            // sphere against the block's sample box independently of the code under test.
            let span = (BLOCK_EDGE_I - 1) as f64 * g.voxel_size();
            let mut expected: Vec<BlockKey3> = Vec::new();
            for bx in -8..8 {
                for by in -8..8 {
                    for bz in -8..8 {
                        let key = [bx, by, bz];
                        let lo = g.corner_position(&[
                            bx * BLOCK_EDGE_I,
                            by * BLOCK_EDGE_I,
                            bz * BLOCK_EDGE_I,
                        ]);
                        let hi = Point3::new(lo.x + span, lo.y + span, lo.z + span);

                        let hit = points.iter().any(|p| {
                            let mut sum = 0.0;
                            for d in 0..3 {
                                let v = p.coords[d];
                                let excess = (lo.coords[d] - v).max(v - hi.coords[d]).max(0.0);
                                sum += excess * excess;
                            }
                            sum.sqrt() <= radius
                        });

                        if hit {
                            expected.push(key);
                        }
                    }
                }
            }
            expected.sort_unstable();

            assert_eq!(
                actual, expected,
                "activation disagreed with brute force at radius {radius}"
            );
        }
    }

    /// The contract the narrow band rests on: no sample within the radius of a point may be left
    /// out of the activated set. Over-inclusion is allowed and costs only a wasted evaluation, but
    /// a sample that should have been evaluated and was not leaves a hole in the band.
    ///
    /// This sweeps every voxel key in a range around the data and checks each one individually,
    /// rather than checking blocks, so a block dropped for being a near miss is caught by the
    /// samples inside it.
    #[test]
    fn activation_omits_no_sample_within_the_radius() {
        let g = grid(0.25);

        let points: Vec<Point3> = (0..40)
            .map(|i| {
                let t = i as f64 * 0.41;
                Point3::new(t.cos() * 3.0, t * 0.15 - 2.0, t.sin() * 3.0)
            })
            .collect();

        for radius in [0.0, 0.25, 0.75, 2.0] {
            let keys = g
                .activate_blocks_near(&points, radius)
                .expect("activation failed");

            let mut checked = 0;
            for kx in -20..20 {
                for ky in -20..20 {
                    for kz in -20..20 {
                        let key = [kx, ky, kz];
                        let sample = g.corner_position(&key);

                        let near = points.iter().any(|p| (p - sample).norm() <= radius + 1e-12);

                        if near {
                            checked += 1;
                            let block = SparseGrid3::<SdfVoxel>::block_key(&key);
                            assert!(
                                keys.binary_search(&block).is_ok(),
                                "sample {key:?} is within {radius} of a point but its block \
                                 {block:?} was not activated"
                            );
                        }
                    }
                }
            }

            if radius > 0.0 {
                assert!(
                    checked > 0,
                    "radius {radius} put no sample in range, so the check proved nothing"
                );
            }
        }
    }

    /// A point sitting in the last voxel slab of a block is farther than zero from every sample in
    /// that block, so a zero radius correctly leaves the block out. This pins the behavior the
    /// documentation describes. This behavior follows from the distinction between samples and
    /// cells.
    #[test]
    fn activation_is_about_samples_rather_than_the_region_a_block_tiles() {
        let g = grid(1.0);

        // Block [0,0,0] holds samples at 0..=7 on each axis. A point at 7.5 is inside the region
        // the block tiles but half a voxel past its last sample.
        let point = [Point3::new(7.5, 0.0, 0.0)];
        assert_eq!(
            SparseGrid3::<SdfVoxel>::block_key(&g.key_of(&point[0])),
            [0, 0, 0]
        );

        let none = g
            .activate_blocks_near(&point, 0.0)
            .expect("activation failed");
        assert!(!none.contains(&[0, 0, 0]));

        // Reaching back to the sample at 7.0 brings it in.
        let some = g
            .activate_blocks_near(&point, 0.5)
            .expect("activation failed");
        assert!(some.contains(&[0, 0, 0]));
    }

    #[test]
    fn activation_of_nothing_is_nothing() {
        let g = grid(1.0);
        let points: Vec<Point3> = Vec::new();
        assert!(
            g.activate_blocks_near(&points, 5.0)
                .expect("activation failed")
                .is_empty()
        );
    }
}

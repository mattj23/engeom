# Choosing a k-d tree backend, and where the tree should live

Measured 2026-08-10 on branch `points`, originally to answer these two questions:

1. **Should a point cloud own its k-d tree, or hand out a borrowed index that is cheap to throw away
   and rebuild?** This decides the shape of the `PointCloud` overhaul.
2. **Is `kiddo` worth returning to?** `engeom` used it until Sept 2025 (commit `236a0c8`, kiddo
   5.2.1) and left it over a correctness bug.

Both ask what a tree build costs relative to the work you would already be doing to the cloud.

## Summary

- **Breakeven is about 0.4N nearest-one queries**, invariant across backend, size and input order.
  A tree pays for itself once you have queried roughly a third of the cloud once.
- **So the build is cheap relative to *using* the tree and expensive relative to *not* using it.** Both
  halves say: never build unless asked. The cloud should not own an eagerly built tree.
- **kiddo 6.0.1 passes the correctness gate** and is 4x to 14x faster to build and 2x to 4x faster
  to query. At sizes where every backend is serial it is still 4.2x to 5.0x on build and 3.7x on
  query, so the margin is not parallelism.
- **`kdtree`'s build degrades badly on monotonically ordered input** (9.4x at 1e6), but real
  mesh-sampled data does not trigger the worst case (6%). A latent hazard, not a live catastrophe.

## Method

Points are surface-distributed rather than uniform in a volume, because every cloud this library
sees is a scan of a surface and a k-d tree over a 2-manifold embedded in 3-space behaves differently
from one over a solid. Two generators: a Fibonacci spiral on a sphere at exact N, and a dense sample
of `tests/data/engine-blade.tcmesh` at 0.25 mm spacing (2,057,975 points) as the realistic
cross-check.

Every set is measured twice, in the generator's natural spatial order and shuffled by a
deterministic Fisher-Yates. **The control is the point of the experiment.** A backend that inserts
points one at a time and splits as it goes can build a badly unbalanced tree when input arrives in
spatial order, and both generators are spatially ordered, as is real scan data arriving along scan
lines. Without the shuffled case beside it there is no way to separate "this backend is slow" from
"this backend is slow on ordered input", and only the second is fixable with a line in our own
wrapper.

Timings are best-of-3, with repetition skipped once a single run exceeds 2 s. Best-of rather than
mean because the noise is one-sided: scheduling, turbo ramp and cache state can only make a run
slower than the work takes. A first pass timed everything once and the O(N) transform pass came out
1.8x apart on identical data between two runs.

Backends: `kdtree` 0.7.0 (what `engeom` ships), `kdtree` 0.8.1, and `kiddo` 6.0.1
`ImmutableKdTree`. Run on a 20-core machine. **kiddo parallelizes construction above 262,144
points**, so the 1e6 and 2e6 rows hand it 20 cores while both `kdtree` builds stay serial; the 1e4
and 1e5 rows are serial for all three.

## Results

Times in milliseconds. Query columns are 10,000 queries in total.

| set | order | backend | build | nearest-one | nearest-7 | within | breakeven q |
|---|---|---|--:|--:|--:|--:|--:|
| sphere 1e4 | ordered | kdtree 0.7 | 12.87 | 32.38 | 58.19 | 38.95 | 3975 |
| sphere 1e4 | ordered | kdtree 0.8.1 | 8.73 | 23.02 | 43.22 | 36.44 | 3791 |
| sphere 1e4 | ordered | kiddo | 1.40 | 4.03 | 10.44 | 8.47 | 3483 |
| sphere 1e4 | shuffled | kdtree 0.7 | 4.80 | 11.99 | 23.64 | 21.04 | 4000 |
| sphere 1e4 | shuffled | kdtree 0.8.1 | 4.80 | 11.92 | 23.54 | 19.49 | 4028 |
| sphere 1e4 | shuffled | kiddo | 1.15 | 3.19 | 11.06 | 9.10 | 3598 |
| sphere 1e5 | ordered | kdtree 0.7 | 164.67 | 43.09 | 69.85 | 61.68 | 38213 |
| sphere 1e5 | ordered | kdtree 0.8.1 | 139.99 | 34.58 | 57.39 | 48.53 | 40482 |
| sphere 1e5 | ordered | kiddo | 11.75 | 2.83 | 8.69 | 6.75 | 41517 |
| sphere 1e5 | shuffled | kdtree 0.7 | 50.64 | 11.81 | 21.86 | 19.78 | 42867 |
| sphere 1e5 | shuffled | kdtree 0.8.1 | 50.35 | 11.89 | 21.44 | 18.38 | 42338 |
| sphere 1e5 | shuffled | kiddo | 10.09 | 3.20 | 9.41 | 7.63 | 31526 |
| sphere 1e6 | ordered | kdtree 0.7 | 11590.36 | 283.92 | 450.65 | 465.41 | 408224 |
| sphere 1e6 | ordered | kdtree 0.8.1 | 10827.21 | 272.85 | 436.67 | 428.34 | 396825 |
| sphere 1e6 | ordered | kiddo | 212.47 | 6.14 | 19.13 | 17.18 | 346263 |
| sphere 1e6 | shuffled | kdtree 0.7 | 1228.01 | 32.35 | 57.40 | 54.58 | 379590 |
| sphere 1e6 | shuffled | kdtree 0.8.1 | 1102.22 | 27.36 | 52.51 | 49.46 | 402799 |
| sphere 1e6 | shuffled | kiddo | 443.47 | 11.87 | 36.23 | 31.45 | 373740 |
| blade 2.06M | ordered | kdtree 0.7 | 3514.83 | 100.43 | 144.99 | 391.83 | 349977 |
| blade 2.06M | ordered | kdtree 0.8.1 | 3159.67 | 104.50 | 150.80 | 356.99 | 302353 |
| blade 2.06M | ordered | kiddo | 257.89 | 21.88 | 49.27 | 138.99 | 117841 |
| blade 2.06M | shuffled | kdtree 0.7 | 3318.33 | 44.64 | 70.22 | 266.05 | 743274 |
| blade 2.06M | shuffled | kdtree 0.8.1 | 3264.89 | 44.28 | 65.70 | 210.38 | 737283 |
| blade 2.06M | shuffled | kiddo | 950.23 | 22.70 | 53.11 | 149.16 | 418544 |

## 1. Breakeven is about 0.4N, and that decides the ownership question

The number of nearest-one queries needed to pay off a build lands between 0.32N and 0.42N in twenty
of the twenty-four rows, and it barely moves across backends that differ by an order of magnitude in
absolute speed. Faster builds come with proportionally faster queries.

That invariance answers the ownership question in both directions at once.

**Cheap relative to using the tree.** Every real consumer is far past 0.4N. Normal estimation
queries every point once. An alignment queries every sample point on every iteration, so a
2000-sample solve over twenty iterations against a 100k cloud is 40,000 queries against a breakeven
of 40,000, and a full robust solve with refinement rounds is several times that. Throwing an index
away and rebuilding it to escape an ownership problem therefore costs a fraction of one
use-episode.

**Expensive relative to not using it.** An eager build on the 2.06M-point blade costs 258 ms with the
best backend and 3.5 s with the one we ship, all of it wasted for a load, downsample and save that
never queries anything. Point clouds in this library come from lptf3 sensor loads and PLY files and
are routinely large.

**Conclusion: never build unless asked.** A `PointCloud3` should not own an eagerly built tree the
way `Mesh3` owns its BVH. Note that `Mesh3` holds its BVH internally because parry forces it, since
`TriMesh::new` unconditionally calls `rebuild_bvh()`, not because engeom chose that shape. The mesh
module already contains the other pattern for accelerators engeom does control: `MeshEdges<'a>` and
`MeshNav<'a>` are borrowed views, built on demand and never cached. The cloud should follow those.

A borrowed index has a second advantage that is not about speed. A type holding `points` and `tree`
as sibling fields needs every `&mut self` method to remember to rebuild, which is the same
forgetting bug as the `tree_uuid` token the `MeshData3` design rejected. A borrow makes staleness a
compile error instead.

## 2. kiddo is faster on both axes, and it is not just parallelism

At 1e4 and 1e5 every backend is serial, and kiddo still builds 4.2x to 5.0x faster and queries 3.7x
faster than `kdtree` 0.7 on shuffled input. The gap widens at 1e6 and on the blade, where kiddo gets
20 cores for construction, so those rows overstate the algorithmic difference.

`kdtree` 0.8.1 is a modest improvement over 0.7.0 at small N and roughly a wash at large N. It is
not an alternative to switching.

## 3. The insertion-order control moderated the headline

On the sphere, `kdtree` builds 9.4x slower from ordered input at 1e6 (11,590 ms against 1228 ms),
and the penalty grows with N: 2.7x at 1e4, 3.3x at 1e5, 9.4x at 1e6. Growth with N is the signature
of a tree that is degenerating rather than merely unlucky.

But the Fibonacci spiral marches monotonically pole to pole, which is the worst case for an
axis-split incremental insert, and **the real mesh-sampled blade barely shows it: 3515 ms ordered
against 3318 ms shuffled, a 6% difference.** Mesh face order is not monotone along any axis.

So the first reading of these numbers, that every tree `engeom` builds today is degenerate, was too
strong. It is a latent hazard on monotone input rather than a live catastrophe, and an lptf3
scan-line load is the most likely candidate to hit it. Query performance does still suffer on the
blade (2.25x on nearest-one), so tree quality is degraded even where build time is not.

A pre-insert shuffle would have been the fix worth considering if we were keeping `kdtree`. A
bulk-loading backend gets order-independence for free, because a median split computed over the
whole set cannot depend on the order it was handed.

## Caveats

- **The `build/xform` ratio reported by the harness is contaminated** and is not reproduced above.
  Best-of-3 hands the O(N) transform pass a warm-cache best case while the build gets no such
  benefit. Breakeven is the sound metric, being a ratio of two quantities measured the same way.
- **Query timings are not comparable between the ordered and shuffled rows of the same set.** The
  harness draws its query sequence by stride from the point array, so reordering the array reorders
  the queries too, mixing tree quality with query locality. Build timings have no such problem, and
  the ordered-versus-shuffled claims above rest on build times.
- **One passing test is not proof of correctness.** The original kiddo bug surfaced in real use, not
  under a targeted test.

## What was kept

`engeom/benches/kd_tree.rs` measures the shipped backend through the `KdTree3` wrapper and carries
no extra dependencies, so it keeps working across a backend swap and is the before-and-after
instrument for one.

The comparison harness itself was not retained. It existed to answer the two questions above, it
depended on both candidate crates being present as dev-dependencies, and a future candidate would
need new code against its own API regardless. The method section above is written to be enough to
rebuild it.

`check_kiddo_bug` in `geom3/mesh/sampling.rs` is the standing correctness gate on whatever backend
is in the tree. It Poisson-samples a sphere at radius `r` and asserts that a radius-`r` query centred
on any sample returns that point and nothing else, which is the shape of the 2025 bug. It should
keep its historical name, since that name is the only record in the source of why the crate changed.

# Benchmark fixtures

Real measured geometry, kept small and normalized to one simple layout so a hand-rolled reader can
load them without pulling in a PLY dependency. Both are `binary_little_endian 1.0` with exactly
three `float` vertex properties (`x`, `y`, `z`) and triangular faces.

These are excluded from the published crate (`exclude` in `Cargo.toml`); they are for developing
the format, not for using it.

| file | verts | faces | units | provenance |
| --- | --: | --: | --- | --- |
| `bunny.ply` | 453 | 948 | m | Stanford bunny, from `engeom/tests/data/bun_zipper_res4.ply`. Re-encoded from ASCII to binary and stripped of its `confidence` and `intensity` properties. |
| `scan-chunk.ply` | 8091 | 15539 | mm | A cropped region of `engeom/tests/data/sample-clip.ply`, real scanner output. Vertices outside the crop box were dropped, faces referencing them removed, and the survivors reindexed. |

Both store coordinates as `f32`, so their own resolution is roughly seven significant digits.
Encoding either at a tolerance below that measures the source file's rounding rather than ours,
which is why the test sweeps stop where they do: 1e-6 m for the bunny, 1e-4 mm for the scan.

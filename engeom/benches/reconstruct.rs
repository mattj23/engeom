//! Where the time goes when a point cloud is turned into a mesh.
//!
//! The reconstruction pipeline runs five distinct stages: two parallel stages over points, two
//! parallel stages over voxel blocks, and a final sequence of serial mesh operations. The code
//! identifies these stages, while this benchmark measures their cost. It determines which stage
//! dominates, how effectively the parallel stages scale, whether serial setup limits that scaling,
//! and how many field evaluations occur at voxels far from the surface.
//!
//! The benchmark does not modify or instrument the library. Every pipeline stage is public, so the
//! benchmark invokes each stage independently without adding instrumentation to the hot path.
//!
//! # Reading this
//!
//! Start with the fixture summary printed before the Criterion groups. Its hit rate is the fraction
//! of field evaluations that produced a value. Activation uses blocks eight voxels across, while
//! the band is approximately two voxels thick. A low hit rate therefore indicates that most of the
//! dominant stage evaluates voxels that cannot contribute to the result.
//!
//! The thread-scaling group measures parallel efficiency. If a stage's eight-thread time is far
//! above one-eighth of its one-thread time, it contains a significant serial section or has a work
//! distribution problem. The shape of the scaling curve distinguishes these causes.
//!
//! # What it found the first time it was run
//!
//! The initial results differed from expectations based on code inspection. With half a million
//! points, field sampling took 3 percent of the total time despite running a tree query for every
//! voxel corner. `repair_buffers` took 80 percent and grew faster than the face count: five times as
//! many faces required six and a half times as much time. The Parry bounding-volume build was the
//! second-largest cost at 13 percent.
//!
//! Both operations are serial and outside the meshing code. The measured band hit rate was 26 to 29
//! percent, close to the prediction, but eliminating all unsuccessful field evaluations would save
//! only approximately 2 percent of the total time. These results are recorded to prevent repeated
//! investigation: converting the triangles into a manifold is more expensive than generating them.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use engeom::geom3::half_edge3::{RepairOpts, repair_buffers};
use engeom::geom3::mesh::{MeshEditor, PatchFilter};
use engeom::geom3::point_cloud::implicit_field::{ImlsField3, compute_sdf_grid};
use engeom::io::read_tc_mesh_from;
use engeom::raster3::{BLOCK_VOXELS, extract_isosurface};
use engeom::{
    Mesh3, MeshData3, NormalOrientation3, NormalSource3, Point3, PointCloud3, ReconstructOpts3,
    UnitVec3, Vector3,
};
use std::collections::HashMap;
use std::hint::black_box;
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

/// Point counts for the scaling group. Reconstruction cost depends on both the point count and the
/// surface area covered by the band, so the voxel size follows the implied point spacing instead of
/// remaining fixed.
const SIZES: [usize; 3] = [20_000, 100_000, 500_000];

/// Build fixtures above this point count only when `ENGEOM_BENCH_BIG` is set.
const BIG_THRESHOLD: usize = 100_000;

/// Thread counts for the scaling group. The original benchmark machine had 20 cores; each run uses
/// only counts supported by the current machine.
const THREADS: [usize; 5] = [1, 2, 4, 8, 20];

const SPHERE_RADIUS: f64 = 100.0;

/// Points on an analytic sphere with their mathematically exact normals.
///
/// A tessellated sphere sampled with `sample_poisson` would also work. This construction is
/// deterministic, inexpensive, and permits the scaling group to select an exact point count.
fn fibonacci_sphere(n: usize, radius: f64) -> (Vec<Point3>, Vec<UnitVec3>) {
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());

    (0..n)
        .map(|i| {
            let y = 1.0 - (i as f64 / (n as f64 - 1.0)) * 2.0;
            let r = (1.0 - y * y).max(0.0).sqrt();
            let theta = golden * i as f64;
            let u = Vector3::new(theta.cos() * r, y, theta.sin() * r);
            (Point3::from(u * radius), UnitVec3::new_unchecked(u))
        })
        .unzip()
}

/// Estimate the spacing of `n` points distributed over a sphere of the specified radius for use as
/// the voxel size.
fn sphere_spacing(n: usize, radius: f64) -> f64 {
    let area = 4.0 * std::f64::consts::PI * radius * radius;
    (area / n as f64).sqrt()
}

/// A cloud carrying oriented normals, ready for the default reconstruction path.
fn sphere_cloud(n: usize) -> PointCloud3 {
    let (points, normals) = fibonacci_sphere(n, SPHERE_RADIUS);
    let mut cloud = PointCloud3::new(points);
    cloud
        .set_point_normals(Some(normals))
        .expect("setting normals failed");
    cloud
}

/// A realistic fixture containing a scan-like dense sample of a real shape, including creases and
/// thin features absent from the sphere.
///
/// The crate's `tests` helpers use `cfg(test)` and are unavailable to a benchmark because it is a
/// separate compilation unit. The benchmark therefore embeds the asset directly.
fn bunny_cloud(spacing: f64) -> PointCloud3 {
    let bytes = include_bytes!("../tests/data/stanford_bun_3.tcmesh");
    let data = read_tc_mesh_from(&mut { bytes.as_slice() }).expect("reading the bunny failed");
    let mesh = Mesh3::from_data(data, false).expect("building the bunny failed");

    mesh.sample_dense(spacing, None)
        .expect("dense sampling failed")
}

struct Fixture {
    name: String,
    cloud: PointCloud3,
    voxel_size: f64,
}

/// Include the largest fixture only when requested because each pipeline run takes approximately 20
/// seconds.
fn want_big() -> bool {
    std::env::var_os("ENGEOM_BENCH_BIG").is_some()
}

/// Fixtures, built once per process.
///
/// Criterion filters the benchmarks that it measures, but it still executes every group body when
/// a filter excludes all benchmarks in that group. Cache the fixtures so a filtered run does not
/// rebuild each fixture several times and spend more time on setup than measurement.
fn fixtures() -> &'static [Fixture] {
    static ALL: OnceLock<Vec<Fixture>> = OnceLock::new();

    ALL.get_or_init(|| {
        let mut out = Vec::new();

        for &n in SIZES.iter() {
            if n > BIG_THRESHOLD && !want_big() {
                continue;
            }
            out.push(Fixture {
                name: format!("sphere/{n}"),
                cloud: sphere_cloud(n),
                voxel_size: sphere_spacing(n, SPHERE_RADIUS),
            });
        }

        let bunny_spacing = 0.002;
        out.push(Fixture {
            name: "bunny".to_string(),
            cloud: bunny_cloud(bunny_spacing),
            voxel_size: bunny_spacing,
        });

        out
    })
}

/// Return the largest sphere fixture that is always present for use by single-size groups.
fn primary() -> &'static Fixture {
    let all = fixtures();
    all.iter()
        .rfind(|f| f.name.starts_with("sphere/"))
        .unwrap_or(&all[0])
}

/// Return prepared data for one fixture, cached by fixture name to avoid repeated setup.
fn prepared(f: &Fixture) -> &'static Prepared {
    static CACHE: OnceLock<Mutex<HashMap<String, &'static Prepared>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));

    let mut guard = cache.lock().expect("prepared cache poisoned");
    if let Some(p) = guard.get(&f.name) {
        return p;
    }

    let leaked: &'static Prepared = Box::leak(Box::new(prepare(f)));
    guard.insert(f.name.clone(), leaked);
    leaked
}

/// One untimed run per fixture, printed as two tables.
///
/// These are single samples suitable only for comparing proportions. The Criterion groups below
/// provide the measurements.
///
/// One pipeline run per fixture supplies both tables. Previously, separate runs for counts and
/// timings formed the largest part of a 67-second fixed cost that every filtered invocation paid
/// before taking a measurement.
fn print_fixture_summary() {
    let runs: Vec<(&Fixture, StageRun)> = fixtures().iter().map(|f| (f, time_stages(f))).collect();

    println!();
    println!("=== fixtures ===");
    println!(
        "{:<14} {:>9} {:>9} {:>8} {:>10} {:>9} {:>9} {:>9} {:>9}",
        "fixture", "points", "voxel", "blocks", "known", "hit", "cells", "skipped", "faces"
    );

    for (f, run) in runs.iter() {
        let attempted = (run.active_blocks * BLOCK_VOXELS) as f64;
        let hit = if attempted > 0.0 {
            run.known_voxels as f64 / attempted
        } else {
            0.0
        };

        println!(
            "{:<14} {:>9} {:>9.4} {:>8} {:>10} {:>8.1}% {:>9} {:>9} {:>9}",
            f.name,
            f.cloud.point_count(),
            f.voxel_size,
            run.active_blocks,
            run.known_voxels,
            hit * 100.0,
            run.cells_visited,
            run.cells_skipped,
            run.faces,
        );
    }

    println!();
    println!("=== stage breakdown, single untimed run ===");
    println!(
        "{:<14} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}",
        "fixture", "normals", "field", "extract", "repair", "(full)", "(bvh)", "total"
    );

    for (f, run) in runs.iter() {
        println!(
            "{:<14} {:>8.0}ms {:>8.0}ms {:>8.0}ms {:>8.0}ms {:>8.0}ms {:>8.0}ms {:>8.0}ms",
            f.name,
            run.normals * 1e3,
            run.field * 1e3,
            run.extract * 1e3,
            run.repair * 1e3,
            run.repair_full * 1e3,
            run.mesh * 1e3,
            run.total * 1e3,
        );
    }

    println!();
    println!(
        "repair is the pipeline default, assuming_oriented_edges; (full) is RepairOpts::default. \
         (bvh) is the caller converting the result to a Mesh3. Neither is counted in the total."
    );
    println!();
}

/// Timings and counts from one hand-walked pipeline run.
struct StageRun {
    normals: f64,
    field: f64,
    extract: f64,
    repair: f64,
    repair_full: f64,
    mesh: f64,
    total: f64,
    active_blocks: usize,
    known_voxels: usize,
    cells_visited: usize,
    cells_skipped: usize,
    faces: usize,
}

/// Walk the pipeline by hand, timing each stage.
///
/// The normals stage uses estimation because the supplied-normal path has no work when the cloud
/// already contains normals and would not provide a useful timing proportion. This choice affects
/// face counts: estimated normals produce a slightly different field from the sphere fixtures'
/// analytic normals, so these counts differ by approximately one percent from the default
/// `reconstruct_surface` result. Both tables use this run and therefore remain consistent.
///
/// Repair is timed twice: once with the pipeline preset and once with the full pass set. This keeps
/// the cost of the two omitted passes visible.
fn time_stages(f: &Fixture) -> StageRun {
    let start = Instant::now();

    let index = f.cloud.compute_index().expect("index failed");
    let radius = f.voxel_size * 2.5;

    let t0 = Instant::now();
    let (estimates, _) = index
        .estimate_normals_oriented(
            radius,
            &NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1.0e6)),
        )
        .expect("normals failed");
    let normals = t0.elapsed().as_secs_f64();

    let band = 2.0 * f.voxel_size;
    let field = ImlsField3::try_new(
        f.cloud.points(),
        &estimates.normals,
        index.tree(),
        band,
        f.voxel_size,
    )
    .expect("field failed");

    let t1 = Instant::now();
    let grid = compute_sdf_grid(&field, f.cloud.points(), f.voxel_size, Point3::origin())
        .expect("grid failed");
    let field_time = t1.elapsed().as_secs_f64();

    let t2 = Instant::now();
    let (raw, stats) = extract_isosurface(&grid).expect("extraction failed");
    let extract = t2.elapsed().as_secs_f64();

    let default_opts = ReconstructOpts3::new(f.voxel_size)
        .repair
        .expect("the pipeline default runs repair");

    let t3 = Instant::now();
    let repaired = repair_buffers(raw.points(), raw.faces(), &default_opts).expect("repair failed");
    let repair = t3.elapsed().as_secs_f64();

    let data = MeshData3::new(repaired.points, repaired.faces).expect("mesh data failed");
    let total = start.elapsed().as_secs_f64();

    // Reconstruction now returns buffers. Measure the caller's optional hierarchy construction
    // separately and exclude it from the pipeline total.
    let t4 = Instant::now();
    let mesh = Mesh3::from_data(data.clone(), false).expect("mesh failed");
    let mesh_time = t4.elapsed().as_secs_f64();
    black_box(&mesh);

    // Measure the full repair outside the total for comparison only.
    let t5 = Instant::now();
    let full =
        repair_buffers(raw.points(), raw.faces(), &RepairOpts::default()).expect("full repair");
    let repair_full = t5.elapsed().as_secs_f64();
    black_box(&full);

    let faces = data.faces().len();

    StageRun {
        normals,
        field: field_time,
        extract,
        repair,
        repair_full,
        mesh: mesh_time,
        total,
        active_blocks: grid.block_count(),
        known_voxels: grid.known_count(),
        cells_visited: stats.cells_visited,
        cells_skipped: stats.cells_skipped_unknown,
        faces,
    }
}

/// Data prepared once so each benchmark measures only the selected stage.
struct Prepared {
    cloud: PointCloud3,
    voxel_size: f64,
    normals: Vec<UnitVec3>,
    raw_points: Vec<Point3>,
    raw_faces: Vec<[u32; 3]>,
    mesh: Mesh3,
}

fn prepare(f: &Fixture) -> Prepared {
    let index = f.cloud.compute_index().expect("index failed");
    let (estimates, _) = index
        .estimate_normals_oriented(
            f.voxel_size * 2.5,
            &NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1.0e6)),
        )
        .expect("normals failed");

    let band = 2.0 * f.voxel_size;
    let field = ImlsField3::try_new(
        f.cloud.points(),
        &estimates.normals,
        index.tree(),
        band,
        f.voxel_size,
    )
    .expect("field failed");

    let grid = compute_sdf_grid(&field, f.cloud.points(), f.voxel_size, Point3::origin())
        .expect("grid failed");
    let (raw, _) = extract_isosurface(&grid).expect("extraction failed");

    let default_opts = ReconstructOpts3::new(f.voxel_size)
        .repair
        .expect("the pipeline default runs repair");
    let repaired = repair_buffers(raw.points(), raw.faces(), &default_opts).expect("repair failed");
    let data = MeshData3::new(repaired.points, repaired.faces).expect("mesh data failed");
    let mesh = Mesh3::from_data(data, false).expect("mesh failed");

    Prepared {
        cloud: f.cloud.clone(),
        voxel_size: f.voxel_size,
        normals: estimates.normals,
        raw_points: raw.points().to_vec(),
        raw_faces: raw.faces().to_vec(),
        mesh,
    }
}

/// Measure each stage independently at one fixture size to identify the dominant stage.
fn stages(c: &mut Criterion) {
    let p = prepared(primary());

    let mut group = c.benchmark_group("reconstruct stages");
    group.sample_size(10);

    group.bench_function("index", |b| {
        b.iter(|| p.cloud.compute_index().expect("index failed"))
    });

    let index = p.cloud.compute_index().expect("index failed");
    let radius = p.voxel_size * 2.5;

    group.bench_function("normals/viewpoint", |b| {
        b.iter(|| {
            index
                .estimate_normals_oriented(
                    radius,
                    &NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1.0e6)),
                )
                .expect("normals failed")
        })
    });

    group.bench_function("normals/propagate", |b| {
        b.iter(|| {
            index
                .estimate_normals_oriented(radius, &NormalOrientation3::Propagate { k: 12 })
                .expect("normals failed")
        })
    });

    group.bench_function("point_spacing", |b| {
        b.iter(|| black_box(index.estimate_point_spacing()))
    });

    let band = 2.0 * p.voxel_size;
    let field = ImlsField3::try_new(
        p.cloud.points(),
        &p.normals,
        index.tree(),
        band,
        p.voxel_size,
    )
    .expect("field failed");

    group.bench_function("field", |b| {
        b.iter(|| {
            compute_sdf_grid(&field, p.cloud.points(), p.voxel_size, Point3::origin())
                .expect("grid failed")
        })
    });

    let grid = compute_sdf_grid(&field, p.cloud.points(), p.voxel_size, Point3::origin())
        .expect("grid failed");

    group.bench_function("extract", |b| {
        b.iter(|| extract_isosurface(&grid).expect("extraction failed"))
    });

    let default_repair = ReconstructOpts3::new(p.voxel_size)
        .repair
        .expect("the pipeline default runs repair");

    group.bench_function("repair/default", |b| {
        b.iter(|| {
            repair_buffers(&p.raw_points, &p.raw_faces, &default_repair).expect("repair failed")
        })
    });

    group.bench_function("repair/full", |b| {
        b.iter(|| {
            repair_buffers(&p.raw_points, &p.raw_faces, &RepairOpts::default())
                .expect("repair failed")
        })
    });

    group.bench_function("mesh_build", |b| {
        b.iter(|| {
            let data = MeshData3::new(p.raw_points.clone(), p.raw_faces.clone())
                .expect("mesh data failed");
            Mesh3::from_data(data, false).expect("mesh failed")
        })
    });

    group.bench_function("patch_filter", |b| {
        b.iter(|| {
            p.mesh
                .remove_small_patches(&PatchFilter::keep_largest())
                .expect("patch filter failed")
        })
    });

    group.bench_function("smooth/1", |b| {
        b.iter(|| {
            let mut editor = MeshEditor::new(&p.mesh).expect("editor failed");
            editor.smooth(1).expect("smoothing failed");
            editor.into_mesh(false, true).expect("mesh failed")
        })
    });

    group.bench_function("full/existing_normals", |b| {
        b.iter(|| {
            p.cloud
                .reconstruct_surface(&ReconstructOpts3::new(p.voxel_size))
                .expect("reconstruction failed")
        })
    });

    group.bench_function("full/estimated_normals", |b| {
        let opts = ReconstructOpts3::new(p.voxel_size).with_normals(NormalSource3::Estimate {
            radius,
            orientation: NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1.0e6)),
        });
        b.iter(|| {
            p.cloud
                .reconstruct_surface(&opts)
                .expect("reconstruction failed")
        })
    });

    group.finish();
}

/// Measure the dominant stages across point counts to determine how each one scales.
fn scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("reconstruct scaling");
    group.sample_size(10);

    for f in fixtures().iter().filter(|f| f.name.starts_with("sphere/")) {
        let p = prepared(f);
        let index = p.cloud.compute_index().expect("index failed");
        let band = 2.0 * p.voxel_size;
        let field = ImlsField3::try_new(
            p.cloud.points(),
            &p.normals,
            index.tree(),
            band,
            p.voxel_size,
        )
        .expect("field failed");
        let grid = compute_sdf_grid(&field, p.cloud.points(), p.voxel_size, Point3::origin())
            .expect("grid failed");

        let n = p.cloud.point_count();

        group.bench_with_input(BenchmarkId::new("field", n), &n, |b, _| {
            b.iter(|| {
                compute_sdf_grid(&field, p.cloud.points(), p.voxel_size, Point3::origin())
                    .expect("grid failed")
            })
        });

        group.bench_with_input(BenchmarkId::new("extract", n), &n, |b, _| {
            b.iter(|| extract_isosurface(&grid).expect("extraction failed"))
        });

        let default_repair = ReconstructOpts3::new(p.voxel_size)
            .repair
            .expect("the pipeline default runs repair");

        group.bench_with_input(BenchmarkId::new("repair/default", n), &n, |b, _| {
            b.iter(|| {
                repair_buffers(&p.raw_points, &p.raw_faces, &default_repair).expect("repair failed")
            })
        });

        group.bench_with_input(BenchmarkId::new("repair/full", n), &n, |b, _| {
            b.iter(|| {
                repair_buffers(&p.raw_points, &p.raw_faces, &RepairOpts::default())
                    .expect("repair failed")
            })
        });

        group.bench_with_input(BenchmarkId::new("mesh_build", n), &n, |b, _| {
            b.iter(|| {
                let data = MeshData3::new(p.raw_points.clone(), p.raw_faces.clone())
                    .expect("mesh data failed");
                Mesh3::from_data(data, false).expect("mesh failed")
            })
        });

        group.bench_with_input(BenchmarkId::new("full", n), &n, |b, _| {
            b.iter(|| {
                p.cloud
                    .reconstruct_surface(&ReconstructOpts3::new(p.voxel_size))
                    .expect("reconstruction failed")
            })
        });
    }

    group.finish();
}

/// Measure each repair pass independently because repair dominates the runtime and the passes have
/// different behavior.
///
/// The passes are interdependent. Removing degenerate and duplicate faces changes the input to
/// later passes, and `drop_isolated_vertices` has no work unless an earlier pass removed a face.
/// Each pass is therefore measured twice: once by itself on the raw buffers and once cumulatively
/// after all preceding passes in the default order. The cumulative measurement shows the pass's
/// cost in its normal position.
fn repair_passes(c: &mut Criterion) {
    let p = prepared(primary());

    let mut group = c.benchmark_group("repair passes");
    group.sample_size(10);

    /// A named repair pass and the builder call that enables it.
    type Pass = (&'static str, fn(RepairOpts) -> RepairOpts);

    let named: [Pass; 6] = [
        ("drop_degenerate", |o| o.with_drop_degenerate(true)),
        ("drop_duplicate_faces", |o| {
            o.with_drop_duplicate_faces(true)
        }),
        ("resolve_nonmanifold_edges", |o| {
            o.with_resolve_nonmanifold_edges(true)
        }),
        ("orient_consistently", |o| o.with_orient_consistently(true)),
        ("split_bowtie_vertices", |o| {
            o.with_split_bowtie_vertices(true)
        }),
        ("drop_isolated_vertices", |o| {
            o.with_drop_isolated_vertices(true)
        }),
    ];

    group.bench_function("none", |b| {
        let opts = RepairOpts::none();
        b.iter(|| repair_buffers(&p.raw_points, &p.raw_faces, &opts).expect("repair failed"))
    });

    for (name, apply) in named.iter() {
        let opts = apply(RepairOpts::none());
        group.bench_with_input(BenchmarkId::new("alone", name), &opts, |b, opts| {
            b.iter(|| repair_buffers(&p.raw_points, &p.raw_faces, opts).expect("repair failed"))
        });
    }

    let mut cumulative = RepairOpts::none();
    for (name, apply) in named.iter() {
        cumulative = apply(cumulative);
        group.bench_with_input(
            BenchmarkId::new("cumulative", name),
            &cumulative,
            |b, opts| {
                b.iter(|| repair_buffers(&p.raw_points, &p.raw_faces, opts).expect("repair failed"))
            },
        );
    }

    group.finish();
}

/// Measure parallel efficiency by running each stage in Rayon pools of different sizes.
///
/// A well-balanced parallel stage should take half as long when the pool size doubles. A curve that
/// flattens early and remains flat indicates that a serial section dominates. A curve that scales
/// well until the thread count approaches the item count indicates that the work is divided too
/// coarsely to occupy additional threads.
fn threads(c: &mut Criterion) {
    let p = prepared(primary());

    let index = p.cloud.compute_index().expect("index failed");
    let radius = p.voxel_size * 2.5;
    let band = 2.0 * p.voxel_size;
    let field = ImlsField3::try_new(
        p.cloud.points(),
        &p.normals,
        index.tree(),
        band,
        p.voxel_size,
    )
    .expect("field failed");
    let grid = compute_sdf_grid(&field, p.cloud.points(), p.voxel_size, Point3::origin())
        .expect("grid failed");

    let available = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);

    let mut group = c.benchmark_group("reconstruct threads");
    group.sample_size(10);

    for &t in THREADS.iter().filter(|&&t| t <= available) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build()
            .expect("thread pool failed");

        group.bench_with_input(BenchmarkId::new("normals", t), &t, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    index
                        .estimate_normals_oriented(
                            radius,
                            &NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1.0e6)),
                        )
                        .expect("normals failed")
                })
            })
        });

        group.bench_with_input(BenchmarkId::new("field", t), &t, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    compute_sdf_grid(&field, p.cloud.points(), p.voxel_size, Point3::origin())
                        .expect("grid failed")
                })
            })
        });

        group.bench_with_input(BenchmarkId::new("extract", t), &t, |b, _| {
            b.iter(|| pool.install(|| extract_isosurface(&grid).expect("extraction failed")))
        });

        group.bench_with_input(BenchmarkId::new("full", t), &t, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    p.cloud
                        .reconstruct_surface(&ReconstructOpts3::new(p.voxel_size))
                        .expect("reconstruction failed")
                })
            })
        });
    }

    group.finish();
}

/// Print the fixture summary once before Criterion takes control of the output.
fn summary(_c: &mut Criterion) {
    print_fixture_summary();
}

criterion_group!(benches, summary, stages, scaling, repair_passes, threads);
criterion_main!(benches);

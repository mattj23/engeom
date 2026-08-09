//! This module is for cutting a point sequence into contiguous runs, where each will end up being
//! written against its own tighter bounding box.
//!
//! A single bounding box charges every point for the full extent of the set. Splitting the set into
//! groups, each with its own box and its own per-axis widths, lowers the width every point in that
//! group pays, at the cost of one partition header per group.
//!
//! This module does that without moving the points. The runs it returns are contiguous spans of the
//! caller's own sequence, so point order survives untouched, which is what lets it be applied to
//! every container in the crate, since at a minimum the polyline containers don't have indices
//! because their order encodes connectivity information in its place.  Also, for now, the mesh
//! vertex order is the first-use order that [`crate::indices`] needs for its high-water coding, so
//! that currently has to be encoded as-is.
//!
//! # Contiguous runs
//!
//! Some real sequences of points will arrive spatially coherent.
//!
//! - Polylines are chains, so any contiguous run of it will be a short piece of curve inside a
//!   small set of bounds, unless it has very odd or degenerate point spacing.
//! - Raw point clouds or meshes have a tendency to arrive in acquisition order, which will often be
//!   something like a raster or a profile sweep.
//!
//! Where a sequence genuinely has no locality the partitioning will stop showing benefit and the
//! boxes will stop shrinking.
//!
//! # The cost model
//!
//! A partition costs `(5 + N) * 8` bits of fixed header before a single point is written, plus two
//! corner codes per axis at the anchor's widths. Unlike the fixed part, that second term is a
//! property of the data rather than of the format: at a typical metrology resolution it comes to
//! about a hundred bits in three dimensions, for a header of around 166 bits, or 20 bytes.
//!
//! Halving a box along one axis drops that axis by one bit, so a cut pays once a run holds a couple
//! of hundred points. The objective is a step function of the extents:
//! shrinking one by 30% is worth nothing and shrinking it by 51% is worth a full bit per point, so
//! candidates are priced with [`crate::points::partition_bytes`], the same arithmetic the writer
//! emits, rather than with a volume or variance proxy.
//!
//! # What is searched
//!
//! An exact dynamic program over every position costs `O(n²)` and is not affordable at the sizes
//! this crate is for, so the search runs over a candidate set instead of over every index.
//!
//! A forward greedy pass proposes a candidate wherever absorbing the next point widens the current
//! run's box, and closes its own run outright once that widening has cost the points already in it
//! more than a fresh header would. Those positions alone are not enough: they cluster where a box
//! is still small and thin out once it has grown, so on a steadily advancing sequence like a
//! polyline the positions a cut actually wants are not among them. Candidates are therefore also
//! placed at a regular [`stride`], sized in the cost model's own terms.
//!
//! A dynamic program then settles each candidate end against every start that can reach it, taking
//! ends in order and starts backwards from the nearest. Runs are priced from precomputed per-gap
//! boxes, and a growing box's widths only ever step upward, so neither is recomputed from scratch.
//! Two rules end a scan early, both of which are what keeps [`LOOKBACK`] affordable: a start that
//! has fallen a header behind the best way to reach that end can never be overtaken by an earlier
//! one, and a run already as wide as the whole set has stopped being a partition in any useful
//! sense.
//!
//! None of that is assumed to be free. `comes_within_a_percent_of_the_exhaustive_optimum` measures
//! the result against an exhaustive `O(n²)` dynamic program over every position. Separated
//! clusters and scattered clouds come out at the optimum; a smooth chain, whose cuts want a cadence
//! that no width jump marks, sits above it, by 0.13% on a loose helix and 0.52% on a tight one.
//!
//! One shape is worse than that and is held separately by
//! `a_raster_sweep_is_the_planners_weakest_shape`. A raster sweep plans 3.40% above the optimum,
//! because a run's cost steps at each row boundary rather than growing smoothly and the candidate
//! positions do not line up with those steps. It is the order a line scanner would deliver, so I
//! think it's worth marking as a place for potential improvement.  For my own use cases, I
//! typically store line scanner data in a form very similar to how it originally comes quantized
//! over the wire (which was the original inspiration for the tools that inspired this crate), so
//! this isn't a big issue for me personally, but if there is need it's something that can be given
//! attention.
//!
//! Whatever comes out is compared against writing the set as one partition and the smaller is
//! returned, so no input is ever planned into something larger than it is today.
//!
//! # Cost
//!
//! Planning a million points takes roughly three times what writing them does, and a sequence with
//! no locality costs less than that, since the dominance rule cuts every scan short almost
//! immediately.

use crate::bounds::Bounds;
use crate::effort::Effort;
use crate::error::Result;
use crate::points::{partition_bytes, widths_for};
use crate::quantize::{Quantizer, width_fits};
use crate::transform::{self, Rotation};

/// The most candidate breakpoints back the dynamic program will look when closing a run.
///
/// This is a ceiling rather than a cost. The two early exits in `search` end almost every scan long
/// before reaching it, so raising it buys reach on coarsely structured data without being paid for
/// on anything else. A run longer still is reachable through the single-partition comparison
/// [`plan`] finishes with, which is the case that wants one.
pub const LOOKBACK: usize = 512;

/// The coarsest spacing of candidate breakpoints, in points, for a set stored at `widths`.
///
/// Width jumps alone are not enough. They cluster where a run's box is still small and thin out
/// once it has grown, so on a steadily advancing sequence like a polyline the positions a cut wants
/// to sit at are not among them at all.
///
/// The spacing is expressed in the cost model's own terms: the number of points whose payload comes
/// to a partition header. Data whose points are expensive is worth placing cuts precisely in, and
/// data whose points are cheap is not, which is what makes this a scale rather than a constant.
/// Halving it again is worth about 0.3% on a smooth chain and costs four times the planning time.
pub fn stride<const N: usize>(widths: &[u8; N], effort: Effort) -> usize {
    let per_point: u64 = widths.iter().map(|&w| u64::from(w)).sum();
    if per_point == 0 {
        // Every axis is degenerate, so no cut can save anything and the greedy pass will propose
        // nothing either. The spacing is immaterial.
        return MAX_STRIDE;
    }

    // [`Effort::Thorough`] places candidates twice as finely, which measured about 0.3% smaller on
    // a smooth chain for about four times the planning time. [`Effort::Quick`] never gets here,
    // since it does not search at all.
    let divisor = match effort {
        Effort::Quick | Effort::Balanced => 1,
        Effort::Thorough => 2,
    };

    ((header_bits(widths) / per_point / divisor) as usize).clamp(1, MAX_STRIDE)
}

/// The widest spacing [`stride`] will return, so that cheap data still gets a look in.
const MAX_STRIDE: usize = 64;

/// The bits a partition spends before a single point of it is written: the count, the transform
/// flag, a width byte per axis, and its two corners as codes on the anchor lattice.
///
/// The corners are why this takes the anchor rather than only the dimension. They are the larger
/// part of the header and their cost is a property of the data, not of the format.
fn header_bits<const N: usize>(anchor_widths: &[u8; N]) -> u64 {
    let corners: u64 = anchor_widths.iter().map(|&w| 2 * u64::from(w)).sum();
    (4 + 1 + N as u64) * 8 + corners
}

/// The same, rounded up to whole bytes, which is what the dominance test in `search` compares.
fn header_bytes<const N: usize>(anchor_widths: &[u8; N]) -> u64 {
    header_bits(anchor_widths).div_ceil(8)
}

/// One contiguous run of the input, with the box and widths it is to be written against.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Run<const N: usize> {
    /// Index of the run's first point.
    pub start: usize,
    /// Index one past the run's last point.
    pub end: usize,
    /// The frame the run is written in, or `None` for world axes.
    ///
    /// When this is set, `bounds` and `widths` describe the points as seen from that frame rather
    /// than from the world, since those are the box and the widths that get written.
    pub transform: Option<Rotation<N>>,
    /// The box enclosing the run's points, in whatever frame `transform` names.
    pub bounds: Bounds<N>,
    /// The per-axis widths meeting the tolerance over `bounds`.
    pub widths: [u8; N],
}

impl<const N: usize> Run<N> {
    /// How many points the run holds.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Whether the run holds no points, which a planned run never does.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// What this run costs to write, header, corner codes and frame included.
    pub fn bytes(&self, anchor_widths: &[u8; N]) -> u64 {
        partition_bytes(
            self.len(),
            &self.widths,
            anchor_widths,
            self.transform.as_ref(),
        )
    }
}

/// A whole point block's worth of plan.
///
/// The anchor is the box enclosing every point, and the widths it would be stored at. Partition
/// corners are stored as codes on that lattice rather than as `f64`, which is most of what a
/// partition header used to cost, so the anchor has to travel with the runs that refer to it.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan<const N: usize> {
    /// The box enclosing every point.
    pub anchor: Bounds<N>,
    /// The widths the set as a whole needs, which is the resolution corners are stored at.
    pub anchor_widths: [u8; N],
    /// The runs, in sequence order, covering every point.
    pub runs: Vec<Run<N>>,
}

impl<const N: usize> Plan<N> {
    /// What the whole block costs to write, anchor included.
    pub fn bytes(&self) -> u64 {
        // The dimension byte, the partition count, and the anchor itself.
        let mut total = 1 + 4;
        if !self.runs.is_empty() {
            total += 16 * N as u64 + N as u64;
        }

        total
            + self
                .runs
                .iter()
                .map(|r| r.bytes(&self.anchor_widths))
                .sum::<u64>()
    }
}

/// What a corner code decodes to.
///
/// The two extreme codes are handled by name rather than by arithmetic. `min + (max - min)` is not
/// bit-identical to `max` in floating point, so a code of `max_int` run through the general formula
/// lands a few ulps below the anchor's own upper corner. That is harmless as a position and fatal as
/// a bound: a partition reaching the edge of the anchor would declare a box that does not quite
/// enclose its own outermost point, and the quantizer would clamp it. Codes 0 and `max_int` mean the
/// anchor's corners, so they return them.
fn corner_value(q: &Quantizer, code: u64, lo: f64, hi: f64) -> f64 {
    if code == 0 {
        lo
    } else if code == q.max_int() {
        hi
    } else {
        q.decode(code)
    }
}

/// The codes a box's corners take on the anchor lattice, rounded outward so the box they decode
/// back to still contains everything the original did.
///
/// Rounding outward rather than to nearest is what keeps the guarantee: a point outside the stored
/// box would be clamped by the quantizer to the boundary, which is an error no width can bound.
pub fn corner_codes<const N: usize>(
    bounds: &Bounds<N>,
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
) -> ([u64; N], [u64; N]) {
    let mut lo = [0u64; N];
    let mut hi = [0u64; N];

    for i in 0..N {
        if anchor_widths[i] == 0 {
            // A degenerate axis has one value and no lattice to place anything on.
            continue;
        }

        let q = Quantizer::new(
            anchor.mins[i],
            anchor.maxs[i] - anchor.mins[i],
            anchor_widths[i],
        );
        let top = q.max_int();

        let (alo, ahi) = (anchor.mins[i], anchor.maxs[i]);
        lo[i] = q.encode(bounds.mins[i]);
        hi[i] = q.encode(bounds.maxs[i]);

        // A millionth of a lattice step, which is nowhere near any error that matters and is still
        // thousands of times larger than the last-place noise a box arrives carrying. Without it,
        // a corner sitting on a lattice point steps outward or not depending on which side of that
        // point the arithmetic happened to land, and a box that moves a whole step moves every one
        // of its points with it.
        let slack = ((ahi - alo) / top as f64) * 1e-6;

        // `encode` rounds to nearest, so step outward until the decoded corners really do enclose.
        // Comparing against the decoded value rather than computing a floor keeps this faithful to
        // the arithmetic the decoder will do.
        while lo[i] > 0 && corner_value(&q, lo[i], alo, ahi) > bounds.mins[i] + slack {
            lo[i] -= 1;
        }
        while hi[i] < top && corner_value(&q, hi[i], alo, ahi) < bounds.maxs[i] - slack {
            hi[i] += 1;
        }
    }

    (lo, hi)
}

/// The box a set of corner codes decodes back to.
pub fn decode_corners<const N: usize>(
    lo: &[u64; N],
    hi: &[u64; N],
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
) -> Bounds<N> {
    let mut mins = anchor.mins;
    let mut maxs = anchor.maxs;

    for i in 0..N {
        if anchor_widths[i] == 0 {
            continue;
        }
        let q = Quantizer::new(
            anchor.mins[i],
            anchor.maxs[i] - anchor.mins[i],
            anchor_widths[i],
        );
        mins[i] = corner_value(&q, lo[i], anchor.mins[i], anchor.maxs[i]);
        maxs[i] = corner_value(&q, hi[i], anchor.mins[i], anchor.maxs[i]);
    }

    Bounds { mins, maxs }
}

/// Plan how to cut `points` into partitions, without reordering anything.
///
/// Returns one run per partition, in sequence order and covering every point. An empty input plans
/// to no runs at all, and [`Effort::Quick`] plans to a single run without searching.
///
/// # Errors
///
/// [`crate::Error::ToleranceNotRepresentable`] if the set as a whole spans a range too wide to meet
/// `tol`, and [`crate::Error::Malformed`] if any coordinate is not finite. The whole set is priced
/// before anything is cut, so a set that could not be written before cannot become writable by
/// being split up: the error behaviour is the same as it has always been.
pub fn plan<const N: usize>(points: &[[f64; N]], tol: f64, effort: Effort) -> Result<Plan<N>> {
    if points.is_empty() {
        return Ok(Plan {
            anchor: Bounds {
                mins: [0.0; N],
                maxs: [0.0; N],
            },
            anchor_widths: [0; N],
            runs: Vec::new(),
        });
    }

    let anchor = Bounds::from_points(points)?;
    let anchor_widths = widths_for(&anchor, tol)?;

    // The whole set as one partition, which is what the format wrote before any of this existed.
    // Its corners are the anchor's own, so they cost their codes and nothing is rounded.
    let single = finish_run(0, points.len(), &anchor, &anchor, &anchor_widths, tol)?;
    let mut best = Plan {
        anchor,
        anchor_widths,
        runs: vec![single],
    };

    // Every sub-box sits inside the whole one, so no run can ever need a wider axis than this. That
    // is what lets the growing widths below step upward without a representability check of their
    // own.
    let axis_tol = tol / (N as f64).sqrt();
    let whole = anchor_widths.iter().map(|&w| u32::from(w)).sum();
    let mut tried: Vec<usize> = Vec::new();

    for level in ladder(effort) {
        let spacing = stride(&anchor_widths, *level);

        // The clamps inside `stride` can land two levels on the same spacing, and running the same
        // search twice would only spend time to reach the same answer.
        if tried.contains(&spacing) {
            continue;
        }
        tried.push(spacing);

        let breaks = candidate_breaks(points, axis_tol, &anchor_widths, spacing);
        let runs = search(
            points,
            tol,
            axis_tol,
            &breaks,
            whole,
            &anchor,
            &anchor_widths,
        )?;
        let candidate = Plan {
            anchor,
            anchor_widths,
            runs,
        };

        if candidate.bytes() < best.bytes() {
            best = candidate;
        }
    }

    Ok(best)
}

/// Turn a span into the run that will be written for it.
///
/// The box is snapped onto the lattice its corners will be coded on and the widths are taken from
/// the snapped box. Snapping here rather than at write time is what keeps the plan's arithmetic and
/// the writer's output the same number. The rounding can widen a box by up to two lattice steps,
/// which very occasionally costs a bit, and a planner pricing the unrounded box would not know.
///
/// This is where a frame estimator belongs, and where one used to sit: fit a rotation over
/// `points[start..end]`, settle the run a second time against it, and keep whichever of the two
/// prices lower. [`crate::transform`] records why there is no estimator at the moment.
fn finish_run<const N: usize>(
    start: usize,
    end: usize,
    bounds: &Bounds<N>,
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
    tol: f64,
) -> Result<Run<N>> {
    // Every run is written on world axes. The format can carry a rotation and both ends of it are
    // implemented, but nothing chooses one; see the `transform` module for why not, and for what a
    // future estimator would slot in here.
    settle_run(start, end, None, bounds, anchor, anchor_widths, tol)
}

/// Snap a box onto the lattice its corners will be coded on and take the widths from what lands.
fn settle_run<const N: usize>(
    start: usize,
    end: usize,
    transform: Option<Rotation<N>>,
    bounds: &Bounds<N>,
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
    tol: f64,
) -> Result<Run<N>> {
    let lattice = match &transform {
        Some(rotation) => transform::rotated_anchor(anchor, rotation),
        None => *anchor,
    };

    let (lo, hi) = corner_codes(bounds, &lattice, anchor_widths);
    let snapped = decode_corners(&lo, &hi, &lattice, anchor_widths);
    let widths = widths_for(&snapped, tol)?;

    Ok(Run {
        start,
        end,
        transform,
        bounds: snapped,
        widths,
    })
}

/// The searches a level runs, which are its own and every one below it.
///
/// This is what makes "a higher level never produces a larger file" structural rather than a
/// property to be hoped for. A finer candidate spacing is not uniformly better: candidates are what
/// [`LOOKBACK`] is counted in, so placing them more finely also shortens the reach of the search,
/// and on noisy data whose runs want to be long that trade came out 149 bytes worse. Running the
/// coarser search as well and keeping the smaller result costs a level nothing but time, which is
/// the currency it was asked for.
fn ladder(effort: Effort) -> &'static [Effort] {
    match effort {
        Effort::Quick => &[],
        Effort::Balanced => &[Effort::Balanced],
        Effort::Thorough => &[Effort::Balanced, Effort::Thorough],
    }
}

/// A box being grown a point at a time, carrying the widths it currently needs.
#[derive(Debug, Clone, Copy)]
struct Growing<const N: usize> {
    mins: [f64; N],
    maxs: [f64; N],
    widths: [u8; N],
    /// The widths summed, which is what a point in this run costs.
    total: u32,
}

impl<const N: usize> Growing<N> {
    /// A run holding one point, which has zero extent on every axis and so costs no bits at all.
    fn seeded(p: &[f64; N]) -> Self {
        Self {
            mins: *p,
            maxs: *p,
            widths: [0; N],
            total: 0,
        }
    }

    /// A run covering one already-computed box.
    fn over(mins: [f64; N], maxs: [f64; N]) -> Self {
        Self {
            mins,
            maxs,
            widths: [0; N],
            total: 0,
        }
    }

    /// Widen to hold `p`, reporting whether any axis had to gain bits to do it.
    fn absorb(&mut self, p: &[f64; N], axis_tol: f64) -> bool {
        let mut widened = false;

        for (i, &c) in p.iter().enumerate() {
            if c < self.mins[i] {
                self.mins[i] = c;
            } else if c > self.maxs[i] {
                self.maxs[i] = c;
            } else {
                continue;
            }
            widened |= self.fit_axis(i, axis_tol);
        }

        widened
    }

    /// Widen to hold another box, reporting whether any axis had to gain bits.
    fn absorb_box(&mut self, mins: &[f64; N], maxs: &[f64; N], axis_tol: f64) -> bool {
        let mut widened = false;

        for i in 0..N {
            let grew = mins[i] < self.mins[i] || maxs[i] > self.maxs[i];
            self.mins[i] = self.mins[i].min(mins[i]);
            self.maxs[i] = self.maxs[i].max(maxs[i]);
            if grew {
                widened |= self.fit_axis(i, axis_tol);
            }
        }

        widened
    }

    /// Step one axis up until it holds its extent, reporting whether it moved.
    ///
    /// Widths only ever grow, since the box only ever grows, so this is amortized over the life of
    /// the run rather than being a search per point.
    fn fit_axis(&mut self, i: usize, axis_tol: f64) -> bool {
        let extent = self.maxs[i] - self.mins[i];
        if width_fits(extent, axis_tol, self.widths[i]) {
            return false;
        }

        while !width_fits(extent, axis_tol, self.widths[i]) {
            self.widths[i] += 1;
            self.total += 1;
        }

        true
    }
}

/// Propose the positions worth considering as breakpoints.
///
/// A cut only changes what a run costs per point where the box would otherwise have widened, so
/// those are the positions collected. The pass also closes a run outright once the width it just
/// took on has cost the points already in it more than a fresh header, which keeps later candidates
/// coming from a box that has not already saturated to the whole set.
fn candidate_breaks<const N: usize>(
    points: &[[f64; N]],
    axis_tol: f64,
    ceiling: &[u8; N],
    spacing: usize,
) -> Vec<usize> {
    let header = header_bits(ceiling);

    let mut breaks = vec![0usize];
    let mut run = Growing::seeded(&points[0]);
    let mut start = 0usize;

    for (j, p) in points.iter().enumerate().skip(1) {
        let before = run.total;
        let widened = run.absorb(p, axis_tol);

        debug_assert!(
            (0..N).all(|i| run.widths[i] <= ceiling[i]),
            "a run's box sits inside the whole set's, so it cannot need a wider axis"
        );

        if !widened {
            if j - breaks[breaks.len() - 1] >= spacing {
                breaks.push(j);
            }
            continue;
        }

        // The box just widened, so every point already in the run pays the difference from here on.
        // Cutting in front of this point is worth considering, and worth taking outright once that
        // retroactive cost has passed what a fresh header would have been.
        if breaks[breaks.len() - 1] != j {
            breaks.push(j);
        }
        if u64::from(run.total - before) * (j - start) as u64 > header {
            run = Growing::seeded(p);
            start = j;
        }
    }

    breaks.push(points.len());
    breaks
}

/// Choose among the candidate breakpoints by dynamic program, looking back [`LOOKBACK`] of them.
#[allow(clippy::too_many_arguments)]
fn search<const N: usize>(
    points: &[[f64; N]],
    tol: f64,
    axis_tol: f64,
    breaks: &[usize],
    whole: u32,
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
) -> Result<Vec<Run<N>>> {
    let gaps = breaks.len() - 1;

    // The box of each gap, so pricing a run is a merge of the gaps it spans and no point is visited
    // more than once no matter how many runs are considered over it.
    let gap_boxes: Vec<([f64; N], [f64; N])> = (0..gaps)
        .map(|t| {
            let span = &points[breaks[t]..breaks[t + 1]];
            let mut mins = span[0];
            let mut maxs = span[0];
            for p in &span[1..] {
                for (i, &c) in p.iter().enumerate() {
                    mins[i] = mins[i].min(c);
                    maxs[i] = maxs[i].max(c);
                }
            }
            (mins, maxs)
        })
        .collect();

    // cost[t] is the fewest bytes the points before breaks[t] can be written in.
    let mut cost = vec![u64::MAX; gaps + 1];
    let mut pred = vec![0usize; gaps + 1];
    cost[0] = 0;

    // Each end is settled against every start that could reach it, walking the starts backwards so
    // that the run's box grows a gap at a time. Taking the ends in order is what makes `cost[s]`
    // final by the time it is read, and taking the starts nearest-first is what gives the dominance
    // test below a good `cost[t]` to measure against straight away.
    for t in 1..=gaps {
        let (mins, maxs) = gap_boxes[t - 1];
        let mut run = Growing::over(mins, maxs);
        for i in 0..N {
            run.fit_axis(i, axis_tol);
        }

        for s in (t.saturating_sub(LOOKBACK)..t).rev() {
            if s + 1 < t {
                let (m, x) = &gap_boxes[s];
                run.absorb_box(m, x, axis_tol);
            }

            let candidate =
                cost[s] + partition_bytes(breaks[t] - breaks[s], &run.widths, anchor_widths, None);
            if candidate < cost[t] {
                cost[t] = candidate;
                pred[t] = s;
            }

            // Dominance. Once starting at `s` has fallen a header behind the best way to reach `t`,
            // no earlier start can beat it either, so the rest of the lookback is wasted work.
            //
            // For any earlier `s'`, reaching `s` costs at most `cost[s'] + header + the points in
            // s'..s`, which bounds `cost[s']` from below. Substituting that into the cost of the run
            // `s' -> t`, whose width is no narrower than this one's, cancels the points before `s`
            // and leaves `cost[s']`-based cost at least this one's less a header. The extra byte
            // covers the payload rounding, which is charged per partition rather than per point.
            if candidate > cost[t] + header_bytes(anchor_widths) + 1 {
                break;
            }

            // And a run as wide as the whole set has stopped being a partition in any useful sense.
            // Carrying it further only reproduces the single-partition option `plan` already
            // compares against.
            if run.total >= whole {
                break;
            }
        }
    }

    // Walk the chosen predecessors back from the end, then turn each span into a run whose bounds
    // and widths are derived the same way the single-partition path derives them.
    let mut chosen = vec![gaps];
    while *chosen.last().unwrap() > 0 {
        chosen.push(pred[*chosen.last().unwrap()]);
    }
    chosen.reverse();

    let mut runs = Vec::with_capacity(chosen.len() - 1);
    for pair in chosen.windows(2) {
        let (start, end) = (breaks[pair[0]], breaks[pair[1]]);
        let bounds = Bounds::from_points(&points[start..end])?;
        runs.push(finish_run(start, end, &bounds, anchor, anchor_widths, tol)?);
    }

    Ok(runs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corpus;
    use crate::testgen::Rng;

    /// The exhaustive optimum: a dynamic program over every position, with no candidate set and no
    /// lookback limit. Quadratic, so only ever run on small inputs, and the only thing that can say
    /// whether the two restrictions in [`plan`] cost anything.
    fn exhaustive<const N: usize>(points: &[[f64; N]], tol: f64) -> u64 {
        let n = points.len();
        let axis_tol = tol / (N as f64).sqrt();
        let anchor = Bounds::from_points(points).unwrap();
        let anchor_widths = widths_for(&anchor, tol).unwrap();

        let mut cost = vec![u64::MAX; n + 1];
        cost[0] = 0;

        for end in 1..=n {
            let mut run = Growing::seeded(&points[end - 1]);
            for i in 0..N {
                run.fit_axis(i, axis_tol);
            }

            for start in (0..end).rev() {
                if start < end - 1 {
                    run.absorb(&points[start], axis_tol);
                }
                if cost[start] == u64::MAX {
                    continue;
                }
                let candidate =
                    cost[start] + partition_bytes(end - start, &run.widths, &anchor_widths, None);
                cost[end] = cost[end].min(candidate);
            }
        }

        // The block's own overhead, so this is comparable with what `Plan::bytes` reports: the
        // dimension byte, the partition count and the anchor.
        cost[n] + 1 + 4 + 16 * N as u64 + N as u64
    }

    /// Runs must tile the input in order, leaving no gap and no overlap, and their bounds and
    /// widths must be the ones the writer would derive for those same points.
    ///
    /// A run's box and widths are stated in whatever frame the run carries, so everything here is
    /// checked in that frame rather than in world axes.
    fn assert_well_formed<const N: usize>(points: &[[f64; N]], plan: &Plan<N>, tol: f64) {
        let centre = transform::centre_of(&plan.anchor);
        let mut at = 0;

        for run in &plan.runs {
            assert_eq!(run.start, at, "runs must be contiguous");
            assert!(run.end > run.start, "a planned run must hold points");

            let span: Vec<[f64; N]> = points[run.start..run.end]
                .iter()
                .map(|p| match &run.transform {
                    Some(rotation) => rotation.to_local(p, &centre),
                    None => *p,
                })
                .collect();
            let lattice = match &run.transform {
                Some(rotation) => transform::rotated_anchor(&plan.anchor, rotation),
                None => plan.anchor,
            };

            // The stored box is the tight one snapped outward onto the corner lattice, and it still
            // has to enclose every point it claims, or the quantizer would clamp them.
            let tight = Bounds::from_points(&span).unwrap();
            for i in 0..N {
                assert!(
                    run.bounds.mins[i] <= tight.mins[i] && run.bounds.maxs[i] >= tight.maxs[i],
                    "run {}..{}: stored box does not enclose its points on axis {i}\n stored  [{:.17e}, {:.17e}]\n tight   [{:.17e}, {:.17e}]\n lattice [{:.17e}, {:.17e}] w={}\n turned  {}",
                    run.start,
                    run.end,
                    run.bounds.mins[i],
                    run.bounds.maxs[i],
                    tight.mins[i],
                    tight.maxs[i],
                    lattice.mins[i],
                    lattice.maxs[i],
                    plan.anchor_widths[i],
                    run.transform.is_some()
                );
            }

            let (lo, hi) = corner_codes(&run.bounds, &lattice, &plan.anchor_widths);
            assert_eq!(
                decode_corners(&lo, &hi, &lattice, &plan.anchor_widths),
                run.bounds,
                "run {}..{}: stored box is not on the corner lattice",
                run.start,
                run.end
            );
            assert_eq!(
                run.widths,
                widths_for(&run.bounds, tol).unwrap(),
                "run {}..{}: widths",
                run.start,
                run.end
            );

            at = run.end;
        }
        assert_eq!(at, points.len(), "runs must cover the whole input");
    }

    /// What the same points cost as a single partition, which is the baseline every plan has to
    /// beat or match.
    fn single_partition<const N: usize>(points: &[[f64; N]], tol: f64) -> u64 {
        let anchor = Bounds::from_points(points).unwrap();
        let widths = widths_for(&anchor, tol).unwrap();
        1 + 4 + 16 * N as u64 + N as u64 + partition_bytes(points.len(), &widths, &widths, None)
    }

    /// Tight clusters visited one after another, which is what an ordered sequence over separated
    /// regions looks like: a scan pass over several features, or a file of stacked cross sections.
    fn clustered(clusters: usize, per: usize, spread: f64, radius: f64) -> Vec<[f64; 3]> {
        let mut rng = Rng::new(7_101);
        let mut points = Vec::with_capacity(clusters * per);

        for _ in 0..clusters {
            let centre = rng.point::<3>(-spread, spread);
            for _ in 0..per {
                points.push(std::array::from_fn(|i| centre[i] + rng.gaussian(radius)));
            }
        }

        points
    }

    /// A raster sweep over a gently curved surface, which is the order a line scanner delivers.
    ///
    /// Its shape is awkward because a run shorter than a row has no extent at all on the row axis,
    /// so the cost of a run jumps up sharply at each row boundary, and where the cheapest cut
    /// actually is has little to do with where the box last widened.
    fn raster(rows: usize, per_row: usize, span: f64) -> Vec<[f64; 3]> {
        let step = span / per_row as f64;
        (0..rows * per_row)
            .map(|i| {
                let u = (i % per_row) as f64 * step;
                let v = (i / per_row) as f64 * step;
                [u, v, 2.0 * (u * 0.05).sin() + (v * 0.03).cos()]
            })
            .collect()
    }

    /// A chain marching along a helix, which is what a polyline or a toolpath is: every contiguous
    /// piece of it sits in a small box even though the whole spans a large one.
    fn chain(count: usize, radius: f64, pitch: f64) -> Vec<[f64; 3]> {
        (0..count)
            .map(|i| {
                let t = i as f64 * 0.01;
                [radius * t.cos(), radius * t.sin(), pitch * t]
            })
            .collect()
    }

    #[test]
    fn an_empty_input_plans_to_nothing() {
        let points: Vec<[f64; 3]> = Vec::new();
        assert!(
            plan(&points, 1e-3, Effort::Balanced)
                .unwrap()
                .runs
                .is_empty()
        );
    }

    #[test]
    fn a_single_point_plans_to_one_run() {
        let plan = plan(&[[1.0, 2.0, 3.0]], 1e-3, Effort::Balanced).unwrap();
        assert_eq!(plan.runs.len(), 1);
        assert_eq!(plan.runs[0].widths, [0, 0, 0]);
    }

    /// The restrictions in [`plan`], the candidate set and the lookback, are where this could
    /// quietly stop finding good cuts. Nothing else in the module would notice.
    ///
    /// Separated clusters come out at the optimum, since every cut wants to sit at a boundary the
    /// candidate set contains anyway. A smoothly advancing chain is the harder case: its cuts want
    /// a steady cadence that no width jump marks, so they land on [`stride`] positions instead and
    /// the result sits a little above the optimum, 0.72% at the worst of these.
    ///
    /// The margin is set to catch a real regression rather than that residue. Dropping the stride
    /// entirely, which is what this module did before it was measured, costs 1.2% here and fails.
    #[test]
    fn comes_within_a_percent_of_the_exhaustive_optimum() {
        let cases: Vec<(&str, Vec<[f64; 3]>, f64)> = vec![
            ("two clusters", clustered(2, 150, 400.0, 0.5), 1e-3),
            ("eight clusters", clustered(8, 200, 500.0, 0.5), 1e-3),
            ("many small clusters", clustered(40, 60, 500.0, 0.2), 1e-3),
            ("helix", chain(2000, 50.0, 4.0), 1e-3),
            ("tight helix", chain(1500, 2.0, 0.1), 1e-4),
            (
                "random cloud",
                Rng::new(7_102).points(1200, -50.0, 50.0),
                1e-3,
            ),
            (
                "short cloud",
                Rng::new(7_103).points(200, -50.0, 50.0),
                1e-3,
            ),
        ];

        for (name, points, tol) in cases {
            let plan = plan(&points, tol, Effort::Balanced).unwrap();
            assert_well_formed(&points, &plan, tol);

            let planned = plan.bytes();
            let optimal = exhaustive(&points, tol);

            assert!(
                planned <= optimal + optimal / 100,
                "{name}: planned {planned} bytes into {} runs against an optimum of {optimal}, \
                 a gap of {:.2}%",
                plan.runs.len(),
                100.0 * (planned as f64 / optimal as f64 - 1.0)
            );
        }
    }

    /// The planner's weakest known shape.
    ///
    /// A raster sweep does not come within a percent of the optimum the way the other shapes do. It
    /// plans 7730 bytes against an optimum of 7476, a gap of 3.40%, where separated clusters land on
    /// the optimum and a helix sits 0.13% above it.
    ///
    /// What makes it hard is the row. Inside a single row the sweep axis has no extent, so a run's
    /// cost steps sharply as it crosses a row boundary rather than growing smoothly, and the greedy
    /// pass proposes candidates from where the box widened, which is not where those steps are. This
    /// is the point order you would find from a panning line scanner, points recovered from a depth
    /// image, or probably the point order from the raw snapshots of a structured light scanner
    /// although I'm not 100% sure about that claim.
    ///
    /// This test is to catch regression so the limits in it are how it currently performs today
    /// while I'm writing this.  I'm not sanctioning these values as the goal, or even saying
    /// they're desireable. It's just here to let me know if I do anything that makes it worse.
    #[test]
    fn a_raster_sweep_is_the_planners_weakest_shape() {
        let points = raster(20, 100, 100.0);
        let tol = 1e-3;

        let plan = plan(&points, tol, Effort::Balanced).unwrap();
        assert_well_formed(&points, &plan, tol);

        let planned = plan.bytes();
        let optimal = exhaustive(&points, tol);
        let gap = 100.0 * (planned as f64 / optimal as f64 - 1.0);

        assert!(
            planned < single_partition(&points, tol),
            "a raster sweep still has to be worth cutting at all"
        );
        assert!(
            gap <= 4.0,
            "planned {planned} bytes into {} runs against an optimum of {optimal}, \
             a gap of {gap:.2}% where 3.40% is the figure on record",
            plan.runs.len()
        );
    }

    /// The guarantee that makes this safe to turn on by default.
    #[test]
    fn never_plans_larger_than_a_single_partition() {
        let mut cases: Vec<(String, Vec<[f64; 3]>, f64)> = corpus::all()
            .into_iter()
            .filter(|c| !c.points.is_empty())
            .map(|c| (c.name.to_string(), c.points, c.tol))
            .collect();

        cases.push(("helix".into(), chain(20_000, 50.0, 4.0), 1e-3));
        cases.push(("clusters".into(), clustered(12, 900, 500.0, 0.5), 1e-3));

        for (name, points, tol) in cases {
            let plan = plan(&points, tol, Effort::Balanced).unwrap();
            assert_well_formed(&points, &plan, tol);

            let single = single_partition(&points, tol);
            assert!(
                plan.bytes() <= single,
                "{name}: planned {} bytes into {} runs against {single} as one",
                plan.bytes(),
                plan.runs.len()
            );
        }
    }

    /// A sequence with no locality has nothing to gain, and must not be cut into pieces that each
    /// pay a header for a box no smaller than the whole.
    #[test]
    fn a_sequence_without_locality_stays_whole() {
        let points: Vec<[f64; 3]> = Rng::new(7_104).points(50_000, -50.0, 50.0);
        let plan = plan(&points, 1e-3, Effort::Balanced).unwrap();

        assert_eq!(
            plan.runs.len(),
            1,
            "expected one run, got {}",
            plan.runs.len()
        );
    }

    /// The shuffled corpus case is the same geometry as the smooth surface with its point order
    /// destroyed, so it is the direct measurement of the thing this module trades on.
    #[test]
    fn locality_is_what_pays_not_geometry() {
        let smooth = corpus::smooth_surface();
        let shuffled = corpus::shuffled();

        let ordered = plan(&smooth.points, smooth.tol, Effort::Balanced).unwrap();
        let scrambled = plan(&shuffled.points, shuffled.tol, Effort::Balanced).unwrap();

        assert!(
            ordered.bytes() < scrambled.bytes(),
            "the same geometry in a coherent order must plan smaller: {} against {}",
            ordered.bytes(),
            scrambled.bytes()
        );
    }

    /// Separated clusters are the case the whole feature exists for, so the saving has to be large
    /// rather than merely present.
    #[test]
    fn separated_clusters_collapse() {
        let points = clustered(8, 2000, 500.0, 0.5);
        let tol = 1e-3;

        let plan = plan(&points, tol, Effort::Balanced).unwrap();
        let single = single_partition(&points, tol);

        // A cluster of this radius needs about 12 bits an axis where the whole spread needs 20, so
        // the floor is around 0.6 and the rest is the eight headers.
        let ratio = plan.bytes() as f64 / single as f64;
        assert!(
            ratio < 0.65,
            "expected clusters to collapse, got {ratio} at {} runs",
            plan.runs.len()
        );
    }

    /// A cut has to earn back its header out of the difference in width across it, so short
    /// sequences stay whole. This falls out of the accounting rather than being a special case.
    ///
    /// There is no universal count below which nothing is ever cut: the floor depends on how much
    /// narrower the pieces are than the whole. This chain clears it at 128 points, where cutting in
    /// two is both what the planner does and what the exhaustive optimum agrees on.
    #[test]
    fn a_short_sequence_is_never_cut() {
        for count in [2usize, 10, 32] {
            let points = chain(count, 50.0, 4.0);
            let plan = plan(&points, 1e-3, Effort::Balanced).unwrap();

            assert_eq!(plan.runs.len(), 1, "{count} points should stay in one run");
        }
    }

    /// Two dimensions carry a smaller header, so they should cut sooner, and the generic code has
    /// to work there at all.
    #[test]
    fn plans_in_two_dimensions() {
        let points: Vec<[f64; 2]> = (0..4000)
            .map(|i| {
                let t = i as f64 * 0.01;
                [50.0 * t.cos(), 50.0 * t.sin()]
            })
            .collect();

        let tol = 1e-3;
        let plan = plan(&points, tol, Effort::Balanced).unwrap();
        assert_well_formed(&points, &plan, tol);
        assert!(plan.bytes() < single_partition(&points, tol));
    }

    /// A set no single box can hold to the tolerance is refused, and is not quietly made writable
    /// by being cut into pieces that each happen to fit.
    #[test]
    fn an_unrepresentable_set_is_still_refused() {
        let points = clustered(4, 100, 1e6, 1e-9);

        assert!(matches!(
            plan(&points, 1e-9, Effort::Balanced),
            Err(crate::Error::ToleranceNotRepresentable { .. })
        ));
    }

    #[test]
    fn non_finite_coordinates_are_rejected() {
        let points = [[0.0, 0.0, 0.0], [f64::NAN, 1.0, 1.0]];

        assert!(matches!(
            plan(&points, 1e-3, Effort::Balanced),
            Err(crate::Error::Malformed(_))
        ));
    }
}

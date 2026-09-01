//! How hard should an encoder look to make a smaller file?
//!
//! Most of what this crate does after basic bit packing is searching: searching for where a point
//! set should be cut into partitions, searching for which index coding is smaller, and whatever
//! else gets added later. Those searches pay encoding time to buy smaller files, and how _much_
//! should be paid is up to the caller.
//!
//! [`Effort`] is the mechanism for the caller to say what it's willing to pay. It's one setting
//! used throughout the library instead of having lots of different options each place that encoding
//! happens. I wanted to have something that would (a) let callers express a general preference for
//! compute cost vs file size rather than forcing them to know the exact best way to combine
//! features to achieve their compression goals, and (b) not require changes to the API if I found
//! and added more compression features.
//!
//! # Limits on what "effort" affects
//!
//! `Effort` has three main rules meant to keep the feature usable:
//!
//! - **It never changes what comes back.** Every level meets the tolerance, and every level answers
//!   with the same points and the same connectivity. Anything that changes what a caller gets back,
//!   such as [`crate::mesh::VertexOrder`] renumbering vertices, stays an explicit setting of its own
//!   and is never folded in here.
//! - **It never reaches the reader.** Effort is not recorded in the file and a decoder cannot tell
//!   which level wrote one. It follows that a technique which moves work onto the decoder does not
//!   belong on this dial, whatever it saves.
//! - **It never costs size.** A higher level never produces a larger file than a lower one, because
//!   each level's search subsumes the one below it and every search ends by comparing what it found
//!   against not having bothered.
//!
//! # Conventions
//!
//! Each container module carries its own `WriteOptions`, holding the file-level metadata, this
//! setting, and whatever semantic settings that container actually has. They are separate types
//! rather than one shared one so that a container never accepts a setting it silently ignores: a
//! polyline's point order is the curve itself and can never be renumbered, so [`crate::polyline`]
//! has no ordering field to set.
//!
//! ```
//! use tol_compress::{Cloud3, Effort, cloud};
//!
//! let scan = Cloud3::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
//! let options = cloud::WriteOptions::new().with_effort(Effort::Thorough);
//!
//! let mut buf = Vec::new();
//! cloud::write_to_with(&mut buf, std::slice::from_ref(&scan), 1e-4, &options)?;
//! # Ok::<(), tol_compress::Error>(())
//! ```

/// How hard an encoder should search for a smaller representation.
///
/// Ordered, so `Quick < Balanced < Thorough`, which is the order of both the time they cost and the
/// sizes they achieve.
///
/// This is marked non-exhaustive because levels may be added. A caller matching on one needs a
/// wildcard arm; constructing and comparing them is unaffected.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum Effort {
    /// Write the representation that falls out of the data directly, searching for nothing.
    ///
    /// The bit widths are still the narrowest that meet the tolerance, since that is the format
    /// rather than an optimization, and the index coding is still whichever of the two measures
    /// smaller, since trying both costs almost nothing. What this level skips is everything that
    /// costs a pass over the data to decide.
    ///
    /// Wanted where the encoder is in a loop with something interactive, or where the file is
    /// about to be thrown away.
    Quick,

    /// Search where the searching is known to pay for itself. The default.
    ///
    /// Measurement data is written once and read many times, often into an archive that outlives
    /// the tooling that made it, so the balance sits nearer size than speed.
    #[default]
    Balanced,

    /// Search harder, for the last few percent.
    ///
    /// The returns here are small and the cost is not. Wanted where the file is a deliverable, or
    /// is going somewhere that charges by the byte, and where encoding time is nobody's problem.
    Thorough,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corpus;
    use crate::{Cloud3, Mesh3, cloud, mesh};

    #[test]
    fn levels_are_ordered_by_the_work_they_do() {
        assert!(Effort::Quick < Effort::Balanced);
        assert!(Effort::Balanced < Effort::Thorough);
        assert_eq!(Effort::default(), Effort::Balanced);
    }

    /// The guarantee that makes a single dial safe to hand a caller: asking for more effort can
    /// cost time, but it can never cost bytes.
    ///
    /// Today the levels all write the same file and this holds trivially. It is here so that it
    /// keeps holding once they stop being the same, which is the moment it could quietly stop being
    /// true and nothing else would notice.
    #[test]
    fn more_effort_never_makes_a_larger_file() {
        for case in corpus::all() {
            let sizes: Vec<usize> = [Effort::Quick, Effort::Balanced, Effort::Thorough]
                .into_iter()
                .map(|effort| {
                    let mut buf = Vec::new();
                    if case.faces.is_empty() {
                        let item = Cloud3::new(case.points.clone());
                        let options = cloud::WriteOptions::new().with_effort(effort);
                        cloud::write_to_with(
                            &mut buf,
                            std::slice::from_ref(&item),
                            case.tol,
                            &options,
                        )
                        .unwrap();
                    } else {
                        let item = Mesh3::new(case.points.clone(), case.faces.clone());
                        let options = mesh::WriteOptions::new().with_effort(effort);
                        mesh::write_to_with(
                            &mut buf,
                            std::slice::from_ref(&item),
                            case.tol,
                            &options,
                        )
                        .unwrap();
                    }
                    buf.len()
                })
                .collect();

            assert!(
                sizes[0] >= sizes[1] && sizes[1] >= sizes[2],
                "{}: sizes went up with effort, {sizes:?}",
                case.name
            );
        }
    }

    /// Effort is an encoder-side choice and leaves no trace in the file, so a reader cannot tell
    /// which level wrote one and never has to care.
    #[test]
    fn the_level_is_not_recorded_in_the_file() {
        let case = corpus::smooth_surface();
        let item = Mesh3::new(case.points.clone(), case.faces.clone()).named("part");

        for effort in [Effort::Quick, Effort::Balanced, Effort::Thorough] {
            let mut buf = Vec::new();
            let options = mesh::WriteOptions::new().with_effort(effort);
            mesh::write_to_with(&mut buf, std::slice::from_ref(&item), case.tol, &options).unwrap();

            let back = mesh::read_one_from(&mut buf.as_slice()).unwrap();
            assert_eq!(back.name.as_deref(), Some("part"), "{effort:?}");
            assert_eq!(back.points.len(), item.points.len(), "{effort:?}");
            assert_eq!(back.faces.len(), item.faces.len(), "{effort:?}");
        }
    }
}

# Releasing

Right now this repo has two crates that it publishes on different versions/schedules, and unfortunately one of them (`engeom`) is dependent on the other (`tol-compress`) in a way that getting the order wrong causes a Github CI failure.

It's not a big deal, except that `tol-compress` publishes infrequently enough that by the time I need to do it I've forgotten the correct order and I make the same mistake every time.  This document exists to be an explicit instruction of how to release and why it has to be done this way.

## `tol-compress` and `engeom`

The `tol-compress` crate ignores the shared workspace version used by `engeom` and its Python bindings `py-engeom` and instead has its own version.  Since the CI/CD pipelines here use `git` tags to mark releases, it has its own tag format to disambiguate it from the `engeom` releases.

That's because `tol-compress` is its own crate with no dependencies inside the workspace.  It was originally part of `engeom` but was useful enough to break it out into a separate crate, so `engeom` requires it but not vice versa.

| | `tol-compress` | `engeom` |
| --- | --- | --- |
| Version source | `tol-compress/Cargo.toml`, stated literally | `[workspace.package] version` in the root `Cargo.toml` |
| Tag | `tol-compress-vX.Y.Z` | `vX.Y.Z` |
| Workflow | `.github/workflows/release-tol-compress.yml` | `.github/workflows/release.yml` |
| Publishes to | crates.io | crates.io and PyPI |

Also note,  `py-engeom` is `publish = false`; it is never published as a crate but instead becomes the PyPI `engeom` wheels built by `wheels.yml`.

Engeom's `Cargo.toml` declares the dependency as such:

```toml
tol-compress = { version = "0.3", path = "../tol-compress" }
```

Local builds and tests use the `path` reference, so changes in the local repo are picked up without having to do a whole release.  

But `cargo package` has to use the `version` reference, since packaging can't preserve a path that points outide of the package directory.  That means that _packaged_ `engeom` needs to depend on `tol-compress` from [crates.io](https://crates.io/crates/tol-compress), and `cargo package` verifies itself by building...so...you can probably see the issue.

If, when `cargo package` runs, the version referenced in `engeom/Cargo.toml`'s `tol-compress.version` field is not _already_ on crates.io, the CI package job fails with the following error:

```
error: failed to prepare local package for uploading

Caused by:
  failed to select a version for the requirement `tol-compress = "^0.3"`
  candidate versions found which didn't match: 0.2.0, 0.1.0
  location searched: crates.io index
```

That will error out on every pull request and every push to `develop` until the actual `tol-compress` crate is published.  Therefore:

> [!IMPORTANT]
>  A `tol-compress` version must be published before any branch that sets `engeom`'s requirement to it will pass CI.

## Releasing `tol-compress`

1. Bump `tol-compress`'s version and `engeom`'s requirement in a single commit. This keeps everything internally self-consistent, but obviously will result in a failing CI if it's done on a pull request or a branch marked to run the CI pipeline.
    - Update `version` in `tol-compress/Cargo.toml`
    - For good hygeine, write the new section in `tol-compress/CHANGELOG.md`
    - Set `version = "X.Y"` in `engeom/Cargo.toml`
2. Check that nothing is broken on `tol-compress` by running its release action as a rehersal.
    - Github Actions > "Release tol-compress" > Run workflow
    - Pick the branch, leave `dry_run` checked. 
    - This should run `check-version`, `verify`, and `msrv`, but won't publish anything.
    - If all of the steps pass, move on.  Otherwise fix it and retry.
3. Push a release git tag at that commit. The tag can be on the feature branch; it doesn't have to be merged first.

   ```bash
   git tag tol-compress-v0.3.0 <commit>
   git push origin tol-compress-v0.3.0
   ```
4. Wait for the run and the deployment to crates.io to finish.

At this point if you were in the middle of a failing pull request you should be able to re-run any failed jobs.  If you did this from a branch that doesn't run CI on every commit and you weren't working on a pull request, you should be good to go already.

Just keep in mind that publishing is not reversible.

Also, note that a re-pushed tag after a partial failure is safe. The `publish-crates` step queries crates.io first and skips the publish if the version is already there.

## Releasing `engeom`

Make sure that the `tol-compress` version that `engeom` requires is already published on crates.io.

1. Bump `[workspace.package] version` in the root `Cargo.toml` and push it.
2. Do a pull request into `main`: all `engeom` releases should be on the `main` branch.
3. Make sure that all CI jobs are green.
4. Merge the pull request and close it.
5. Apply the version tag to the resulting commit on `main` and push it:

   ```bash
   git tag v0.6.0 <commit>
   git push origin v0.6.0
   ```

The CI will run again, with the additional release job.  You can dry-run this the same way that you do for `tol-compress` but I haven't ever needed to so far.

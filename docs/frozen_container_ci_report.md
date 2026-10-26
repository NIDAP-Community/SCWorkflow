# SCWorkflow Frozen-Container CI Strategy And Setup

Date: 2026-10-26
Repository: NIDAP-Community/SCWorkflow
Task branch: `DEV-frozen-container-ci`
Branch base: `DEV` commit `87634480f6194e70ec38d254ce16c4a158aee694`
Status: implementation prepared for PR review; container execution and GitHub results pending.

## Purpose

Test the current SCWorkflow source against an existing, preserved dependency environment instead of reconstructing legacy R dependencies from today's CRAN or Bioconductor repositories. Older dependency versions may no longer resolve through the current repository endpoints even when they remain usable in the archived container.

This work changes SCWorkflow package validation workflows and makes `ggalluvial` optional at package load. It does not change capsule workflows, dependency declarations, Seurat, package versions, or the container build. Capsule CI troubleshooting remains deferred.

## Strategy

1. Check out the source revision selected by the GitHub event.
2. Run package validation inside a prebuilt container pinned to an immutable Linux image digest.
3. Use the R runtime and dependency packages already present in that image. Do not run dependency installation, dependency resolution, or package-update actions in these validation jobs.
4. Build and install the checked-out SCWorkflow package as required by package checking and coverage. Do not validate only the older SCWorkflow installation that may already be in the image.
5. Fail clearly if required validation tools or namespace-import dependencies are unavailable. Do not automatically download replacements or silently skip validation.
6. Preserve environment evidence and genuine check failures for review.

This freezes the R environment, not the entire GitHub runner. GitHub Actions, image pulling, artifact uploads, Codecov, and tests that access external data can still require network access. No claim of completely offline execution is made.

## Selected Image

The candidate was taken from the image already configured in the package's coverage job; a new image was not invented or built.

| Item | Value |
| --- | --- |
| Existing tag | `nciccbr/scworkflow:v1.0.2_79e5d37` |
| Registry index digest observed | `sha256:8d441f51fac5d3be23ebc3ef8c19cfa014202d91364b5109b6fe60094d910495` |
| Pinned Linux image | `nciccbr/scworkflow@sha256:1bf587adfb1f4703e1e653ee1628dfeb93937fc16a2cc706aa47064a6a710b6a` |
| Selected platform | `linux/amd64` |
| Runtime metadata | `R_VERSION=4.1.3`, `R_HOME=/usr/local/lib/R` |
| Declared image base | `docker.io/rocker/r-ver:4.1.3` |

The remote manifest and image configuration were readable on 2026-10-26 without pulling the image layers. Runtime metadata is not a substitute for executing R in the image. Installed Seurat, Matrix, and CI-tool versions have not yet been verified inside this image.

Collaborator confirmation is required that this is the intended archived dependency image. The digest identifies the exact artifact inspected; it does not by itself establish that its contents are complete, correct, secure, or the intended release environment.

### Ggalluvial Optional Dependency

The user reports that the existing container does not include `ggalluvial`. This has not yet been independently verified by executing the pinned image.

The original `DEV` namespace imported `geom_flow` unconditionally despite `ggalluvial` being listed in `Suggests`. This branch removes the corresponding roxygen import and regenerates [NAMESPACE](../NAMESPACE). [compareCellPopulations()](../R/Compare_Cell_Populations.R) now checks the optional dependency when called and uses `ggalluvial::geom_flow()` explicitly. Missing `ggalluvial` therefore no longer creates that unconditional package-load dependency; calling the function without it gives a clear error.

The preflight checks the checked-out namespace rather than hard-coding `ggalluvial` as required. Other unavailable namespace imports or tools still fail clearly. No dependency is installed automatically and no comparison tests are skipped by this change. Existing comparison tests that call the function can still fail when the optional package is absent, and package checks can separately report unavailable suggested packages. Passing full CI is not guaranteed by removing the namespace import.

Confirm the actual image inventory and GitHub results before deciding on any compatible image extension. No `ggalluvial` version has been selected and no image extension has been made. `ggalluvial` remains in `Suggests`; moving it to `Imports` is not part of this optional-dependency approach.

## Changes Made

| File | Change |
| --- | --- |
| [R-CMD-check workflow](../.github/workflows/R-CMD-check.yaml) | Run package checking, lint, and coverage in the same pinned image. Remove `install-r-pak` and `setup-r-dependencies`. Retain existing package-check, goodpractice, lint, coverage, and upload actions. |
| [Separate coverage workflow](../.github/workflows/test-coverage.yaml) | Use the same image and preflight instead of constructing dependencies with pak. Preserve its existing triggers and coverage upload behavior. |
| [Shared preflight](../.github/scripts/ci-preflight.R) | Verify the package root, record the source SHA and installed environment, and check mode-specific tools and namespace-import loadability without installing anything. |
| [Build exclusions](../.Rbuildignore) and [Git exclusions](../.gitignore) | Exclude generated `ci-artifacts/` from the package build and version control. |

The package-check job is now named `Frozen container check` rather than a matrix label that assumes a separately provisioned R installation. The aggregate job ID `check` is retained. Its success condition now requires all three validation jobs to succeed, including rejecting skipped or cancelled jobs. The aggregate job itself runs on the host because it executes no R code.

The two existing coverage runs have not been consolidated as part of this change. Deployment, documentation-site, and image-building workflows are outside scope. Existing warning-only lint behavior is unchanged; the setup does not introduce a stricter lint policy or suppress package-check failures.

## Preflight And Evidence

| Mode | Required preinstalled tools |
| --- | --- |
| `check` | `rcmdcheck`, `testthat` |
| `lint` | `goodpractice`, `lintr` |
| `coverage` | `covr`, `testthat` |

Each job records these files before raising a missing-tool or missing-import error:

- `ci-artifacts/session-info.txt`: source SHA from `GITHUB_SHA`, run ID, validation mode, and R session information.
- `ci-artifacts/installed-packages.csv`: package names, versions, library paths, and build versions.
- `ci-artifacts/namespace-dependencies.csv`: packages imported by the checked-out `NAMESPACE` and whether each namespace is loadable in the environment.

The workflows attempt to upload these files even if validation fails, as `frozen-check-environment`, `frozen-lint-environment`, or `frozen-coverage-environment`. Existing package-check and failed-test artifacts remain available through their original steps. If checkout or container startup fails, the preflight cannot produce an inventory.

The preflight checks validation tooling and namespace imports, not every package-qualified or conditional runtime dependency of every SCWorkflow feature. Actual package checks and tests must still expose other missing feature dependencies or compatibility failures. Unavailable namespaces can also reflect incompatible installed dependencies, not only completely absent packages. Image gaps require a separately reviewed image decision, not a package-install fallback in CI.

## Verification Performed

| Verification | Result |
| --- | --- |
| Branch created from fetched current `DEV` | Confirmed; base commit recorded above. |
| Existing image manifest and runtime configuration | Read successfully from the registry; Linux image digest and R metadata recorded above. |
| YAML parsing of both changed workflows | Passed using the installed R `yaml` package. |
| All four R-running jobs use the same immutable digest and shared preflight | Passed static checks. |
| Dependency installer actions removed from the two changed workflows | Passed static checks. |
| Embedded workflow R code and preflight syntax | Parsed successfully. |
| Isolated preflight behavior | Eleven scenarios cover three success modes, three missing-tool failures, optional `ggalluvial` absent in all three modes, invalid-mode rejection, and wrong-package rejection. Namespace and tool availability are simulated, not container results. |
| Optional-dependency regression | The new missing-optional-package test and existing non-Seurat input test passed. Namespace regeneration removed only the `ggalluvial` import, and `Suggests` was retained. |
| Inventory files in isolated scenarios | All three inventories are checked on successful, missing-tool, and optional-`ggalluvial`-absent paths. |
| Container execution | Not performed: the local Docker daemon was not running. |
| Actions-specific lint | Not performed: `actionlint` was not installed; no installation was attempted. |
| Full package check, tests, lint, or coverage inside the frozen image | Not performed. |
| GitHub workflow execution for this branch | Not performed; no successful remote result is claimed. |

The isolated harness tests validate preflight control flow, not the contents of the frozen image. Local R results must not be presented as R 4.1.3 container results. No installed R packages or package-version constraints were changed during setup.

## Collaborator Verification Procedure

### Review The Revision

After the branch is committed and pushed, fetch `DEV-frozen-container-ci` and record the exact revision under review:

```bash
git fetch origin DEV-frozen-container-ci
git switch DEV-frozen-container-ci
git rev-parse HEAD
git diff origin/DEV -- .github/workflows .github/scripts .Rbuildignore .gitignore
```

Confirm that function changes are limited to the optional dependency guard, qualified `geom_flow` call, and corresponding roxygen/namespace import removal. Dependency declarations and `.github/package-versions.txt` are unchanged. The historical pins remain Seurat 4.1.1 and Matrix 1.5.1; compare the actual image inventory with the intended frozen environment rather than assuming that the tag guarantees these versions.

### Inspect And Run The Image

From the package root, with Docker running:

```bash
image='nciccbr/scworkflow@sha256:1bf587adfb1f4703e1e653ee1628dfeb93937fc16a2cc706aa47064a6a710b6a'
docker buildx imagetools inspect "$image"
docker pull --platform linux/amd64 "$image"
docker run --rm --platform linux/amd64 \
  -e CI=true -e NOT_CRAN=true \
  -e GITHUB_SHA="$(git rev-parse HEAD)" \
  -v "$PWD:/workspace" -w /workspace "$image" \
  Rscript .github/scripts/ci-preflight.R check
```

Repeat the final command with mode `lint` and then `coverage`. Review all three inventories, especially `ggalluvial`, Seurat, Matrix, testthat, rcmdcheck, covr, lintr, and goodpractice. With the optional dependency change, `ggalluvial` is no longer a namespace requirement, but other missing imports or tools can still block the preflight. Apple Silicon hosts may require amd64 emulation; GitHub's selected Ubuntu runners are amd64.

These preflight commands do not run the package test suite. They can create root-owned local artifacts because the container writes to a mounted directory; use a disposable checkout for verification if needed.

### Execute GitHub Validation

After publishing the branch, open a PR targeting `DEV` to exercise both workflows, or use the existing manual dispatch for `R-CMD-check`:

```bash
gh workflow run R-CMD-check.yaml \
  --repo NIDAP-Community/SCWorkflow \
  --ref DEV-frozen-container-ci
```

A push to the task branch alone does not trigger the current push filters. The separate coverage workflow has no manual-dispatch trigger and is exercised by a PR to a supported target branch.

Review the run's source SHA, selected digest, preflight inventories, actual check/test output, and aggregate result. Verify that the checked-out source is built and tested, not just the package already installed in the image. Review any required-status configuration affected by the package-check display-name change. Do not mark this strategy verified based on a manifest lookup or preflight alone.

## Acceptance Checklist

- [ ] Collaborator confirms the selected image is the intended archived environment.
- [ ] Runtime inventory confirms the intended R, Seurat, Matrix, and validation-tool versions.
- [ ] The reported `ggalluvial` gap is confirmed in the image; package loading and the optional-function error are reviewed separately from full test-suite results.
- [ ] All preflight modes succeed in the pinned image without dependency installation.
- [ ] GitHub runs against the reviewed revision and produces environment artifacts.
- [ ] Package checks and coverage execute against the checked-out source, with failures investigated rather than bypassed.
- [ ] Aggregate status and any required-status configuration are verified.
- [ ] Any remaining runtime, metadata, data-access, or Codecov failures are recorded separately.

## Limits And Next Steps

Passing in a frozen container demonstrates compatibility with that environment. It does not prove a fresh installation from `DESCRIPTION`, correctness of all dependency declarations, support for newer R/Seurat, or ongoing security support for the legacy stack. Dependency-metadata gaps identified in the release audit remain open and must not be hidden by the preinstalled libraries.

Confirm the image contents and execute the branch on GitHub before merging this setup. Then return to the small dependency-metadata patch. If the image is incomplete, agree on a replacement built from preserved package sources and pin its digest in a separate, reviewable change; do not silently fetch current replacements during CI.
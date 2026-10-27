# SCWorkflow Frozen-Container CI Strategy And Setup

Date: 2026-10-26
Repository: NIDAP-Community/SCWorkflow
Task branch: `DEV-frozen-container-ci`
Branch base: `DEV` commit `87634480f6194e70ec38d254ce16c4a158aee694`
PR: [#107](https://github.com/NIDAP-Community/SCWorkflow/pull/107), draft, targeting `DEV`; no merge authorized.
Verified published frozen-runtime/workflow revision: `fd4dd6f58b0caa295cc70a23bfd7b729ef79ec14`; the later Gitflow alignment below has local validation, with GitHub verification pending.
Status: At published revision `93903bc`, frozen preparation, lint, and both coverage jobs pass. The package check reports examples OK and 0 errors, 6 warnings, 6 notes; the action fails because its existing threshold is warnings. The legacy Gitflow failure is addressed by the validation alignment below, verified locally in frozen R 4.1.3; its GitHub result remains pending. All-green CI is not claimed.

## Local Gitflow Validation Alignment

The user authorized aligning Gitflow with the preserved environment. A direct image substitution in the shared workflow is insufficient: its push check calls `conda env list`, and its scripts assume an old renv library, call invalid `libPaths()`, and parse log text to detect failures. The shared parser has no setup hook for restoring the isolated library.

The [Gitflow caller](../.github/workflows/gitflow-R-action.yml) now runs validation locally using the same preparation workflow, rs413s digest, artifact restoration, and namespace preflight as the other frozen validation jobs. The [local checker](../.github/scripts/ci-gitflow-check.R) retains the existing event-specific intent:

- Pushes lint changed R files and run changed unit tests plus tests matching changed R sources, retaining fixture/snapshot/helper and legacy integration-test exclusions for unit-test selection.
- Release PRs run all tests matching package function filenames, including suffixed test variants rather than only changed files.
- Main PRs regenerate documentation and run a package check, with an explicit warning failure threshold consistent with the separate package-check workflow.

Each selected unit-test file runs in a fresh R process. Testthat failures and process errors fail the job directly; no test-log parsing, hard-coded library path, Conda command, or package installation is needed. `.libPaths()` reports the actual restored library. Changed-file lint and unit-test warnings remain visible, matching the existing advisory lint and test-warning behavior. New-branch pushes use the merge base with the default branch when the prior SHA is all zeroes; unresolvable comparison commits fail instead of silently selecting nothing.

The original event/branch filters are unchanged. The shared parser runs only for closed, merged PRs, retaining its original deployment image, artifact RIDs, package/R settings, and production/development secret mappings. Deployment runtime migration is not part of this validation change. DEV and main report no required status-check contexts, so the new validation-job name has no observed conflict; repository protection settings were not changed.

Actionlint passes, and executable checks cover event routing, changed-test/source matching, case-insensitive names, exclusions, deduplication, documentation-only changes, release variants, push comparison commits, and new-branch fallback. An offline/read-only frozen R 4.1.3 run loads the package and passes all three existing trajectory test cases using the restored library. Its six test warnings and 21 namespace-import warnings remain visible. A deliberately failing scratch test is correctly rejected, proving child-process failure propagation. The main-PR tools devtools, roxygen2, and rcmdcheck already load without installation.

Evidence is under `/tmp/scworkflow-gitflow-check.JS3luM`: `check.log`, `check-isolated.log`, and `check-final.log`. Original package sources, tests, base image, and tool library were read-only inputs; execution used a separate scratch package copy with networking disabled. No SIF was accessed or modified, and no dependency or Seurat version changed. At local verification, full main/release validation and GitHub execution of this follow-up had not been run. The user subsequently authorized commit and push; the six package-check warnings are not resolved or suppressed.

## Latest Published Results: Two Separate Blockers

The user's continuing error was checked against published commit `93903bc195b01ef7b4bd6f36224cea776ef3771f`, not an older run. [R-CMD-check run 40439487233, package-check job 121976870667](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40439487233/job/121976870667) confirms:

```text
* checking examples ... OK
Status: 6 WARNINGs, 6 NOTEs
Error: R CMD check found WARNINGs
```

There are no R CMD check errors in this result. The action logs show `error-on: "warning"`, so the six warnings still fail the job. They are:

- Nonportable file names and paths in vignettes, templates, and fixtures.
- Installation warnings from 21 namespace-import replacements.
- Undeclared code dependencies (Matrix, SingleCellExperiment, SummarizedExperiment, nidapFunctions, scater, scuttle, shiny) and package-code library/require usage.
- Missing documentation for exported `build_modscore_plots`.
- Documentation/usage mismatches in compareCellPopulations and filterQC.
- Vignette files without the expected inst/doc output or VignetteBuilder configuration.

Preparation, lint, and main coverage have completed successfully; [separate coverage run 40439487286](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40439487286) also succeeds. The aggregate job was still finishing at evidence capture. These results resolve the trajectory example failure but do not validate its raw input or establish that every unit test ran.

[Gitflow run 40439484735](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40439484735) fails separately. The [Gitflow caller](../.github/workflows/gitflow-R-action.yml) delegates to `fnlcr-bids-sdsi/gitflow-R/.github/workflows/parser.yml@master` with `ghcr.io/nidap-community/single-cell-rna-r4:latest`, not the pinned rs413s image and isolated tool library. Its test output includes `The package "EnhancedVolcano" is required` while loading imports, and its diagnostic command calls nonexistent `libPaths()` instead of `.libPaths()`. These are environment/workflow problems, not evidence that a trajectory unit-test assertion failed.

The earlier four frozen R jobs refer only to R-CMD-check and test-coverage, not this additional legacy Gitflow path. Gitflow needed a separate alignment of its runtime and required tools; merely fixing package documentation warnings would not repair that environment. At that diagnostic stage, no warning gate, workflow, runtime function, unit test, dependency version, deployment setting, or original SIF was changed. That result update was then local pending publication authorization; the subsequent alignment is described above.

## Trajectory Example Policy Follow-Up

The user explicitly requested that trajectory follow the same example strategy as the other 21 functions. Those source examples use `\dontrun{}` with illustrative inputs; trajectory alone used `\donttest{}`, which runs under the existing `--as-cran` check. Changed only that wrapper in the trajectory roxygen source and regenerated its help page. The example body and all parsed runtime expressions are unchanged.

Local and offline/read-only frozen R 4.1.3 `tools::Rd2ex()` checks confirm that normal extraction skips the example even when `donttest` examples are enabled, while explicit `dontrun` extraction retains it. This is a user-approved example execution policy change, not a repair of missing reductions or evidence that the full trajectory workflow succeeds. Unit tests, global check flags, warning-level gates, and dependency versions are unchanged. The required parameter inventory and Graphify refresh succeeded.

At local verification, this follow-up was uncommitted; full package check and fresh GitHub CI had not been run. Existing package warnings and the example's unsuitable raw input remain separate findings. The earlier sections below preserve their original revision and verification context. The original SIF was neither accessed nor modified.

## Patchwork Declaration Follow-Up

The user's report of a continuing error was confirmed in [run 40430987972](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40430987972), job `121949158014`, at published revision `41b3897`. R CMD check still stops at the missing DESCRIPTION declaration; passing environment and coverage checks do not remove that package error.

The follow-up adds `patchwork` to DESCRIPTION Imports, matching its existing unconditional NAMESPACE import. It does not install or upgrade patchwork, change version constraints, modify function code/tests, or touch the original SIF. The unchanged R 4.1.3 image already provides patchwork 1.2.0.

A focused offline/read-only container check parses NAMESPACE as R expressions and DESCRIPTION as DCF, verifies that every non-base imported package is declared, and confirms patchwork loads. That check passes. Standard R base packages are exempted using installed-package Priority metadata. At that local verification point, GitHub results were pending; the subsequent published result is recorded next, and earlier findings below retain their original revision context.

The correction was published as `2a8f64c85b8a4b20ef0cc64061482f4a1ae2714b`. [R-CMD-check run 40432161914](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40432161914) confirms the patchwork error is gone and the check advances beyond dependency validation. Preparation, lint, and both coverage jobs pass; [separate coverage run 40432162075](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40432162075) succeeds. The package check remains failed with 1 error, 6 warnings, and 6 notes.

## Trajectory Declarations And Remaining Example Failure

The next error in that GitHub run is `Missing required packages: TrajectoryUtils, TSCAN, mclust` while executing the unchanged trajectory example. Direct read-only inspection shows that all three already load in the frozen image: TrajectoryUtils 1.2.0, TSCAN 1.32.0, and mclust 6.0.0. Although the user authorized compatible missing additions to the isolated library, no package additions, rebuilds, or upgrades are needed for these three.

The follow-up declares these existing optional dependencies in DESCRIPTION Suggests. A focused offline R 4.1.3 check then advances the unchanged example past namespace loading. It stops later with:

```text
'PCA' not in 'reducedDimNames(<SingleCellExperiment>)'
```

The example constructs a small Seurat object without a PCA reduction, while the function expects that reduction. This remaining example/input issue must be handled separately from environment preparation; no function/example/test code or validation gates were changed to bypass it. Package warnings also remain, including namespace-import replacement, other undeclared code dependencies, documentation mismatches, vignette layout, and file portability. The warning-level CI gate is retained, so fixing the example alone would not prove all-green CI.

Local evidence is stored under `/tmp/scworkflow-trajectory-check.8qW76S`. The focused check uses `--no-tests --no-vignettes --no-codoc` but still runs the original examples; it is not presented as full CI. An initial build exhausted a 1 GiB temporary filesystem while copying fixtures; a disk-backed TMPDIR fixed only the validation harness. The image and prepared tool library remained read-only with networking disabled. Fresh GitHub verification of the optional-dependency declarations remains pending.

## Purpose

Test the current SCWorkflow source against an existing, preserved dependency environment instead of reconstructing legacy R dependencies from today's CRAN or Bioconductor repositories. Older dependency versions may no longer resolve through the current repository endpoints even when they remain usable in the archived container.

This work changes SCWorkflow package validation workflows, makes `ggalluvial` optional at package load, corrects existing patchwork and optional trajectory dependency declarations, and aligns the trajectory documentation example's execution policy with the other examples. It does not change capsule workflows, Seurat, installed dependency versions, or the container build. Capsule CI troubleshooting remains deferred.

## Strategy

1. Check out the source revision selected by the GitHub event.
2. Run package validation inside a prebuilt container pinned to an immutable Linux image digest.
3. Reuse the R runtime and dependency packages already present in that image. A separate preparation job builds only missing, checksum-pinned tools into an isolated library; no base package is replaced. Validation jobs restore that artifact without dependency installation, dependency resolution, or package updates.
4. Build and install the checked-out SCWorkflow package as required by package checking and coverage. Do not validate only the older SCWorkflow installation that may already be in the image.
5. Fail clearly if required validation tools or namespace-import dependencies are unavailable. Do not automatically download replacements or silently skip validation.
6. Preserve environment evidence and genuine check failures for review.

This freezes the R environment, not the entire GitHub runner. GitHub Actions, image pulling, artifact uploads, Codecov, and tests that access external data can still require network access. No claim of completely offline execution is made.

## Current Verified Environment

All four R validation jobs now select `rfinney/rs413s@sha256:7cb32eaac452ad32abb34e9a9e3edb42590e4dfcf696caeb72bb375a5b3c83cf`, the public Docker source named in the supplied SIF. Direct Docker execution confirms R 4.1.3, Seurat/SeuratObject 4.1.1, Matrix 1.5-1, and ggplot2 3.3.6, matching the SIF's inspected versions. All 43 namespace imports load. This verifies the relevant runtime, not byte-for-byte equivalence of the Docker and SIF filesystems. The original SIF is neither used nor modified by the workflows.

The shared [environment-preparation workflow](../.github/workflows/frozen-environment.yaml) downloads only the exact tool archives in [the source lock](../.github/scripts/ci-tool-sources.json). Host HTTPS verification stays enabled, and every source is checked by SHA-256 both on the host and before installation. This avoids the legacy container's observed certificate-trust failure without disabling TLS validation or changing its installed packages.

The [preparation script](../.github/scripts/ci-prepare-tools.R) runs with the base filesystem read-only and networking disabled. It installs the following additions into `/ci-tools/library`, with dependency resolution disabled:

| Addition | Version |
| --- | --- |
| covr | 3.6.1 |
| lintr | 3.0.2 |
| goodpractice | 1.0.4 |
| ggalluvial | 0.12.3 |
| rex | 1.2.1 |
| cyclocomp | 1.1.0 |
| xmlparsedata | 1.0.5 |
| clisymbols | 1.2.0 |
| whoami | 1.3.0 |

The script refuses to shadow any base package and verifies the complete base package/version/library-path inventory is unchanged. Existing rcmdcheck 1.4.0, testthat 3.1.6, and all other installed prerequisites are reused. No dependency already installed in the preserved image is rebuilt or upgraded.

Preparation verifies all three preflight modes, then uploads `frozen-ci-tools` containing `library.tar.gz`, `base-packages.csv`, and `tools-manifest.csv`. Both caller workflows wait for preparation; each validation job restores the library and sets `R_LIBS_USER` through `GITHUB_ENV`. The aggregate requires preparation and all validation jobs to succeed. Failed setup is not skipped or treated as a passing code check.

Local verification restored the archive into a fresh executable temporary filesystem in a separate offline container. `check`, `lint`, and `coverage` each exited 0, and the base inventory still matched exactly. A deliberately corrupted source was rejected before installation. The current source also completed `R CMD INSTALL` and loaded from its new temporary installation, not a preinstalled SCWorkflow copy. Its 21 import-replacement warnings are recorded separately; no warning was suppressed to obtain a pass. Actionlint 1.7.7 passes all three workflows.

These are environment and installation checks, not full package tests, coverage, lint results, or a successful GitHub run. Source downloads and image pulls still require network access. The tool library must be on an executable filesystem: an initial local `noexec` tmpfs prevented covr's shared library from loading, and the verification harness was corrected to match an executable library mount.

## GitHub Environment Readiness Verified

At revision `fd4dd6f58b0caa295cc70a23bfd7b729ef79ec14`, [R-CMD-check run 40430026101](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40430026101) and [coverage run 40430026047](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40430026047) confirm the following environment gates:

| Job | Verified environment result |
| --- | --- |
| Main preparation (`121945155067`) | Source preparation, all three offline preflight modes, and artifact upload pass. |
| Frozen container check (`121945906519`) | Artifact restoration and check-mode preflight pass. |
| Lint (`121945907007`) | Artifact restoration and lint-mode preflight pass. |
| Main coverage (`121945906572`) | Artifact restoration and coverage-mode preflight pass. |
| Separate coverage preparation (`121945154888`) | Preparation and all-mode verification pass. |
| Separate coverage (`121945986184`) | Artifact restoration and coverage-mode preflight pass. |

This establishes that the real runner can prepare and restore the preserved environment before testing checked-out package code. Their package-validation outcomes remain separate from the successful environment gates; no tests, package warnings, or aggregate failure conditions were bypassed.

## Package Validation Findings

Both coverage jobs complete successfully at the verified revision, including the complete separate coverage workflow. The package-check job reaches R CMD check and stops at dependency validation with **1 error, 0 warnings, 0 notes**:

```text
Namespace dependency missing from DESCRIPTION Imports/Depends entries: patchwork
```

This is a package metadata error, not a missing installation: `patchwork` loads in the preserved environment, but its namespace import is not declared in DESCRIPTION. The check stops before its later checks/tests, so the zero warning/note counts are not an exhaustive clean result. No package metadata, code, or tests were altered to bypass it. Reviewing and correcting that declaration is the next separate package-side task.

At this evidence capture, the goodpractice step has completed and the lint step is still running. No lint-success claim is made. The aggregate must remain non-successful when package check fails, even though environment preparation and coverage pass. The original SIF and all existing dependency versions remain unchanged.

## First GitHub Preparation Results And Path Correction

At revision `7ab36a1a2e8bbce654d9375157b641db7b688fb9`, the preparation jobs in [R-CMD-check run 40429160106](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40429160106) and [coverage run 40429159946](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40429159946) both passed. The real GitHub runners downloaded and verified the sources, prepared the isolated library, passed all three preflight modes, and uploaded the artifact. The downstream package jobs did not reach code validation.

Their download action correctly mapped the host temp directory to `/__w/_temp/frozen-ci-tools` inside the job container, but the shell restore commands still referenced `/home/runner/work/_temp/frozen-ci-tools`. Shell paths now use `$RUNNER_TEMP`; action inputs retain the runner context because the actions translate it. The same correction applies to coverage diagnostic output paths. No package or dependency change is involved.

A focused local check extracted the actual four restore blocks from the workflow YAML, required them to be identical, executed them under `RUNNER_TEMP=/__w/_temp`, restored the artifact on an executable temporary filesystem, and passed check/lint/coverage preflight again with networking disabled. The subsequent GitHub results above confirm restoration and preflight pass in every downstream validation job.

## Original CI Image (Rejected)

The original candidate was taken from the image already configured in the package's coverage job. The evidence below explains why it was replaced; no new Docker image was built.

| Item | Value |
| --- | --- |
| Existing tag | `nciccbr/scworkflow:v1.0.2_79e5d37` |
| Registry index digest observed | `sha256:8d441f51fac5d3be23ebc3ef8c19cfa014202d91364b5109b6fe60094d910495` |
| Pinned Linux image | `nciccbr/scworkflow@sha256:1bf587adfb1f4703e1e653ee1628dfeb93937fc16a2cc706aa47064a6a710b6a` |
| Selected platform | `linux/amd64` |
| Runtime metadata | `R_VERSION=4.1.3`, `R_HOME=/usr/local/lib/R` |
| Declared image base | `docker.io/rocker/r-ver:4.1.3` |

The tag was rechecked on 2026-10-26 and still resolves to the recorded index and Linux image digests. The first GitHub runs execute R 4.1.3 on Ubuntu 20.04.6 LTS. Their inventories record Matrix 1.5-1 in the site library and Matrix 1.4-0 in the base library, but do not list Seurat or the required CI tools. Subsequent local execution reproduces that inventory under both default and GitHub-style `HOME`; a full filesystem scan finds no installed directories for the blocking packages. See the execution evidence below.

Collaborator confirmation is required that this is the intended archived dependency image. The digest identifies the exact artifact inspected; it does not by itself establish that its contents are complete, correct, secure, or the intended release environment.

### Ggalluvial Optional Dependency

The GitHub inventory and both local R environments do not list `ggalluvial` in their active libraries. The full filesystem investigation finds no alternate complete dependency library that could resolve the blocking imports and tools.

The original `DEV` namespace imported `geom_flow` unconditionally despite `ggalluvial` being listed in `Suggests`. This branch removes the corresponding roxygen import and regenerates [NAMESPACE](../NAMESPACE). [compareCellPopulations()](../R/Compare_Cell_Populations.R) now checks the optional dependency when called and uses `ggalluvial::geom_flow()` explicitly. Missing `ggalluvial` therefore no longer creates that unconditional package-load dependency; calling the function without it gives a clear error.

The preflight checks the checked-out namespace rather than hard-coding `ggalluvial` as required. Other unavailable namespace imports or tools still fail clearly. No dependency is installed automatically and no comparison tests are skipped by this change. Existing comparison tests that call the function can still fail when the optional package is absent, and package checks can separately report unavailable suggested packages. Passing full CI is not guaranteed by removing the namespace import.

The current isolated library supplies checksum-pinned `ggalluvial` 0.12.3 for comparison tests. It is not installed into the original image or SIF. `ggalluvial` remains in `Suggests`; moving it to `Imports` is not part of this optional-dependency approach.

## Changes Made

| File | Change |
| --- | --- |
| [R-CMD-check workflow](../.github/workflows/R-CMD-check.yaml) | Run package checking, lint, and coverage in the same pinned image. Remove `install-r-pak` and `setup-r-dependencies`. Retain existing package-check, goodpractice, lint, coverage, and upload actions. |
| [Separate coverage workflow](../.github/workflows/test-coverage.yaml) | Use the same image and preflight instead of constructing dependencies with pak. Preserve its existing triggers and coverage upload behavior. |
| [Shared environment preparation](../.github/workflows/frozen-environment.yaml) | Download locked sources on the host, prepare only missing tools offline in a read-only base, verify all modes, and upload the isolated library. |
| [Tool preparation script](../.github/scripts/ci-prepare-tools.R) and [source lock](../.github/scripts/ci-tool-sources.json) | Pin versions and source checksums, disable dependency resolution, and reject changes to existing base packages. |
| [Shared preflight](../.github/scripts/ci-preflight.R) | Verify the package root, record the source SHA and installed environment, and check mode-specific tools and namespace-import loadability without installing anything. |
| [Build exclusions](../.Rbuildignore) and [Git exclusions](../.gitignore) | Exclude generated `ci-artifacts/` from the package build and version control. |

The package-check job is now named `Frozen container check` rather than a matrix label that assumes a separately provisioned R installation. The aggregate job ID `check` is retained. Its success condition requires preparation and all three validation jobs to succeed, including rejecting skipped or cancelled jobs. The aggregate job itself runs on the host because it executes no R code.

The two existing coverage runs have not been consolidated as part of this change. Deployment, documentation-site, and image-building workflows are outside scope. Existing warning-only lint behavior is unchanged; the setup does not introduce a stricter lint policy or suppress package-check failures.

## Preflight And Evidence

| Mode | Required tools from the base or prepared library |
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
| Container execution | GitHub and local Docker both execute the pinned image's R 4.1.3. Local read-only filesystem scans and default/GitHub-style `HOME` comparisons are complete. |
| Actions-specific lint | Actionlint 1.7.7 passes all three updated workflows using a tooling container; no host installation. |
| Current environment preflights | Check, lint, and coverage all exit 0 after archive restoration in an offline read-only container. Full base inventory is unchanged. |
| Current source installation/load | Passed in the new environment; 21 namespace-import replacement warnings remain. |
| GitHub environment restoration/preflights | Pass in all four validation jobs at `fd4dd6f`, after both preparation jobs succeed. |
| Actual coverage | Both GitHub coverage jobs pass; the separate coverage workflow succeeds. |
| Full package check and lint | R CMD check stops with one `patchwork` dependency-declaration error; lint is still in progress at evidence capture. |
| GitHub workflow execution for this branch | Both validation workflows ran against `a637050` and failed in preflight; no successful package-validation result is claimed. |

The early isolated harness tests validate preflight control flow, not image contents. The current environment checks above execute the real R 4.1.3 Docker runtime. Existing base packages and package-version constraints were not changed; only the nine explicitly pinned missing additions were installed into a separate temporary library.

### First GitHub Execution

- [R-CMD-check run 40390514732](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40390514732): package-check, lint, and coverage jobs all fail during preflight. The aggregate job correctly reports failure.
- [Separate coverage run 40390514752](https://github.com/NIDAP-Community/SCWorkflow/actions/runs/40390514752): the same runtime namespace failures and unavailable `covr` stop preflight.
- Both runs test commit `a6370507ef928aaf6c9d470ed16a097fbdd5bcd6`. Image pulling, container startup, checkout, R execution, and environment artifact upload work. No package tests or coverage calculations start.

The `frozen-check-environment` artifact records these gaps in the active R libraries:

| Category | Packages not listed as installed |
| --- | --- |
| Failed namespace imports | `Seurat`, `celldex`, `quantmod`, `EnhancedVolcano`, `scDblFinder` |
| Validation tools | `rcmdcheck`, `covr`, `goodpractice`, `lintr` |
| Other inspected packages | `SeuratObject`, `ggalluvial` |

The inventory also records `ggplot2` 4.0.2 and `testthat` 3.3.2. The image tag is not evidence of an entirely historical dependency stack. `ggalluvial` is no longer a namespace-preflight blocker after the optional-dependency change.

The initial GitHub inventory alone did not distinguish an incomplete image from libraries outside `.libPaths()`. The local checks below now rule out the proposed `HOME` explanation and find no hidden installation of the blocking packages. Confirm the intended complete archived image before considering any separately approved additions; do not rebuild dependencies merely to diagnose them.

The Dockerfile at tag-source commit `79e5d37` uses `rocker/rstudio:4.1.3` and attempts versioned installs with CRAN/Bioconductor fallbacks. Its package lists describe intended installs, not proof that those packages were successfully installed in the published artifact.

### Local Filesystem And Library Verification

The exact pinned image was pulled and inspected on macOS using `linux/amd64`. Diagnostic containers used `--read-only --network none`; R's temporary directories were RAM-backed. No host source was mounted and no dependency installation or image rebuild was performed.

| Setting | Default environment | GitHub-style environment |
| --- | --- | --- |
| `HOME` | `/root` | `/github/home` |
| R executable | `/usr/local/bin/Rscript` | `/usr/local/bin/Rscript` |
| R home | `/usr/local/lib/R` | `/usr/local/lib/R` |
| `.libPaths()` | `/usr/local/lib/R/site-library`, `/usr/local/lib/R/library` | Same two paths |
| `R_LIBS` | `/usr/local/lib/R/site-library:/usr/local/lib/R/library` | Same value |
| Packages absent from the active libraries | All eleven packages listed above | Same eleven packages |

A full `/` scan for R libraries and the five failing namespace imports, `SeuratObject`, and four CI-tool package directories found only the standard R libraries, plus unrelated directory-name matches in Boost headers and lifecycle configuration. None of the blocking package installations was found elsewhere. A separate scan for `.tar.gz`, `.tar.xz`, `.tgz`, and `.zip` files found only installed packages' test/example data, not a preserved dependency source or library bundle.

Conclusion: changing `HOME` or adding a library path cannot expose packages that are not installed in this image. Leave the workflow's library settings unchanged. The next decision requires the image name/tag/digest for the intended complete archived environment, not a fresh dependency build. The separate `Dockerfiles/Plain_R_container/Dockerfile` recipe uses `/renv_cache`, but does not identify a published image to select; the inspected pinned image has no corresponding complete library. No replacement image has been selected or built.

### User-Supplied Preserved SIF

The user supplied `/Volumes/CCBR/NIDAP/SIFS/rs413s_latest.sif` for read-only inspection. Original SIF files are immutable: never modify, overwrite, rebuild in place, delete, move, or change their permissions or ownership. If writable image access is required, first create a distinct scratch copy and verify its SHA-256 against the original; operate only on that copy. All extraction, conversion, build, and overlay outputs must remain separate from the original.

| Identity | Observed value |
| --- | --- |
| Original SHA-256 before and after inspection | `37c4797e39de97680ae46cf9282d654322088ba3a97657a3e76ac010648b3899` (unchanged) |
| SIF identifier | `8d04c93e-624e-46c1-ba8f-52f39e90d073` |
| SIF creation/modification timestamp | `2025-12-19 13:41:07 UTC` |
| Filesystem | SquashFS, amd64, partition 4 at byte offset 40960 |
| Embedded Docker source | `rfinney/rs413s` |
| Current public source-image digest | `rfinney/rs413s@sha256:7cb32eaac452ad32abb34e9a9e3edb42590e4dfcf696caeb72bb375a5b3c83cf` |

The Docker source reference was subsequently executed directly. Its key runtime versions match the SIF and all 43 namespace imports load. Full filesystem equivalence to the SIF is not asserted; the selected immutable Docker digest is recorded separately.

Inspection used Apptainer 1.5.5 tooling pinned to `ghcr.io/apptainer/apptainer@sha256:a2016f1606a4f79d24a0578424cd63d490141f86bbec384dfdcdc19c33dd0d14`. Docker mounted the original SIF read-only and disabled container networking. The filesystem was extracted to a separate Docker volume, omitting SELinux extended attributes after Docker rejected their preservation. Runtime inspection then used that copy mounted read-only, RAM-backed temporary/device paths, and a read-only source checkout. The runtime inspection did not mount the original SIF at all. No R dependency was compiled, installed, or updated.

| Package/runtime | Version in the SIF copy |
| --- | --- |
| R | 4.1.3 |
| Seurat / SeuratObject | 4.1.1 / 4.1.1 |
| Matrix | 1.5-1 in the site library; 1.4-0 in the base library |
| ggplot2 | 3.3.6 |
| celldex | 1.4.0 |
| quantmod | 0.4.20 |
| EnhancedVolcano | 1.12.0 |
| scDblFinder | 1.8.0 |
| rcmdcheck / testthat | 1.4.0 / 3.1.6 |

R used `/usr/local/lib/R/site-library` and `/usr/local/lib/R/library`. All 43 packages imported by the current checked-out `NAMESPACE` loaded successfully through `loadNamespace()`, as did `rcmdcheck` and `testthat`. The active libraries do not contain `covr`, `goodpractice`, `lintr`, or `ggalluvial`. Coverage and lint therefore still lack their required tools; comparison-function tests can still require optional `ggalluvial`.

This is a copied-filesystem R runtime check, not a full Apptainer launch or package-validation run. The local chroot diagnostic emitted `/dev/null` permission messages and a missing `/proc/stat` warning; the import and tool checks nevertheless completed. No full package check, tests, lint, coverage, or GitHub run against this SIF is claimed.

The subsequent Docker verification and isolated-tool preparation are described above. CI now selects that public Docker digest, not the local SIF; no SIF upload or Apptainer runner access is required. The original SIF remains immutable and unchanged.

## Collaborator Verification Procedure

### Review The Revision

The branch is published as draft PR #107. Fetch `DEV-frozen-container-ci` and record the exact revision under review:

```bash
git fetch origin DEV-frozen-container-ci
git switch DEV-frozen-container-ci
git rev-parse HEAD
git diff origin/DEV -- .github/workflows .github/scripts .Rbuildignore .gitignore
```

Confirm that function changes are limited to the optional dependency guard, qualified `geom_flow` call, and corresponding roxygen/namespace import removal. Dependency declarations and `.github/package-versions.txt` are unchanged. The historical pins remain Seurat 4.1.1 and Matrix 1.5.1; compare the actual image inventory with the intended frozen environment rather than assuming that the tag guarantees these versions.

### Inspect And Run The Image

From the package root, with Docker running, first obtain the `frozen-ci-tools` artifact from a successful preparation job, or run the source-download and preparation steps in the shared workflow. Set `tools_directory` to the absolute directory containing its `library.tar.gz`. On Apple Silicon, use `--platform linux/amd64` for preparation as well as verification.

```bash
image='rfinney/rs413s@sha256:7cb32eaac452ad32abb34e9a9e3edb42590e4dfcf696caeb72bb375a5b3c83cf'
tools_directory='/absolute/path/to/frozen-ci-tools'
docker buildx imagetools inspect "$image"
docker pull --platform linux/amd64 "$image"
docker run --rm --platform linux/amd64 --read-only --network none \
  --tmpfs /tmp:rw,exec,size=1g \
  --mount "type=bind,source=$PWD,target=/source,readonly" \
  --mount "type=bind,source=$tools_directory,target=/ci-tools,readonly" \
  --env R_LIBS_USER=/tmp/restored/library --entrypoint /bin/sh "$image" -ec '
    mkdir /tmp/restored /tmp/preflight
    tar -xzf /ci-tools/library.tar.gz -C /tmp/restored
    cp /source/DESCRIPTION /source/NAMESPACE /tmp/preflight/
    cd /tmp/preflight
    for mode in check lint coverage; do
      Rscript /source/.github/scripts/ci-preflight.R "$mode"
    done
  '
```

The final command checks all three modes. Review the base inventory and source manifest, especially Seurat, Matrix, testthat, rcmdcheck, covr, lintr, goodpractice, and optional ggalluvial. The original SIF is not accessed by these commands. Apple Silicon hosts require amd64 emulation; GitHub's selected Ubuntu runners are amd64.

These preflight commands do not run the package test suite. Source and prepared-library mounts are read-only; verification artifacts are temporary and removed with the container.

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
- [x] Local runtime inventory confirms the recorded R, Seurat, Matrix, and validation-tool versions.
- [ ] The reported `ggalluvial` gap is confirmed in the image; package loading and the optional-function error are reviewed separately from full test-suite results.
- [x] All preflight modes succeed locally after prepared-library restoration, without installation during validation.
- [x] GitHub runs against the recorded revision and produces the prepared artifact; restoration and preflight pass in all four validation jobs.
- [x] Package checks and coverage execute against the checked-out source; the package metadata failure is investigated and reported rather than bypassed.
- [ ] Aggregate status and any required-status configuration are verified.
- [ ] Any remaining runtime, metadata, data-access, or Codecov failures are recorded separately.

## Limits And Next Steps

Passing in a frozen container demonstrates compatibility with that environment. It does not prove a fresh installation from `DESCRIPTION`, correctness of all dependency declarations, support for newer R/Seurat, or ongoing security support for the legacy stack. Dependency-metadata gaps identified in the release audit remain open and must not be hidden by the preinstalled libraries.

The incomplete original CI image is replaced in the published task branch with the verified public Docker source and isolated pinned tool preparation. Environment checks, lint, and both coverage jobs pass on GitHub, and the published patchwork correction resolves the original error. Optional trajectory dependency declarations pass a focused local loading check; the unchanged example then fails on missing PCA input. Fresh GitHub verification of those declarations, the example/input issue, and remaining package warnings must be reviewed separately. No tests were changed or skipped to force success. The SIF and every pre-existing base package remain untouched. No image was rebuilt or published. The draft PR must not be merged without authorization; broader package cleanup remains separate.
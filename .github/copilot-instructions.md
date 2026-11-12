# SCWorkflow AI Coding Instructions

## Project Overview
SCWorkflow is an R package for single-cell RNA-seq analysis built on the Seurat framework. It's designed for analyzing multimodal 10x Genomics data, with support for CITE-Seq, cell hashing, and TCR-seq data. The package is deployed as both an R package and Docker container for use in NIDAP (Palantir Foundry) and Biowulf HPC environments.

## Core Architecture

### Sequential Workflow Pattern
Functions follow a numbered workflow sequence:
1. **processRawData()** - Process H5 files into Seurat objects
2. **filterQC()** - Quality control and filtering  
3. **combineNormalize()** - Merge samples, normalize, dimension reduction
4. **Harmony integration** (optional) - Batch correction
5. **annotateCellTypes()** - Automatic cell type annotation via SingleR
6. **Analysis functions** - DEG, visualization, clustering

### Seurat Object as Central Data Structure
- All functions expect/return Seurat objects as primary data containers
- Functions often modify objects in-place and return modified versions
- Metadata is heavily used for sample tracking and analysis parameters
- Multiple assays supported: RNA, SCT (SCTransform), ADT (CITE-seq), etc.

### Function Design Patterns
- **Comprehensive parameter lists**: Most functions have 15-30+ parameters with sensible defaults
- **Conditional workflows**: Functions check input types (H5 vs Seurat) and adapt behavior
- **Multi-sample handling**: Functions can process lists of Seurat objects or file paths
- **Plotting integration**: Most analysis functions return both data and visualization outputs

## Key Conventions

### File Organization
- `R/` contains one function per file, named descriptively (e.g., `Process_Raw_Data.R`)
- Function names use camelCase: `processRawData()`, `annotateCellTypes()`
- File names use Snake_Case with capitalization: `Process_Raw_Data.R`

### Parameter Naming Patterns
- `object` - Primary Seurat object input
- `samples.to.include` - Character vector for sample subsetting
- `reduction.type` - Visualization method ("umap", "tsne", "pca")
- `organism` - Species specification ("Human" or "Mouse")
- Boolean parameters use `.` separator: `do.normalize.data`, `draw.umap`

### Documentation Standards
- Extensive roxygen2 documentation with `@details` sections explaining workflow step numbers
- `@importFrom` statements for specific function imports
- `@export` for all user-facing functions
- Parameter descriptions include defaults and valid options

## Development Workflows

## Output standards
- Generate R CMD check friendly code.
- Use roxygen2 for exported functions.
- Add utils::globalVariables() for NSE variables.
- Prefer explicit namespaces (dplyr::, ggplot2::).
- Avoid side effects in tests.

### Testing Structure
- Comprehensive `tests/testthat/` with both unit tests and integration tests
- `fixtures/` directory contains real Seurat objects for testing
- Helper functions in `helper-*.R` files for test setup
- Tests follow naming convention: `test-Function_Name.R`

#### Raw-data integration tests

The raw-data integration cases in `test-Process_Raw_Data.R` run sequentially,
each in a fresh R subprocess, through `testRawDataIsolated()` in the matching
helper. This releases Seurat objects and process allocations between cases.
The same isolation applies to local tests, R CMD check, GitHub coverage, and
Gitflow validation; no separate case list or opt-in flag is needed. Original
assertions and skip rules remain in the test file, and child expectations are
reported by the parent testthat runner. Subprocess logs are temporary files
under `file.path(tempdir(), "Process_Raw_Data-isolated")`.

To run only these cases from the package root in a configured R environment:

```r
devtools::test(filter = "Process_Raw_Data")
```

### Test fixture releases

All files from `tests/testthat/fixtures/`, including previously untracked RDS
and metadata files, are hosted in the versioned
[test-fixtures-v1 release](https://github.com/NIDAP-Community/SCWorkflow/releases/tag/test-fixtures-v1).
They are not tracked by Git or included in package source archives.
`tests/fixture-bundle.json` pins the release tag, asset name, archive size and
SHA-256, and the relative path, size, and SHA-256 of every fixture file.

Normal testthat runs restore missing fixtures automatically through
`helper-00-fixtures.R`. This applies to `devtools::test()`, individual
`testthat::test_file()` calls with helpers enabled, R CMD check, and coverage.
The first restore uses `piggyback::pb_download()` and needs internet access;
public downloads need no GitHub token. Install the `piggyback` test dependency
before downloading fixtures locally.
Existing fixtures are verified and reused without downloading. Local files
that differ from the pin cause an error rather than being overwritten.

To prepare explicitly from the package root:

```sh
Rscript tests/setup-fixtures.R
```

To prepare the verified archive before an offline test run:

```sh
Rscript tests/setup-fixtures.R --download-only
```

The archive cache defaults to `tools::R_user_dir("SCWorkflow", "cache")`.
Set `SCWORKFLOW_FIXTURE_CACHE` to choose another cache directory, or
`SCWORKFLOW_FIXTURE_ARCHIVE` to use an already-downloaded archive. The latter
must match the pinned size and checksum. Test outputs remain temporary files;
the ignored fixture directory is only a local input-data copy.

GitHub check, coverage, and Gitflow jobs prefetch the archive through the shared
NIDAP setup action, with an Actions cache keyed by the manifest. The local
Podman check runner prefetches into its run-artifact directory. Both use the
checksum-pinned piggyback 0.1.5 source in the isolated CI tool library; no
base-image packages are updated. The fixture wrapper retains checksum and
archive-path validation, offline reuse, and active test-directory restoration.

Fixture versions are independent of package releases. To change the data:

1. Prepare and validate the intended local fixture files on a task branch.
2. Run `.github/scripts/bundle-test-fixtures.R` with `Rscript`, an output
  directory outside the package repository, and a new tag such as
  `test-fixtures-v2`. It includes every fixture file, including ignored files,
  and generates an archive, `fixture-bundle.json`, and `SHA256SUMS` externally.
3. Publish those three assets as a data-only GitHub release at an already
  published commit, without making it the latest package release. Do not
  overwrite published assets. GitHub CLI supports this repository's existing
  draft releases; piggyback 0.1.5 currently fails while listing those drafts.
4. Verify a fresh public download against the generated checksums, then copy
  only the generated manifest to `tests/fixture-bundle.json` and test both
  restoration and offline reuse. Never force-add the fixture directory or
  release archive to Git. Review the small manifest change with the code.

Moving fixtures out of tracked files does not rewrite existing Git history.

### Local Podman resources

The frozen image fixes R and package versions, not resource allocation. Local
runner resource settings do not change GitHub CI, test selection, assertions,
or fixture data. Unset settings preserve the existing runner behavior.

The VM and the check container have separate limits. On the current 12-core,
18 GiB Mac, a starting VM allocation is 8 CPUs and 12 GiB RAM, leaving host
headroom. This is not a measured speedup. Change the VM only after every active
container has finished; restarting it interrupts running validation:

```sh
podman ps
podman machine stop scworkflow-ci
podman machine set --cpus 8 --memory 12288 scworkflow-ci
podman machine start scworkflow-ci
podman machine ssh scworkflow-ci free -m
```

Check available swap after restarting. Existing temporary VM swap may need to
be re-enabled separately; do not recreate or overwrite an existing swap file.
`SCWORKFLOW_PODMAN_MEMORY_SWAP` is the total RAM-plus-swap container limit, not
additional swap, and it does not create swap for the VM.

After the VM allocation and swap are ready, an optional check profile is:

```sh
SCWORKFLOW_PODMAN_CPUS=8 \
SCWORKFLOW_PODMAN_MEMORY=10g \
SCWORKFLOW_PODMAN_MEMORY_SWAP=16g \
SCWORKFLOW_PODMAN_THREADS=2 \
bash .github/scripts/check-r-package-podman.sh
```

The resource controls apply only to the final check container; the small
preparation/preflight containers retain their existing limits. The thread
setting controls OpenBLAS, OpenMP, MKL, and RcppParallel for that process. It
does not register parallel R/Bioconductor workers or parallelize test files.
In particular, QC doublet detection explicitly uses `SerialParam`, and many
test cases run sequentially. More RAM alone does not speed up CPU-bound work.

Native Mac R uses Apple Accelerate while this AMD64 image uses OpenBLAS under
emulation. Resource tuning does not make those implementations identical.
Benchmark thread counts and any future isolated parallel test runner before
claiming an improvement; higher counts can be slower under emulation.

### Maintaining pinned NIDAP CI packages

GitHub Actions uses a digest-pinned NIDAP image with R 4.1.3. The
[source manifest](scripts/ci-tool-sources.json) pins additional packages
installed into an isolated CI library, including required runtime packages
missing from the image. It is not a lockfile for the entire base environment.
Existing base packages, including Seurat, must not be upgraded or shadowed.

Pins are maintainer-reviewed, not automatically refreshed to the latest
versions. Change them only for a documented need, such as a new required tool,
a compatibility fix, or a security fix. Use this process:

1. Explain the reason and proposed versions in a pull request. Check each
  package's R and dependency requirements against R 4.1.3 and the packages
  already in the pinned image. If an update requires changing a base package
  or the image, stop and discuss that as a separate environment change.
2. Obtain the exact source archive from a trusted HTTPS source, such as the
  CRAN archive. Compute its SHA-256 with `shasum -a 256 package_version.tar.gz`.
  Update `Version`, `URL`, and `SHA256` together in the same manifest change;
  the URL and checksum must identify that exact package version. Verify these
  values from the actual archive rather than trusting generated metadata.
3. Add any missing dependencies explicitly, before the packages that need
  them, with compatible versions and checksums. The
  [preparer](scripts/ci-prepare-tools.R) installs with `dependencies = FALSE`;
  it does not resolve dependencies or replace packages already in the base
  image. Do not use a manifest entry to override a base-package version.
4. Push the change for validation on GitHub Actions. Require the
  [NIDAP preparation workflow](workflows/NIDAP-environment.yaml) to pass source
  checksum verification, offline installation, and all preflight modes. Review
  the downstream package-check, lint, coverage, and Gitflow results under their
  existing failure policies. Successful preparation alone is not a successful
  package check; do not bypass checksum checks or relax validation gates.
5. Inspect the `NIDAP-ci-tools` artifact. Confirm that `tools-manifest.csv`
  contains only the intended pin changes and that `base-packages.csv` matches
  the last known-good run using the same image digest. The preparer also
  compares the full base inventory before and after installation and rejects
  additions that shadow base packages. Record the comparison and GitHub run
  links in the pull request before merging.

Keep the maintenance rationale and validation evidence in the pull request,
not in user-facing release notes. Do not add a NEWS entry for CI-only changes.

### Deliverables
- Update NEWS.md only for user-facing changes; keep CI details and internal investigations out of release notes.
- Keep AI scratch artifacts outside the package repository, following the rules below.
- Add or update testthat tests for changed package behavior.

### Documentation placement
- Before editing documentation, identify its intended audience and choose the existing file that serves that audience; do not append details merely because a related heading is nearby.
- Keep `.github/CONTRIBUTING.md` focused on new human contributors: setup, ordinary test commands, contribution workflow, and brief descriptions and links for required data assets.
- Put advanced agent and maintainer guidance in `.github/copilot-instructions.md`, including raw-data subprocess isolation, fixture release and offline-cache procedures, local Podman resource tuning, and pinned CI package maintenance.
- Keep change-specific rationale, investigation results, benchmarks, and validation evidence in the pull request discussion, not in newcomer documentation or user-facing release notes. Keep personal scratch notes in the external location defined below.
- Prefer one detailed source of guidance with a short summary and link elsewhere. Before finishing, verify that each edited section fits its destination's audience and that moved links still resolve.

### AI Scratch Artifacts
- Store personal plans, scratch prompts, session notes, transcripts, research summaries, completion reports, and decision logs at `.agents/SCWorkflow/` under the multi-repository workspace root, outside the `SCWorkflow/` repository. For a standalone checkout, use an external local scratch directory instead.
- Create scratch files only when needed; prefer chat for status updates and task summaries. Never force-add scratch artifacts or copy them into tracked README.md, NEWS.md, or docs/ files.
- Classify files by purpose, not AI authorship or directory name. Code, tests, and requested project documentation are deliverables, not scratch artifacts.
- Keep intentional shared instructions, skills (including `.agents/skills/`), reusable contributor prompts, and workflows tracked. Do not blanket-ignore `.agents/`; ignore only specific local scratch paths if needed.
- Add new shared agent configuration or promote scratch material into project documentation only when explicitly requested.
- Before staging or committing, review changed and untracked files and the staged file list. Stage explicit task files, not `git add .` or `git add -A`, and exclude scratch artifacts.
- Do not delete or relocate existing user files merely because they appear AI-generated; get approval unless that cleanup was explicitly requested.

### Docker Development
- Base image: `nciccbr/ccbr_ubuntu_22.04:v4`
- Conda-based R environment (R 4.3.2)
- Multi-stage build supporting both development and production
- Container designed for HPC environments (Biowulf) and cloud platforms (NIDAP)

### CI/CD Integration  
- GitFlow-based R package workflow via `gitflow-R-action.yml`
- Automatic NIDAP deployment on successful builds
- pkgdown documentation generation
- Docker image building and deployment

## Critical Dependencies and Integration Points

### Bioconductor/Seurat Ecosystem
- **Seurat 4.1.1+**: Core framework for single-cell analysis
- **SingleR**: Automated cell type annotation
- **MAST/limma/edgeR**: Differential expression analysis
- **Harmony**: Batch effect correction and integration

### External API Integration
- `palantir_api_call.R`: Custom integration with Palantir Foundry APIs
- Designed for NIDAP cloud environment deployment
- Handles authentication and data transfer patterns specific to NIH infrastructure

### Memory Management Patterns
- `cell.count.limit` parameters (default: 35000) trigger memory conservation
- `only.var.genes` option reduces memory footprint for large datasets
- SCTransform normalization levels: sample-wise vs merged strategies

## Common Debugging Areas

### Data Type Handling
- Functions check for H5 vs Seurat object inputs: `is(input, "Seurat")`
- Cell filtering can drastically reduce object sizes - verify cell counts
- Assay availability: functions may require specific assays (RNA, SCT, ADT)

### Visualization Rendering
- Plot functions return both Seurat objects AND plot objects
- `reduction.type` must match available reductions in object
- Color palettes automatically generated but can be overridden

### Environment-Specific Issues
- Package designed for both local R and containerized environments
- File path handling differs between NIDAP, Biowulf, and local development
- Some functions have environment-specific conditional logic

## Quick Start for AI Agents
When working with this codebase:
1. Always check if input is a Seurat object: `class(object)`
2. Verify required reductions exist: `object@reductions`
3. Check available assays: `names(object@assays)`
4. Use `str(object@meta.data)` to understand metadata structure
5. Test functions start with small datasets from `tests/testthat/fixtures/`

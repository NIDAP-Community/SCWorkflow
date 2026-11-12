# Contributing to SCWorkflow

## Proposing changes with issues

If you want to make a change, it's a good idea to first
[open an issue](https://code-review.tidyverse.org/issues/)
and make sure someone from the team agrees that it’s needed.

If you've decided to work on an issue,
[assign yourself to the issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/assigning-issues-and-pull-requests-to-other-github-users#assigning-an-individual-issue-or-pull-request)
so others will know you're working on it.

## Pull request process

We use [GitHub Flow](https://docs.github.com/en/get-started/using-github/github-flow)
as our collaboration process.
Follow the steps below for detailed instructions on contributing changes to
SCWorkflow.

![GitHub Flow diagram](https://raw.githubusercontent.com/CCBR/CCBR_NextflowTemplate/main/.github/img/GitHub-Flow_bg-white.png)


### Clone the repo

If you are a member of [NIDAP-Community](https://github.com/NIDAP-Community),
you can clone this repository to your computer or development environment.
Otherwise, you will first need to
[fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo)
the repo and clone your fork. You only need to do this step once.

```sh
git clone https://github.com/NIDAP-Community/SCWorkflow
```

> Cloning into 'SCWorkflow'... <br>
> remote: Enumerating objects: 1136, done. <br>
> remote: Counting objects: 100% (463/463), done. <br>
> remote: Compressing objects: 100% (357/357), done. <br>
> remote: Total 1136 (delta 149), reused 332 (delta 103), pack-reused 673 <br>
> Receiving objects: 100% (1136/1136), 11.01 MiB | 9.76 MiB/s, done. <br>
> Resolving deltas: 100% (530/530), done. <br>

```sh
cd SCWorkflow
```

### If this is your first time cloning the repo, install dependencies

- In an R console, install the R development dependencies with
  `devtools::install_dev_deps()`, and then make sure the package passes R CMD
  check by running `devtools::check()`. If R CMD check doesn't pass cleanly,
  it's a good idea to ask for help before continuing.

- Install [`pre-commit`](https://pre-commit.com/#install) if you don't already
  have it. Then from the repo's root directory, run

  ```sh
  pre-commit install
  ```

  This will install the repo's pre-commit hooks.
  You'll only need to do this step the first time you clone the repo.

### Create a branch

  Create a Git branch for your pull request (PR). Give the branch a descriptive
  name for the changes you will make, such as `iss-10` if it is for a specific
  issue.

  ```sh
  # create a new branch and switch to it
  git branch iss-10
  git switch iss-10
  ```

  > Switched to a new branch 'iss-10'

### Make your changes

Edit the code, write unit tests, and update the documentation as needed.

#### style

New code should follow the [tidyverse style guide](https://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with
your PR.

A brief overview of conventions according to the tidyverse style guide:

- most object names (variables and functions) should be in [snake_case](https://style.tidyverse.org/syntax.html#sec-objectnames)
- function names should use [verbs](https://style.tidyverse.org/functions.html#naming) where possible
- use `<-` for assignment
- use [pipes](https://style.tidyverse.org/pipes.html) to chain operations on a single object

Please see the [tidyverse style guide](https://style.tidyverse.org) for more details.

#### test

Most changes to the code will also need unit tests to demonstrate that the
changes work as intended.
Use [`testthat`](https://testthat.r-lib.org/) to create your unit tests and test
the code.
Test files are organized as described in
<https://style.tidyverse.org/tests.html>.
Take a look at the existing code in this package for examples.

To run the raw-data integration tests from the package root in a configured
R environment:

```r
devtools::test(filter = "Process_Raw_Data")
```

#### document

If you have written a new function or changed the API of an existing function,
you will need to update the function's documentation using
[roxygen2](https://cran.r-project.org/package=roxygen2) with
[Markdown syntax](https://roxygen2.r-lib.org/articles/rd-formatting.html).
See instructions on writing roxygen2 comments here:
<https://r-pkgs.org/man.html>.
If the function is used in a vignette, you may also need to update the vignette.

#### check

After making your changes, run `devtools::check()` from an R console to make
sure the package still passes R CMD check.

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

### Commit and push your changes

If you're not sure how often you should commit or what your commits should
consist of, we recommend following the "atomic commits" principle where each
commit contains one new feature, fix, or task.
Learn more about atomic commits here:
<https://www.freshconsulting.com/insights/blog/atomic-commits/>

First, add the files that you changed to the staging area:

```sh
git add path/to/changed/files/
```

Then make the commit.
Your commit message should follow the
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/)
specification.
Briefly, each commit should start with one of the approved types such as
`feat`, `fix`, `docs`, etc. followed by a description of the commit.
Take a look at the [Conventional Commits specification](https://www.conventionalcommits.org/en/v1.0.0/#summary)
for more detailed information about how to write commit messages.


```sh
git commit -m 'feat: create function for awesome feature'
```

pre-commit will enforce that your commit message and the code changes are
styled correctly and will attempt to make corrections if needed.

> Check for added large files..............................................Passed <br>
> Fix End of Files.........................................................Passed <br>
> Trim Trailing Whitespace.................................................Failed <br>
> - hook id: trailing-whitespace <br>
> - exit code: 1 <br>
> - files were modified by this hook <br>
> <br>
> Fixing path/to/changed/files/file.txt <br>
> <br>
> codespell................................................................Passed <br>
> style-files..........................................(no files to check)Skipped <br>
> readme-rmd-rendered..................................(no files to check)Skipped <br>
> use-tidy-description.................................(no files to check)Skipped <br>

In the example above, one of the hooks modified a file in the proposed commit,
so the pre-commit check failed. You can run `git diff` to see the changes that
pre-commit made and `git status` to see which files were modified. To proceed
with the commit, re-add the modified file(s) and re-run the commit command:

```sh
git add path/to/changed/files/file.txt
git commit -m 'feat: create function for awesome feature'
```

This time, all the hooks either passed or were skipped
(e.g. hooks that only run on R code will not run if no R files were
committed).
When the pre-commit check is successful, the usual commit success message
will appear after the pre-commit messages showing that the commit was created.

> Check for added large files..............................................Passed <br>
> Fix End of Files.........................................................Passed <br>
> Trim Trailing Whitespace.................................................Passed <br>
> codespell................................................................Passed <br>
> style-files..........................................(no files to check)Skipped <br>
> readme-rmd-rendered..................................(no files to check)Skipped <br>
> use-tidy-description.................................(no files to check)Skipped <br>
> Conventional Commit......................................................Passed <br>
> [iss-10 9ff256e] feat: create function for awesome feature <br>
> 1 file changed, 22 insertions(+), 3 deletions(-) <br>

Finally, push your changes to GitHub:

```sh
git push
```

If this is the first time you are pushing this branch, you may have to
explicitly set the upstream branch:

```sh
git push --set-upstream origin iss-10
```

> Enumerating objects: 7, done. <br>
> Counting objects: 100% (7/7), done. <br>
> Delta compression using up to 10 threads <br>
> Compressing objects: 100% (4/4), done. <br>
> Writing objects: 100% (4/4), 648 bytes | 648.00 KiB/s, done. <br>
> Total 4 (delta 3), reused 0 (delta 0), pack-reused 0 <br>
> remote: Resolving deltas: 100% (3/3), completed with 3 local objects. <br>
> remote: <br>
> remote: Create a pull request for 'iss-10' on GitHub by visiting: <br>
> remote: https://github.com/NIDAP-Community/SCWorkflow/pull/new/iss-10 <br>
> remote: <br>
> To https://github.com/NIDAP-Community/SCWorkflow <br>
> <br>
>   [new branch] iss-10 -> iss-10 <br>
>   branch 'iss-10' set up to track 'origin/iss-10'. <br>

We recommend pushing your commits often so they will be backed up on GitHub.
You can view the files in your branch on GitHub at
`https://github.com/NIDAP-Community/SCWorkflow/tree/<your-branch-name>`
(replace `<your-branch-name>` with the actual name of your branch).

### Create the PR

Once your branch is ready, create a PR on GitHub:
<https://github.com/NIDAP-Community/SCWorkflow/pull/new/>

Select the branch you just pushed:

![Create a new PR from your branch](https://raw.githubusercontent.com/CCBR/CCBR_NextflowTemplate/main/.github/img/new-PR.png)

Edit the PR title and description.
The title should briefly describe the change.
Follow the comments in the template to fill out the body of the PR, and
you can delete the comments (everything between `<!--` and `-->`) as you go.
When you're ready, click 'Create pull request' to open it.

![Open the PR after editing the title and description](https://raw.githubusercontent.com/CCBR/CCBR_NextflowTemplate/main/.github/img/create-PR.png)

Optionally, you can mark the PR as a draft if you're not yet ready for it to
be reviewed, then change it later when you're ready.

### Wait for a maintainer to review your PR

We will do our best to follow the tidyverse code review principles:
<https://code-review.tidyverse.org/>.
The reviewer may suggest that you make changes before accepting your PR in
order to improve the code quality or style.
If that's the case, continue to make changes in your branch and push them to
GitHub, and they will appear in the PR.

Once the PR is approved, the maintainer will merge it and the issue(s) the PR
links will close automatically.
Congratulations and thank you for your contribution!

### After your PR has been merged

After your PR has been merged, update your local clone of the repo by
switching to the main branch and pulling the latest changes:

```sh
git checkout main
git pull
```

It's a good idea to run `git pull` before creating a new branch so it will
start from the most recent commits in main.

## Helpful links for more information

- This contributing guide was adapted from the [tidyverse contributing guide](https://github.com/tidyverse/tidyverse/blob/main/.github/CONTRIBUTING.md)
- [GitHub Flow](https://docs.github.com/en/get-started/using-github/github-flow)
- [tidyverse style guide](https://style.tidyverse.org)
- [tidyverse code review principles](https://code-review.tidyverse.org)
- [reproducible examples](https://www.tidyverse.org/help/#reprex)
- [R packages book](https://r-pkgs.org/)
- packages:
  - [usethis](https://usethis.r-lib.org/)
  - [devtools](https://devtools.r-lib.org/)
  - [testthat](https://testthat.r-lib.org/)
  - [styler](https://styler.r-lib.org/)
  - [roxygen2](https://roxygen2.r-lib.org)

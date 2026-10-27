#!/usr/bin/env bash
# Run SCWorkflow R CMD check locally in the same pinned environment as CI.
# Usage: .github/scripts/check-r-package-podman.sh [--help]
# Optional: set SCWORKFLOW_PODMAN_CHECK_DIR to a new directory for the run artifacts.

set -Eeuo pipefail

usage() {
  cat <<'USAGE'
Usage: .github/scripts/check-r-package-podman.sh [--help]

Prepares the checksum-pinned CI tool library and runs the check-mode preflight
and R CMD check in the same digest-pinned R 4.1.3 image used by GitHub Actions.

Requirements: Podman with a running machine and sufficient VM disk space on
macOS, jq, curl, and shasum.

By default, sources, the prepared tools library, preflight evidence, and check
logs are retained under a unique ci-artifacts/podman-r-cmd-check.* directory.
Set SCWORKFLOW_PODMAN_CHECK_DIR to choose a new output directory explicitly.
USAGE
}

fail() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--help" ) ]]; then
  usage >&2
  exit 2
fi
if [[ $# -eq 1 ]]; then
  usage
  exit 0
fi

for required_command in podman jq curl shasum; do
  command -v "$required_command" >/dev/null 2>&1 || fail "Required command not found: $required_command"
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "$script_dir/../.." && pwd -P)"
source_manifest="$repo_root/.github/scripts/ci-tool-sources.json"
ci_image="rfinney/rs413s@sha256:7cb32eaac452ad32abb34e9a9e3edb42590e4dfcf696caeb72bb375a5b3c83cf"
podman_platform="${PODMAN_PLATFORM:-linux/amd64}"

[[ -f "$source_manifest" ]] || fail "Missing CI tool manifest: $source_manifest"
podman info >/dev/null 2>&1 || fail "Podman is unavailable. Start Podman Desktop or run 'podman machine start', then retry."

if [[ -n "${SCWORKFLOW_PODMAN_CHECK_DIR:-}" ]]; then
  work_dir="$SCWORKFLOW_PODMAN_CHECK_DIR"
  [[ "$work_dir" = /* ]] || work_dir="$repo_root/$work_dir"
  [[ ! -e "$work_dir" ]] || fail "Output directory already exists; choose a new SCWORKFLOW_PODMAN_CHECK_DIR: $work_dir"
  mkdir -p "$work_dir"
else
  mkdir -p "$repo_root/ci-artifacts"
  work_dir="$(mktemp -d "$repo_root/ci-artifacts/podman-r-cmd-check.XXXXXXXX")"
fi

sources_dir="$work_dir/sources"
output_dir="$work_dir/output"
mkdir -p "$sources_dir" "$output_dir"

printf 'Pulling CI image for %s: %s\n' "$podman_platform" "$ci_image"
podman pull --platform "$podman_platform" "$ci_image"

printf 'Downloading and verifying pinned CI tool sources...\n'
while IFS=$'\t' read -r package version checksum url; do
  [[ -n "$package" && -n "$version" && -n "$checksum" && -n "$url" ]] ||
    fail "Invalid row in CI tool manifest: $source_manifest"
  archive="$sources_dir/${package}_${version}.tar.gz"
  curl --fail --location --retry 3 --proto '=https' --proto-redir '=https' \
    "$url" --output "$archive"
  actual_checksum="$(shasum -a 256 "$archive" | awk '{print $1}')"
  [[ "$actual_checksum" == "$checksum" ]] ||
    fail "SHA-256 mismatch for $package $version (expected $checksum; got $actual_checksum)"
done < <(jq -r '.[] | [.Package, .Version, .SHA256, .URL] | @tsv' "$source_manifest")

printf 'Preparing isolated tools offline...\n'
podman run --rm --platform "$podman_platform" \
  --read-only --network none --tmpfs /tmp:rw,size=1g \
  --volume "$repo_root:/workspace:ro" \
  --volume "$work_dir:/ci-tools:rw" \
  --env R_LIBS_USER=/ci-tools/library \
  --entrypoint Rscript "$ci_image" \
  /workspace/.github/scripts/ci-prepare-tools.R /ci-tools

printf 'Verifying check, lint, and coverage preflights offline...\n'
podman run --rm --platform "$podman_platform" \
  --read-only --network none --tmpfs /tmp:rw,size=1g \
  --volume "$repo_root:/workspace:ro" \
  --volume "$work_dir:/ci-tools:ro" \
  --volume "$output_dir:/ci-output:rw" \
  --env R_LIBS_USER=/ci-tools/library \
  --workdir /tmp --entrypoint /bin/sh "$ci_image" -ec '
    mkdir -p /ci-output/preflight
    cp /workspace/DESCRIPTION /workspace/NAMESPACE /ci-output/preflight/
    cd /ci-output/preflight
    for mode in check lint coverage; do
      Rscript /workspace/.github/scripts/ci-preflight.R "$mode"
    done
  '

printf 'Running R CMD check in the CI image...\n'
podman run --rm --platform "$podman_platform" \
  --tmpfs /tmp:rw,size=1g \
  --volume "$repo_root:/workspace:ro" \
  --volume "$work_dir:/ci-tools:ro" \
  --volume "$output_dir:/ci-output:rw" \
  --env R_LIBS_USER=/ci-tools/library \
  --env R_KEEP_PKG_SOURCE=yes \
  --workdir /workspace --entrypoint /bin/sh "$ci_image" -ec '
    test "$(Rscript -e '\''cat(as.character(getRversion()))'\'')" = "4.1.3"
    Rscript -e '\''rcmdcheck::rcmdcheck(
      path = "/workspace",
      args = c("--no-manual", "--as-cran"),
      build_args = "--no-manual",
      error_on = "warning",
      check_dir = "/ci-output/check"
    )'\''
  '

printf 'R CMD check passed. Run artifacts: %s\n' "$work_dir"

#!/usr/bin/env bash
set -euo pipefail

# add scworkflow executable to the path
scworkflow=$(R -s -e "cat(system.file('exec','scworkflow', package='SCWorkflow'))")
export PATH="$PATH:$(dirname $scworkflow)"

# run functions from SCWorkflow
scworkflow processRawData --json=json_args/processRawData.json
scworkflow filterQC --json=json_args/filterQC.json

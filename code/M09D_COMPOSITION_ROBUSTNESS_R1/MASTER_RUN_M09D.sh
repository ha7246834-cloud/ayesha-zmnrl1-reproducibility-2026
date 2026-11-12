#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 PROJECT_ROOT [--resume /absolute/path/to/.STAGING_token]" >&2
  exit 2
fi

PROJECT_ROOT=$(realpath "$1")
shift
TOOL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON=/home/hammadse/miniconda3/envs/wrky_q1/bin/python
STAGING=""

if [[ $# -gt 0 ]]; then
  if [[ $# -ne 2 || "$1" != "--resume" ]]; then
    echo "Only --resume STAGING_DIR is accepted after PROJECT_ROOT" >&2
    exit 2
  fi
  STAGING=$(realpath "$2")
  "$PYTHON" "$TOOL_DIR/m09d_workflow.py" "$PROJECT_ROOT" prepare --staging "$STAGING"
else
  PREPARE_LOG=$(mktemp /tmp/m09d_prepare.XXXXXX.log)
  "$PYTHON" "$TOOL_DIR/m09d_workflow.py" "$PROJECT_ROOT" prepare | tee "$PREPARE_LOG"
  STAGING=$(sed -n 's/^M09D_STAGING=//p' "$PREPARE_LOG" | tail -n 1)
  rm -f "$PREPARE_LOG"
  if [[ -z "$STAGING" ]]; then
    echo "Could not determine M09D staging directory" >&2
    exit 1
  fi
fi

for ANALYSIS in R1 R2 R3 R4; do
  "$PYTHON" "$TOOL_DIR/m09d_workflow.py" "$PROJECT_ROOT" run \
    --staging "$STAGING" --analysis "$ANALYSIS"
done

export MPLCONFIGDIR=${TMPDIR:-/tmp}/m09d_matplotlib_cache
mkdir -p "$MPLCONFIGDIR"
"$PYTHON" "$TOOL_DIR/m09d_finalize.py" "$PROJECT_ROOT" build --staging "$STAGING"

echo "M09D_BUILD_COMPLETE_PENDING_VISUAL_QA=$STAGING"
echo "After independent visual inspection and PASS attestation, publish with:"
echo "$PYTHON $TOOL_DIR/m09d_finalize.py $PROJECT_ROOT publish --staging $STAGING"

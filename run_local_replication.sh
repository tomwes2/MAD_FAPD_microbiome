#!/usr/bin/env bash
# Replicate Paper 1 (paired fecal vs cecal) analyses for the integrated draft.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export MAD_DATA_DIR="${MAD_DATA_DIR:-$ROOT/MAD_analysis_replication}"
export MAD_OUTPUT_DIR="${MAD_OUTPUT_DIR:-$ROOT}"
cd "$MAD_DATA_DIR"
echo "MAD_DATA_DIR=$MAD_DATA_DIR"
echo "MAD_OUTPUT_DIR=$MAD_OUTPUT_DIR"
exec Rscript run_replication.R

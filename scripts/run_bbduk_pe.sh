#!/usr/bin/env bash
# Paired-end adapter/quality trimming with BBTools bbduk for MAD shotgun libraries.
#
# This study used bbduk-trimmed FASTQs for SingleM and HUMAnN profiling.
# Assembly and read-mapping were not used for the site-comparison analyses.
#
# Requirements: bbduk.sh on PATH (BBMap / BBTools; manuscript used 39.33).
#
# Usage:
#   ./scripts/run_bbduk_pe.sh /path/to/raw_fastq [/path/to/trimmed_out] [threads]
#
# Recognized paired-end name patterns (R1 required; R2 inferred):
#   *_R1.fastq.gz / *_R2.fastq.gz
#   *_r1.fq.gz    / *_r2.fq.gz
#   *_1.fastq.gz  / *_2.fastq.gz
#   (also works with .fastq / .fq uncompressed)
#
set -euo pipefail

IN_DIR="${1:-}"
OUT_DIR="${2:-./bbduk_trimmed}"
THREADS="${3:-8}"

if [[ -z "${IN_DIR}" || ! -d "${IN_DIR}" ]]; then
  echo "Usage: $0 <input_fastq_dir> [output_dir] [threads]" >&2
  exit 1
fi

if ! command -v bbduk.sh >/dev/null 2>&1; then
  echo "ERROR: bbduk.sh not found on PATH. Install BBMap/BBTools and retry." >&2
  exit 1
fi

# Prefer the adapters file shipped with BBMap; override with BBDUK_ADAPTERS if needed.
ADAPTERS="${BBDUK_ADAPTERS:-}"
if [[ -z "${ADAPTERS}" ]]; then
  BBDUK_BIN="$(command -v bbduk.sh)"
  BBMAP_HOME="$(cd "$(dirname "${BBDUK_BIN}")/.." && pwd)"
  for cand in \
      "${BBMAP_HOME}/resources/adapters.fa" \
      "$(dirname "${BBDUK_BIN}")/resources/adapters.fa" \
      "${CONDA_PREFIX:-}/opt/bbmap-*/resources/adapters.fa"
  do
    # shellcheck disable=SC2086
    for hit in ${cand}; do
      if [[ -f "${hit}" ]]; then
        ADAPTERS="${hit}"
        break 2
      fi
    done
  done
fi

if [[ -z "${ADAPTERS}" || ! -f "${ADAPTERS}" ]]; then
  echo "ERROR: could not find adapters.fa. Set BBDUK_ADAPTERS=/path/to/adapters.fa" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}" "${OUT_DIR}/logs"

echo "Input:     ${IN_DIR}"
echo "Output:    ${OUT_DIR}"
echo "Threads:   ${THREADS}"
echo "Adapters:  ${ADAPTERS}"
echo

shopt -s nullglob
R1_FILES=(
  "${IN_DIR}"/*_R1.fastq.gz
  "${IN_DIR}"/*_R1.fastq
  "${IN_DIR}"/*_r1.fq.gz
  "${IN_DIR}"/*_r1.fq
  "${IN_DIR}"/*_1.fastq.gz
  "${IN_DIR}"/*_1.fastq
  "${IN_DIR}"/*_1.fq.gz
  "${IN_DIR}"/*_1.fq
)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "ERROR: no R1 FASTQ files found in ${IN_DIR}" >&2
  exit 1
fi

pair_r2() {
  local r1="$1"
  case "${r1}" in
    *_R1.fastq.gz) echo "${r1/_R1.fastq.gz/_R2.fastq.gz}" ;;
    *_R1.fastq)    echo "${r1/_R1.fastq/_R2.fastq}" ;;
    *_r1.fq.gz)    echo "${r1/_r1.fq.gz/_r2.fq.gz}" ;;
    *_r1.fq)       echo "${r1/_r1.fq/_r2.fq}" ;;
    *_1.fastq.gz)  echo "${r1/_1.fastq.gz/_2.fastq.gz}" ;;
    *_1.fastq)     echo "${r1/_1.fastq/_2.fastq}" ;;
    *_1.fq.gz)     echo "${r1/_1.fq.gz/_2.fq.gz}" ;;
    *_1.fq)        echo "${r1/_1.fq/_2.fq}" ;;
    *) return 1 ;;
  esac
}

out_stem() {
  local r1="$1" base
  base="$(basename "${r1}")"
  base="${base%_R1.fastq.gz}"
  base="${base%_R1.fastq}"
  base="${base%_r1.fq.gz}"
  base="${base%_r1.fq}"
  base="${base%_1.fastq.gz}"
  base="${base%_1.fastq}"
  base="${base%_1.fq.gz}"
  base="${base%_1.fq}"
  echo "${base}"
}

n_ok=0
n_skip=0
for r1 in "${R1_FILES[@]}"; do
  r2="$(pair_r2 "${r1}")"
  stem="$(out_stem "${r1}")"
  if [[ ! -f "${r2}" ]]; then
    echo "SKIP (missing R2): ${r1}" >&2
    n_skip=$((n_skip + 1))
    continue
  fi

  out1="${OUT_DIR}/${stem}_R1.fastq.gz"
  out2="${OUT_DIR}/${stem}_R2.fastq.gz"
  stats="${OUT_DIR}/logs/${stem}.bbduk.stats.txt"
  log="${OUT_DIR}/logs/${stem}.bbduk.log"

  if [[ -s "${out1}" && -s "${out2}" ]]; then
    echo "SKIP (exists): ${stem}"
    n_skip=$((n_skip + 1))
    continue
  fi

  echo "bbduk: ${stem}"
  # Adapter trim (left) + quality trim + length filter.
  # Parameters match common BBTools PE shotgun defaults used with BBMap 39.x.
  bbduk.sh \
    in1="${r1}" in2="${r2}" \
    out1="${out1}" out2="${out2}" \
    ref="${ADAPTERS}" \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=rl trimq=20 \
    minlen=50 \
    threads="${THREADS}" \
    stats="${stats}" \
    2>"${log}"

  n_ok=$((n_ok + 1))
done

echo
echo "Finished. Trimmed pairs written: ${n_ok}; skipped: ${n_skip}"
echo "Logs/stats: ${OUT_DIR}/logs/"

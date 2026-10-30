#!/usr/bin/env bash
# =============================================================================
# F26-DATA-002: Robust download & validation pipeline for F26 cross-species
#               FASTQ datasets (fastBact Phase F26, Stage B).
#
# Slurm batch script. Submit from the worktree root with:
#     mkdir -p logs && sbatch scripts/download_f26_datasets.sh
#
# Per-accession workflow:
#   1. Resume check: skip accession iff both R1/R2 exist, pass `pigz -t`,
#      line counts divide by 4, and both mate read counts are > 0 and equal.
#   2. Download with `fasterq-dump --split-3 -e 8` (retried on network
#      dropouts: MAX_ATTEMPTS attempts with exponential backoff).
#   3. Compress with `pigz -p 8`; verify archives with `pigz -t`.
#   4. Count exact reads per mate (line_count / 4); require R1 == R2 > 0.
#   5. Append SHA-256, byte size, and read counts for every compressed FASTQ
#      to the persistent checksum ledger ${OUTDIR}/f26_sha256_checksums.tsv.
#
# Downstream exact validation (ID synchronization, JSON report):
#     python3 scripts/validate_fastq_pairs.py --fastq-dir "${OUTDIR}" --manifest F26_DATASET_MANIFEST.tsv
# ==============================================================================
#SBATCH -p qfree
#SBATCH -w vnode1
#SBATCH -J f26_download
#SBATCH -c 8
#SBATCH --mem=14G
#SBATCH -t 08:00:00
#SBATCH -o /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree/logs/download_f26_%j.log
#SBATCH -e /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree/logs/download_f26_%j.err

set -euo pipefail

# --- Tooling (env-overridable solely for deterministic stub testing) ---------
FASTERQ_BIN="${FASTERQ_BIN:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/tools/sra-tools/bin/fasterq-dump}"
PIGZ_BIN="${PIGZ_BIN:-/usr/bin/pigz}"

# --- Locations (env-overridable for sandboxed testing; defaults are the
#     F26 production targets) ---------------------------------------------------
OUTDIR="${F26_OUTDIR:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/fastq}"
TMP_BASE="${F26_TMP_BASE:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/tmp_fasterq}"
CHECKSUM_TSV="${F26_CHECKSUM_TSV:-${OUTDIR}/f26_sha256_checksums.tsv}"

# --- Retry policy (network dropouts) ------------------------------------------
MAX_ATTEMPTS="${F26_MAX_ATTEMPTS:-3}"
RETRY_BACKOFF_S="${F26_RETRY_BACKOFF_S:-60}"   # base backoff; doubles per retry

# F26 cross-species datasets (see F26_DATASET_MANIFEST.tsv, role
# "F26 Cross-Species Evaluation"; historical baselines live in R1-*).
ACCESSIONS=(
    "CS-SA01:SRR33682404"   # Staphylococcus aureus MW2           (Gram+)
    "CS-BS01:SRR40739147"   # Bacillus subtilis str. 168          (Gram+)
    "CS-PA01:SRR41244734"   # Pseudomonas aeruginosa PAO1         (Gram-)
    "CS-SE01:SRR41010285"   # Salmonella enterica LT2             (Gram-)
)

log() { echo "[f26-download $(date -u '+%Y-%m-%dT%H:%M:%SZ')] $*"; }

mkdir -p "${OUTDIR}" "${TMP_BASE}"
touch "${CHECKSUM_TSV}"

# Count lines of a decompressed FASTQ. Returns the count via stdout.
# Fails (non-zero) if the gzip stream is corrupt, per `set -o pipefail`.
count_lines() {
    "${PIGZ_BIN}" -dc "$1" | wc -l
}

# Exact read count for a well-formed FASTQ: line count must divide by 4.
# Echoes the read count; returns non-zero on divisibility violation.
reads_of_file() {
    local file="$1" lines
    lines="$(count_lines "$1")"
    if (( lines % 4 != 0 )); then
        log "ERROR: ${file} has ${lines} lines, not divisible by 4"
        return 1
    fi
    echo "$(( lines / 4 ))"
}

# Resume verification: succeeds iff both mates exist, both gzip streams are
# intact, line counts divide by 4, and read counts are non-zero and equal.
pair_verified() {
    local r1="$1" r2="$2" n1 n2
    [[ -f "${r1}" && -f "${r2}" ]] || return 1
    "${PIGZ_BIN}" -t "${r1}" "${r2}" || return 1
    n1="$(reads_of_file "${r1}")" || return 1
    n2="$(reads_of_file "${r2}")" || return 1
    (( n1 > 0 && n1 == n2 )) || return 1
    return 0
}

# Record one compressed FASTQ in the persistent checksum ledger + stdout.
# The caller must have verified gzip integrity of ${file} beforehand.
record_checksum() {
    local acc="$1" file="$2" reads="$3"
    local sha size
    sha="$(sha256sum "${file}" | cut -d' ' -f1)"
    size="$(stat -c '%s' "${file}")"
    log "SHA256 ${file}: ${sha} (bytes=${size}, reads=${reads})"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${acc}" "$(basename "${file}")" "${size}" "${sha}" "${reads}" \
        "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "${CHECKSUM_TSV}"
}

# Download one accession with retries. Any non-retryable integrity failure
# aborts the whole run (data corruption must never be papered over).
process_accession() {
    local dsid="$1" acc="$2"
    local r1="${OUTDIR}/${acc}_1.fastq.gz"
    local r2="${OUTDIR}/${acc}_2.fastq.gz"
    local sing_uncompressed="${OUTDIR}/${acc}.fastq"
    local sing_gz="${OUTDIR}/${acc}.fastq.gz"
    local acc_tmp="${TMP_BASE}/${acc}_tmp"
    local attempt reads1 reads2 sing_lines backoff="${RETRY_BACKOFF_S}"

    log "--- ${dsid} (${acc}) ---"

    if pair_verified "${r1}" "${r2}"; then
        reads1="$(reads_of_file "${r1}")"
        log "${acc}: RESUME - already verified (R1=R2=${reads1} read pairs); skipping download"
        return 0
    fi
    log "${acc}: no verified pair on disk - downloading"

    # Stale artifacts from any previous interrupted attempt must not leak into
    # the new attempt.
    rm -rf "${acc_tmp}" "${sing_gz}"
    rm -f "${r1}" "${r2}" "${sing_uncompressed}" \
          "${OUTDIR}/${acc}_1.fastq" "${OUTDIR}/${acc}_2.fastq"

    for (( attempt = 1; attempt <= MAX_ATTEMPTS; attempt++ )); do
        log "${acc}: fasterq-dump attempt ${attempt}/${MAX_ATTEMPTS}"
        mkdir -p "${acc_tmp}"
        if "${FASTERQ_BIN}" --split-3 -e 8 -O "${OUTDIR}" -t "${acc_tmp}" "${acc}"; then
            break
        fi
        rm -rf "${acc_tmp}"
        rm -f "${OUTDIR}/${acc}_1.fastq" "${OUTDIR}/${acc}_2.fastq" "${sing_uncompressed}"
        if (( attempt == MAX_ATTEMPTS )); then
            log "ERROR: ${acc}: fasterq-dump failed after ${MAX_ATTEMPTS} attempts"
            return 1
        fi
        log "${acc}: attempt failed (likely network dropout); retrying in ${backoff}s"
        sleep "${backoff}"
        backoff=$(( backoff * 2 ))
    done
    rm -rf "${acc_tmp}"

    # Singleton file from --split-3: compress if non-empty, drop if empty.
    if [[ -f "${sing_uncompressed}" ]]; then
        if [[ -s "${sing_uncompressed}" ]]; then
            sing_lines="$(count_lines "${sing_uncompressed}")"
            log "${acc}: singleton file present (lines=${sing_lines}, reads=$(( sing_lines / 4 ))) - compressing"
            "${PIGZ_BIN}" -p 8 -f "${sing_uncompressed}"
        else
            log "${acc}: empty singleton file - removing"
            rm -f "${sing_uncompressed}"
        fi
    fi

    log "${acc}: compressing mates with pigz -p 8"
    "${PIGZ_BIN}" -p 8 -f "${OUTDIR}/${acc}_1.fastq"
    "${PIGZ_BIN}" -p 8 -f "${OUTDIR}/${acc}_2.fastq"

    # Gzip integrity of the freshly produced archives.
    "${PIGZ_BIN}" -t "${r1}" "${r2}"
    log "${acc}: gzip integrity OK"

    reads1="$(reads_of_file "${r1}")"
    reads2="$(reads_of_file "${r2}")"
    log "${acc}: exact read counts R1=${reads1} R2=${reads2}"
    if (( reads1 == 0 )); then
        log "ERROR: ${acc}: zero reads in R1"
        return 1
    fi
    if (( reads1 != reads2 )); then
        log "ERROR: ${acc}: read count mismatch R1=${reads1} != R2=${reads2}"
        return 1
    fi

    record_checksum "${acc}" "${r1}" "${reads1}"
    record_checksum "${acc}" "${r2}" "${reads2}"
    if [[ -f "${sing_gz}" ]]; then
        record_checksum "${acc}" "${sing_gz}" "$(( $(count_lines "${sing_gz}") / 4 ))"
    fi
    log "${acc}: OK"
}

log "=== F26 cross-species dataset download started (node=$(hostname), OUTDIR=${OUTDIR}) ==="
failures=()
for entry in "${ACCESSIONS[@]}"; do
    dsid="${entry%%:*}"
    acc="${entry##*:}"
    if ! process_accession "${dsid}" "${acc}"; then
        failures+=( "${acc}" )
        log "ERROR: ${acc}: FAILED; continuing with remaining accessions"
    fi
done

if (( ${#failures[@]} > 0 )); then
    log "=== F26 download FAILED for: ${failures[*]} ==="
    exit 1
fi
log "=== F26 cross-species dataset download complete; checksums in ${CHECKSUM_TSV} ==="

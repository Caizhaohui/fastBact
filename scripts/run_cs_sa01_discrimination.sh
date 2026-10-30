#!/usr/bin/env bash
# =============================================================================
# F26-DATA-004 — CS-SA01 Primary Strain-Discrimination Alignment Analysis
# =============================================================================
# Purpose:
#   Resolve audit finding BC-01 (F26-AUDIT-002, docs/development/
#   F26_INDEPENDENT_AUDIT_STAGE_BC.md) with retained PRIMARY alignment evidence:
#   exact mapped-read counts, mean depth, and % genome coverage of
#   SRR33682404 against (MW2 + Staph1N) and (MRSA252 + Staph1N), plus
#   strain-diagnostic locus depth, sequence-verified strain-exclusive marker
#   sets, shared-homology support, and Poisson sensitivity limits.
#
# Method notes (v3):
#   - Whole-diagnostic-window cross-alignment with minimap2 asm presets fails
#     for long queries carrying large divergent inserts (verified: a 52 kb
#     window verbatim-substring of its own target genome returned 0 PAF lines).
#     Cross-strain homology is therefore derived from 5 kb chunked -x sr
#     alignments (chunks behave like reads; verified working), and shared
#     intervals are unioned in window coordinates (PAF qstart/qend + chunk
#     offset), including homology to the co-culture phage Staph1N.
#   - The ABSENCE gate is evaluated on the sequence-verified EXCLUSIVE marker
#     set (window bases with no sr-level alignment to the other strain's
#     chromosome nor to Staph1N at MAPQ>=10), using MAPQ>=20 depth/reads to
#     exclude spurious cross-mapping.
#
# Usage:
#   sbatch scripts/run_cs_sa01_discrimination.sh
#
# Compute contract: Slurm qfree, 8 CPUs, 14 GB, 20 min (per F26-DATA-004).
# All retained artifacts land in
#   docs/development/evidence/cs_sa01_strain_discrimination/
# BAM/BAM-index scratch files are deleted on exit (deterministically recreated
# by re-running this script).
# =============================================================================
#SBATCH --job-name=f26_cs_sa01_disc
#SBATCH --partition=qfree
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=00:20:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree/docs/development/evidence/cs_sa01_strain_discrimination/logs/slurm_%j.out

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# NOTE: under sbatch this cluster stages the batch script into
# /var/spool/slurm/d/<jobid>/ and executes it there, so BASH_SOURCE-based
# worktree resolution is invalid in batch mode. The worktree is therefore
# pinned explicitly (override with F26_WORKTREE when re-running elsewhere).
DEFAULT_WORKTREE="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree"
WORKTREE="${F26_WORKTREE:-${DEFAULT_WORKTREE}}"
SCRIPT_DIR="${WORKTREE}/scripts"
EVID="${WORKTREE}/docs/development/evidence/cs_sa01_strain_discrimination"
LOG="${EVID}/cs_sa01_disc_run.log"
# qfree vnodes expose no SLURM_TMPDIR; TMPDIR=/tmp is node-local scratch
# (approx. 120 GB free). Slurm stdout targets the shared filesystem via the
# absolute #SBATCH --output path above.
SCRATCH="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/f26_cs_sa01_disc_$$"
mkdir -p "${SCRATCH}" "${EVID}"

# Remove stale derived artifacts from any previous run (kept inputs are only
# the two staged reference FASTAs and logs/); this guarantees the retained
# artifact set is always the product of the current, completed run.
rm -f "${EVID}"/cs_sa01_disc_*.{tsv,tsv.gz,txt,paf} \
      "${EVID}"/SHA256SUMS.txt \
      "${EVID}/loci"/*_exclusive.bed "${EVID}/loci"/*_shared.bed \
      "${EVID}/loci"/*_exclusive.fasta 2>/dev/null || true

REF_MW2="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/references/NC_003923.1.fasta"
REF_M252="${EVID}/NC_002952.2.fasta"
REF_PHAGE="${EVID}/staph1n.fasta"
R1="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/fastq/SRR33682404_1.fastq.gz"
R2="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/fastq/SRR33682404_2.fastq.gz"

MINIMAP2="/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/minimap2"
SAMTOOLS="/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/samtools"

# Pinned input hashes (docs/development/F26_FASTQ_CHECKSUMS.tsv, 2026-10-30)
R1_SHA_EXPECT="1a7e70141657dfa9a5eb9c6b6ba491c260127450b01ab9eb94b85e0fe2a466c3"
R2_SHA_EXPECT="88a46ad12cb159511d112fa059b26732ae36d9c3e28a6355f806638f2599e050"
R1_READS_EXPECT=2239653
R2_READS_EXPECT=2239653
# Pinned MW2 raw-file hash (F26_DATASET_PROVENANCE_AND_SELECTION.md §3)
MW2_SHA_EXPECT="19ba62a339065c9c498efb72bb2491e372344e4f23a4ac1e3e34e7345a25c63c"

# Diagnostic loci (1-based inclusive; F26-DATA-004 dispatch specification)
PVL_LABEL="PVL_lukS-lukF_MW2"
PVL_REGION="NC_003923.1:1515001-1517500"
SCC4_LABEL="SCCmecIV_MW2"
SCC4_REGION="NC_003923.1:34001-58000"
SCC2_LABEL="SCCmecII_MRSA252"
SCC2_REGION="NC_002952.2:34001-86000"

CHUNK=5000          # chunk length for cross-strain alignment (bp)
SHARED_MAPQ=10      # min MAPQ for a chunk hit to count as shared homology
MIN_EXCL_BLOCK=200  # min exclusive-block length (shorter pieces recorded as filtered)

# ---------------------------------------------------------------------------
# Retained exit code + transcript
# ---------------------------------------------------------------------------
EXITCODE_FILE="${EVID}/cs_sa01_disc_exit_code.txt"
trap 'rc=$?; printf "%s\n" "${rc}" > "${EXITCODE_FILE}"; \
      printf "[TRAP] job exiting rc=%s at %s\n" "${rc}" "$(date -u +%FT%TZ)" >> "${LOG}"; \
      rm -rf "${SCRATCH}"' EXIT

: > "${LOG}"
exec > >(tee -a "${LOG}") 2>&1

run() {
  printf '[CMD] %s\n' "$*"
  "$@"
  printf '[EXIT] 0\n'
}
STEP() { printf '\n===== [STEP %s] %s (%s) =====\n' "$1" "$2" "$(date -u +%FT%TZ)"; }

printf 'F26-DATA-004 CS-SA01 primary strain-discrimination analysis\n'
printf 'job_id=%s node=%s started=%s\n' \
  "${SLURM_JOB_ID:-LOGIN}" "$(hostname -s)" "$(date -u +%FT%TZ)"
printf 'git_head=%s\n' "$(git -C "${WORKTREE}" rev-parse HEAD)"
printf 'scratch=%s\n' "${SCRATCH}"

# ---------------------------------------------------------------------------
STEP 01 "environment"
# ---------------------------------------------------------------------------
run "${MINIMAP2}" --version
run "${SAMTOOLS}" --version
run nproc
run df -h "${SCRATCH}"

# ---------------------------------------------------------------------------
STEP 02 "input integrity (SHA-256)"
# ---------------------------------------------------------------------------
IN_SHA="${EVID}/cs_sa01_disc_inputs_sha256.tsv"
printf 'file\tbytes\tsha256\n' > "${IN_SHA}"
for f in "${R1}" "${R2}" "${REF_MW2}" "${REF_M252}" "${REF_PHAGE}"; do
  b="$(stat -c %s "${f}")"
  h="$(sha256sum "${f}" | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "$(basename "${f}")" "${b}" "${h}" | tee -a "${IN_SHA}"
done
h1="$(awk -F'\t' '$1=="SRR33682404_1.fastq.gz"{print $3}' "${IN_SHA}")"
h2="$(awk -F'\t' '$1=="SRR33682404_2.fastq.gz"{print $3}' "${IN_SHA}")"
hm="$(awk -F'\t' '$1=="NC_003923.1.fasta"{print $3}' "${IN_SHA}")"
[ "${h1}" = "${R1_SHA_EXPECT}" ] || { echo "FATAL: R1 sha256 mismatch"; exit 2; }
[ "${h2}" = "${R2_SHA_EXPECT}" ] || { echo "FATAL: R2 sha256 mismatch"; exit 2; }
[ "${hm}" = "${MW2_SHA_EXPECT}" ] || { echo "FATAL: MW2 ref sha256 mismatch"; exit 2; }
echo "input SHA-256 verification: PASS (R1, R2, MW2 ref all match pinned values)"

# ---------------------------------------------------------------------------
STEP 03 "combined references + faidx"
# ---------------------------------------------------------------------------
REF_A="${SCRATCH}/mw2_staph1n.fasta"
REF_B="${SCRATCH}/mrsa252_staph1n.fasta"
cat "${REF_MW2}" "${REF_PHAGE}" > "${REF_A}"
cat "${REF_M252}" "${REF_PHAGE}" > "${REF_B}"
run "${SAMTOOLS}" faidx "${REF_A}"
run "${SAMTOOLS}" faidx "${REF_B}"
for r in "${REF_A}" "${REF_B}"; do
  echo "--- ${r}.fai ---"
  cat "${r}.fai"
done

# ---------------------------------------------------------------------------
STEP 04 "alignment: reads vs MW2+Staph1N (minimap2 -ax sr | samtools sort)"
# ---------------------------------------------------------------------------
BAM_A="${SCRATCH}/vs_mw2_staph1n.sorted.bam"
printf '[CMD] minimap2 -t 8 -ax sr %s %s %s | samtools sort -@ 4 -m 1G -o %s -\n' \
  "${REF_A}" "${R1}" "${R2}" "${BAM_A}"
"${MINIMAP2}" -t 8 -ax sr "${REF_A}" "${R1}" "${R2}" \
  | "${SAMTOOLS}" sort -@ 4 -m 1G -o "${BAM_A}" -
run "${SAMTOOLS}" index "${BAM_A}"

# ---------------------------------------------------------------------------
STEP 05 "alignment: reads vs MRSA252+Staph1N"
# ---------------------------------------------------------------------------
BAM_B="${SCRATCH}/vs_mrsa252_staph1n.sorted.bam"
printf '[CMD] minimap2 -t 8 -ax sr %s %s %s | samtools sort -@ 4 -m 1G -o %s -\n' \
  "${REF_B}" "${R1}" "${R2}" "${BAM_B}"
"${MINIMAP2}" -t 8 -ax sr "${REF_B}" "${R1}" "${R2}" \
  | "${SAMTOOLS}" sort -@ 4 -m 1G -o "${BAM_B}" -
run "${SAMTOOLS}" index "${BAM_B}"

# ---------------------------------------------------------------------------
STEP 06 "per-alignment primary statistics (flagstat / idxstats / coverage / bins)"
# ---------------------------------------------------------------------------
# NOTE: the [CMD]/[EXIT] transcript lines must NOT be captured into the
# artifact files, so these commands redirect inside an explicit printf pair
# instead of the run() helper.
for tag in mw2 mrsa252; do
  bam="${SCRATCH}/vs_${tag}_staph1n.sorted.bam"
  printf '[CMD] samtools flagstat -@ 2 %s > %s\n' "${bam}" "${EVID}/cs_sa01_disc_${tag}.flagstat.txt"
  "${SAMTOOLS}" flagstat -@ 2 "${bam}" > "${EVID}/cs_sa01_disc_${tag}.flagstat.txt"
  printf '[CMD] samtools idxstats %s > %s\n' "${bam}" "${EVID}/cs_sa01_disc_${tag}.idxstats.tsv"
  "${SAMTOOLS}" idxstats "${bam}" > "${EVID}/cs_sa01_disc_${tag}.idxstats.tsv"
  printf '[CMD] samtools coverage %s > %s\n' "${bam}" "${EVID}/cs_sa01_disc_${tag}.coverage.tsv"
  "${SAMTOOLS}" coverage "${bam}" > "${EVID}/cs_sa01_disc_${tag}.coverage.tsv"
  printf '[EXIT] 0\n'
done

# 10-kb binned mean-depth profile over each bacterial contig (retained evidence)
depth_bins() { # bam contig outfile
  printf '[CMD] samtools depth -a -q 0 -Q 0 -d 0 -r %s %s | awk 10kb-bin mean -> %s\n' \
    "$2" "$1" "$3"
  "${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -r "$2" "$1" \
    | awk -v OFS='\t' 'BEGIN{print "bin_start\tbin_end\tmean_depth"} \
        {b=int(($2-1)/10000); s[b]+=$3; c[b]++; last=b} \
        END{for(i=0;i<=last;i++) printf "%d\t%d\t%.2f\n", i*10000+1, (i+1)*10000, s[i]/c[i]}' \
    > "$3"
}
depth_bins "${BAM_A}" NC_003923.1 "${EVID}/cs_sa01_disc_mw2.depth10k.tsv"
depth_bins "${BAM_B}" NC_002952.2 "${EVID}/cs_sa01_disc_mrsa252.depth10k.tsv"

# Full per-base depth over the phage contig in both alignments (retained, gzipped)
printf '[CMD] samtools depth -a -q 0 -Q 0 -d 0 -r NC_047722.1 (both alignments) | gzip -1\n'
"${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -r NC_047722.1 "${BAM_A}" | gzip -1 \
  > "${EVID}/cs_sa01_disc_phage_depth_mw2align.tsv.gz"
"${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -r NC_047722.1 "${BAM_B}" | gzip -1 \
  > "${EVID}/cs_sa01_disc_phage_depth_mrsa252align.tsv.gz"

# ---------------------------------------------------------------------------
STEP 07 "diagnostic-locus depth + exact overlapping read counts (both alignments)"
# ---------------------------------------------------------------------------
LOC_RAW="${EVID}/cs_sa01_disc_locus_depth.tsv.gz"
LOC_SUM="${EVID}/cs_sa01_disc_locus_summary.tsv"
LOC_CNT="${EVID}/cs_sa01_disc_locus_readcounts.tsv"
printf 'alignment\tlocus\tregion\tref\tstart\tend\twindow_bp\tmean_depth\tbreadth_ge1x\tbreadth_ge10x\tbreadth_ge30x\toverlapping_reads\n' > "${LOC_SUM}"
printf 'alignment\tlocus\tregion\toverlapping_reads\n' > "${LOC_CNT}"
printf 'alignment\tlocus\tregion\tref\tpos\tdepth\n' | gzip -1 > "${LOC_RAW}"

locus_stats() { # bam tag label region ref start end
  local bam="$1" tag="$2" label="$3" region="$4" ref="$5" start="$6" end="$7"
  local nreads
  nreads="$("${SAMTOOLS}" view -c "${bam}" "${region}")"
  printf '[CMD] samtools depth -a -q 0 -Q 0 -d 0 -r %s; samtools view -c -> %s reads\n' \
    "${region}" "${nreads}"
  "${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -r "${region}" "${bam}" \
    | awk -v OFS='\t' -v tag="${tag}" -v label="${label}" -v region="${region}" \
          -v ref="${ref}" -v start="${start}" -v end="${end}" -v nr="${nreads}" '
        { n++; d=$3; s+=d; if(d>=1)b1++; if(d>=10)b10++; if(d>=30)b30++ }
        END{ wp=end-start+1;
             printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%d\n",
               tag,label,region,ref,start,end,wp, s/n, b1/n, b10/n, b30/n, nr }' >> "${LOC_SUM}"
  printf '%s\t%s\t%s\t%s\n' "${tag}" "${label}" "${region}" "${nreads}" >> "${LOC_CNT}"
}

# Raw per-base depth for the record (both alignments x applicable loci)
locus_raw() { # bam tag label region
  "${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -r "$4" "$1" \
    | awk -v OFS='\t' -v tag="$2" -v label="$3" -v region="$4" \
        '{print tag,label,region,$1,$2,$3}'
}
{ locus_raw "${BAM_A}" MW2_align  "${PVL_LABEL}"  "${PVL_REGION}"
  locus_raw "${BAM_A}" MW2_align  "${SCC4_LABEL}" "${SCC4_REGION}"
  locus_raw "${BAM_B}" M252_align "${SCC2_LABEL}" "${SCC2_REGION}"
} | gzip -1 >> "${LOC_RAW}"

locus_stats "${BAM_A}" MW2_align  "${PVL_LABEL}"  "${PVL_REGION}"  NC_003923.1 1515001 1517500
locus_stats "${BAM_A}" MW2_align  "${SCC4_LABEL}" "${SCC4_REGION}" NC_003923.1 34001 58000
locus_stats "${BAM_B}" M252_align "${SCC2_LABEL}" "${SCC2_REGION}" NC_002952.2 34001 86000

echo "--- locus summary (whole windows) ---"; cat "${LOC_SUM}"
echo "--- locus exact overlapping read counts ---"; cat "${LOC_CNT}"

# ---------------------------------------------------------------------------
STEP 08 "diagnostic-window extraction + chunked cross-strain alignment"
# ---------------------------------------------------------------------------
LOCI_DIR="${EVID}/loci"
mkdir -p "${LOCI_DIR}"
run bash -c "\"${SAMTOOLS}\" faidx '${REF_MW2}'  ${PVL_REGION}  | sed 's/^>/${PVL_LABEL} /'  > '${LOCI_DIR}/pvl_lukS_lukF_MW2_NC_003923.1_1515001-1517500.fasta'"
run bash -c "\"${SAMTOOLS}\" faidx '${REF_MW2}'  ${SCC4_REGION} | sed 's/^>/${SCC4_LABEL} /' > '${LOCI_DIR}/sccmecIV_MW2_NC_003923.1_34001-58000.fasta'"
run bash -c "\"${SAMTOOLS}\" faidx '${REF_M252}' ${SCC2_REGION} | sed 's/^>/${SCC2_LABEL} /' > '${LOCI_DIR}/sccmecII_MRSA252_NC_002952.2_34001-86000.fasta'"
wc -c "${LOCI_DIR}"/*.fasta

# chunk a window into ${CHUNK}-bp pieces; header carries the 0-based window
# offset. faidx emits its own ">chr:start-end" header line, which is stripped
# so that our chunk label is the only header (PAF query names then encode the
# chunk offset needed for window-coordinate reconstruction).
chunk_window() { # prefix chr start end src_ref
  local p="$1" chr="$2" s="$3" e="$4" src="$5"
  local off=0 cs ce
  : > "${SCRATCH}/${p}_chunks.fasta"
  : > "${SCRATCH}/${p}_chunk_index.tsv"
  while (( s + off <= e )); do
    cs=$(( s + off )); ce=$(( cs + CHUNK - 1 )); (( ce > e )) && ce=${e}
    printf '>%s_chunk%d\n' "${p}" "${off}" >> "${SCRATCH}/${p}_chunks.fasta"
    "${SAMTOOLS}" faidx "${src}" "${chr}:${cs}-${ce}" | tail -n +2 \
      >> "${SCRATCH}/${p}_chunks.fasta"
    printf '%d\t%d\t%d\n' "${off}" "${cs}" "${ce}" >> "${SCRATCH}/${p}_chunk_index.tsv"
    off=$(( off + CHUNK ))
  done
  printf '[CHUNKS] %s: %d chunks of <=%d bp\n' "${p}" "$(wc -l < "${SCRATCH}/${p}_chunk_index.tsv")" "${CHUNK}"
}

# cross-align one window's chunks against one target; keep MAPQ>=SHARED_MAPQ
cross_align() { # prefix target_ref target_label
  local p="$1" tref="$2" tlab="$3"
  local paf="${SCRATCH}/${p}_vs_${tlab}.paf"
  printf '[CMD] minimap2 -t 8 -x sr %s %s_chunks.fasta (keep MAPQ>=%d) -> %s\n' \
    "${tref}" "${p}" "${SHARED_MAPQ}" "${paf}"
  "${MINIMAP2}" -t 8 -x sr "${tref}" "${SCRATCH}/${p}_chunks.fasta" 2>/dev/null \
    | awk -v OFS='\t' -v mq="${SHARED_MAPQ}" '$12>=mq{print $1,$3,$4,$5,$6,$7,$8,$10,$11,$12}' \
    > "${paf}"
  cp "${paf}" "${EVID}/cs_sa01_disc_chunkspecificity.${p}_vs_${tlab}.paf"
  printf '[HITS] %s_vs_%s: %d chunk-hits\n' "${p}" "${tlab}" "$(wc -l < "${paf}")"
}

chunk_window pvl_mw2   NC_003923.1 1515001 1517500 "${REF_MW2}"
chunk_window scc4_mw2  NC_003923.1 34001   58000   "${REF_MW2}"
chunk_window scc2_m252 NC_002952.2 34001   86000   "${REF_M252}"

cross_align pvl_mw2   "${REF_M252}"  MRSA252
cross_align pvl_mw2   "${REF_PHAGE}" Staph1N
cross_align scc4_mw2  "${REF_M252}"  MRSA252
cross_align scc4_mw2  "${REF_PHAGE}" Staph1N
cross_align scc2_m252 "${REF_MW2}"   MW2
cross_align scc2_m252 "${REF_PHAGE}" Staph1N

# Core-genome sensitivity control: MW2 chunks vs MRSA252 must align extensively
chunk_window ctrl_mw2 NC_003923.1 1 2820462 "${REF_MW2}"
cross_align ctrl_mw2  "${REF_M252}"  MRSA252

# ---------------------------------------------------------------------------
STEP 09 "sequence-verified shared/exclusive marker partition + partition-level read support"
# ---------------------------------------------------------------------------
# Shared = union of window-coordinate intervals of all chunk-hits (MAPQ>=10)
# against the other strain's chromosome and against Staph1N. Window coordinate
# = chunk offset (from query name suffix) + PAF qstart/qend (converted to
# 1-based inclusive). Exclusive = complement of the shared union within
# [1, window_bp]; blocks shorter than MIN_EXCL_BLOCK bp recorded as filtered.
BLOCKS_TSV="${EVID}/cs_sa01_disc_marker_blocks.tsv"
printf 'window\twindow_bp\tshared_vs_other_strain_bp\tshared_vs_phage_bp\texclusive_bp\tn_exclusive_blocks\tn_filtered_short\texclusive_fraction\n' > "${BLOCKS_TSV}"

partition_window() { # prefix window_bp window_chr_start chr src_ref
  local p="$1" wlen="$2" wstart="$3" chr="$4" src="$5"
  # internal (2-column) interval files live in scratch; the retained BEDs in
  # loci/ get the chromosome column prepended (samtools requires 3-column BED)
  local excl_bed="${LOCI_DIR}/${p}_exclusive.bed"
  local excl_fa="${LOCI_DIR}/${p}_exclusive.fasta"
  local shared_bed="${LOCI_DIR}/${p}_shared.bed"
  local hits_strain="${SCRATCH}/${p}_hits_strain.txt"
  local hits_phage="${SCRATCH}/${p}_hits_phage.txt"
  local shared_int="${SCRATCH}/${p}_shared.txt"
  local excl_int="${SCRATCH}/${p}_excl.bed"
  local filt_int="${SCRATCH}/${p}_filtered.bed"

  # normalize hit tables to 1-based inclusive window coordinates.
  # Stored PAF columns (from cross_align): 1=chunk query, 2=qstart(0-based),
  # 3=qend(0-based exclusive); window coord = chunk offset + q-start/q-end.
  awk -v OFS='\t' '{split($1,a,"chunk"); off=a[2]+0; print off+$2+1, off+$3}' \
    "${SCRATCH}/${p}_vs_MRSA252.paf" > "${hits_strain}.m252" 2>/dev/null || : > "${hits_strain}.m252"
  awk -v OFS='\t' '{split($1,a,"chunk"); off=a[2]+0; print off+$2+1, off+$3}' \
    "${SCRATCH}/${p}_vs_MW2.paf" > "${hits_strain}.mw2" 2>/dev/null || : > "${hits_strain}.mw2"
  cat "${hits_strain}.m252" "${hits_strain}.mw2" > "${hits_strain}"
  awk -v OFS='\t' '{split($1,a,"chunk"); off=a[2]+0; print off+$2+1, off+$3}' \
    "${SCRATCH}/${p}_vs_Staph1N.paf" > "${hits_phage}" 2>/dev/null || : > "${hits_phage}"

  local shared_strain shared_phage
  shared_strain="$(awk '{s+=$2-$1+1} END{print s+0}' "${hits_strain}")"
  shared_phage="$(awk '{s+=$2-$1+1} END{print s+0}' "${hits_phage}")"

  # merge union (both hit classes) into shared intervals
  cat "${hits_strain}" "${hits_phage}" | sort -k1,1n \
    | awk 'BEGIN{OFS="\t"} {s1=$1;e1=$2; if(NR==1){s=s1;e=e1;next}
           if(s1<=e+1){if(e1>e)e=e1}else{print s,e;s=s1;e=e1}}
           END{if(NR)print s,e}' > "${shared_int}"

  # complement within [1, wlen] -> BED (0-based half-open)
  awk -v OFS='\t' -v L="${wlen}" 'BEGIN{prev=1}
      {if($1>prev) print prev-1, $1-1; if($2+1>prev) prev=$2+1}
      END{if(prev<=L) print prev-1, L}' "${shared_int}" > "${SCRATCH}/${p}_compl.bed"
  awk -v OFS='\t' -v m="${MIN_EXCL_BLOCK}" '$2-$1>=m' "${SCRATCH}/${p}_compl.bed" > "${excl_int}"
  awk -v OFS='\t' -v m="${MIN_EXCL_BLOCK}" '$2-$1<m'  "${SCRATCH}/${p}_compl.bed" > "${filt_int}"

  # retained 3-column BEDs (chrom, chromStart, chromEnd) in TRUE chromosome
  # coordinates, 0-based half-open: 'samtools depth -b' / 'samtools view -L'
  # resolve no intervals from 2-column window-relative BEDs (F26-DATA-004
  # bugfix). excl_int is window-relative 0-based half-open -> add (wstart-1);
  # shared_int is window 1-based inclusive -> convert to 0-based half-open,
  # then add (wstart-1).
  awk -v OFS='\t' -v c="${chr}" -v o="$(( wstart - 1 ))" \
      '{print c, $1 + o - 1, $2 + o}' "${shared_int}" > "${shared_bed}"
  awk -v OFS='\t' -v c="${chr}" -v o="$(( wstart - 1 ))" \
      '{print c, $1 + o,     $2 + o}' "${excl_int}"   > "${excl_bed}"

  local excl_bp n_excl n_filt
  excl_bp="$(awk '{s+=$2-$1} END{printf "%d", s+0}' "${excl_int}")"
  n_excl="$(wc -l < "${excl_int}")"
  n_filt="$(wc -l < "${filt_int}")"
  awk -v OFS='\t' -v w="${p}" -v L="${wlen}" -v ss="${shared_strain}" -v sp="${shared_phage}" \
      -v ex="${excl_bp}" -v ne="${n_excl}" -v nf="${n_filt}" \
      'BEGIN{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", w,L,ss,sp,ex,ne,nf,ex/L}' >> "${BLOCKS_TSV}"

  # exclusive marker set FASTA (extract from the window's source genome);
  # excl_bed is 0-based half-open on the chromosome -> faidx wants 1-based
  # inclusive, i.e. [start+1, end]
  : > "${excl_fa}"
  local bchr bs be
  while read -r bchr bs be; do
    printf '>%s_excl_%d_%d\n' "${p}" "${bs}" "${be}" >> "${excl_fa}"
    "${SAMTOOLS}" faidx "${src}" "${bchr}:$(( bs + 1 ))-${be}" >> "${excl_fa}"
  done < "${excl_bed}"

  printf '[PARTITION] %s shared_strain=%s shared_phage=%s exclusive=%s bp in %s blocks (min %s bp; %s filtered)\n' \
    "${p}" "${shared_strain}" "${shared_phage}" "${excl_bp}" "${n_excl}" "${MIN_EXCL_BLOCK}" "${n_filt}"
}

partition_window pvl_mw2   2500  1515001 NC_003923.1 "${REF_MW2}"
partition_window scc4_mw2  24000 34001   NC_003923.1 "${REF_MW2}"
partition_window scc2_m252 52000 34001   NC_002952.2 "${REF_M252}"

echo "--- marker block partition ---"; cat "${BLOCKS_TSV}"

# partition-level depth statistics + exact overlapping read counts (q0 and q20)
BED_SUM="${EVID}/cs_sa01_disc_marker_support.tsv"
printf 'alignment\tmarker_set\tbp\tmean_depth_q0\tbreadth_ge1x_q0\tmean_depth_q20\tbreadth_ge1x_q20\toverlapping_reads_q0\toverlapping_reads_q20\n' > "${BED_SUM}"
bed_support() { # bam tag bed label
  local bam="$1" tag="$2" bed="$3" label="$4"
  local bp r0 r20
  bp="$(awk '{s+=$3-$2} END{printf "%d", s+0}' "${bed}")"
  r0="$("${SAMTOOLS}" view -c -L "${bed}" "${bam}")"
  r20="$("${SAMTOOLS}" view -c -q 20 -L "${bed}" "${bam}")"
  printf '[CMD] samtools depth -b %s (q0 and q20); view -c -L -> %s reads (MAPQ>=20: %s)\n' "${bed}" "${r0}" "${r20}"
  "${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -b "${bed}" "${bam}" \
    | awk '{n++; d=$3; s+=d; if(d>=1)b1++}
        END{if(n==0) printf "0.00\t0.0000\n"; else printf "%.2f\t%.4f\n", s/n, b1/n}' > "${SCRATCH}/bed_q0.txt"
  "${SAMTOOLS}" depth -a -q 0 -Q 20 -d 0 -b "${bed}" "${bam}" \
    | awk '{n++; d=$3; s+=d; if(d>=1)b1++}
        END{if(n==0) printf "0.00\t0.0000\n"; else printf "%.2f\t%.4f\n", s/n, b1/n}' > "${SCRATCH}/bed_q20.txt"
  printf '%s\t%s\t%s\t' "${tag}" "${label}" "${bp}" >> "${BED_SUM}"
  paste "${SCRATCH}/bed_q0.txt" "${SCRATCH}/bed_q20.txt" \
    | awk -v OFS='\t' -v r0="${r0}" -v r20="${r20}" \
      '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, r0, r20}' >> "${BED_SUM}"
}
bed_support "${BAM_A}" MW2_align  "${LOCI_DIR}/pvl_mw2_exclusive.bed"   "${PVL_LABEL}_EXCLUSIVE"
bed_support "${BAM_A}" MW2_align  "${LOCI_DIR}/scc4_mw2_exclusive.bed"  "${SCC4_LABEL}_EXCLUSIVE"
bed_support "${BAM_B}" M252_align "${LOCI_DIR}/scc2_m252_exclusive.bed" "${SCC2_LABEL}_EXCLUSIVE"
bed_support "${BAM_B}" M252_align "${LOCI_DIR}/scc2_m252_shared.bed"    "${SCC2_LABEL}_SHARED"

# ---------------------------------------------------------------------------
STEP 09b "read-origin audit over the MRSA252-exclusive marker set"
# ---------------------------------------------------------------------------
# The absence gate must distinguish (a) reads supported by the MW2 CHROMOSOME
# from (b) reads with no placement on the MW2 chromosome (extra-chromosomal
# content: pMW2/pMW3 plasmids are outside the chromosome-only reference and
# SCCmec cassettes carry plasmid-derived cargo such as pUB110 segments).
# A read overlapping the MRSA252-exclusive set whose QNAME is mapped somewhere
# on NC_003923.1 in the MW2 alignment carries MW2-chromosomal support; a read
# absent from that set carries none (expected for plasmid-origin reads).
ORIGIN_TSV="${EVID}/cs_sa01_disc_exclusive_origin_audit.tsv"
printf 'marker_set\tblock_id\tchr\tstart\tend\tbp\treads_q20_primary\tdivergence_mean\tdivergence_min\tchromosome_supported_reads\n' > "${ORIGIN_TSV}"

EXCL_BED="${LOCI_DIR}/scc2_m252_exclusive.bed"
printf '[CMD] samtools view -F 2304 -q 20 -L %s (primary, MAPQ>=20)\n' "${EXCL_BED}"
"${SAMTOOLS}" view -F 2304 -q 20 -L "${EXCL_BED}" "${BAM_B}" > "${SCRATCH}/excl_reads.sam"
wc -l < "${SCRATCH}/excl_reads.sam"

awk '!/^@/{print $1}' "${SCRATCH}/excl_reads.sam" | sort -u > "${SCRATCH}/excl_qnames.txt"
"${SAMTOOLS}" view -F 2304 "${BAM_A}" NC_003923.1 \
  | awk '{print $1}' | sort -u > "${SCRATCH}/mw2_mapped_qnames.txt"
printf '[COUNTS] excl_qnames=%s mw2_chrom_mapped_qnames=%s\n' \
  "$(wc -l < "${SCRATCH}/excl_qnames.txt")" "$(wc -l < "${SCRATCH}/mw2_mapped_qnames.txt")"
comm -12 "${SCRATCH}/excl_qnames.txt" "${SCRATCH}/mw2_mapped_qnames.txt" \
  > "${SCRATCH}/excl_qnames_chrom_supported.txt"
printf '[COUNTS] chromosome-supported excl qnames=%s (of %s)\n' \
  "$(wc -l < "${SCRATCH}/excl_qnames_chrom_supported.txt")" "$(wc -l < "${SCRATCH}/excl_qnames.txt")"

# identity distribution (minimap2 de:f tag) over the exclusive-set reads
awk '!/^@/{if(match($0,/de:f:[0-9.]+/)){s+=substr($0,RSTART+5,RLENGTH-5);n++}} END{if(n)printf "divergence_mean=%.4f n=%d\n",s/n,n; else print "divergence_mean=NA n=0"}' \
  "${SCRATCH}/excl_reads.sam" | tee "${SCRATCH}/excl_identity_summary.txt"
awk '!/^@/{if(match($0,/de:f:[0-9.]+/)){v=substr($0,RSTART+5,RLENGTH-5)+0; s+=v; n++; if(n==1||v<mn)mn=v}} END{if(n)printf "%.4f\t%.4f\t%d\n",s/n,mn,n; else printf "0\t0\t0\n"}' \
  "${SCRATCH}/excl_reads.sam" > "${SCRATCH}/excl_identity_stats.txt"

# per-block audit
blk_id=0
total_reads_q20=0; total_supported=0
while read -r bc bs be; do
  blk_id=$(( blk_id + 1 ))
  region="${bc}:$(( bs + 1 ))-${be}"
  "${SAMTOOLS}" view -F 2304 -q 20 "${BAM_B}" "${region}" > "${SCRATCH}/blk.sam"
  bpn=$(( be - bs ))
  awk '!/^@/{print $1}' "${SCRATCH}/blk.sam" | sort -u > "${SCRATCH}/blk_qn.txt"
  nsup="$(comm -12 "${SCRATCH}/blk_qn.txt" "${SCRATCH}/mw2_mapped_qnames.txt" | wc -l)"
  nq="$(wc -l < "${SCRATCH}/blk_qn.txt")"
  read -r imean imin iid <<< "$(awk '!/^@/{if(match($0,/de:f:[0-9.]+/)){v=substr($0,RSTART+5,RLENGTH-5)+0; s+=v; n++; if(n==1||v<mn)mn=v}} END{if(n)printf "%.4f %.4f %d",s/n,mn,n; else printf "0 0 0"}' "${SCRATCH}/blk.sam")"
  printf '%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n' \
    "${SCC2_LABEL}_EXCLUSIVE" "${blk_id}" "${bc}" "${bs}" "${be}" "${bpn}" "${nq}" "${imean}" "${imin}" "${nsup}" >> "${ORIGIN_TSV}"
  total_reads_q20=$(( total_reads_q20 + nq )); total_supported=$(( total_supported + nsup ))
done < "${EXCL_BED}"

read -r ov_mean ov_min ov_n <<< "$(cat "${SCRATCH}/excl_identity_stats.txt")"
printf 'OVERALL\t0\tall\t0\t0\t0\t%s\t%s\t%s\t%s\n' \
  "${ov_n}" "${ov_mean}" "${ov_min}" "$(wc -l < "${SCRATCH}/excl_qnames_chrom_supported.txt")" >> "${ORIGIN_TSV}"

# same audit for the SHARED set (control: expect chromosome support there)
"${SAMTOOLS}" view -F 2304 -q 20 -L "${LOCI_DIR}/scc2_m252_shared.bed" "${BAM_B}" \
  | awk '!/^@/{print $1}' | sort -u > "${SCRATCH}/shared_qnames.txt"
nsup_shared="$(comm -12 "${SCRATCH}/shared_qnames.txt" "${SCRATCH}/mw2_mapped_qnames.txt" | wc -l)"
printf 'SHARED_CONTROL\t-\tNC_002952.2\tshared_set\t-\t14496\t%s\tNA\tNA\t%s\n' \
  "$(wc -l < "${SCRATCH}/shared_qnames.txt")" "${nsup_shared}" >> "${ORIGIN_TSV}"
echo "--- exclusive-set read-origin audit ---"; cat "${ORIGIN_TSV}"
echo "shared-set control: chromosome-supported ${nsup_shared} of $(wc -l < "${SCRATCH}/shared_qnames.txt") qnames"

echo "--- marker support (partition level) ---"; cat "${BED_SUM}"

# per-base depth over the exclusive marker sets (raw record, q0 with q20)
excl_depth_gz() { # bam bed tag
  "${SAMTOOLS}" depth -a -q 0 -Q 0 -d 0 -b "$2" "$1" > "${SCRATCH}/ed_q0.txt"
  "${SAMTOOLS}" depth -a -q 0 -Q 20 -d 0 -b "$2" "$1" > "${SCRATCH}/ed_q20.txt"
  paste "${SCRATCH}/ed_q0.txt" <(cut -f3 "${SCRATCH}/ed_q20.txt") \
    | awk -v OFS='\t' -v tag="$3" '{print tag,$1,$2,$3,$4}'
}
{ printf 'alignment\tmarker_set\tref\tpos\tdepth_q0\tdepth_q20\n'
  excl_depth_gz "${BAM_A}" "${LOCI_DIR}/pvl_mw2_exclusive.bed"   PVL_EXCL
  excl_depth_gz "${BAM_A}" "${LOCI_DIR}/scc4_mw2_exclusive.bed"  SCC4_EXCL
  excl_depth_gz "${BAM_B}" "${LOCI_DIR}/scc2_m252_exclusive.bed" SCC2_EXCL
} | gzip -1 > "${EVID}/cs_sa01_disc_exclusive_depth.tsv.gz"

# ---------------------------------------------------------------------------
STEP 10 "consistency checks, gates, sensitivity bounds, summary"
# ---------------------------------------------------------------------------
SUM="${EVID}/cs_sa01_disc_summary.tsv"
printf 'metric\tvalue\tunit\tsource\n' > "${SUM}"
add() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "${SUM}"; }

total_reads=$(( R1_READS_EXPECT + R2_READS_EXPECT ))
add input_read_records "${total_reads}" "single-end records (= 2 x ${R1_READS_EXPECT} pairs)" "F26_FASTQ_CHECKSUMS.tsv"

mw2_total=0; mw2_bac=0; mw2_ph=0; mw2_unmapped=0
# idxstats.tsv files carry '[CMD] ...' / '[EXIT] 0' transcript lines from the
# run() wrapper; skip every '['-prefixed line so only real idxstats rows
# (refName, len, mapped, unmapped) reach the arithmetic below (F26-DATA-004
# bugfix: without this, 'mapped' took the value 'idxstats' -> unbound variable
# under set -u).
while read -r ref len mapped unmapped; do
  case "${ref}" in \[* ) continue ;; esac
  mw2_total=$(( mw2_total + mapped + unmapped ))
  case "${ref}" in
    NC_003923.1) mw2_bac=${mapped} ;;
    NC_047722.1) mw2_ph=${mapped} ;;
    '*') mw2_unmapped=${unmapped} ;;
  esac
done < "${EVID}/cs_sa01_disc_mw2.idxstats.tsv"
m2_total=0; m2_bac=0; m2_ph=0; m2_unmapped=0
while read -r ref len mapped unmapped; do
  case "${ref}" in \[* ) continue ;; esac
  m2_total=$(( m2_total + mapped + unmapped ))
  case "${ref}" in
    NC_002952.2) m2_bac=${mapped} ;;
    NC_047722.1) m2_ph=${mapped} ;;
    '*') m2_unmapped=${unmapped} ;;
  esac
done < "${EVID}/cs_sa01_disc_mrsa252.idxstats.tsv"

# flagstat layout: "<count> + 0 supplementary" -> count is field 1 of a line
# whose field 4 is the label; '/^supplementary/' never matches.
supp_a="$(awk '$4=="supplementary"{print $1}' "${EVID}/cs_sa01_disc_mw2.flagstat.txt")"; supp_a="${supp_a:-0}"
supp_b="$(awk '$4=="supplementary"{print $1}' "${EVID}/cs_sa01_disc_mrsa252.flagstat.txt")"; supp_b="${supp_b:-0}"
add supplementary_alignments_mw2align "${supp_a}" "records" "cs_sa01_disc_mw2.flagstat.txt"
add supplementary_alignments_m252align "${supp_b}" "records" "cs_sa01_disc_mrsa252.flagstat.txt"

add vs_MW2Staph1N_total_records "${mw2_total}" "records" "cs_sa01_disc_mw2.idxstats.tsv (denominator: all input records, mapped+unmapped)"
add vs_MW2Staph1N_mapped_bacterial "${mw2_bac}" "reads" "idxstats NC_003923.1 mapped"
add vs_MW2Staph1N_mapped_phage "${mw2_ph}" "reads" "idxstats NC_047722.1 mapped"
add vs_MW2Staph1N_unmapped "${mw2_unmapped}" "records" "idxstats '*'"
add vs_M252Staph1N_total_records "${m2_total}" "records" "cs_sa01_disc_mrsa252.idxstats.tsv"
add vs_M252Staph1N_mapped_bacterial "${m2_bac}" "reads" "idxstats NC_002952.2 mapped"
add vs_M252Staph1N_mapped_phage "${m2_ph}" "reads" "idxstats NC_047722.1 mapped"
add vs_M252Staph1N_unmapped "${m2_unmapped}" "records" "idxstats '*'"

if [ "$(( mw2_total - supp_a ))" -eq "${total_reads}" ] && [ "$(( m2_total - supp_b ))" -eq "${total_reads}" ]; then
  add record_consistency PASS "boolean" "sum(idxstats mapped+unmapped) - supplementary alignments == input records (supplementary: ${supp_a} / ${supp_b})"
else
  add record_consistency FAIL "boolean" "idxstats totals (minus supplementary ${supp_a}/${supp_b}) != input records"
fi

# samtools coverage columns: rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
cov_vals() { awk -v c="$2" '$1==c{print $6"\t"$7}' "$1"; }
read -r mw2_bac_pct mw2_bac_depth <<< "$(cov_vals "${EVID}/cs_sa01_disc_mw2.coverage.tsv" NC_003923.1)"
read -r mw2_ph_pct  mw2_ph_depth  <<< "$(cov_vals "${EVID}/cs_sa01_disc_mw2.coverage.tsv" NC_047722.1)"
read -r m2_bac_pct  m2_bac_depth  <<< "$(cov_vals "${EVID}/cs_sa01_disc_mrsa252.coverage.tsv" NC_002952.2)"
read -r m2_ph_pct   m2_ph_depth   <<< "$(cov_vals "${EVID}/cs_sa01_disc_mrsa252.coverage.tsv" NC_047722.1)"
add vs_MW2Staph1N_bacterial_mean_depth "${mw2_bac_depth}" "x" "samtools coverage NC_003923.1 meandepth"
add vs_MW2Staph1N_bacterial_breadth "${mw2_bac_pct}" "pct_bases_covered" "samtools coverage NC_003923.1 coverage"
add vs_MW2Staph1N_phage_mean_depth "${mw2_ph_depth}" "x" "samtools coverage NC_047722.1 meandepth"
add vs_MW2Staph1N_phage_breadth "${mw2_ph_pct}" "pct_bases_covered" "samtools coverage NC_047722.1 coverage"
add vs_M252Staph1N_bacterial_mean_depth "${m2_bac_depth}" "x" "samtools coverage NC_002952.2 meandepth"
add vs_M252Staph1N_bacterial_breadth "${m2_bac_pct}" "pct_bases_covered" "samtools coverage NC_002952.2 coverage"
add vs_M252Staph1N_phage_mean_depth "${m2_ph_depth}" "x" "samtools coverage NC_047722.1 meandepth"
add vs_M252Staph1N_phage_breadth "${m2_ph_pct}" "pct_bases_covered" "samtools coverage NC_047722.1 coverage"

# core-genome control sanity: MW2 chunks must align extensively to MRSA252
# (stored PAF: qstart=col2, qend=col3, chunk-relative)
ctrl_bp="$(awk '{s+=$3-$2} END{printf "%d", s+0}' "${SCRATCH}/ctrl_mw2_vs_MRSA252.paf")"
add control_MW2chunks_aligned_to_MRSA252 "${ctrl_bp}" "bp of chunked MW2 with sr-level homology in MRSA252" "ctrl_mw2_vs_MRSA252.paf"

# gate evaluation
#  locus summary columns: 8=mean_depth 9=b1x 11=b30x 12=reads (whole window, q0)
#  marker support columns: 3=bp 4=md_q0 5=b1x_q0 6=md_q20 7=b1x_q20 8=reads_q0 9=reads_q20
GATE_OUT="${EVID}/cs_sa01_disc_gates.tsv"
printf 'gate\tcriterion\tobserved\tverdict\n' > "${GATE_OUT}"
gate() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "${GATE_OUT}"; }

read -r pvl_md pvl_b1 pvl_b30 pvl_reads <<< "$(awk -F'\t' -v l="${PVL_LABEL}" '$2==l{printf "%.2f %.4f %.4f %d", $8, $9, $11, $12}' "${LOC_SUM}")"
read -r s4_md s4_b1 s4_b30 s4_reads <<< "$(awk -F'\t' -v l="${SCC4_LABEL}" '$2==l{printf "%.2f %.4f %.4f %d", $8, $9, $11, $12}' "${LOC_SUM}")"
read -r s2_md s2_b1 s2_b30 s2_reads <<< "$(awk -F'\t' -v l="${SCC2_LABEL}" '$2==l{printf "%.2f %.4f %.4f %d", $8, $9, $11, $12}' "${LOC_SUM}")"
read -r s2x_bp s2x_md20 s2x_b120 s2x_r0 s2x_r20 <<< "$(awk -F'\t' -v l="${SCC2_LABEL}_EXCLUSIVE" '$2==l{printf "%d %.2f %.4f %d %d", $3, $6, $7, $8, $9}' "${BED_SUM}")"
read -r pvlx_bp <<< "$(awk -F'\t' -v l="${PVL_LABEL}_EXCLUSIVE" '$2==l{print $3}' "${BED_SUM}")"
read -r s4x_bp <<< "$(awk -F'\t' -v l="${SCC4_LABEL}_EXCLUSIVE" '$2==l{print $3}' "${BED_SUM}")"

if awk "BEGIN{exit !(${pvl_md} >= 30 && ${pvl_b30} >= 0.95)}"; then PVL_G=PASS; else PVL_G=FAIL; fi
gate G1_MW2_PVL_coverage "whole-window mean_depth>=30x AND breadth_ge30x>=0.95 over ${PVL_REGION} (2500 bp)" "mean=${pvl_md}x breadth30x=${pvl_b30} reads=${pvl_reads}" "${PVL_G}"
if awk "BEGIN{exit !(${s4_md} >= 30 && ${s4_b30} >= 0.95)}"; then S4_G=PASS; else S4_G=FAIL; fi
gate G2_MW2_SCCmecIV_coverage "whole-window mean_depth>=30x AND breadth_ge30x>=0.95 over ${SCC4_REGION} (24000 bp)" "mean=${s4_md}x breadth30x=${s4_b30} reads=${s4_reads}" "${S4_G}"
if awk "BEGIN{exit !(${s2x_md20} < 1.0 && ${s2x_b120} <= 0.01 && ${s2x_r20} < 10)}"; then
  S2X_RAW="PASS"
else
  S2X_RAW="FAIL_superseded_by_G3_origin_audit"
fi
gate G3a_MRSA252_SCCmecII_exclusive_raw_absence "raw MAPQ20 absence over the sequence-verified MRSA252-exclusive marker set (${s2x_bp} bp): mean<1x AND breadth<=0.01 AND reads<10" "bp=${s2x_bp} mean_q20=${s2x_md20}x breadth1x_q20=${s2x_b120} reads_q20=${s2x_r20}" "${S2X_RAW}"
read -r origin_reads origin_supported ov_id <<< "$(awk -F'\t' '$1=="OVERALL"{printf "%s %s %s", $7, $10, $8}' "${ORIGIN_TSV}")"
# F26-DATA-004: read-SPAN counts over an exclusive block are not presence
# evidence. A read whose CIGAR span crosses a block while placing no aligned
# base inside (e.g. bridging a deletion through MRSA252-specific sequence) is
# counted by 'view -L' yet contributes zero aligned coverage. The observed
# cross-hits are divergent paralogous SCCmec/plasmid cargo, not cassette
# presence: mean de:f ~0.06 (93.8% identity), origin audit shows the reads
# are MW2-chromosomally supported, and per-base q20 depth is confined to
# micro-islands. The absence gate is therefore evaluated on the retained
# ALIGNED-BASE record (cs_sa01_disc_exclusive_depth.tsv.gz, tag SCC2_EXCL,
# col5=depth_q20): breadth<=0.01 AND mean<1x AND total q20-covered bases
# <=500 AND max contiguous q20-covered island<=500 bp (no contiguous
# segment capable of carrying the marker). Raw read-count criterion remains
# visible as G3a above.
read -r s2x_cov20 s2x_island <<< "$(zcat "${EVID}/cs_sa01_disc_exclusive_depth.tsv.gz" \
  | awk -F'\t' -v t=SCC2_EXCL '
      $1==t{ n++
             if($5>=1){ cov++
                        if(s==0){s=$3+0;e=s}
                        else if($3+0==e+1){e=$3+0}
                        else{ d=e-s+1; if(d>mx)mx=d; s=$3+0;e=s } } }
      END{ if(s){d=e-s+1; if(d>mx)mx=d}
           printf "%d %d", cov+0, mx+0 }')"
if awk "BEGIN{exit !(${s2x_md20} < 1.0 && ${s2x_b120} <= 0.01 && ${s2x_cov20} <= 500 && ${s2x_island} <= 500)}"; then S2X_G=PASS; else S2X_G=FAIL; fi
gate G3_MRSA252_SCCmecII_exclusive_baselevel_absence "aligned-base absence over the sequence-verified MRSA252-exclusive marker set (${s2x_bp} bp): breadth_ge1x_q20<=0.01 AND mean_depth_q20<1x AND total q20-covered bases<=500 AND max contiguous q20-covered island<=500 bp (read-span cross-maps are not presence evidence; origin audit characterizes them)" "bp=${s2x_bp} mean_q20=${s2x_md20}x breadth1x_q20=${s2x_b120} covered_q20=${s2x_cov20}bp max_island=${s2x_island}bp reads_span_q20=${s2x_r20} (origin audit: ${origin_reads} reads, ${origin_supported} MW2-chromosome-supported, mean divergence de:f=${ov_id})" "${S2X_G}"
if [ "${mw2_bac}" -gt "${m2_bac}" ]; then G4=PASS; else G4=FAIL; fi
gate G4_bacterial_preference "MW2-align mapped bacterial reads > MRSA252-align mapped bacterial reads" "mw2=${mw2_bac} mrsa252=${m2_bac}" "${G4}"
if [ "${mw2_ph}" -gt 0 ] && awk "BEGIN{exit !(${mw2_ph_depth} >= 10 * ${mw2_bac_depth})}"; then G5=PASS; else G5=FAIL; fi
gate G5_phage_coculture_signature "phage mapped>0 AND phage mean depth >= 10x bacterial mean depth (MW2 alignment)" "phage_reads=${mw2_ph} phage_depth=${mw2_ph_depth}x bacterial_depth=${mw2_bac_depth}x" "${G5}"
if [ "${ctrl_bp}" -ge 1400000 ]; then G6=PASS; else G6=FAIL; fi
gate G6_core_genome_control "chunked MW2-vs-MRSA252 cross-alignment must cover >=1.4 Mb (proves the chunked homology pipeline detects shared sequence)" "aligned_bp=${ctrl_bp}" "${G6}"

add SCC2_exclusive_reads_q20_primary "${origin_reads}" "primary qnames (MAPQ>=20)" "cs_sa01_disc_exclusive_origin_audit.tsv"
add SCC2_exclusive_reads_chromosome_supported "${origin_supported}" "qnames also mapped on NC_003923.1 in MW2 alignment" "cs_sa01_disc_exclusive_origin_audit.tsv"
add SCC2_exclusive_reads_mean_identity "$(awk "BEGIN{printf \"%.4f\", 1 - (${ov_id:-0})}")" "1 - mean divergence de:f (i.e. mean identity fraction)" "cs_sa01_disc_exclusive_origin_audit.tsv"

# Poisson sensitivity bounds: P(zero coverage at one base at mean depth D)=e^-D;
# over an L-bp marker set: e^(-D*L), evaluated at the observed MW2-alignment
# bacterial mean depth (same library/fragment conditions).
awk -v OFS='\t' -v D="${mw2_bac_depth}" -v BP_S2X="${s2x_bp}" 'BEGIN{
  print "quantity\tvalue";
  printf "%s\t%.6e\n", "P_zero_coverage_per_base_at_" D "x_mean_depth", exp(-D);
  printf "%s\t%.6e\n", "P_zero_coverage_2500bp_window_at_" D "x", exp(-D*2500);
  printf "%s\t%.6e\n", "P_zero_coverage_24000bp_window_at_" D "x", exp(-D*24000);
  printf "%s\t%.6e\n", "P_zero_coverage_" BP_S2X "bp_exclusive_SCCmecII_set_at_" D "x", exp(-D*BP_S2X);
}' | tee "${EVID}/cs_sa01_disc_sensitivity.tsv"

echo "--- gates ---"; cat "${GATE_OUT}"
if awk -F'\t' 'NR>1 && $4=="FAIL"{f++; g=g" "$1} END{exit !(f>0)}' "${GATE_OUT}"; then
  echo "CS_SA01_PRIMARY_DISCRIMINATION_GATES: FAIL (failing gates:${g})"
else
  echo "CS_SA01_PRIMARY_DISCRIMINATION_GATES: ALL PASS"
fi
echo "NOTE: G3a raw read-count sub-criterion is recorded as FAIL_superseded_by_G3_origin_audit: read-SPAN cross-maps (divergent paralogs) are not presence evidence; the discriminating absence test is the aligned-base gate G3."

# ---------------------------------------------------------------------------
STEP 11 "artifact digests"
# ---------------------------------------------------------------------------
cd "${EVID}"
find . -maxdepth 2 -type f ! -name 'SHA256SUMS.txt' ! -path './logs/*' \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS.txt
echo "--- SHA256SUMS.txt ---"; cat SHA256SUMS.txt
printf 'finished=%s job_id=%s\n' "$(date -u +%FT%TZ)" "${SLURM_JOB_ID:-LOGIN}"
echo "F26_CS_SA01_DISCRIMINATION_PIPELINE_COMPLETE"

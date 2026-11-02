# F26 Stage H Audit Closure Addendum

**Task:** F26-AUDIT-005 / DISP-F26-026
**Date:** 2026-11-02
**Auditor:** Codex Sol, independent scientific auditor under Policy V3.2
**Audited worktree HEAD:** `444f95ba460a914ced2eb6b2058e7b9e47f7b6e2`
**Frozen comparator:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528`
**Verdict:** **PASS: `F26_SCIENTIFICALLY_VALIDATED / READY_FOR_USER_ACCEPTANCE`**

## Determination and scope

F26-H01 and F26-H02 are **CLOSED by measured remedy evidence**, without a waiver or change to either frozen acceptance contract. This addendum supersedes the acceptance-withheld determination and open H01/H02 dispositions in the 2026-10-30 [final audit](F26_FINAL_SCIENTIFIC_ACCEPTANCE_REPORT.md). That audit's qualified assessment of Stages B–F, scientific tradeoffs, and F26-H03/H04 interpretation limitations remains applicable. Scientific validation here means completion of the agreed, bounded experimental panel; it does not certify universal assembly accuracy, superiority, or invariants from successful process exits.

The Controller's HANDOFF and ACTIVE_TASK still describe the pre-remedy withheld state. They were read but are outside this dispatch's two-file write ownership and require Controller synchronization. This certification does not adjudicate the separately recorded F23/C0-02R status conflict or authorize scientific-contract changes, integration, push, or release. Final user acceptance remains outstanding.

## Independent checks performed

Read the raw directories and reconstructed collation in memory, without regenerating stored benchmark outputs. The present collection contains **44 records**, not 56: 12 replacement cross-species full-depth runs, 12 retained competitor runs, 16 retained gradient runs, and four new historical oracle replays. Its records and determinism reports equal the current collator's reconstruction; all 44 TSV rows equal their flattened JSON records. Independently calculated SHA-256 matches **44/44 contig FASTAs**. Direct Decimal comparisons against QUAST `report.tsv` match **572/572 metric cells** (13 per run). Per-job summaries agree with their loaded raw metrics for assembly hashes, metrics, and complete identity. Older competitor/gradient metrics use the documented in-memory backfill from per-job summaries; all **16 remedy metrics files themselves contain all 12 identity fields**, with no backfill needed.

Freshly rehashed **22 distinct remedy FASTQ/reference files** (16 mates plus six references) and the release executable. Every digest matches the recorded provenance. Remedy commands use the full registered inputs, eight threads, fastBact `--min-len 500 --kmers auto --depth 120`, and QUAST v5.2.0 with the registered reference and `--min-contig 200`. These checks read retained artifacts; no assembly, QUAST execution, Rust compilation or release test was run on the login node.

Read-only scheduler accounting (`sacct -X -n -P --format=JobID,JobName%40,Partition,State,ExitCode,AllocCPUS,ReqMem,NodeList`) confirms all 16 listed jobs **COMPLETED, exit 0:0, qfree, 8 CPUs, 14G**, on vnode1/vnode2. Each corresponding stdout log identifies the job and dataset/replicate output directory and records harness PASS. All raw assembler and QUAST exit codes are zero.

| Dataset | Remedy jobs, ordered by replicate | Recorded execution Git HEAD |
| --- | --- | --- |
| CS-SA01 | 3095492, 3095493, 3095494 | `640fd0b916121d4fff42c8dd3953f14ca2cb19a4` |
| CS-BS01 | 3095495, 3095496, 3095497 | same |
| CS-PA01 | 3095498, 3095499, 3095500 | same |
| CS-SE01 | 3095501, 3095502, 3095503 | same |
| R1-001 | 3095417 | `aeb1bb648bcf68da19f8ff2b8402762fa748ffc4` |
| R1-003 | 3095418 | same |
| R1-001-Val | 3095419 | same |
| R1-003-Val | 3095420 | same |

Execution HEADs are recorded run provenance, not the later documentation audit HEAD. Their `src` trees and current HEAD's source tree are all `fe98d216f11b2b579cd358de2418773d9e4ff683`. `git diff --exit-code 02872f81811402d2bd603ff7c0a2ecf8b3dbe528 -- src/` exits 0 with no output. Release binary SHA-256 is `517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1` in all remedy records and in the current executable. The retained independent release-test log for qfree job 3076232 keeps its original digest and records 213 passed, zero failed/ignored and exit 0; this audit does not claim a new release-suite execution.

## F26-H01: CLOSED

For each species, independently verified replicate indices exactly {1,2,3}, identical complete identity tuples (dataset/tool, read paths and SHA-256, reference path and SHA-256, binary SHA-256, Git HEAD, threads, memory), PASS statuses, identical FASTA bytes, and exact equality of all 13 QUAST metrics. Current aggregate certification agrees with these checks.

| Dataset | Replicates | Common assembly SHA-256 | Stored variance_max |
| --- | --- | --- | --- |
| CS-BS01 | 3 | `2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861` | 2.0194839173657902e-28 |
| CS-PA01 | 3 | `b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff` | 0.0 |
| CS-SA01 | 3 | `86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962` | 0.0 |
| CS-SE01 | 3 | `4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb` | 0.0 |

BS01's `2.0194839173657902e-28` is floating-point averaging residue from three exactly equal decimal values; independently compared metric values are exactly equal, so their exact population variance is zero. No tolerance was used to hide a metric difference. Runtime/RSS are excluded by contract §10.5.

Inspected identity stamping and grouping code, and ran both existing safe self-tests: `python3 -B scripts/run_competitor_benchmark.py --self-test` and `python3 -B scripts/collate_benchmark_results.py --self-test`; both exit 0 with SELF_TEST_PASS. They exercise disk identity propagation, missing/empty identity rejection, changed-identity separation, legacy in-memory backfill, and orphan rejection. Identity differences produce separate groups that cannot collectively claim N=3; the collator need not reject the entire collection when distinct complete experiments exist. The successful present groups therefore satisfy §10.7 without pooling the old mismatched HEADs.

## F26-H02: CLOSED

Compared each replay directly with the frozen F26 ledger and with retained F25 Job 3072539 artifacts in `../fastBact-phase-f25-correctness-worktree/benchmark-results/f25_validation_20261030/`. All seven ledger quality metrics have exact zero deltas. Additionally all four replay FASTAs equal their F25 FASTAs byte for byte, and all 13 comparable QUAST metrics agree (F25 `# contigs` equals current `# contigs` and `# contigs (>= 200 bp)` under the same 200-bp cutoff). This provides direct paired baseline/output identity, beyond source preservation or rounded summary claims.

| Dataset | N50 | NGA50 | GF % | Extensive | Local | Mismatch/100kb | Indel/100kb |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R1-001 | 48115 | 45092 | 98.115 | 5 | 0 | 0.85 | 0.53 |
| R1-001-Val | 54554 | 54554 | 98.233 | 0 | 0 | 0.81 | 0.24 |
| R1-003 | 44896 | 44702 | 98.06 | 0 | 1 | 0.13 | 0.0 |
| R1-003-Val | 40664 | 39314 | 98.097 | 0 | 1 | 0.26 | 0.0 |

Each historical dataset has one new replay, not a fresh N=3 F26 determinism certification. Its bitwise equality to the retained F25 assembly and exact paired metric equality demonstrate the frozen non-regression predicate. These runs are not pooled with F25 under a shared Git/binary identity; the F25 binary differs. No execution-environment identity across phases is asserted, and no runtime/RSS equality is required. Raw R1-001 5/0 and BW25113 0/1 discrepancy counts remain visible; cross-species events receive no new biological exemption.

| Dataset | F25 and F26 common FASTA SHA-256 | F25 QUAST report SHA-256 |
| --- | --- | --- |
| R1-001 | `62d0c77f0daaec2aef1cc5f8ba92c387ffdae8cb9aa359071b87dc2080be1564` | `4bdcd8f96243c5c586452ebdbf9b0ebc770f1f06fd406edb5fabfce914faccb0` |
| R1-001-Val | `108679e2c060fcfbbc1bd94121af5f56f10dd75be4a203294bfd48ced7b9307c` | `0ad296b2dd56b1084b553dc9ec401d7e3d010ceb1d4f95f8b8b3f5dd3e6388a1` |
| R1-003 | `bf262998da3bd0fb1ddc3f3720cc80d3a5e28da00e9e450bda09113e57ec215e` | `d55d11eefc0eaaa753930de79d5ef681cc04c3afa6a40de4494c72f136fc9e42` |
| R1-003-Val | `550a7002aea1e542c6f5a75364d548db5984ce7834c6942a1033c914d1ff0d5d` | `40445d750e1858e7618611ec0b9f8f83396c49a44492839ff73a4977fbe1afa7` |

## Remaining qualifications and authorization boundary

F26-H03 and H04 remain interpretation limitations: MW2 phage dilution and nominal depth, PAO1's clamped top-depth input, lineage differences and unresolved raw cross-species events. Retain original competitor measurements when quoting their historical speed ratios; replacement fastBact replicate runtimes in the current summary are new measurements. Single-run competitor comparisons, poorer contiguity, taxon-dependent accuracy and LT2 memory cost remain material limits. This closure grants no universal determinism claim and supplies no new biological attribution.

Only the two dispatched reports were changed. Protected source/data/scripts and benchmark artifacts were read only. Local main remains `02872f81811402d2bd603ff7c0a2ecf8b3dbe528`; no commit, merge, push, tag or release was performed. Remote state was not queried in this closure audit. The old report's remote observations remain dated historical evidence.

## Evidence digests and reproduction

SHA-256 of the snapshots inspected:

| Path | SHA-256 |
| --- | --- |
| `F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` | `f2913504f8d0c533a7012b9ae22d5b69842f5522ca247a4ce32be3e40e4d1f0d` |
| `F26_COMPETITOR_BENCHMARK_CONTRACT.md` | `bfe46ff5e4caacd06399fc0fcabd65627b54a45271dc8be25740b121407be80a` |
| `F26_DATASET_MANIFEST.tsv` | `e23dd244b42d7d565b46b68483c2b68ab1cedfd024b7f28228b3f8e31ac70dc9` |
| `F26_BENCHMARK_LEDGER.md` | `188377319f64ee776b4d6d8d67631c2195cb4f572c242217d7e2751ddd77572d` |
| `scripts/run_competitor_benchmark.py` | `861042f6fb49c74321aaa2a4076c716c8e5533b18b1af2dfcaf55bb73aff68c4` |
| `scripts/collate_benchmark_results.py` | `cbfac1d6c00d1f21613f382e7321819f79da129049cfd7cb2b17c2d5c499bc83` |
| `scripts/submit_f26_remedy_benchmarks.sh` | `dc0b52b2c0b73b98d669360bfc3af6bf4f26428e5a1b454c6641c92cebbee7a7` |
| `benchmarks/f26_results/summary.json` | `f43c73552428430de90b96077ef1a88b50c878a865fa9873a52713674c02fa56` |
| `benchmarks/f26_results/summary.tsv` | `d90fcf3169c6b224d18c964e4ed5daf36408b2b40a8103b43ad3ae22008e7605` |
| `logs/f26_sol_release_audit_3076232.log` | `9eb653edca0e843f1f37f5870a7dab1bfdecea54d40151b63a02c3df618386a7` |

The following 96-file manifest binds the remedy metrics, per-job summaries, FASTAs, QUAST tables and 32 Slurm logs. Its SHA-256 is **`8e23809c3f796cb2d2f014f52cd6d48c53928b6114eb4132550df4946686450e`**. To reproduce, sort the relative paths lexicographically and concatenate UTF-8 lines `SHA256`, two spaces, `relative/path`, LF, then SHA-256 that byte string. Hash files in binary mode. Input checks can be repeated from the 16 metrics records' r1/r2/reference fields, deduplicating paths; compare hashes with their corresponding `_sha256` fields. Reconstruct collation in memory or a temporary output directory, never over the retained evidence.

```text
2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861  benchmarks/f26_results/raw/CS-BS01_fastbact_r1/CS-BS01/fastbact/rep_1/contigs.fa
92bd3b529cbd21e6a46761d48625f53f722e6cc027ee7f6fe117ad83f12fd88b  benchmarks/f26_results/raw/CS-BS01_fastbact_r1/CS-BS01/fastbact/rep_1/metrics.json
7eec676dbe2cb623b5d732830e6592439d1a61dc382400cd5cfed89213f9d38a  benchmarks/f26_results/raw/CS-BS01_fastbact_r1/CS-BS01/fastbact/rep_1/quast/report.tsv
3e1c863a7770572a55c200fcc30bb8275d41d26d5add6d624b355f23810c85e3  benchmarks/f26_results/raw/CS-BS01_fastbact_r1/summary.json
2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861  benchmarks/f26_results/raw/CS-BS01_fastbact_r2/CS-BS01/fastbact/rep_2/contigs.fa
61d53170bdc9d33fe034d85732934fbcf746d0cd0f4ed54e1d1681c2ac95d346  benchmarks/f26_results/raw/CS-BS01_fastbact_r2/CS-BS01/fastbact/rep_2/metrics.json
f7a86c96a6cff9e89637afa963d99f726aee1e5e932a6b2bace532c2927deff1  benchmarks/f26_results/raw/CS-BS01_fastbact_r2/CS-BS01/fastbact/rep_2/quast/report.tsv
b5882a1358f19c85b6ecf707d20aa79073ef0eb7e40f346cdeb918b85231156d  benchmarks/f26_results/raw/CS-BS01_fastbact_r2/summary.json
2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861  benchmarks/f26_results/raw/CS-BS01_fastbact_r3/CS-BS01/fastbact/rep_3/contigs.fa
d7b87f0a555a73e58264fb48341127ea4ff67e3f7d20ccc4a6bc8ce1ac6775ae  benchmarks/f26_results/raw/CS-BS01_fastbact_r3/CS-BS01/fastbact/rep_3/metrics.json
ffdf9d8090cbf09ebeb0fe712f32c1adcd2b9f8206fbbf6f2d508b069e8ecee9  benchmarks/f26_results/raw/CS-BS01_fastbact_r3/CS-BS01/fastbact/rep_3/quast/report.tsv
b5a8097ce7d1b0267a884494ec2bc9a961937de19ee9d0f943173551aa80298b  benchmarks/f26_results/raw/CS-BS01_fastbact_r3/summary.json
b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff  benchmarks/f26_results/raw/CS-PA01_fastbact_r1/CS-PA01/fastbact/rep_1/contigs.fa
4f2bcb1c9d6de200c8c14400cbb4672e6f7aea818e9d42f2efe1c032a7c07865  benchmarks/f26_results/raw/CS-PA01_fastbact_r1/CS-PA01/fastbact/rep_1/metrics.json
86686d79a645987d16a4910b1aa93d092f736e40965467f800831977b7de58be  benchmarks/f26_results/raw/CS-PA01_fastbact_r1/CS-PA01/fastbact/rep_1/quast/report.tsv
189f14a2d4b21774eb9f214a6feabeecf233eaf62b142d5bd4bfc00cd2bc8394  benchmarks/f26_results/raw/CS-PA01_fastbact_r1/summary.json
b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff  benchmarks/f26_results/raw/CS-PA01_fastbact_r2/CS-PA01/fastbact/rep_2/contigs.fa
36f520c176a0da63662ef1c9698cdc05b2a428fb770009cfd63412b5d5703f89  benchmarks/f26_results/raw/CS-PA01_fastbact_r2/CS-PA01/fastbact/rep_2/metrics.json
80276dee5d448ac250691a42c9ebee6b00cafe6ae83e0c1e9e0ddb28fede37cc  benchmarks/f26_results/raw/CS-PA01_fastbact_r2/CS-PA01/fastbact/rep_2/quast/report.tsv
eeaedc5742ac3127022856e204e8c8c55930d038bd4b61df61e628dace54cc9d  benchmarks/f26_results/raw/CS-PA01_fastbact_r2/summary.json
b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff  benchmarks/f26_results/raw/CS-PA01_fastbact_r3/CS-PA01/fastbact/rep_3/contigs.fa
aa3598e6acec711bbfee4dc5d04a9b98e4b6478b4ac8b23c3d9246eb3b8fbb50  benchmarks/f26_results/raw/CS-PA01_fastbact_r3/CS-PA01/fastbact/rep_3/metrics.json
bf85dcaafa6f992e320db6b54b56857b7fddefbbcda2d888975c0a2c2b31ba96  benchmarks/f26_results/raw/CS-PA01_fastbact_r3/CS-PA01/fastbact/rep_3/quast/report.tsv
8aa6b7bb1641751f79049a7c2649c2e88a4c607cb470c646307eab14d1fdaa77  benchmarks/f26_results/raw/CS-PA01_fastbact_r3/summary.json
86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962  benchmarks/f26_results/raw/CS-SA01_fastbact_r1/CS-SA01/fastbact/rep_1/contigs.fa
58b4f64bc5bb0a6da59c636e6b0d29768f6fafc896d2c503880d1867d1d5f577  benchmarks/f26_results/raw/CS-SA01_fastbact_r1/CS-SA01/fastbact/rep_1/metrics.json
9810dc689b8f2c40033b8e2ab004571397f70a165d2f7be31e566ea0249eeb26  benchmarks/f26_results/raw/CS-SA01_fastbact_r1/CS-SA01/fastbact/rep_1/quast/report.tsv
d4812675b95d32eeb7920efb5e1d70d387b712a4f9dfeb3b39a9efe6768485ab  benchmarks/f26_results/raw/CS-SA01_fastbact_r1/summary.json
86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962  benchmarks/f26_results/raw/CS-SA01_fastbact_r2/CS-SA01/fastbact/rep_2/contigs.fa
4b08b5e335ba054554149c1ac86bdccf1d7a809549af29ff6ff7adffb219975f  benchmarks/f26_results/raw/CS-SA01_fastbact_r2/CS-SA01/fastbact/rep_2/metrics.json
20afdfdd20bbb421bb6bc391ef9063259764916a7163c54e6a577610a9cbdc86  benchmarks/f26_results/raw/CS-SA01_fastbact_r2/CS-SA01/fastbact/rep_2/quast/report.tsv
134c69360bad8fd62740220ca751ac0e2dde1c7e790fb5844ddfac75ed80cb9f  benchmarks/f26_results/raw/CS-SA01_fastbact_r2/summary.json
86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962  benchmarks/f26_results/raw/CS-SA01_fastbact_r3/CS-SA01/fastbact/rep_3/contigs.fa
e242bb45cf946044a81350c6e050d6bf514156611826b101ada1b1681b405f87  benchmarks/f26_results/raw/CS-SA01_fastbact_r3/CS-SA01/fastbact/rep_3/metrics.json
bd5bc0e149f5f28ec4641a320a2ced5d6d5a34d7dc64305e93cc1b579a9398fb  benchmarks/f26_results/raw/CS-SA01_fastbact_r3/CS-SA01/fastbact/rep_3/quast/report.tsv
cc4ace5d1928ca8d4911e05189a044bc1800e6c95c3ac94e5f2292ffa2e28b64  benchmarks/f26_results/raw/CS-SA01_fastbact_r3/summary.json
4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb  benchmarks/f26_results/raw/CS-SE01_fastbact_r1/CS-SE01/fastbact/rep_1/contigs.fa
90958f4c55491fa2e0ca3668eceef4f918fc440eceb040f8f8b9bfd445075d99  benchmarks/f26_results/raw/CS-SE01_fastbact_r1/CS-SE01/fastbact/rep_1/metrics.json
13517e787c804748da76c6b7cb7a36529a2fb199adfe6893afdc6beab75625d2  benchmarks/f26_results/raw/CS-SE01_fastbact_r1/CS-SE01/fastbact/rep_1/quast/report.tsv
4fb4ef02cfb31b2a2081bb62f30fd7050ac4b0bc5d43624fb9511bf77680c329  benchmarks/f26_results/raw/CS-SE01_fastbact_r1/summary.json
4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb  benchmarks/f26_results/raw/CS-SE01_fastbact_r2/CS-SE01/fastbact/rep_2/contigs.fa
5010894f939a0cf7670707d64e096b5a0a79ea1fe17c553327133f222ecf1505  benchmarks/f26_results/raw/CS-SE01_fastbact_r2/CS-SE01/fastbact/rep_2/metrics.json
8cd67a7bfe1922e20f4a6953a143c22f109ded1ae009bc9a001dfcdafe40d25c  benchmarks/f26_results/raw/CS-SE01_fastbact_r2/CS-SE01/fastbact/rep_2/quast/report.tsv
90b5dbadaadb74dd23e9e21aab79bfb078e00561d5a3c1d3cba3670cf8e5d5d6  benchmarks/f26_results/raw/CS-SE01_fastbact_r2/summary.json
4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb  benchmarks/f26_results/raw/CS-SE01_fastbact_r3/CS-SE01/fastbact/rep_3/contigs.fa
6d68d7a1ba7cee1e4f6bf42472e6ac324bbaab97bb0ebb9fbb0c695245da47fc  benchmarks/f26_results/raw/CS-SE01_fastbact_r3/CS-SE01/fastbact/rep_3/metrics.json
5c16dc795b4f679f723b8a92b6ea3e5eae48e64e26d365493275dd57e0d96a6d  benchmarks/f26_results/raw/CS-SE01_fastbact_r3/CS-SE01/fastbact/rep_3/quast/report.tsv
f0bffd4fdfa10e231c6b589da58369cc66fb6d603dcf81acd748a3ed607f9ee7  benchmarks/f26_results/raw/CS-SE01_fastbact_r3/summary.json
108679e2c060fcfbbc1bd94121af5f56f10dd75be4a203294bfd48ced7b9307c  benchmarks/f26_results/raw/R1-001-Val_fastbact_r1/R1-001-Val/fastbact/rep_1/contigs.fa
5680fb000f3eb49927dadffef783eaa71651d4f88c67974af2bc69f429db3a74  benchmarks/f26_results/raw/R1-001-Val_fastbact_r1/R1-001-Val/fastbact/rep_1/metrics.json
ae63b003c274ffed90acbdf2480d85f16e29785cac316f2b71b3e0a2c040f6bf  benchmarks/f26_results/raw/R1-001-Val_fastbact_r1/R1-001-Val/fastbact/rep_1/quast/report.tsv
9edbee93abb512c50ceb864c0a139b2972efa98673b73269f5319c3ddf095487  benchmarks/f26_results/raw/R1-001-Val_fastbact_r1/summary.json
62d0c77f0daaec2aef1cc5f8ba92c387ffdae8cb9aa359071b87dc2080be1564  benchmarks/f26_results/raw/R1-001_fastbact_r1/R1-001/fastbact/rep_1/contigs.fa
56a89c01d88ce22009fedeb581f7d440ac7e7327961daa20c989a2d1bc7ca5e8  benchmarks/f26_results/raw/R1-001_fastbact_r1/R1-001/fastbact/rep_1/metrics.json
9de890f2d848ec3dbef714f629817c9af3a577de1421c109c0d6ed0ee9d8a567  benchmarks/f26_results/raw/R1-001_fastbact_r1/R1-001/fastbact/rep_1/quast/report.tsv
bb0a37c2d275122212d6da38bfe2660885d6574c6dbce593c91ee5592c5453dd  benchmarks/f26_results/raw/R1-001_fastbact_r1/summary.json
550a7002aea1e542c6f5a75364d548db5984ce7834c6942a1033c914d1ff0d5d  benchmarks/f26_results/raw/R1-003-Val_fastbact_r1/R1-003-Val/fastbact/rep_1/contigs.fa
2a48311a89fbad703c3e5fc694e7bebf8f9b7fb3d332ca93eb6ec3ee9440b7b7  benchmarks/f26_results/raw/R1-003-Val_fastbact_r1/R1-003-Val/fastbact/rep_1/metrics.json
9df73daa0f3af4d0367d79f9d99988584f556ead10abdebfbfeaa8cae5ffd9d1  benchmarks/f26_results/raw/R1-003-Val_fastbact_r1/R1-003-Val/fastbact/rep_1/quast/report.tsv
da22f151e941697984cc76fd4cba3871905078afdbb655a169f875c37ef06771  benchmarks/f26_results/raw/R1-003-Val_fastbact_r1/summary.json
bf262998da3bd0fb1ddc3f3720cc80d3a5e28da00e9e450bda09113e57ec215e  benchmarks/f26_results/raw/R1-003_fastbact_r1/R1-003/fastbact/rep_1/contigs.fa
8b29bac74cb109da68eed93caa5099a27620de8aa59ff4b78c9c9525f2fb59d6  benchmarks/f26_results/raw/R1-003_fastbact_r1/R1-003/fastbact/rep_1/metrics.json
2c62457413e0a423f4fb688fabe58f621e78b66b61a10e22ae9c047783057cae  benchmarks/f26_results/raw/R1-003_fastbact_r1/R1-003/fastbact/rep_1/quast/report.tsv
f22d6377636a583a9e468c7fbfdb5d757031c57f8690724013725491433a0268  benchmarks/f26_results/raw/R1-003_fastbact_r1/summary.json
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095417.err
c9795a105f826a4119e17bae46f6c66b250b772fe3e2ccb9fbefe0a9f08b0a3d  logs/f26_benchmark_3095417.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095418.err
ce37834434b624bf00a3f20c4b20d2f874a5c9fcd4467e28c9a25edb12d63f66  logs/f26_benchmark_3095418.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095419.err
b6af34d011301a75a9b68452c285e91da6fc698d7216693b6e1c7737f11025c2  logs/f26_benchmark_3095419.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095420.err
633729cb9579badf5ec859703f7defa30ac6a9004bc7ebf1c470a36fd9779acc  logs/f26_benchmark_3095420.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095492.err
739bec44f3f0f90d8396b23691155692b23dacf8e0eea8f7287d9853f4bd6e6b  logs/f26_benchmark_3095492.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095493.err
8c7ecc7f7701ac62cba3663a5cdf0092bf7372807d882e9000f37e5d95114ae6  logs/f26_benchmark_3095493.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095494.err
69d5a3266706e8499875c76c3a57f11c9fb775a83b4e6e9a1acb7690bc6134c7  logs/f26_benchmark_3095494.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095495.err
fee428400c0c65fc1778d4258738c664c592e8dcce65bfd90120abad9643856f  logs/f26_benchmark_3095495.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095496.err
750881aa044fed116ba2ac87c0d13eadf1fa20762f0bbe0852964597d8eb10f9  logs/f26_benchmark_3095496.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095497.err
c6dd5204f620d7ee4275d41be7f1df3967abd2483fea9df18eef71b6bb4268de  logs/f26_benchmark_3095497.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095498.err
69ef60a0ae5f7e4c24805426e65145337440b9c626abcf61b3e6fba7caab7d70  logs/f26_benchmark_3095498.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095499.err
673f97f572f976c1c713526ed1bfd2970af31c6a4324a29b3521e16a92fa55cc  logs/f26_benchmark_3095499.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095500.err
f7db8eb453c518aed07fb8c6a0b66a5e1ae772d6ba8d5eb3864ba71ca37d9185  logs/f26_benchmark_3095500.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095501.err
f5772f9ef0e4a6caeb484b7b89c3fe51849c72a532dce991d537d48502cccf52  logs/f26_benchmark_3095501.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095502.err
b0a938e7379c8c648658404a76433a338b451a99e977749bc4e0161e0d857ac3  logs/f26_benchmark_3095502.out
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  logs/f26_benchmark_3095503.err
985e73bf17d63b8f252ad4fb26211d8b42740e6e040690628e5e4b0948cf9147  logs/f26_benchmark_3095503.out
```

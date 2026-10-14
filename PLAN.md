# Project Specification: `fastBact` (Rust Native Shovill + SPAdes Bacterial Assembler)

## 0. Agent Execution Directives & HPC Cluster Protocol (CRITICAL RULES)

You are an expert Rust bioinformatics systems engineer implementing **`fastBact`**, an ultra-fast, low-memory, short-read bacterial genome assembler in pure Rust that fuses **Shovill's preprocessing/postprocessing engineering pipeline** with **SPAdes' iterative multi-k Compacted De Bruijn Graph (cDBG) and `exSPAnder` repeat resolution algorithms** in a single zero-disk-IO streaming binary.

### 0.1 HPC Cluster & Login Node Ban (STRICTLY ENFORCED)
1. **NEVER Run Compute or Compilation on the Login Node**:
   - You are strictly FORBIDDEN from running `cargo build`, `cargo check`, `cargo test`, `cargo clippy`, `rustc`, or any sequence simulation/assembly command directly on the login node.
   - The login node may ONLY be used for:
     - Editing source files and configuration files.
     - Lightweight environment inspection (`conda info --envs`, `which cargo`, `sinfo`, `squeue`, `ls`, `cat`).
     - Submitting jobs to the cluster queue via the wrapper script `./scripts/hpc_run.sh`.
2. **Use an Existing Server Conda Environment**:
   - Do NOT create a new Conda environment (`conda create` is forbidden).
   - In **Phase 0**, inspect existing Conda environments on the server using `conda info --envs` (or `conda env list`), identify an existing environment suitable for development (and locate `~/.cargo/bin` if Rust is installed via `rustup`), and hardcode its activation into `./scripts/hpc_run.sh`.
3. **Synchronous Queue Submission (`sbatch --wait`)**:
   - All `cargo` commands (`check`, `test`, `build --release`) MUST be executed through `./scripts/hpc_run.sh <command...>`, which submits the task to a compute node via `sbatch --wait`, limits CPU threads (`CARGO_BUILD_JOBS` and `RAYON_NUM_THREADS`) to the allocated cores, and streams the compute node's stdout/stderr log back to you upon completion.

### 0.2 Core Architectural & Rust Coding Constraints
1. **NO `String` or `Vec<u8>` per-read heap allocations in hot loops**: All retained reads and unitigs MUST be stored in contiguous 2-bit packed arenas (`Vec<u64>`) indexed by lightweight handles (`ReadId`, `UnitigId`).
2. **NO pointer-based graph nodes (`Rc`, `Arc`, `RefCell`, `Mutex<Node>`)**: The assembly graph MUST use flat Arena allocation (`Vec<UnitigNode>`) with inline small-array edges (`SmallVec` or fixed `[Edge; 4]`) and `u32` node indices.
3. **Compile-time K-mer Polymorphism**: Implement a `KmerInt` trait supporting `u64` ($k \le 31$), `u128` ($33 \le k \le 63$), and `[u64; 4]` ($65 \le k \le 127$). All $k$-mers must be odd integers to avoid palindromic self-reverse-complements.
4. **Zero External Binary Dependencies**: Do NOT shell out to `bwa`, `flash`, `lighter`, `spades.py`, `samtools`, or `pilon`. Avoid C-library build scripts (`cmake`) that might fail on compute nodes by using pure-Rust crates (`flate2` with `rust_backend`).
5. **Incremental Execution**: Execute strictly **Phase 0 -> Phase 6**. Never leave `todo!()` or `unimplemented!()` stubs in completed phases.

### 0.3 Reference Repositories Usage (`references/`)
Local reference files are provided in `./references/` to guide exact heuristic formulas and thresholds:
1. **When implementing Phase 2 & Phase 6 (Ingest & Curation)**:
   - Read `./references/shovill/bin/shovill` (especially the `estimate_genome_size`, dynamic k-mer selection logic, FLASH overlap parameters, and contig header/filtering rules). Port these exact empirical formulas natively into Rust.
2. **When implementing Phase 4 & Phase 5 (Graph Simplification & exSPAnder)**:
   - Read `./references/spades/configs/debruijn/isolate_mode.info` (or under `assembler/configs/debruijn/`), `simplification.info`, and `pe_params.info`.
   - Extract the exact numerical thresholds for tip clipping (`tc`), bulge removal (`br`), low-coverage edge removal (`ecr`), and paired-end extension (`pe_params`) to parameterize `fastBact`'s Rust structs.
   - **DO NOT** copy C++ OOP patterns; implement the mathematical rules using flat Rust `BitSeqArena` and `UnitigNode` structures.
---

## 1. Repository Layout & `Cargo.toml`

Initialize the workspace with the following directory structure:

```text
fastBact/
├── Cargo.toml
├── PLAN.md
├── logs/                       # Compute queue job stdout/stderr logs
├── scripts/
│   └── hpc_run.sh              # Queue submission wrapper (created in Phase 0)
├── src/
│   ├── main.rs                 # CLI entry point & pipeline orchestrator
│   ├── lib.rs                  # Public library exports
│   ├── cli.rs                  # clap derive CLI definitions
│   ├── common/
│   │   ├── mod.rs
│   │   ├── bitseq.rs           # 2-bit DNA packing, arena storage, reverse complement
│   │   ├── kmer.rs             # KmerInt trait (u64, u128, [u64; 4]) & rolling hash
│   │   └── sketch.rs           # HyperLogLog & Count-Min Sketch
│   ├── ingest/
│   │   ├── mod.rs
│   │   ├── estimator.rs        # Genome size, coverage estimation & dynamic k-mer list
│   │   └── trim_merge.rs       # Adapter trimming, subsampling & FLASH-style PE overlap merge
│   ├── graph/
│   │   ├── mod.rs
│   │   ├── cdbg.rs             # Arena-based Compacted De Bruijn Graph data structures
│   │   ├── build.rs            # Solid k-mer filtering + Minimizer-bucketed unitig compaction
│   │   ├── simplify.rs         # Tip clipping, Myers bit-parallel bubble popping, chimeric removal
│   │   └── local_multik.rs     # Local multi-k subgraph bridging around branching nodes
│   ├── resolve/
│   │   ├── mod.rs
│   │   ├── mapper.rs           # O(1) unique k-mer pseudo-alignment & insert size estimation
│   │   └── exspander.rs        # SPAdes exSPAnder paired-end path extension engine
│   └── polish/
│       ├── mod.rs
│       ├── consensus.rs        # In-memory pileup base polishing (replaces BWA + Pilon)
│       └── output.rs           # Terminal overlap trimming, plasmid circularity & FASTA/GFA writer
└── tests/
    └── integration_test.rs     # Synthetic bacterial genome + plasmid end-to-end assembly benchmark
```

### `Cargo.toml` Specification

```toml
[package]
name = "fastBact"
version = "0.1.0"
edition = "2021"
description = "Ultra-fast Rust bacterial short-read genome assembler fusing Shovill and SPAdes"

[[bin]]
name = "fastBact"
path = "src/main.rs"

[dependencies]
clap = { version = "4.5", features = ["derive"] }
needletail = { version = "0.5", default-features = false, features = ["flate2"] }
flate2 = { version = "1.0", default-features = false, features = ["rust_backend"] }
rayon = "1.10"
crossbeam-channel = "0.5"
hashbrown = { version = "0.14", features = ["rayon"] }
rustc-hash = "2.0"
smallvec = { version = "1.13", features = ["union", "const_generics"] }
tracing = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
anyhow = "1.0"

[dev-dependencies]
rand = "0.8"
tempfile = "3.10"

[profile.release]
opt-level = 3
lto = "fat"
codegen-units = 1
panic = "abort"
```

---

## 2. Core Data Contracts (Exact Type Signatures)

Implement these exact foundational types in `src/common/` and `src/graph/`:

```rust
/// 2-bit DNA encoding: A = 0b00, C = 0b01, G = 0b10, T = 0b11
#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Base { A = 0, C = 1, G = 2, T = 3 }

/// Contiguous 2-bit sequence storage arena avoiding heap fragmentation
pub struct BitSeqArena {
    pub data: Vec<u64>,       // 32 bases per u64 word
    pub offsets: Vec<u64>,    // Base start index in `data`
    pub lengths: Vec<u32>,    // Sequence length in bp
}

/// Packed Read Pool separating merged single-end long reads from unmerged paired-end reads
pub struct ReadBank {
    pub merged_reads: BitSeqArena,   // Overlapping PE reads merged into 180-280bp reads
    pub pe_forward: BitSeqArena,     // Unmerged R1 reads
    pub pe_reverse: BitSeqArena,     // Unmerged R2 reads
}

/// Directed stranded node pointer (Unitig ID + Orientation)
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DirUnitigId(pub u32); // Lowest bit = strand (0: forward, 1: reverse complement), upper 31 bits = UnitigId

#[derive(Copy, Clone, Debug)]
pub struct GraphEdge {
    pub target: DirUnitigId,
    pub overlap_k: u8,        // Overlap length (typically k - 1)
    pub support: u16,         // Read transition support count
}

pub struct UnitigNode {
    pub id: u32,
    pub seq_offset: u64,      // Offset inside graph's BitSeqArena
    pub seq_len: u32,         // Length in bases (>= k)
    pub kmer_cov: f32,        // Mean k-mer coverage depth
    pub edges_fwd: smallvec::SmallVec<[GraphEdge; 4]>, // Outgoing from 3' end of forward strand
    pub edges_rev: smallvec::SmallVec<[GraphEdge; 4]>, // Outgoing from 3' end of reverse strand
    pub is_deleted: bool,
    pub is_circular: bool,
}
```

---

## 3. Phase-by-Phase Execution Plan

### Phase 0: HPC Queue & Conda Environment Bootstrap (RUN FIRST)
**Goal**: Discover the server's existing Conda environment and queue partition, and create `./scripts/hpc_run.sh` so ZERO builds or tests ever run on the login node.

#### Tasks:
1. Run lightweight inspection commands on the login node:
   - `conda info --envs` to list existing Conda environments and select an appropriate existing environment (check which environment or `~/.cargo/bin` provides `cargo` and `rustc`).
   - `sinfo -s` (or `qstat` if PBS/Torque) to identify an available idle compute partition/queue (e.g., `cpu`, `compute`, `normal`, or default partition).
2. Create `mkdir -p logs scripts` and write an executable `./scripts/hpc_run.sh`:
   ```bash
   #!/usr/bin/env bash
   set -euo pipefail
   mkdir -p logs

   # 1. Configure discovered Conda base & existing environment name
   CONDA_BASE=$(conda info --base)
   TARGET_ENV="<INSERT_DISCOVERED_CONDA_ENV_NAME>"
   PARTITION="<INSERT_DISCOVERED_SLURM_PARTITION>"
   CPUS=8

   JOB_LOG="logs/job_$(date +%Y%m%d_%H%M%S)_$$.log"

   # 2. Submit synchronously to compute node via sbatch --wait (0 CPU used on login node)
   sbatch --wait \
     --job-name="fastBact_ci" \
     --partition="${PARTITION}" \
     --nodes=1 \
     --ntasks=1 \
     --cpus-per-task="${CPUS}" \
     --mem=16G \
     --output="${JOB_LOG}" \
     --error="${JOB_LOG}" \
     --wrap="
       source \"${CONDA_BASE}/etc/profile.d/conda.sh\"
       conda activate \"${TARGET_ENV}\"
       export PATH=\"\$HOME/.cargo/bin:\$PATH\"
       export CARGO_BUILD_JOBS=${CPUS}
       export RAYON_NUM_THREADS=${CPUS}
       echo \"[HPC Node: \$(hostname)] Running: $*\"
       $*
     "
   EXIT_CODE=$?
   cat "${JOB_LOG}"
   exit ${EXIT_CODE}
   ```
   *(Note: If the default partition requires no `--partition` flag, omit `--partition`. Before first offline compute run, if compute nodes lack internet access, run `cargo fetch` once or ensure vendor cache is accessible).*

#### Phase 0 Verification Gate:
- Run `chmod +x ./scripts/hpc_run.sh`
- Test queue submission with: `./scripts/hpc_run.sh cargo --version`
- Verify that the output log prints `[HPC Node: <compute-node-hostname>]` (confirming it did NOT run on the login node) and returns a valid `cargo` version.

---

### Phase 1: Foundation — 2-Bit Sequence Arena, Generic `KmerInt` & Rolling Sketches
**Goal**: Build zero-allocation sequence encoding, canonical $k$-mer bit operations, HyperLogLog, and Count-Min Sketch.

#### Tasks:
1. **`src/common/bitseq.rs`**:
   - Implement ASCII-to-2bit lookup table (`A/a->0, C/c->1, G/g->2, T/t->3`, non-ACGT bases deterministically mapped).
   - Implement `BitSeqArena` with `push_seq(&mut self, ascii: &[u8]) -> u32`, `get_base(&self, id: u32, pos: usize) -> u8`, `extract_ascii(&self, id: u32) -> Vec<u8>`, and fast bit-level reverse complement.
2. **`src/common/kmer.rs`**:
   - Define `pub trait KmerInt: Copy + Eq + Ord + std::hash::Hash + Send + Sync + 'static` with implementations for `u64` ($k \le 31$), `u128` ($33 \le k \le 63$), and `[u64; 4]` ($65 \le k \le 127$).
   - Implement O(1) rolling forward/reverse-complement updates (`roll_next(base: u8, k: usize)`) and `canonical(k: usize) -> (Self, bool)`.
   - Implement `ntHash` rolling 64-bit hash and sliding window minimizer extractor ($m = 15$).
3. **`src/common/sketch.rs`**:
   - Implement `HyperLogLog` (14-bit precision, 16,384 registers) for unique $k$-mer cardinality estimation.
   - Implement a 4-row `CountMinSketch` using 4-bit saturating counters (nibbles packed in `Vec<AtomicU8>`) for lock-free multi-threaded $k$-mer frequency counting up to 15.

#### Phase 1 Verification Gate:
- Write unit tests in `src/common/` verifying:
  1. Packing and unpacking random DNA sequences of lengths $1\dots 500\text{ bp}$ matches original strings and their reverse complements.
  2. Rolling canonical $k$-mer extraction for $k=21, 55, 99, 127$ yields identical canonical representations for a sequence and its reverse complement.
  3. `HyperLogLog` estimates unique 31-mers within $\pm 3\%$ error on $100,000$ synthetic $k$-mers.
- Execute via queue: `./scripts/hpc_run.sh cargo test common::`

---

### Phase 2: Shovill-Style Streaming Ingest — Estimation, Subsampling & Overlap Merging
**Goal**: Replace `KMC/Mash + seqtk + Trimmomatic + FLASH` with a single streaming pass over paired-end FASTQ files.

#### Tasks:
1. **`src/ingest/estimator.rs`**:
   - Stream the first $250,000$ read pairs using `needletail`.
   - Feed $k=21$ canonical $k$-mers into `HyperLogLog` and a frequency histogram to estimate:
     - Genome size $\hat{G}$ (excluding singleton error $k$-mers).
     - Total sequencing depth $\hat{D} = \frac{N_{\text{total\_bases}}}{\hat{G}}$.
     - Subsampling keep probability $p_{\text{keep}} = \min\left(1.0, \frac{\text{target\_cov}}{\hat{D}}\right)$ where default `target_cov = 120.0`.
   - Dynamically compute the odd multi-$k$ schedule based on median read length $L_{\text{med}}$:
     - For standard $150\text{ bp}$ PE reads (which reach $220\text{–}260\text{ bp}$ after overlap merging): select $[31, 55, 79, 103, 127]$ capped at $\text{odd}(\lfloor 0.85 \times L_{\text{effective}} \rfloor)$.
2. **`src/ingest/trim_merge.rs`**:
   - **3' Adapter & Low-Quality Trimming**: Trim trailing bases with Phred quality $< 20$ via sliding window (window size 4), and clip standard Illumina TruSeq/Nextera adapter prefixes at the 3' end.
   - **FLASH-style Paired-End Overlap Merger**:
     - Test overlaps $d \in [15, L_1]$ between the 3' suffix of R1 and the 5' prefix of $\text{RC}(\text{R2})$.
     - Compute mismatch density $m / d$ using XOR + `popcnt` bit-parallel comparisons.
     - If the best overlap has $m / d \le 0.10$, merge R1 and $\text{RC}(\text{R2})$ into a single long read in `ReadBank::merged_reads` (at mismatch positions, pick the base with higher Phred score).
     - Otherwise, store the pair in `ReadBank::pe_forward` and `ReadBank::pe_reverse`.

#### Phase 2 Verification Gate:
- Write unit tests verifying:
  1. Synthetic $2 \times 150\text{ bp}$ read pairs simulated from a $230\text{ bp}$ fragment with 2 sequencing errors are accurately merged into a single $230\text{ bp}$ consensus read.
  2. Read pairs from a $450\text{ bp}$ fragment remain unmerged in `pe_forward` / `pe_reverse`.
- Execute via queue: `./scripts/hpc_run.sh cargo test ingest::`

---

### Phase 3: Online Error Filtering & Minimizer-Bucketed Compacted DBG
**Goal**: Replace `Lighter` + `BayesHammer` and build the compacted De Bruijn graph (cDBG) directly into maximal unitigs.

#### Tasks:
1. **`src/graph/build.rs`**:
   - **Pass 1 (Parallel Count-Min Sketch Population)**:
     - Iterate over `ReadBank` in parallel via `rayon`, populating the lock-free `CountMinSketch` with canonical $k$-mers.
     - Compute the adaptive solid $k$-mer threshold $T_{\text{solid}} = \max\left(3, \lfloor 0.05 \times \text{median\_cov} \rfloor\right)$.
   - **Pass 2 (256-Shard Minimizer Partitioning & Unitig Compaction)**:
     - Extract all solid $k$-mers (count $\ge T_{\text{solid}}$) along with their left/right 4-bit extension flags (`u8`: upper 4 bits = incoming `A,C,G,T`, lower 4 bits = outgoing `A,C,G,T`).
     - Partition solid $k$-mers into 256 shards based on the top 8 bits of their minimizer hash (`FxHashMap<Kmer, (u16, u8)>` per shard: count + adjacency bitmask).
     - **Compaction**: Identify all branching/terminal $k$-mers (in-degree $\ne 1$ or out-degree $\ne 1$) as unitig endpoints. Walk unbranched chains of $1\text{-in}, 1\text{-out}$ $k$-mers in parallel to emit maximal `UnitigNode`s and populate `edges_fwd` / `edges_rev` with `overlap_k = k - 1`.
     - Detect isolated circular loops ($1\text{-in}, 1\text{-out}$ components with no branching nodes, e.g., small plasmids) so they are not skipped during traversal.

#### Phase 3 Verification Gate:
- Write a test that constructs a cDBG from simulated reads (with $1\%$ random sequencing errors) of a $10\text{ kb}$ circular sequence containing one $200\text{ bp}$ exact repeat:
  - Verify that singleton error $k$-mers are filtered out and the graph compacts into exactly 3 unitigs (the repeat unitig with in-degree 2 / out-degree 2, plus 2 flanking unique unitigs).
- Execute via queue: `./scripts/hpc_run.sh cargo test graph::build`

---

### Phase 4: SPAdes-Style Graph Simplification & Local Multi-$k$ Bridging
**Goal**: Clean graph artifacts while preserving high-copy plasmids, and iteratively bridge local repeats using increasing $k$-mer sizes without rebuilding the entire genome graph.

#### Tasks:
1. **`src/graph/simplify.rs`**:
   - Implement iterative priority-queue graph cleaning (repeat until convergence):
     1. **Tip Clipping**: Identify dead-end unitigs (out-degree $0$ on one end, length $< 2k$). Remove if an alternative outgoing edge exists at the branching node with coverage $\ge 4 \times \text{cov}(\text{tip})$. Re-compact unbranched paths immediately after removal.
     2. **Bubble / Bulge Popping**: For any branching node, perform bounded Dijkstra/BFS (up to length $\max(3k, 500\text{ bp})$) to find alternative paths $P_1, P_2$ converging at the same sink node. Compute sequence edit distance via Myers bit-parallel algorithm. If similarity $\ge 90\%$, remove the lower-coverage path, transfer its coverage to the higher-coverage path, and re-compact.
     3. **Local Relative Chimeric Edge Removal (Plasmid-Safe)**: Never delete an entire connected component solely because its coverage is $10\times\text{–}100\times$ higher than the chromosome. Only sever an internal edge $e=(u, v)$ if both $u$ and $v$ have alternative strong edges whose coverage is $> 10\times \text{cov}(e)$.
2. **`src/graph/local_multik.rs`**:
   - Given the simplified graph $G(k_{i-1})$ and next $k_i \in [55, 79, 103, 127]$:
     - **Freeze Long Linear Regions**: Unitigs longer than $1,000\text{ bp}$ get their internal core frozen.
     - **Collect Junction Reads & Flanking Seeds**: Index the $500\text{ bp}$ terminal regions of all branching unitigs in $G(k_{i-1})$. Collect all `merged_reads` and paired reads from `ReadBank` that share a minimizer with any branching terminal, plus the terminal unitig sequences themselves.
     - **Local Subgraph Upgrade**: Build a local $k_i$-cDBG solely from this subset of reads and flanking unitig ends, run `simplify`, and splice the resolved $k_i$ bridges back into the global graph.

#### Phase 4 Verification Gate:
- Write a test with a synthetic $20\text{ kb}$ sequence containing a $65\text{ bp}$ repeat (unresolvable at $k=31$, resolvable at $k=79$) and a $5\text{ kb}$ high-copy plasmid ($25\times$ relative coverage):
  - Verify that after running `simplify` + `local_multik` ($31 \to 79$), the $20\text{ kb}$ sequence assembles into a single contig and the $5\text{ kb}$ plasmid is preserved intact.
- Execute via queue: `./scripts/hpc_run.sh cargo test graph::`

---

### Phase 5: `exSPAnder` Paired-End Repeat Resolution
**Goal**: Resolve repeats longer than $k_{\max}$ (such as $800\text{–}1,500\text{ bp}$ bacterial IS elements) using paired-end insert size constraints.

#### Tasks:
1. **`src/resolve/mapper.rs`**:
   - Index unique $k_1$-mers within the first/last $1,000\text{ bp}$ of every unitig into a fast lookup table `FxHashMap<u64, (DirUnitigId, u32)>`.
   - Map `pe_forward` and `pe_reverse` in parallel via $O(1)$ terminal $k$-mer lookups:
     - From pairs mapping to the same long unitig ($\ge 2,000\text{ bp}$), estimate the empirical insert size distribution $\Psi(d)$ (median $\mu$, MAD-based standard deviation $\sigma$, valid support window $[d_{\min}, d_{\max}] = [\mu - 3\sigma, \mu + 3\sigma]$).
     - From pairs mapping to different unitigs, build a sparse `PairLinkMatrix` recording observed read-pair links $(u_1, u_2, \text{dist\_offset}) \to \text{count}$.
2. **`src/resolve/exspander.rs`**:
   - Sort unique long seed unitigs (length $\ge 500\text{ bp}$, coverage within $[0.5\times, 1.5\times]$ local median) by length descending.
   - Maintain active path $P = (e_1, e_2, \dots, e_m)$. At the 3' end of $e_m$, evaluate all candidate outgoing edges $c_j \in \text{Out}(e_m)$:
     - Compute confidence support score:
       $$\text{Score}_P(c_j) = \frac{\sum_{e_i \in P \text{ in window } [d_{\min}, d_{\max}]} \text{ObservedPairs}(e_i, c_j)}{\sum_{e_i \in P \text{ in window } [d_{\min}, d_{\max}]} \text{ExpectedPairs}(e_i, c_j \mid \Psi)}$$
     - **Winner-Take-All Decision Rule**: Extend path $P$ by $c_{\text{best}}$ if and only if:
       1. Total raw supporting read pairs $\ge \Theta_{\min}$ (default $5$), AND
       2. $\text{Score}_P(c_{\text{best}}) \ge C \cdot \text{Score}_P(c_{\text{second}})$ (default $C = 2.5$).
     - Stop extension if multiple candidates tie or loop count exceeds estimated copy number $\text{round}(\text{cov}(e) / \text{cov}_{\text{seed}})$.

#### Phase 5 Verification Gate:
- Write a test simulating a $30\text{ kb}$ chromosome containing two copies of a $350\text{ bp}$ IS-like repeat sequenced with $2 \times 150\text{ bp}$ reads at $450\text{ bp}$ insert size:
  - Verify `exspander` accurately traverses across both copies of the $350\text{ bp}$ repeat without misassembly, producing a single $30\text{ kb}$ contig.
- Execute via queue: `./scripts/hpc_run.sh cargo test resolve::`

---

### Phase 6: In-Memory Consensus Polishing, Curation & End-to-End CLI
**Goal**: Eliminate the need for `BWA + Pilon` via fast in-memory pileup polishing, trim terminal artifacts, and output production-grade FASTA and GFA files.

#### Tasks:
1. **`src/polish/consensus.rs`**:
   - Anchor all reads from `ReadBank` onto final assembled contigs using exact $31$-mer seeds spaced every $50\text{ bp}$.
   - For reads with a single unambiguous anchor (or consistent paired-end anchor), perform fast ungapped/banded alignment against the target contig window and accumulate per-position base counts `[u16; 4]` (`A, C, G, T`).
   - Replace any reference base where an alternative allele has support $> 75\%$ and depth $\ge 5$.
2. **`src/polish/output.rs`**:
   - **Shovill-style Post-Curation**:
     - Trim terminal $k$-mer self-overlaps on circular contigs and mark `circular=true` in the FASTA header (`>contig_1 len=4639675 cov=85.4 circular=true`).
     - Filter out spurious contigs with length $< \text{min\_len}$ (default $500\text{ bp}$) or coverage $< \text{min\_cov}$ (default $2.0\times$).
   - Write final contigs to `<outdir>/contigs.fa` and the simplified assembly graph to `<outdir>/assembly_graph.gfa` (GFA v1.0 `S` and `L` records).
3. **`src/cli.rs` & `src/main.rs`**:
   - Expose clean CLI flags:
     ```bash
     fastBact -1 R1.fq.gz -2 R2.fq.gz -o out_dir -t 8 --depth 120 --min-len 500 --min-cov 2.0
     ```
   - Log stage timings and peak RSS memory usage (reading `/proc/self/status` `VmHWM` on Linux) via `tracing`.

#### Phase 6 Final Integration Verification Gate:
- Implement `tests/integration_test.rs`:
  1. Generate a synthetic $100\text{ kb}$ bacterial micro-chromosome (with 3 repeats of $60\text{ bp}$, $90\text{ bp}$, and $320\text{ bp}$) plus a $4\text{ kb}$ circular plasmid at $8\times$ copy number.
  2. Simulate $80\times$ depth of $2 \times 150\text{ bp}$ paired-end reads (insert size $350 \pm 30\text{ bp}$, $0.5\%$ base error rate) into temporary `.fq.gz` files.
  3. Run the full `fastBact` pipeline end-to-end.
  4. Assert that:
     - Output contains the complete $100\text{ kb}$ chromosome contig and the $4\text{ kb}$ circular plasmid contig (`circular=true`).
     - Total execution time is $< 2$ seconds in `--release` mode.
     - Base identity against ground truth is $> 99.99\%$ with $0$ misassemblies.
- Execute via queue: `./scripts/hpc_run.sh cargo test --release`
use fastbact::graph::cdbg::DirUnitigId;
use fastbact::resolve::mapper::InsertSizeDist;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;
use std::time::Instant;

/// 16-byte revised ReadPairEvidence using u32 coordinates
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ReadPairEvidence {
    pub pair_id: u32,
    pub dist_u1_exit: u32,
    pub dist_u2_entry: u32,
    pub r1_len: u16,
    pub r2_len: u16,
}

/// Baseline scalar link matrix
pub struct BaselineLinkMatrix {
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), u32>,
}

/// Candidate A: SmallVec inline 4
pub struct EvidenceMatrixSmallVec4 {
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), SmallVec<[ReadPairEvidence; 4]>>,
}

/// Candidate B: SmallVec inline 8
pub struct EvidenceMatrixSmallVec8 {
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), SmallVec<[ReadPairEvidence; 8]>>,
}

/// Candidate C: Standard Vec (pure heap per entry)
pub struct EvidenceMatrixVec {
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), Vec<ReadPairEvidence>>,
}

fn run_microbenchmark(num_pairs: usize, avg_depth: usize, label: &str) {
    println!("\n=======================================================");
    println!(" WORKLOAD: {} ({} unique unitig pairs, depth {}x)", label, num_pairs, avg_depth);
    println!("=======================================================");

    let dist = InsertSizeDist {
        median: 450.0,
        mad: 25.0,
        min_insert: 360,
        max_insert: 540,
    };

    // 1. Baseline Benchmark
    let t0 = Instant::now();
    let mut baseline = BaselineLinkMatrix {
        links: FxHashMap::default(),
    };
    for i in 0..num_pairs {
        let u1 = DirUnitigId::new((i % 1000) as u32, false);
        let u2 = DirUnitigId::new(((i + 1) % 1000) as u32, false);
        for _ in 0..avg_depth {
            *baseline.links.entry((u1, u2)).or_insert(0) += 1;
        }
    }
    let baseline_insert_time = t0.elapsed();
    let baseline_entry_size = std::mem::size_of::<(DirUnitigId, DirUnitigId, u32)>();
    let baseline_approx_mem = baseline.links.capacity() * (baseline_entry_size + 1); // +1 control byte

    // 2. Candidate A: SmallVec 4
    let t0 = Instant::now();
    let mut matrix_sv4 = EvidenceMatrixSmallVec4 {
        links: FxHashMap::default(),
    };
    let mut sv4_heap_spills = 0usize;
    for i in 0..num_pairs {
        let u1 = DirUnitigId::new((i % 1000) as u32, false);
        let u2 = DirUnitigId::new(((i + 1) % 1000) as u32, false);
        let entry = matrix_sv4.links.entry((u1, u2)).or_insert_with(SmallVec::new);
        for d in 0..avg_depth {
            entry.push(ReadPairEvidence {
                pair_id: (i * avg_depth + d) as u32,
                dist_u1_exit: (100 + d * 5) as u32,
                dist_u2_entry: (120 + d * 5) as u32,
                r1_len: 150,
                r2_len: 150,
            });
        }
        if entry.spilled() {
            sv4_heap_spills += 1;
        }
    }
    let sv4_insert_time = t0.elapsed();
    let sv4_entry_size = std::mem::size_of::<((DirUnitigId, DirUnitigId), SmallVec<[ReadPairEvidence; 4]>)>();
    let sv4_heap_bytes: usize = matrix_sv4
        .links
        .values()
        .filter(|v| v.spilled())
        .map(|v| v.capacity() * std::mem::size_of::<ReadPairEvidence>())
        .sum();
    let sv4_total_mem = matrix_sv4.links.capacity() * (sv4_entry_size + 1) + sv4_heap_bytes;

    // 3. Candidate B: SmallVec 8
    let t0 = Instant::now();
    let mut matrix_sv8 = EvidenceMatrixSmallVec8 {
        links: FxHashMap::default(),
    };
    let mut sv8_heap_spills = 0usize;
    for i in 0..num_pairs {
        let u1 = DirUnitigId::new((i % 1000) as u32, false);
        let u2 = DirUnitigId::new(((i + 1) % 1000) as u32, false);
        let entry = matrix_sv8.links.entry((u1, u2)).or_insert_with(SmallVec::new);
        for d in 0..avg_depth {
            entry.push(ReadPairEvidence {
                pair_id: (i * avg_depth + d) as u32,
                dist_u1_exit: (100 + d * 5) as u32,
                dist_u2_entry: (120 + d * 5) as u32,
                r1_len: 150,
                r2_len: 150,
            });
        }
        if entry.spilled() {
            sv8_heap_spills += 1;
        }
    }
    let sv8_insert_time = t0.elapsed();
    let sv8_entry_size = std::mem::size_of::<((DirUnitigId, DirUnitigId), SmallVec<[ReadPairEvidence; 8]>)>();
    let sv8_heap_bytes: usize = matrix_sv8
        .links
        .values()
        .filter(|v| v.spilled())
        .map(|v| v.capacity() * std::mem::size_of::<ReadPairEvidence>())
        .sum();
    let sv8_total_mem = matrix_sv8.links.capacity() * (sv8_entry_size + 1) + sv8_heap_bytes;

    // 4. Candidate C: Standard Vec
    let t0 = Instant::now();
    let mut matrix_vec = EvidenceMatrixVec {
        links: FxHashMap::default(),
    };
    for i in 0..num_pairs {
        let u1 = DirUnitigId::new((i % 1000) as u32, false);
        let u2 = DirUnitigId::new(((i + 1) % 1000) as u32, false);
        let entry = matrix_vec.links.entry((u1, u2)).or_insert_with(Vec::new);
        for d in 0..avg_depth {
            entry.push(ReadPairEvidence {
                pair_id: (i * avg_depth + d) as u32,
                dist_u1_exit: (100 + d * 5) as u32,
                dist_u2_entry: (120 + d * 5) as u32,
                r1_len: 150,
                r2_len: 150,
            });
        }
    }
    let vec_insert_time = t0.elapsed();
    let vec_entry_size = std::mem::size_of::<((DirUnitigId, DirUnitigId), Vec<ReadPairEvidence>)>();
    let vec_heap_bytes: usize = matrix_vec
        .links
        .values()
        .map(|v| v.capacity() * std::mem::size_of::<ReadPairEvidence>())
        .sum();
    let vec_total_mem = matrix_vec.links.capacity() * (vec_entry_size + 1) + vec_heap_bytes;

    // 5. Query Throughput Benchmark: evaluate 100,000 path span queries
    let num_queries = 100_000usize;
    let path_sep = 200i64; // repeat gap

    let t0 = Instant::now();
    let mut valid_count_sv4 = 0u32;
    for q in 0..num_queries {
        let u1 = DirUnitigId::new((q % 1000) as u32, false);
        let u2 = DirUnitigId::new(((q + 1) % 1000) as u32, false);
        if let Some(list) = matrix_sv4.links.get(&(u1, u2)) {
            for ev in list {
                let span = ev.dist_u1_exit as i64 + path_sep + ev.dist_u2_entry as i64;
                if span >= dist.min_insert as i64 && span <= dist.max_insert as i64 {
                    valid_count_sv4 += 1;
                }
            }
        }
    }
    let query_time_sv4 = t0.elapsed();

    println!("--- Memory Model & Consumption ---");
    println!("  Baseline (u32 count):      {:>8.2} KB  (entry size: {} B)", baseline_approx_mem as f64 / 1024.0, baseline_entry_size);
    println!("  SmallVec<[Ev; 4]>:         {:>8.2} KB  (entry size: {} B, heap spilled pairs: {}/{})", sv4_total_mem as f64 / 1024.0, sv4_entry_size, sv4_heap_spills, num_pairs);
    println!("  SmallVec<[Ev; 8]>:         {:>8.2} KB  (entry size: {} B, heap spilled pairs: {}/{})", sv8_total_mem as f64 / 1024.0, sv8_entry_size, sv8_heap_spills, num_pairs);
    println!("  Standard Vec<Ev>:          {:>8.2} KB  (entry size: {} B, heap buffers: {})", vec_total_mem as f64 / 1024.0, vec_entry_size, num_pairs);

    println!("--- Execution Runtime ---");
    println!("  Baseline Insert Time:      {:>8.3} ms", baseline_insert_time.as_secs_f64() * 1000.0);
    println!("  SmallVec<[Ev; 4]> Insert:  {:>8.3} ms", sv4_insert_time.as_secs_f64() * 1000.0);
    println!("  SmallVec<[Ev; 8]> Insert:  {:>8.3} ms", sv8_insert_time.as_secs_f64() * 1000.0);
    println!("  Vec<Ev> Insert:            {:>8.3} ms", vec_insert_time.as_secs_f64() * 1000.0);

    println!("--- Query Performance ---");
    println!("  SmallVec<[Ev; 4]> Queries: {} queries in {:.3} ms ({:.1} ns/query, valid_hits: {})",
        num_queries,
        query_time_sv4.as_secs_f64() * 1000.0,
        query_time_sv4.as_nanos() as f64 / num_queries as f64,
        valid_count_sv4
    );
}

fn main() {
    println!("fastBact P1 ReadPairEvidence Microbenchmark on Compute Node");
    run_microbenchmark(500, 8, "Small Isolate (500 links, 8x depth)");
    run_microbenchmark(3000, 20, "Standard Bacterial (3,000 links, 20x depth)");
    run_microbenchmark(20000, 30, "High-Complexity Metagenome (20,000 links, 30x depth)");
}

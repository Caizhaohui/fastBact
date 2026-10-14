//! Fast generator for 2.5 Mb synthetic bacterial genome with research-grade repeats & plasmid.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};

fn splitmix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

fn rand_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(len);
    let mut s = seed;
    for _ in 0..len {
        s = splitmix64(s);
        seq.push(bases[(s & 3) as usize]);
    }
    seq
}

fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'A',
        })
        .collect()
}

fn main() {
    let outdir = "benchmark";
    create_dir_all(outdir).unwrap();

    eprintln!("[1/3] Generating 2.5 Mb chromosome with 1.2 kb IS and 5 kb rRNA repeats...");
    let chrom_len = 2_500_000usize;
    let mut chrom = rand_dna(chrom_len, 42);

    // 1.2 kb IS element repeat (6 identical copies)
    let is_elem = rand_dna(1_200, 101);
    let is_positions = [200_000, 600_000, 1_100_000, 1_500_000, 1_900_000, 2_300_000];
    for &pos in &is_positions {
        chrom[pos..pos + is_elem.len()].copy_from_slice(&is_elem);
    }

    // 5.0 kb rRNA operon repeat (4 identical copies)
    let rrna_operon = rand_dna(5_000, 202);
    let rrna_positions = [400_000, 900_000, 1_600_000, 2_100_000];
    for &pos in &rrna_positions {
        chrom[pos..pos + rrna_operon.len()].copy_from_slice(&rrna_operon);
    }

    // 20 kb high-copy circular plasmid
    let plasmid_len = 20_000usize;
    let plasmid = rand_dna(plasmid_len, 303);

    // Write reference FASTA
    let ref_path = format!("{}/ref_genome.fa", outdir);
    let mut ref_file = BufWriter::new(File::create(&ref_path).unwrap());
    writeln!(ref_file, ">chromosome len={}", chrom_len).unwrap();
    for chunk in chrom.chunks(80) {
        ref_file.write_all(chunk).unwrap();
        ref_file.write_all(b"\n").unwrap();
    }
    writeln!(ref_file, ">plasmid_20kb len={} circular=true", plasmid_len).unwrap();
    for chunk in plasmid.chunks(80) {
        ref_file.write_all(chunk).unwrap();
        ref_file.write_all(b"\n").unwrap();
    }
    ref_file.flush().unwrap();
    eprintln!("  Reference saved to: {}", ref_path);

    // 2. Generate Paired-End 150 bp Reads
    eprintln!("[2/3] Simulating 50x depth on chromosome and 150x depth on plasmid...");
    let read_len = 150usize;
    let chrom_pairs = (chrom_len * 50) / (read_len * 2); // 416,666 pairs
    let plasmid_pairs = (plasmid_len * 150) / (read_len * 2); // 10,000 pairs
    let total_pairs = chrom_pairs + plasmid_pairs;

    let r1_path = format!("{}/reads_R1.fq.gz", outdir);
    let r2_path = format!("{}/reads_R2.fq.gz", outdir);
    let r1_file = File::create(&r1_path).unwrap();
    let r2_file = File::create(&r2_path).unwrap();

    let mut r1_out = BufWriter::new(GzEncoder::new(r1_file, Compression::fast()));
    let mut r2_out = BufWriter::new(GzEncoder::new(r2_file, Compression::fast()));

    let qual_str = vec![b'I'; read_len]; // Q40 high quality
    let mut rng_state = 123456789u64;

    for i in 0..total_pairs {
        let (is_plasmid, ref_seq, ref_len) = if i < chrom_pairs {
            (false, &chrom, chrom_len)
        } else {
            (true, &plasmid, plasmid_len)
        };

        rng_state = splitmix64(rng_state);
        let insert_var = ((rng_state % 81) as i64) - 40; // 360..440
        let insert_size = (400 + insert_var) as usize;

        rng_state = splitmix64(rng_state);
        let max_start = if is_plasmid { ref_len } else { ref_len.saturating_sub(insert_size) };
        if max_start == 0 {
            continue;
        }
        let start = (rng_state as usize) % max_start;

        // Extract fragment
        let mut fragment = Vec::with_capacity(insert_size);
        if is_plasmid {
            // Circular wrapping
            for offset in 0..insert_size {
                fragment.push(ref_seq[(start + offset) % ref_len]);
            }
        } else {
            fragment.extend_from_slice(&ref_seq[start..start + insert_size]);
        }

        // R1: 5' prefix of fragment
        let mut r1_seq = fragment[..read_len].to_vec();
        // R2: RC of 3' suffix of fragment
        let mut r2_seq = rev_comp(&fragment[insert_size - read_len..]);

        // 0.5% substitution error
        for pos in 0..read_len {
            rng_state = splitmix64(rng_state);
            if (rng_state % 200) == 0 {
                let bases = [b'A', b'C', b'G', b'T'];
                r1_seq[pos] = bases[(rng_state >> 8) as usize & 3];
            }
            rng_state = splitmix64(rng_state);
            if (rng_state % 200) == 0 {
                let bases = [b'A', b'C', b'G', b'T'];
                r2_seq[pos] = bases[(rng_state >> 8) as usize & 3];
            }
        }

        // Write FASTQ records
        writeln!(r1_out, "@read_{}/1", i + 1).unwrap();
        r1_out.write_all(&r1_seq).unwrap();
        r1_out.write_all(b"\n+\n").unwrap();
        r1_out.write_all(&qual_str).unwrap();
        r1_out.write_all(b"\n").unwrap();

        writeln!(r2_out, "@read_{}/2", i + 1).unwrap();
        r2_out.write_all(&r2_seq).unwrap();
        r2_out.write_all(b"\n+\n").unwrap();
        r2_out.write_all(&qual_str).unwrap();
        r2_out.write_all(b"\n").unwrap();
    }

    r1_out.flush().unwrap();
    r2_out.flush().unwrap();

    eprintln!("[3/3] Successfully generated {} read pairs in '{}'", total_pairs, outdir);
}

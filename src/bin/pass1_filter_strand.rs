use std::collections::HashSet;
use clap::Parser;
use rayon::prelude::*;

use alphasplitter::monomer::{revcomp, build_anchor_kmers, count_anchor_hits};
use alphasplitter::io::{read_monomers_tsv, write_monomers_tsv};

#[derive(Parser)]
#[command(name = "pass1_filter_strand", about = "E1 Pass 1: Filter alpha monomers and strand-normalize per array")]
struct Args {
    /// Input ArraySplitter .monomers.tsv
    input: String,

    /// Output TSV path
    #[arg(short, long, default_value = "pass1_strand_normalized.tsv")]
    output: String,

    /// Output stats JSON path
    #[arg(long, default_value = "pass1_stats.json")]
    stats: String,

    /// Min monomer length
    #[arg(long, default_value_t = 140)]
    min_length: u32,

    /// Max monomer length
    #[arg(long, default_value_t = 210)]
    max_length: u32,

    /// K-mer size for strand anchors
    #[arg(long, default_value_t = 6)]
    anchor_k: usize,

    /// Number of anchor k-mers
    #[arg(long, default_value_t = 50)]
    anchor_n: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    // --- Step 1: Read and filter ---
    eprintln!("Reading {}...", args.input);
    let (mut arrays, mut stats) = read_monomers_tsv(&args.input, args.min_length, args.max_length);
    eprintln!("Filtered: {} alpha monomers in {} arrays", stats.n_monomers, stats.n_arrays);
    eprintln!(
        "  (removed: {} non-base_monomer, {} non-alpha period, {} length outliers)",
        stats.filtered_type, stats.filtered_period, stats.filtered_length
    );

    if stats.n_monomers == 0 {
        eprintln!("ERROR: No monomers passed filters!");
        std::process::exit(1);
    }

    // --- Step 2: Find seed consensus (medoid by k-mer overlap) ---
    eprintln!("Finding seed consensus...");
    let seed = find_seed_medoid(&arrays, 2000);
    eprintln!("  Seed length: {} bp", seed.len());
    stats.seed_length = seed.len();

    // Build anchor k-mers from seed
    let anchor_fwd = build_anchor_kmers(seed.as_bytes(), args.anchor_k, args.anchor_n);
    let anchor_rc: HashSet<u64> = anchor_fwd.iter().map(|&h| {
        // Compute revcomp of k-mer hash
        let mut rc_hash: u64 = 0;
        let mut fwd = h;
        for _ in 0..args.anchor_k {
            let base = fwd & 3;
            let comp = 3 - base; // A<->T, C<->G
            rc_hash = (rc_hash << 2) | comp;
            fwd >>= 2;
        }
        rc_hash
    }).collect();

    eprintln!("  Anchor k-mers: {} fwd, {} rc", anchor_fwd.len(), anchor_rc.len());

    // --- Step 3: Strand normalize per array (parallel) ---
    eprintln!("Strand normalizing per array...");

    let array_ids: Vec<String> = arrays.keys().cloned().collect();
    let flip_decisions: Vec<(String, bool)> = array_ids.par_iter().map(|array_id| {
        let monomers = &arrays[array_id];
        let mut fwd_score: u64 = 0;
        let mut rc_score: u64 = 0;

        for m in monomers {
            let seq = m.sequence.as_bytes();
            fwd_score += count_anchor_hits(seq, &anchor_fwd, args.anchor_k) as u64;
            rc_score += count_anchor_hits(seq, &anchor_rc, args.anchor_k) as u64;
        }

        (array_id.clone(), rc_score > fwd_score)
    }).collect();

    let mut n_flipped: u64 = 0;
    let mut n_kept: u64 = 0;

    for (array_id, flip) in &flip_decisions {
        if *flip {
            n_flipped += 1;
            let monomers = arrays.get_mut(array_id).unwrap();
            for m in monomers.iter_mut() {
                m.sequence = String::from_utf8(revcomp(m.sequence.as_bytes())).unwrap();
                m.strand_flipped = true;
            }
        } else {
            n_kept += 1;
        }
    }

    stats.arrays_kept_fwd = n_kept;
    stats.arrays_flipped = n_flipped;
    eprintln!("  Arrays kept forward: {}", n_kept);
    eprintln!("  Arrays flipped to revcomp: {}", n_flipped);

    // --- Step 4: Write output ---
    eprintln!("Writing {}...", args.output);
    write_monomers_tsv(&args.output, &arrays);

    // --- Step 5: Write stats ---
    let stats_json = serde_json::to_string_pretty(&stats).unwrap();
    std::fs::write(&args.stats, stats_json).unwrap();

    eprintln!("Done. Stats written to {}", args.stats);
}

fn find_seed_medoid(arrays: &std::collections::HashMap<String, Vec<alphasplitter::monomer::Monomer>>, sample_size: usize) -> String {
    

    // Collect all sequences
    let all_seqs: Vec<&str> = arrays.values()
        .flat_map(|v| v.iter().map(|m| m.sequence.as_str()))
        .collect();

    // Sample
    let sample: Vec<&str> = if all_seqs.len() <= sample_size {
        all_seqs
    } else {
        
        
        let mut selected = Vec::with_capacity(sample_size);
        let step = all_seqs.len() / sample_size;
        for i in (0..all_seqs.len()).step_by(step.max(1)).take(sample_size) {
            selected.push(all_seqs[i]);
        }
        selected
    };

    let k = 6usize;

    // Build k-mer sets for each sampled monomer
    let kmer_sets: Vec<HashSet<u64>> = sample.par_iter().map(|seq| {
        let bytes = seq.as_bytes();
        let mut set = HashSet::new();
        for i in 0..bytes.len().saturating_sub(k - 1) {
            if let Some(h) = alphasplitter::monomer::kmer_hash(&bytes[i..i + k], k) {
                set.insert(h);
            }
        }
        set
    }).collect();

    // Find medoid: highest total Jaccard
    let check_size = sample.len().min(300);
    let check_step = sample.len() / check_size.max(1);

    let best_idx = (0..sample.len()).into_par_iter().map(|i| {
        let mut total: f64 = 0.0;
        for j in (0..sample.len()).step_by(check_step.max(1)).take(check_size) {
            if i == j { continue; }
            let inter = kmer_sets[i].intersection(&kmer_sets[j]).count();
            let union = kmer_sets[i].len() + kmer_sets[j].len() - inter;
            if union > 0 {
                total += inter as f64 / union as f64;
            }
        }
        (i, total)
    }).max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
    .map(|(i, _)| i)
    .unwrap_or(0);

    sample[best_idx].to_string()
}

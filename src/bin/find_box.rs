use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;

#[derive(Parser)]
#[command(name = "find_box", about = "Search for TIGD4/CENP-B box pattern in satellite arrays")]
struct Args {
    input: String,
    #[arg(long, default_value = "TTNNGGNNANGGNNGG")]
    pattern: String,
    #[arg(long, default_value_t = 7)]
    min_score: u32,
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

fn main() {
    let args = Args::parse();
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    }

    // Parse pattern: uppercase = fixed, N = variable
    let pattern: Vec<(usize, u8)> = args.pattern.bytes().enumerate()
        .filter(|(_, b)| *b != b'N' && *b != b'n')
        .map(|(i, b)| (i, b.to_ascii_uppercase()))
        .collect();
    let plen = args.pattern.len();
    let n_fixed = pattern.len();
    eprintln!("Pattern: {} ({}bp, {} fixed positions)", args.pattern, plen, n_fixed);

    // Read ALL arrays
    let arrays = read_all_arrays(&args.input);
    let total_bp: usize = arrays.iter().map(|(_, s)| s.len()).sum();
    eprintln!("{} arrays, {:.1}Mb", arrays.len(), total_bp as f64 / 1e6);

    // Search on both strands
    let results: Vec<(String, usize, Vec<(usize, u32, char)>)> = arrays.par_iter().map(|(name, seq)| {
        let mut hits = Vec::new();
        let rc = revcomp(seq);

        for (strand_seq, strand) in [(&seq[..], '+'), (&rc[..], '-')] {
            if strand_seq.len() < plen { continue; }
            for i in 0..=(strand_seq.len() - plen) {
                let score: u32 = pattern.iter()
                    .filter(|&&(pos, base)| strand_seq.as_bytes()[i + pos].to_ascii_uppercase() == base)
                    .count() as u32;
                if score >= args.min_score {
                    let real_pos = if strand == '+' { i } else { seq.len() - i - plen };
                    hits.push((real_pos, score, strand));
                }
            }
        }
        hits.sort_by_key(|h| h.0);
        // dedup
        let mut deduped = Vec::new();
        for h in &hits {
            if deduped.last().map_or(true, |last: &(usize, u32, char)| h.0.abs_diff(last.0) > 10) {
                deduped.push(*h);
            }
        }
        (name.clone(), seq.len(), deduped)
    }).collect();

    // Aggregate
    let total_hits: usize = results.iter().map(|(_, _, h)| h.len()).sum();
    let total_monomers = total_bp as f64 / 717.0;
    eprintln!("\nTotal hits (score>={}): {}", args.min_score, total_hits);
    eprintln!("Per monomer: {:.2}", total_hits as f64 / total_monomers);

    // Spacing in largest array
    let (_, _, ref biggest_hits) = results.iter().max_by_key(|(_, l, _)| *l).unwrap();
    if biggest_hits.len() > 5 {
        let mut dist_counts: HashMap<usize, usize> = HashMap::new();
        for i in 0..biggest_hits.len() - 1 {
            let d = biggest_hits[i + 1].0.abs_diff(biggest_hits[i].0);
            let bin = d / 50 * 50;
            *dist_counts.entry(bin).or_insert(0) += 1;
        }
        eprintln!("\nSpacing distribution (largest array, {} hits):", biggest_hits.len());
        let mut sorted: Vec<_> = dist_counts.into_iter().collect();
        sorted.sort_by_key(|&(d, _)| d);
        for (d, cnt) in &sorted {
            if *cnt >= 5 {
                let bar: String = "#".repeat((*cnt).min(50));
                eprintln!("  {}bp: {} {}", d, cnt, bar);
            }
        }
    }

    // Score distribution
    eprintln!("\nScore distribution:");
    let mut score_counts = vec![0usize; n_fixed + 1];
    for (_, seq, _) in &results {
        // Recalculate for score dist
    }

    // Strand distribution
    let fwd_hits: usize = results.iter().flat_map(|(_, _, h)| h.iter()).filter(|h| h.2 == '+').count();
    let rc_hits: usize = results.iter().flat_map(|(_, _, h)| h.iter()).filter(|h| h.2 == '-').count();
    eprintln!("\nStrand: fwd={} rc={}", fwd_hits, rc_hits);
}

fn revcomp(seq: &str) -> String {
    seq.bytes().rev().map(|b| match b {
        b'A' | b'a' => b'T', b'T' | b't' => b'A',
        b'C' | b'c' => b'G', b'G' | b'g' => b'C', o => o,
    } as char).collect()
}

fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T', b'T' | b't' => b'A',
        b'C' | b'c' => b'G', b'G' | b'g' => b'C', o => o,
    }).collect()
}

fn read_all_arrays(path: &str) -> Vec<(String, String)> {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(path).unwrap();
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);
    let mut arrays = Vec::new();
    let mut name = String::new();
    let mut seq = String::new();
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            if !name.is_empty() && !seq.is_empty() {
                arrays.push((name.clone(), seq.clone()));
            }
            name = line[1..].trim().to_string();
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() && !seq.is_empty() {
        arrays.push((name, seq));
    }
    arrays
}

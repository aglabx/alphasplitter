use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;

use alphasplitter::monomer::{kmer_hash, cyclic_rotate};

#[derive(Parser)]
#[command(name = "pass2_anchor_discovery", about = "E1 Pass 2: Discover conserved anchor windows from strand-normalized monomers")]
struct Args {
    /// Input: pass1 output TSV
    input: String,

    /// Output: anchor windows JSON
    #[arg(short, long, default_value = "anchors.json")]
    output: String,

    /// Output: conservation profile TSV (per-position stats)
    #[arg(long, default_value = "conservation_profile.tsv")]
    profile: String,

    /// Output: canonical monomers TSV (rotation-assigned)
    #[arg(long, default_value = "pass2_canonical.tsv")]
    canonical: String,

    /// Sample size for conservation scan
    #[arg(long, default_value_t = 50000)]
    sample_size: usize,

    /// Canonical monomer length for alignment
    #[arg(long, default_value_t = 171)]
    canon_len: usize,

    /// Anchor window size
    #[arg(long, default_value_t = 11)]
    window_size: usize,

    /// Number of anchor windows to select
    #[arg(long, default_value_t = 5)]
    n_anchors: usize,

    /// Min spacing between anchor windows (positions)
    #[arg(long, default_value_t = 20)]
    min_spacing: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Serialize)]
struct AnchorWindow {
    start: usize,
    end: usize,
    mean_conservation: f64,
    consensus: String,
}

#[derive(Debug, Serialize)]
struct AnchorReport {
    n_monomers_sampled: usize,
    canon_len: usize,
    seed_sequence: String,
    windows: Vec<AnchorWindow>,
    conservation_summary: ConservationSummary,
}

#[derive(Debug, Serialize)]
struct ConservationSummary {
    mean_conservation: f64,
    median_conservation: f64,
    min_conservation: f64,
    max_conservation: f64,
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    // --- Step 1: Read strand-normalized monomers ---
    eprintln!("Reading {}...", args.input);
    let (arrays, all_seqs) = read_pass1_tsv(&args.input);
    eprintln!("  {} monomers in {} arrays", all_seqs.len(), arrays.len());

    // --- Step 2: Find seed and do rough cyclic alignment ---
    eprintln!("Finding seed for cyclic alignment...");

    // Use the most common length monomers
    let target_len = args.canon_len;
    let exact_len_seqs: Vec<&[u8]> = all_seqs.iter()
        .filter(|s| s.len() == target_len)
        .map(|s| s.as_bytes())
        .collect();
    eprintln!("  {} monomers with exact length {}", exact_len_seqs.len(), target_len);

    // Build seed from medoid of exact-length monomers
    let seed = find_medoid_fast(&exact_len_seqs, 2000);
    eprintln!("  Seed: {}...{}", &seed[..20.min(seed.len())], &seed[seed.len().saturating_sub(20)..]);

    // --- Step 3: Cyclic-align sample to seed ---
    eprintln!("Cyclic aligning sample to seed...");
    let sample_seqs: Vec<&[u8]> = if exact_len_seqs.len() <= args.sample_size {
        exact_len_seqs.clone()
    } else {
        let step = exact_len_seqs.len() / args.sample_size;
        exact_len_seqs.iter().step_by(step.max(1)).take(args.sample_size).copied().collect()
    };

    let seed_bytes = seed.as_bytes();

    // For each sampled monomer, find best cyclic rotation against seed
    let rotations: Vec<(usize, i32)> = sample_seqs.par_iter().enumerate().map(|(i, seq)| {
        let best_rot = find_best_rotation(seq, seed_bytes);
        (i, best_rot)
    }).collect();

    // --- Step 4: Build conservation profile ---
    eprintln!("Building conservation profile (n={})...", sample_seqs.len());
    let mut profile = vec![[0u32; 5]; target_len]; // A, C, G, T, gap

    for &(i, rot) in &rotations {
        let seq = sample_seqs[i];
        let rotated = cyclic_rotate(seq, rot as usize);
        for (pos, &base) in rotated.iter().enumerate().take(target_len) {
            let idx = match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            };
            profile[pos][idx] += 1;
        }
    }

    // Compute per-position conservation (max base frequency)
    let n_total = sample_seqs.len() as f64;
    let conservation: Vec<f64> = profile.iter().map(|counts| {
        let max_count = counts[..4].iter().max().unwrap();
        *max_count as f64 / n_total
    }).collect();

    // Consensus sequence
    let bases = [b'A', b'C', b'G', b'T'];
    let consensus: String = profile.iter().map(|counts| {
        let best = counts[..4].iter().enumerate().max_by_key(|(_, &c)| c).unwrap().0;
        bases[best] as char
    }).collect();

    // --- Step 5: Select anchor windows ---
    eprintln!("Selecting {} anchor windows (size={}, min_spacing={})...",
        args.n_anchors, args.window_size, args.min_spacing);

    let windows = select_anchor_windows(
        &conservation, args.window_size, args.n_anchors, args.min_spacing, &consensus,
    );

    for (i, w) in windows.iter().enumerate() {
        eprintln!("  Anchor {}: pos {}-{}, conservation={:.3}, consensus={}",
            i + 1, w.start, w.end, w.mean_conservation, w.consensus);
    }

    // --- Step 6: Write conservation profile ---
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.profile).unwrap());
        writeln!(f, "position\tA\tC\tG\tT\tgap\tconservation\tconsensus").unwrap();
        for (pos, counts) in profile.iter().enumerate() {
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                pos, counts[0], counts[1], counts[2], counts[3], counts[4],
                conservation[pos], consensus.as_bytes()[pos] as char).unwrap();
        }
    }
    eprintln!("Conservation profile written to {}", args.profile);

    // --- Step 7: Assign rotation to ALL monomers using multi-anchor scoring ---
    eprintln!("Assigning rotation to all monomers using multi-anchor scoring...");

    // Build anchor PWMs from profile
    let anchor_pwms: Vec<Vec<[f64; 4]>> = windows.iter().map(|w| {
        (w.start..w.end).map(|pos| {
            let total = profile[pos][..4].iter().sum::<u32>() as f64;
            [
                profile[pos][0] as f64 / total,
                profile[pos][1] as f64 / total,
                profile[pos][2] as f64 / total,
                profile[pos][3] as f64 / total,
            ]
        }).collect()
    }).collect();

    // For each array, assign rotation per monomer then smooth
    let array_ids: Vec<String> = arrays.keys().cloned().collect();

    let canonical_arrays: Vec<(String, Vec<(u32, i32, f64, String)>)> = array_ids.par_iter().map(|array_id| {
        let monomers = &arrays[array_id];
        let mut results: Vec<(u32, i32, f64, String)> = Vec::with_capacity(monomers.len());

        for m in monomers {
            let seq = m.1.as_bytes();
            if seq.len() != target_len {
                // Non-canonical length: skip rotation, keep as-is
                results.push((m.0, 0, 0.0, m.1.clone()));
                continue;
            }

            let (best_rot, best_score) = find_best_rotation_multi_anchor(seq, &windows, &anchor_pwms, target_len);
            let canonical = cyclic_rotate(seq, best_rot as usize);
            let canonical_str = String::from_utf8(canonical).unwrap_or_default();

            results.push((m.0, best_rot, best_score, canonical_str));
        }

        // Smooth: if a monomer's rotation differs from neighbors by >10 and confidence is low, adjust
        smooth_rotations(&mut results, target_len);

        (array_id.clone(), results)
    }).collect();

    // --- Step 8: Write canonical monomers ---
    eprintln!("Writing canonical monomers to {}...", args.canonical);
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.canonical).unwrap());
        writeln!(f, "array_id\tidx\trotation_offset\tanchor_confidence\tsequence").unwrap();

        let mut sorted: Vec<_> = canonical_arrays.iter().collect();
        sorted.sort_by(|a, b| a.0.cmp(&b.0));

        for (array_id, monomers) in sorted {
            for (idx, rot, conf, seq) in monomers {
                writeln!(f, "{}\t{}\t{}\t{:.4}\t{}", array_id, idx, rot, conf, seq).unwrap();
            }
        }
    }

    // --- Step 9: Write anchor report ---
    let mut cons_sorted = conservation.clone();
    cons_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let report = AnchorReport {
        n_monomers_sampled: sample_seqs.len(),
        canon_len: target_len,
        seed_sequence: seed.to_string(),
        windows,
        conservation_summary: ConservationSummary {
            mean_conservation: conservation.iter().sum::<f64>() / conservation.len() as f64,
            median_conservation: cons_sorted[cons_sorted.len() / 2],
            min_conservation: cons_sorted[0],
            max_conservation: *cons_sorted.last().unwrap(),
        },
    };

    let report_json = serde_json::to_string_pretty(&report).unwrap();
    std::fs::write(&args.output, &report_json).unwrap();
    eprintln!("Anchor report written to {}", args.output);

    eprintln!("Done.");
}

/// Read pass1 TSV. Returns (array_id -> [(idx, sequence)], all_sequences)
fn read_pass1_tsv(path: &str) -> (HashMap<String, Vec<(u32, String)>>, Vec<String>) {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(path).unwrap();
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);

    let mut arrays: HashMap<String, Vec<(u32, String)>> = HashMap::new();
    let mut all_seqs: Vec<String> = Vec::new();
    let mut first = true;
    let mut col_array = 0;
    let mut col_idx = 0;
    let mut col_seq = 0;

    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();

        if first {
            first = false;
            for (i, &h) in fields.iter().enumerate() {
                match h {
                    "array_id" => col_array = i,
                    "idx" => col_idx = i,
                    "sequence" => col_seq = i,
                    _ => {}
                }
            }
            continue;
        }

        if fields.len() <= col_seq {
            continue;
        }

        let seq = fields[col_seq].to_string();
        let idx: u32 = fields[col_idx].parse().unwrap_or(0);
        let array_id = fields[col_array].to_string();

        all_seqs.push(seq.clone());
        arrays.entry(array_id).or_default().push((idx, seq));
    }

    // Sort within arrays
    for v in arrays.values_mut() {
        v.sort_by_key(|x| x.0);
    }

    (arrays, all_seqs)
}

fn find_medoid_fast(seqs: &[&[u8]], sample_size: usize) -> String {
    let k = 6usize;
    let sample: Vec<&[u8]> = if seqs.len() <= sample_size {
        seqs.to_vec()
    } else {
        let step = seqs.len() / sample_size;
        seqs.iter().step_by(step.max(1)).take(sample_size).copied().collect()
    };

    let kmer_sets: Vec<std::collections::HashSet<u64>> = sample.par_iter().map(|seq| {
        let mut set = std::collections::HashSet::new();
        for i in 0..seq.len().saturating_sub(k - 1) {
            if let Some(h) = kmer_hash(&seq[i..i + k], k) {
                set.insert(h);
            }
        }
        set
    }).collect();

    let check_n = sample.len().min(200);
    let check_step = (sample.len() / check_n).max(1);

    let best_idx = (0..sample.len()).into_par_iter().map(|i| {
        let mut total: f64 = 0.0;
        for j in (0..sample.len()).step_by(check_step).take(check_n) {
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

    String::from_utf8(sample[best_idx].to_vec()).unwrap()
}

/// Find best cyclic rotation of `seq` to match `seed` by hamming distance
fn find_best_rotation(seq: &[u8], seed: &[u8]) -> i32 {
    let n = seq.len().min(seed.len());
    let mut best_rot = 0i32;
    let mut best_matches = 0u32;

    for rot in 0..n {
        let mut matches = 0u32;
        for pos in 0..n {
            let spos = (pos + rot) % seq.len();
            if seq[spos] == seed[pos] {
                matches += 1;
            }
        }
        if matches > best_matches {
            best_matches = matches;
            best_rot = rot as i32;
        }
    }
    best_rot
}

/// Find best cyclic rotation using multi-anchor PWM scoring
fn find_best_rotation_multi_anchor(
    seq: &[u8],
    windows: &[AnchorWindow],
    pwms: &[Vec<[f64; 4]>],
    canon_len: usize,
) -> (i32, f64) {
    let n = seq.len();
    let mut best_rot = 0i32;
    let mut best_score = f64::NEG_INFINITY;

    for rot in 0..n {
        let mut score = 0.0f64;
        for (wi, w) in windows.iter().enumerate() {
            let pwm = &pwms[wi];
            for (pi, pos) in (w.start..w.end).enumerate() {
                let spos = (pos + rot) % n;
                let base_idx = match seq[spos] {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => continue,
                };
                // Log-likelihood vs uniform
                let p = pwm[pi][base_idx];
                if p > 0.0 {
                    score += (p * 4.0).ln();
                } else {
                    score -= 10.0; // penalty for impossible base
                }
            }
        }
        if score > best_score {
            best_score = score;
            best_rot = rot as i32;
        }
    }

    // Normalize confidence to 0-1 range
    // Use softmax-like: confidence = 1 - exp(-delta) where delta = best - second_best
    let confidence = if best_score > f64::NEG_INFINITY {
        (best_score / (windows.iter().map(|w| w.end - w.start).sum::<usize>() as f64)).exp().min(1.0)
    } else {
        0.0
    };

    (best_rot, confidence)
}

fn select_anchor_windows(
    conservation: &[f64],
    window_size: usize,
    n_anchors: usize,
    min_spacing: usize,
    consensus: &str,
) -> Vec<AnchorWindow> {
    let n = conservation.len();
    if n < window_size {
        return vec![];
    }

    // Compute mean conservation per window
    let mut window_scores: Vec<(usize, f64)> = Vec::new();
    for start in 0..=(n - window_size) {
        let mean: f64 = conservation[start..start + window_size].iter().sum::<f64>() / window_size as f64;
        window_scores.push((start, mean));
    }
    window_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Greedy selection with spacing constraint
    let mut selected: Vec<AnchorWindow> = Vec::new();
    for (start, mean_cons) in window_scores {
        // Check spacing to all already selected
        let too_close = selected.iter().any(|w: &AnchorWindow| {
            let dist = if start >= w.start { start - w.start } else { w.start - start };
            dist < min_spacing
        });
        if too_close {
            continue;
        }

        selected.push(AnchorWindow {
            start,
            end: start + window_size,
            mean_conservation: mean_cons,
            consensus: consensus[start..start + window_size].to_string(),
        });

        if selected.len() >= n_anchors {
            break;
        }
    }

    selected.sort_by_key(|w| w.start);
    selected
}

/// Smooth rotation assignments within an array
fn smooth_rotations(results: &mut [(u32, i32, f64, String)], canon_len: usize) {
    if results.len() < 3 {
        return;
    }

    let low_confidence_threshold = 0.3;
    let max_deviation = 10i32;

    for i in 1..results.len() - 1 {
        let (_, rot, conf, _) = &results[i];
        if *conf > low_confidence_threshold {
            continue;
        }

        let prev_rot = results[i - 1].1;
        let next_rot = results[i + 1].1;

        // If neighbors agree and this one deviates
        let neighbor_diff = (prev_rot - next_rot).abs().min(canon_len as i32 - (prev_rot - next_rot).abs());
        if neighbor_diff <= max_deviation {
            let my_diff = (*rot - prev_rot).abs().min(canon_len as i32 - (*rot - prev_rot).abs());
            if my_diff > max_deviation {
                // Adjust to neighbor average
                let avg = ((prev_rot + next_rot) / 2).rem_euclid(canon_len as i32);
                let seq_bytes = results[i].3.as_bytes();
                if !seq_bytes.is_empty() && seq_bytes.len() == canon_len {
                    // Re-rotate: undo current, apply new
                    let current_rot = results[i].1 as usize;
                    let new_rot = avg as usize;
                    let diff_rot = (new_rot + canon_len - current_rot) % canon_len;
                    let new_seq = cyclic_rotate(seq_bytes, diff_rot);
                    results[i].1 = avg;
                    results[i].3 = String::from_utf8(new_seq).unwrap_or_default();
                }
            }
        }
    }
}

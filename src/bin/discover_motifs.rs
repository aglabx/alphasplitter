use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;

#[derive(Parser)]
#[command(name = "discover_motifs", about = "Ab initio motif discovery from satellite arrays: find conserved k-mers with periodic spacing")]
struct Args {
    /// Input: .10kb.fasta file
    input: String,

    /// Output: discovered motifs JSON
    #[arg(short, long, default_value = "motifs.json")]
    output: String,

    /// Target period (monomer length). 0 = auto-detect from most common period in fasta headers
    #[arg(long, default_value_t = 0)]
    period: usize,

    /// Period tolerance (accept periods within this range)
    #[arg(long, default_value_t = 30)]
    period_tol: usize,

    /// Motif length (k-mer size)
    #[arg(long, default_value_t = 11)]
    motif_len: usize,

    /// Number of motifs to discover
    #[arg(long, default_value_t = 5)]
    n_motifs: usize,

    /// Min spacing between motif positions within monomer
    #[arg(long, default_value_t = 20)]
    min_spacing: usize,

    /// Max arrays to sample for discovery
    #[arg(long, default_value_t = 50)]
    max_arrays: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Serialize)]
struct DiscoveredMotif {
    name: String,
    sequence: String,
    position_in_monomer: usize,
    conservation: f64,
    count_in_arrays: usize,
    periodic_score: f64,
}

#[derive(Debug, Serialize)]
struct EnrichedArray {
    name: String,
    length: usize,
    period: usize,
    n_motif_hits: usize,
    distinct_motifs: usize,
    motif_fingerprint: String,
    family: String,
    chr: String,
    strand: char,
}

#[derive(Debug, Serialize)]
struct Family {
    name: String,
    fingerprint_pattern: String,
    n_arrays: usize,
    total_bp: usize,
    periods: Vec<(usize, usize)>,
    n_chromosomes: usize,
}

#[derive(Debug, Serialize)]
struct DiscoveryReport {
    target_period: usize,
    n_arrays_used: usize,
    total_bp: usize,
    motifs: Vec<DiscoveredMotif>,
    enrichment: EnrichmentReport,
    rust_const: String,
}

#[derive(Debug, Serialize)]
struct EnrichmentReport {
    total_arrays_scanned: usize,
    arrays_with_3plus_motifs: usize,
    families: Vec<Family>,
    arrays: Vec<EnrichedArray>,
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    // --- Read arrays ---
    eprintln!("Reading {}...", args.input);
    let all_arrays = read_all_arrays(&args.input);
    eprintln!("  {} total arrays", all_arrays.len());

    // --- Detect or use target period ---
    let target_period = if args.period > 0 {
        args.period
    } else {
        // Auto-detect: most common period from headers
        let mut period_counts: HashMap<usize, usize> = HashMap::new();
        for (name, _) in &all_arrays {
            let fields: Vec<&str> = name.split('_').collect();
            if let Some(p) = fields.last().and_then(|s| s.parse::<usize>().ok()) {
                *period_counts.entry(p).or_insert(0) += 1;
            }
        }
        let mut sorted: Vec<_> = period_counts.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));
        eprintln!("  Top periods: {:?}", &sorted[..sorted.len().min(5)]);
        sorted[0].0
    };
    eprintln!("  Target period: {} bp", target_period);

    // Filter arrays matching target period
    let period_arrays: Vec<&(String, Vec<u8>)> = all_arrays.iter().filter(|(name, _)| {
        let fields: Vec<&str> = name.split('_').collect();
        if let Some(p) = fields.last().and_then(|s| s.parse::<usize>().ok()) {
            p.abs_diff(target_period) <= args.period_tol
        } else {
            false
        }
    }).collect();
    eprintln!("  {} arrays with period ~{}", period_arrays.len(), target_period);

    // Sample arrays
    let sample: Vec<&(String, Vec<u8>)> = if period_arrays.len() <= args.max_arrays {
        period_arrays
    } else {
        let step = period_arrays.len() / args.max_arrays;
        period_arrays.into_iter().step_by(step.max(1)).take(args.max_arrays).collect()
    };
    let total_bp: usize = sample.iter().map(|(_, s)| s.len()).sum();
    eprintln!("  Using {} arrays ({} bp) for motif discovery", sample.len(), total_bp);

    // --- Step 1: Extract monomers by period slicing ---
    // Simple approach: slice each array into chunks of target_period length
    eprintln!("Extracting monomer chunks (period={})...", target_period);

    let mut all_monomers: Vec<Vec<u8>> = Vec::new();
    for (_, seq) in &sample {
        let n_monomers = seq.len() / target_period;
        for i in 0..n_monomers {
            let start = i * target_period;
            let end = start + target_period;
            if end <= seq.len() {
                all_monomers.push(seq[start..end].to_vec());
            }
        }
    }
    eprintln!("  {} monomer chunks", all_monomers.len());

    if all_monomers.len() < 10 {
        eprintln!("ERROR: Too few monomers for discovery");
        std::process::exit(1);
    }

    // --- Step 2: Find seed (medoid) and rough-align all to it ---
    eprintln!("Finding seed monomer...");
    let seed_idx = find_medoid(&all_monomers, 500);
    let seed = &all_monomers[seed_idx];
    eprintln!("  Seed: {}...{} (len={})",
        String::from_utf8_lossy(&seed[..20.min(seed.len())]),
        String::from_utf8_lossy(&seed[seed.len().saturating_sub(20)..]),
        seed.len());

    // Cyclic-align all monomers to seed
    eprintln!("Cyclic-aligning {} monomers to seed...", all_monomers.len());
    let rotations: Vec<usize> = all_monomers.par_iter().map(|mono| {
        find_best_rotation(mono, seed)
    }).collect();

    // --- Step 3: Build conservation profile ---
    eprintln!("Building conservation profile...");
    let n_mono = all_monomers.len();
    let mut profile = vec![[0u32; 5]; target_period]; // A,C,G,T,other

    for (i, mono) in all_monomers.iter().enumerate() {
        let rot = rotations[i];
        for pos in 0..target_period {
            let spos = (pos + rot) % mono.len();
            if spos >= mono.len() { continue; }
            let idx = match mono[spos] {
                b'A' | b'a' => 0, b'C' | b'c' => 1,
                b'G' | b'g' => 2, b'T' | b't' => 3, _ => 4,
            };
            profile[pos][idx] += 1;
        }
    }

    let conservation: Vec<f64> = profile.iter().map(|c| {
        let total = c[..4].iter().sum::<u32>() as f64;
        if total == 0.0 { return 0.0; }
        *c[..4].iter().max().unwrap() as f64 / total
    }).collect();

    let bases = [b'A', b'C', b'G', b'T'];
    let consensus: Vec<u8> = profile.iter().map(|c| {
        bases[c[..4].iter().enumerate().max_by_key(|(_, &v)| v).unwrap().0]
    }).collect();

    let mean_cons = conservation.iter().sum::<f64>() / conservation.len() as f64;
    eprintln!("  Mean conservation: {:.3}", mean_cons);

    // --- Step 4: Select best motif windows ---
    eprintln!("Selecting {} anchor motifs (len={}, min_spacing={})...",
        args.n_motifs, args.motif_len, args.min_spacing);

    let k = args.motif_len;
    let mut window_scores: Vec<(usize, f64)> = Vec::new();
    for start in 0..=(target_period.saturating_sub(k)) {
        let mean: f64 = conservation[start..start + k].iter().sum::<f64>() / k as f64;
        window_scores.push((start, mean));
    }
    window_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let mut motifs: Vec<DiscoveredMotif> = Vec::new();
    for (start, mean_cons) in &window_scores {
        let too_close = motifs.iter().any(|m| {
            m.position_in_monomer.abs_diff(*start) < args.min_spacing
        });
        if too_close { continue; }

        let motif_seq = String::from_utf8_lossy(&consensus[*start..*start + k]).to_string();

        // Count how many arrays contain this motif (with <=2 mismatches)
        let motif_bytes = motif_seq.as_bytes();
        let arrays_with: usize = sample.iter().filter(|(_, seq)| {
            for i in 0..seq.len().saturating_sub(k) {
                let mm = hamming(&seq[i..i+k], motif_bytes);
                if mm <= 2 { return true; }
            }
            false
        }).count();

        // Periodic score: how often does this motif appear at multiples of target_period
        let periodic_score = compute_periodic_score(&sample, motif_bytes, target_period);

        motifs.push(DiscoveredMotif {
            name: format!("M{}", motifs.len() + 1),
            sequence: motif_seq,
            position_in_monomer: *start,
            conservation: *mean_cons,
            count_in_arrays: arrays_with,
            periodic_score,
        });

        if motifs.len() >= args.n_motifs { break; }
    }

    motifs.sort_by_key(|m| m.position_in_monomer);
    // Rename after sorting
    for (i, m) in motifs.iter_mut().enumerate() {
        m.name = format!("M{}", i + 1);
    }

    eprintln!("\n=== DISCOVERED MOTIFS ===");
    for m in &motifs {
        eprintln!("  {}: pos={:>4}, cons={:.3}, arrays={}/{}, periodic={:.3}, seq={}",
            m.name, m.position_in_monomer, m.conservation,
            m.count_in_arrays, sample.len(), m.periodic_score, m.sequence);
    }

    // --- Build Rust const for motif_cut ---
    let rust_lines: Vec<String> = motifs.iter().map(|m| {
        format!("    (\"{}\", \"{}\"),", m.name, m.sequence)
    }).collect();
    let rust_const = format!(
        "const ANCHORS: [(&str, &str); {}] = [\n{}\n];",
        motifs.len(), rust_lines.join("\n")
    );

    eprintln!("\n{}", rust_const);

    // --- Step 5: Enrichment — scan ALL arrays for discovered motifs ---
    eprintln!("\n=== ENRICHMENT: scanning ALL {} arrays for discovered motifs ===", all_arrays.len());

    let motif_seqs: Vec<&[u8]> = motifs.iter().map(|m| m.sequence.as_bytes()).collect();
    let motif_rcs: Vec<Vec<u8>> = motif_seqs.iter().map(|s| revcomp(s)).collect();
    let n_motifs = motifs.len();
    let k = args.motif_len;

    let enriched: Vec<EnrichedArray> = all_arrays.par_iter().map(|(name, seq)| {
        // Count hits on forward strand
        let mut fwd_hits = [0u32; 10];
        // Count hits on reverse strand
        let mut rc_hits = [0u32; 10];
        for (mi, (fwd_motif, rc_motif)) in motif_seqs.iter().zip(motif_rcs.iter()).enumerate() {
            if mi >= 10 { break; }
            for i in 0..seq.len().saturating_sub(k) {
                let window = &seq[i..i + k];
                if hamming(window, fwd_motif) <= 2 { fwd_hits[mi] += 1; }
                if hamming(window, rc_motif) <= 2 { rc_hits[mi] += 1; }
            }
        }

        // Choose best strand: use whichever has more total hits
        let fwd_total: u32 = fwd_hits[..n_motifs].iter().sum();
        let rc_total: u32 = rc_hits[..n_motifs].iter().sum();
        let (hits, strand) = if fwd_total >= rc_total {
            (fwd_hits, '+')
        } else {
            (rc_hits, '-')
        };

        let distinct = (0..n_motifs).filter(|&i| hits[i] > 0).count();
        let total_hits: u32 = hits[..n_motifs].iter().sum();
        let fingerprint: String = (0..n_motifs).map(|i| if hits[i] > 0 { '1' } else { '0' }).collect();

        // Extract chr and period from header: chr_start_end_length_period
        let fields: Vec<&str> = name.split('_').collect();
        let period: usize = fields.last().and_then(|s| s.parse().ok()).unwrap_or(0);
        let chr = if fields.len() >= 5 {
            // header is like NC_060925.1_start_end_len_period => chr = NC_060925.1
            fields[..fields.len() - 4].join("_")
        } else {
            name.clone()
        };

        EnrichedArray {
            name: name.clone(),
            length: seq.len(),
            period,
            n_motif_hits: total_hits as usize,
            distinct_motifs: distinct,
            motif_fingerprint: fingerprint.clone(),
            family: if distinct >= 3 { fingerprint } else { "unassigned".to_string() },
            chr,
            strand,
        }
    }).collect();

    let n_with_3plus = enriched.iter().filter(|a| a.distinct_motifs >= 3).count();
    eprintln!("  Arrays with >= 3 motifs: {} / {}", n_with_3plus, all_arrays.len());

    // Group into families by fingerprint
    let mut family_map: HashMap<String, Vec<&EnrichedArray>> = HashMap::new();
    for a in &enriched {
        if a.distinct_motifs >= 3 {
            family_map.entry(a.motif_fingerprint.clone()).or_default().push(a);
        }
    }

    let mut families: Vec<Family> = family_map.iter().map(|(fp, members)| {
        let total_bp: usize = members.iter().map(|a| a.length).sum();

        // Count periods
        let mut period_counts: HashMap<usize, usize> = HashMap::new();
        for a in members {
            *period_counts.entry(a.period).or_insert(0) += 1;
        }
        let mut periods: Vec<(usize, usize)> = period_counts.into_iter().collect();
        periods.sort_by(|a, b| b.1.cmp(&a.1));

        // Count chromosomes
        let mut chr_counts: HashMap<String, usize> = HashMap::new();
        for a in members {
            *chr_counts.entry(a.chr.clone()).or_insert(0) += 1;
        }
        let n_chr = chr_counts.len();

        Family {
            name: format!("F_{}", fp),
            fingerprint_pattern: fp.clone(),
            n_arrays: members.len(),
            total_bp,
            periods: periods.into_iter().take(10).collect(),
            n_chromosomes: n_chr,
        }
    }).collect();
    families.sort_by(|a, b| b.n_arrays.cmp(&a.n_arrays));

    eprintln!("\n=== FAMILIES ===");
    eprintln!("{:<20} {:>8} {:>12} {:>5}  top periods", "fingerprint", "arrays", "total_bp", "chrs");
    for f in &families {
        let periods_str: String = f.periods.iter().take(5)
            .map(|(p, c)| format!("{}x{}", p, c))
            .collect::<Vec<_>>().join(", ");
        eprintln!("{:<20} {:>8} {:>12} {:>5}  {}", f.fingerprint_pattern, f.n_arrays, f.total_bp, f.n_chromosomes, periods_str);
    }

    let enrichment = EnrichmentReport {
        total_arrays_scanned: all_arrays.len(),
        arrays_with_3plus_motifs: n_with_3plus,
        families,
        arrays: enriched,
    };

    // --- Write report ---
    let report = DiscoveryReport {
        target_period,
        n_arrays_used: sample.len(),
        total_bp,
        motifs,
        enrichment,
        rust_const,
    };
    std::fs::write(&args.output, serde_json::to_string_pretty(&report).unwrap()).unwrap();
    eprintln!("\nWritten to {}", args.output);
}

fn compute_periodic_score(arrays: &[&(String, Vec<u8>)], motif: &[u8], period: usize) -> f64 {
    let k = motif.len();
    let mut periodic_hits = 0u64;
    let mut total_hits = 0u64;

    for (_, seq) in arrays {
        let mut hit_positions: Vec<usize> = Vec::new();
        for i in 0..seq.len().saturating_sub(k) {
            if hamming(&seq[i..i+k], motif) <= 2 {
                hit_positions.push(i);
            }
        }
        // Dedup nearby hits
        let mut deduped: Vec<usize> = Vec::new();
        for &p in &hit_positions {
            if deduped.is_empty() || p - *deduped.last().unwrap() > k {
                deduped.push(p);
            }
        }

        total_hits += deduped.len() as u64;
        for i in 0..deduped.len().saturating_sub(1) {
            let spacing = deduped[i + 1] - deduped[i];
            // Is spacing close to a multiple of period?
            let remainder = spacing % period;
            let deviation = remainder.min(period - remainder);
            if deviation <= period / 10 { // within 10%
                periodic_hits += 1;
            }
        }
    }

    if total_hits <= 1 { return 0.0; }
    periodic_hits as f64 / (total_hits - 1) as f64
}

fn find_medoid(seqs: &[Vec<u8>], sample_size: usize) -> usize {
    let n = seqs.len().min(sample_size);
    let step = (seqs.len() / n).max(1);
    let sample: Vec<usize> = (0..seqs.len()).step_by(step).take(n).collect();

    let check_n = n.min(100);
    let check_step = (n / check_n).max(1);

    let best = sample.par_iter().map(|&i| {
        let mut total = 0u64;
        for &j in sample.iter().step_by(check_step).take(check_n) {
            if i == j { continue; }
            total += hamming(&seqs[i], &seqs[j]) as u64;
        }
        (i, total)
    }).min_by_key(|&(_, t)| t).unwrap().0;

    best
}

fn find_best_rotation(seq: &[u8], seed: &[u8]) -> usize {
    let n = seq.len().min(seed.len());
    let mut best_rot = 0;
    let mut best_matches = 0u32;
    for rot in 0..n {
        let mut matches = 0u32;
        for pos in 0..n {
            let spos = (pos + rot) % seq.len();
            if seq[spos] == seed[pos] { matches += 1; }
        }
        if matches > best_matches { best_matches = matches; best_rot = rot; }
    }
    best_rot
}

fn hamming(a: &[u8], b: &[u8]) -> u32 {
    a.iter().zip(b.iter()).filter(|(x, y)| x.to_ascii_uppercase() != y.to_ascii_uppercase()).count() as u32
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T', b'T' | b't' => b'A',
        b'C' | b'c' => b'G', b'G' | b'g' => b'C', o => o,
    }).collect()
}

fn read_all_arrays(path: &str) -> Vec<(String, Vec<u8>)> {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(path).unwrap();
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);
    let mut arrays = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            if !name.is_empty() && !seq.is_empty() {
                arrays.push((name.clone(), seq.clone()));
            }
            name = line[1..].trim().to_string();
            seq.clear();
        } else {
            seq.extend(line.trim().as_bytes());
        }
    }
    if !name.is_empty() && !seq.is_empty() {
        arrays.push((name, seq));
    }
    arrays
}

use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;

#[derive(Parser)]
#[command(name = "pass3_cluster_alphabet", about = "E1 Pass 3: Cluster canonical monomers into alphabet letters")]
struct Args {
    /// Input: pass2 canonical monomers TSV
    input: String,

    /// Output: alphabet assignment TSV
    #[arg(short, long, default_value = "alphabet.tsv")]
    output: String,

    /// Output: cluster report JSON
    #[arg(long, default_value = "alphabet_report.json")]
    report: String,

    /// Output: consensus sequences FASTA
    #[arg(long, default_value = "alphabet_consensus.fa")]
    consensus_fa: String,

    /// Canonical monomer length to cluster
    #[arg(long, default_value_t = 171)]
    canon_len: usize,

    /// Subsample size for distance matrix clustering
    #[arg(long, default_value_t = 10000)]
    subsample: usize,

    /// Target number of clusters (0 = auto)
    #[arg(long, default_value_t = 0)]
    target_k: usize,

    /// Max k to try for auto selection
    #[arg(long, default_value_t = 60)]
    max_k: usize,

    /// Min k to try
    #[arg(long, default_value_t = 10)]
    min_k: usize,

    /// Number of k-medoids iterations
    #[arg(long, default_value_t = 50)]
    iterations: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Clone)]
struct MonomerEntry {
    array_id: String,
    idx: u32,
    sequence: Vec<u8>,
}

#[derive(Debug, Serialize)]
struct ClusterInfo {
    id: usize,
    letter: String,
    size: usize,
    consensus: String,
    mean_intra_hamming: f64,
    max_intra_hamming: u32,
    fraction_of_total: f64,
    cenp_b_box_score: f64,
}

#[derive(Debug, Serialize)]
struct AlphabetReport {
    n_monomers_total: usize,
    n_monomers_clustered: usize,
    n_clusters: usize,
    silhouette_score: f64,
    clusters: Vec<ClusterInfo>,
    alphabet_string: String,
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    // --- Step 1: Read canonical monomers ---
    eprintln!("Reading {}...", args.input);
    let monomers = read_canonical_tsv(&args.input);
    eprintln!("  {} total monomers", monomers.len());

    let exact: Vec<usize> = (0..monomers.len())
        .filter(|&i| monomers[i].sequence.len() == args.canon_len)
        .collect();
    eprintln!("  {} monomers with length={}", exact.len(), args.canon_len);

    // --- Step 2: Subsample and compute distance matrix ---
    let sub_indices: Vec<usize> = if exact.len() <= args.subsample {
        exact.clone()
    } else {
        let step = exact.len() / args.subsample;
        exact.iter().step_by(step.max(1)).take(args.subsample).copied().collect()
    };
    let n_sub = sub_indices.len();
    eprintln!("Computing {}x{} hamming distance matrix...", n_sub, n_sub);

    let sub_seqs: Vec<&[u8]> = sub_indices.iter().map(|&i| monomers[i].sequence.as_slice()).collect();

    // Compute distance matrix (upper triangle stored as flat vec)
    // dist[i][j] = dist_flat[i * n_sub + j]
    let dist_matrix: Vec<u16> = (0..n_sub).into_par_iter().flat_map(|i| {
        let mut row = vec![0u16; n_sub];
        for j in 0..n_sub {
            if i != j {
                row[j] = hamming(sub_seqs[i], sub_seqs[j]) as u16;
            }
        }
        row
    }).collect();

    eprintln!("  Distance matrix computed ({:.1} MB)", (n_sub * n_sub * 2) as f64 / 1e6);

    // --- Step 3: Find optimal k using silhouette score ---
    let (_best_k, best_assignments, best_silhouette) = if args.target_k > 0 {
        eprintln!("Running k-medoids with k={}...", args.target_k);
        let (assign, _medoids) = kmedoids(&dist_matrix, n_sub, args.target_k, args.iterations);
        let sil = silhouette_score(&dist_matrix, n_sub, &assign);
        eprintln!("  k={}, silhouette={:.4}", args.target_k, sil);
        (args.target_k, assign, sil)
    } else {
        eprintln!("Scanning k={}..{} for optimal silhouette...", args.min_k, args.max_k);
        let mut best_k = args.min_k;
        let mut best_sil = f64::NEG_INFINITY;
        let mut best_assign = vec![0usize; n_sub];

        // Scan in steps for speed, then refine
        let step = 5;
        let candidates: Vec<usize> = (args.min_k..=args.max_k).step_by(step).collect();

        for &k in &candidates {
            let (assign, _) = kmedoids(&dist_matrix, n_sub, k, args.iterations);
            let sil = silhouette_score(&dist_matrix, n_sub, &assign);
            eprintln!("  k={:3}, silhouette={:.4}", k, sil);
            if sil > best_sil {
                best_sil = sil;
                best_k = k;
                best_assign = assign;
            }
        }

        // Refine around best
        let refine_lo = best_k.saturating_sub(step).max(args.min_k);
        let refine_hi = (best_k + step).min(args.max_k);
        for k in refine_lo..=refine_hi {
            if candidates.contains(&k) { continue; }
            let (assign, _) = kmedoids(&dist_matrix, n_sub, k, args.iterations);
            let sil = silhouette_score(&dist_matrix, n_sub, &assign);
            eprintln!("  k={:3}, silhouette={:.4} (refine)", k, sil);
            if sil > best_sil {
                best_sil = sil;
                best_k = k;
                best_assign = assign;
            }
        }

        eprintln!("Best k={}, silhouette={:.4}", best_k, best_sil);
        (best_k, best_assign, best_sil)
    };

    // --- Step 4: Build consensus for each cluster from subsample ---
    let all_seqs: Vec<&[u8]> = exact.iter().map(|&i| monomers[i].sequence.as_slice()).collect();

    let mut cluster_members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, &c) in best_assignments.iter().enumerate() {
        cluster_members.entry(c).or_default().push(i);
    }

    let mut consensuses: Vec<(usize, Vec<u8>)> = cluster_members.iter().map(|(&cid, indices)| {
        let seqs_for_cons: Vec<&[u8]> = indices.iter().map(|&i| sub_seqs[i]).collect();
        let cons = build_consensus_from_refs(&seqs_for_cons, args.canon_len);
        (cid, cons)
    }).collect();
    consensuses.sort_by_key(|c| c.0);

    // --- Step 5: Assign ALL exact-length monomers to nearest consensus ---
    eprintln!("Assigning all {} monomers to nearest consensus...", all_seqs.len());

    let consensus_refs: Vec<&[u8]> = consensuses.iter().map(|(_, c)| c.as_slice()).collect();
    let cluster_ids: Vec<usize> = consensuses.iter().map(|(id, _)| *id).collect();

    let all_assignments: Vec<(usize, u32)> = all_seqs.par_iter().map(|seq| {
        let mut best_cid = 0;
        let mut best_dist = u32::MAX;
        for (i, cons) in consensus_refs.iter().enumerate() {
            let d = hamming(seq, cons);
            if d < best_dist {
                best_dist = d;
                best_cid = cluster_ids[i];
            }
        }
        (best_cid, best_dist)
    }).collect();

    // --- Step 6: Build final stats ---
    let mut final_sizes: HashMap<usize, usize> = HashMap::new();
    for &(cid, _) in &all_assignments {
        *final_sizes.entry(cid).or_insert(0) += 1;
    }

    // Sort by size, assign letters
    let mut size_sorted: Vec<(usize, usize)> = final_sizes.into_iter().collect();
    size_sorted.sort_by(|a, b| b.1.cmp(&a.1));

    let mut id_to_letter: HashMap<usize, String> = HashMap::new();
    for (i, (cid, _)) in size_sorted.iter().enumerate() {
        let letter = if i < 26 {
            format!("{}", (b'A' + i as u8) as char)
        } else {
            format!("Z{}", i - 25)
        };
        id_to_letter.insert(*cid, letter);
    }

    // Rebuild consensus from ALL assigned monomers
    let n_total = all_seqs.len();
    let mut cluster_infos: Vec<ClusterInfo> = size_sorted.par_iter().map(|(cid, size)| {
        let indices: Vec<usize> = all_assignments.iter().enumerate()
            .filter(|(_, (c, _))| *c == *cid)
            .map(|(i, _)| i)
            .collect();

        let seqs_for_cons: Vec<&[u8]> = indices.iter().map(|&i| all_seqs[i]).collect();
        let cons = build_consensus_from_refs(&seqs_for_cons, args.canon_len);
        let cons_str = String::from_utf8(cons.clone()).unwrap_or_default();

        // Sample stats
        let sample_n = indices.len().min(2000);
        let step = (indices.len() / sample_n).max(1);
        let sample_idx: Vec<usize> = indices.iter().step_by(step).take(sample_n).copied().collect();

        let dists: Vec<u32> = sample_idx.iter().map(|&i| hamming(all_seqs[i], &cons)).collect();
        let mean_dist = dists.iter().sum::<u32>() as f64 / dists.len() as f64;
        let max_dist = *dists.iter().max().unwrap_or(&0);

        let cenp_b_score = score_cenp_b_box(&cons);

        ClusterInfo {
            id: *cid,
            letter: id_to_letter[cid].clone(),
            size: *size,
            consensus: cons_str,
            mean_intra_hamming: mean_dist,
            max_intra_hamming: max_dist,
            fraction_of_total: *size as f64 / n_total as f64,
            cenp_b_box_score: cenp_b_score,
        }
    }).collect();

    cluster_infos.sort_by(|a, b| a.letter.cmp(&b.letter));

    let n_clusters_final = cluster_infos.len();
    eprintln!("\n=== ALPHABET: {} letters (silhouette={:.4}) ===", n_clusters_final, best_silhouette);
    for ci in &cluster_infos {
        let b_marker = if ci.cenp_b_box_score > 0.5 { "B+" } else { "B-" };
        eprintln!("  {}: n={:>5} ({:>4.1}%), mean_d={:>5.1}, max_d={:>3}, {} {}...{}",
            ci.letter, ci.size, ci.fraction_of_total * 100.0,
            ci.mean_intra_hamming, ci.max_intra_hamming, b_marker,
            &ci.consensus[..20.min(ci.consensus.len())],
            &ci.consensus[ci.consensus.len().saturating_sub(20)..]);
    }

    // --- Step 7: Write outputs ---
    eprintln!("\nWriting outputs...");

    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.output).unwrap());
        writeln!(f, "array_id\tidx\tletter\tcluster_id\tdist_to_consensus\tsequence").unwrap();
        for (i, &orig_idx) in exact.iter().enumerate() {
            let entry = &monomers[orig_idx];
            let (cid, dist) = all_assignments[i];
            let letter = &id_to_letter[&cid];
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}",
                entry.array_id, entry.idx, letter, cid, dist,
                std::str::from_utf8(&entry.sequence).unwrap_or("")
            ).unwrap();
        }
    }

    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.consensus_fa).unwrap());
        for ci in &cluster_infos {
            writeln!(f, ">letter_{} n={} mean_dist={:.1} cenp_b={:.2}",
                ci.letter, ci.size, ci.mean_intra_hamming, ci.cenp_b_box_score).unwrap();
            writeln!(f, "{}", ci.consensus).unwrap();
        }
    }

    let alphabet_str: String = cluster_infos.iter().map(|c| c.letter.as_str()).collect::<Vec<_>>().join("");
    let report = AlphabetReport {
        n_monomers_total: monomers.len(),
        n_monomers_clustered: exact.len(),
        n_clusters: n_clusters_final,
        silhouette_score: best_silhouette,
        clusters: cluster_infos,
        alphabet_string: alphabet_str,
    };
    let json = serde_json::to_string_pretty(&report).unwrap();
    std::fs::write(&args.report, json).unwrap();

    eprintln!("Done.");
}

// ============ K-MEDOIDS ============

fn kmedoids(dist: &[u16], n: usize, k: usize, max_iter: usize) -> (Vec<usize>, Vec<usize>) {
    // PAM-like k-medoids
    // Initialize medoids: first = most central point, then farthest-point seeding

    // Find most central point
    let centrality: Vec<u64> = (0..n).into_par_iter().map(|i| {
        let mut total: u64 = 0;
        for j in 0..n {
            total += dist[i * n + j] as u64;
        }
        total
    }).collect();

    let mut medoids: Vec<usize> = Vec::with_capacity(k);
    let first = centrality.iter().enumerate().min_by_key(|(_, &c)| c).unwrap().0;
    medoids.push(first);

    // Farthest-point seeding for remaining medoids
    let mut min_dist_to_medoid = vec![u16::MAX; n];
    for _ in 1..k {
        // Update min distances
        let last_med = *medoids.last().unwrap();
        for i in 0..n {
            let d = dist[i * n + last_med];
            if d < min_dist_to_medoid[i] {
                min_dist_to_medoid[i] = d;
            }
        }
        // Pick farthest
        let next = min_dist_to_medoid.iter().enumerate()
            .filter(|(i, _)| !medoids.contains(i))
            .max_by_key(|(_, &d)| d)
            .unwrap().0;
        medoids.push(next);
    }

    let mut assignments = vec![0usize; n];

    for _iter in 0..max_iter {
        // Assign each point to nearest medoid
        let new_assignments: Vec<usize> = (0..n).into_par_iter().map(|i| {
            let mut best_m = 0;
            let mut best_d = u16::MAX;
            for (mi, &med) in medoids.iter().enumerate() {
                let d = dist[i * n + med];
                if d < best_d {
                    best_d = d;
                    best_m = mi;
                }
            }
            best_m
        }).collect();

        // Update medoids: for each cluster, find the point that minimizes total distance
        let new_medoids: Vec<usize> = (0..k).into_par_iter().map(|mi| {
            let members: Vec<usize> = new_assignments.iter().enumerate()
                .filter(|(_, &a)| a == mi)
                .map(|(i, _)| i)
                .collect();

            if members.is_empty() {
                return medoids[mi]; // keep old
            }

            let best = members.iter().min_by_key(|&&i| {
                members.iter().map(|&j| dist[i * n + j] as u64).sum::<u64>()
            }).unwrap();
            *best
        }).collect();

        let converged = new_medoids == medoids;
        assignments = new_assignments;
        medoids = new_medoids;

        if converged {
            break;
        }
    }

    (assignments, medoids)
}

fn silhouette_score(dist: &[u16], n: usize, assignments: &[usize]) -> f64 {
    // Sample for speed
    let sample_n = n.min(3000);
    let step = (n / sample_n).max(1);

    let k = *assignments.iter().max().unwrap_or(&0) + 1;

    let scores: Vec<f64> = (0..n).step_by(step).take(sample_n).collect::<Vec<usize>>()
        .par_iter().map(|&i| {
        let my_cluster = assignments[i];

        // a(i) = mean distance to same cluster
        let mut same_total: u64 = 0;
        let mut same_count: u64 = 0;
        for j in 0..n {
            if j != i && assignments[j] == my_cluster {
                same_total += dist[i * n + j] as u64;
                same_count += 1;
            }
        }
        if same_count == 0 { return 0.0; }
        let a = same_total as f64 / same_count as f64;

        // b(i) = min mean distance to other clusters
        let mut b = f64::MAX;
        for c in 0..k {
            if c == my_cluster { continue; }
            let mut other_total: u64 = 0;
            let mut other_count: u64 = 0;
            for j in 0..n {
                if assignments[j] == c {
                    other_total += dist[i * n + j] as u64;
                    other_count += 1;
                }
            }
            if other_count > 0 {
                let mean = other_total as f64 / other_count as f64;
                if mean < b { b = mean; }
            }
        }

        if a.max(b) == 0.0 { 0.0 } else { (b - a) / a.max(b) }
    }).collect();

    scores.iter().sum::<f64>() / scores.len() as f64
}

// ============ UTILITIES ============

fn read_canonical_tsv(path: &str) -> Vec<MonomerEntry> {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(path).unwrap();
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);
    let mut entries = Vec::new();
    let mut first = true;
    let mut col_array = 0usize;
    let mut col_idx = 0usize;
    let mut col_seq = 0usize;

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
        if fields.len() <= col_seq { continue; }
        let seq = fields[col_seq].as_bytes().to_vec();
        if seq.is_empty() { continue; }

        entries.push(MonomerEntry {
            array_id: fields[col_array].to_string(),
            idx: fields[col_idx].parse().unwrap_or(0),
            sequence: seq,
        });
    }
    entries
}

fn hamming(a: &[u8], b: &[u8]) -> u32 {
    let n = a.len().min(b.len());
    let mut dist = 0u32;
    for i in 0..n {
        if a[i] != b[i] { dist += 1; }
    }
    dist += (a.len() as i32 - b.len() as i32).unsigned_abs();
    dist
}

fn build_consensus_from_refs(seqs: &[&[u8]], len: usize) -> Vec<u8> {
    let mut counts = vec![[0u32; 4]; len];
    for seq in seqs {
        for (pos, &b) in seq.iter().enumerate().take(len) {
            let idx = match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => continue,
            };
            counts[pos][idx] += 1;
        }
    }
    let bases = [b'A', b'C', b'G', b'T'];
    counts.iter().map(|c| {
        let best = c.iter().enumerate().max_by_key(|(_, &v)| v).unwrap().0;
        bases[best]
    }).collect()
}

fn score_cenp_b_box(seq: &[u8]) -> f64 {
    let pattern: [(usize, u8); 9] = [
        (1, b'T'), (2, b'T'), (3, b'C'), (4, b'G'),
        (9, b'A'), (12, b'C'), (13, b'G'), (14, b'G'), (15, b'G'),
    ];
    if seq.len() < 17 { return 0.0; }

    let mut best_score = 0.0f64;
    for start in 0..=(seq.len() - 17) {
        let mut matches = 0;
        for &(offset, expected) in &pattern {
            if seq[start + offset].to_ascii_uppercase() == expected {
                matches += 1;
            }
        }
        let score = matches as f64 / pattern.len() as f64;
        if score > best_score { best_score = score; }
    }
    best_score
}

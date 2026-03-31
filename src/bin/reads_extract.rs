use std::io::{BufRead, BufReader, Write, Read as IoRead};
use std::process::{Command, Stdio};

/// Pass 1: Extract reads containing satellite DNA anchor sites.
/// Streams FASTQ(.gz), checks for chain anchor sites, outputs matching reads as FASTA.
///
/// Usage: reads_extract <input.fastq.gz> <chains.json> [min_sites] [threads]
///
/// A read passes if it contains >= min_sites distinct anchor sites (default: 3).

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: reads_extract <input.fastq.gz> <chains.json> [min_sites=3]");
        std::process::exit(1);
    }
    let input = &args[1];
    let chains_path = &args[2];
    let min_sites: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3);

    // Load chain sites from JSON
    let sites = load_chain_sites(chains_path);
    eprintln!("Loaded {} anchor sites from {}", sites.len(), chains_path);
    for (i, s) in sites.iter().enumerate() {
        eprintln!("  site{}: {} ({}bp)", i, String::from_utf8_lossy(s), s.len());
    }
    eprintln!("Min sites per read: {}", min_sites);

    // Also prepare homopolymer-compressed versions of sites
    let sites_hpc: Vec<Vec<u8>> = sites.iter().map(|s| hpc(s)).collect();
    eprintln!("HPC sites:");
    for (i, s) in sites_hpc.iter().enumerate() {
        eprintln!("  site{}_hpc: {} ({}bp)", i, String::from_utf8_lossy(s), s.len());
    }

    // Open input
    let is_gzipped = input.ends_with(".gz");
    let is_fastq = input.contains(".fastq") || input.contains(".fq");

    let reader: Box<dyn IoRead + Send> = if is_gzipped {
        let decompress = if which_exists("pigz") { "pigz" } else { "gzip" };
        eprintln!("Streaming: {} -dc {}", decompress, input);
        let child = Command::new(decompress)
            .args(["-dc", input])
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("Failed to start decompressor");
        Box::new(child.stdout.unwrap())
    } else {
        Box::new(std::fs::File::open(input).expect("Cannot open file"))
    };

    let buf_reader = BufReader::with_capacity(32 * 1024 * 1024, reader);
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::with_capacity(16 * 1024 * 1024, stdout.lock());

    let mut lines = buf_reader.lines();
    let mut total_reads = 0usize;
    let mut passed_reads = 0usize;
    let mut total_bp = 0usize;
    let mut passed_bp = 0usize;

    loop {
        // Read one FASTQ record
        let header = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        if is_fastq {
            let _ = lines.next(); // +
            let _ = lines.next(); // qual
        }

        total_reads += 1;
        total_bp += seq.len();

        // Check: how many distinct sites found in this read?
        let seq_bytes = seq.as_bytes();
        let seq_upper: Vec<u8> = seq_bytes.iter().map(|b| b.to_ascii_uppercase()).collect();
        let seq_hpc = hpc(&seq_upper);
        let rc = revcomp(&seq_upper);
        let rc_hpc = hpc(&rc);

        let mut found_sites = 0usize;
        for (site, site_hpc) in sites.iter().zip(sites_hpc.iter()) {
            // Try exact match first
            let found = contains_subseq(&seq_upper, site)
                || contains_subseq(&rc, site)
                // Then HPC match
                || contains_subseq(&seq_hpc, site_hpc)
                || contains_subseq(&rc_hpc, site_hpc);
            if found {
                found_sites += 1;
            }
        }

        if found_sites >= min_sites {
            passed_reads += 1;
            passed_bp += seq.len();
            // Output as FASTA
            let name = if header.starts_with('@') { &header[1..] } else { &header[1..] };
            writeln!(out, ">{} sites={}", name.split_whitespace().next().unwrap_or(name), found_sites).unwrap();
            writeln!(out, "{}", seq).unwrap();
        }

        if total_reads % 100000 == 0 {
            eprintln!("\r  {} reads, {} passed ({:.1}%), {:.1}Gb",
                total_reads, passed_reads,
                passed_reads as f64 / total_reads as f64 * 100.0,
                total_bp as f64 / 1e9);
        }
    }

    out.flush().unwrap();

    eprintln!("\n\nDone:");
    eprintln!("  Total reads: {}", total_reads);
    eprintln!("  Passed reads: {} ({:.2}%)", passed_reads, passed_reads as f64 / total_reads as f64 * 100.0);
    eprintln!("  Total bp: {:.2}Gb", total_bp as f64 / 1e9);
    eprintln!("  Passed bp: {:.2}Gb ({:.2}%)", passed_bp as f64 / 1e9, passed_bp as f64 / total_bp as f64 * 100.0);
}

/// Homopolymer compression: AAACCCGT → ACGT
fn hpc(seq: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(seq.len());
    for &b in seq {
        if result.last().copied() != Some(b) {
            result.push(b);
        }
    }
    result
}

fn contains_subseq(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.len() > haystack.len() { return false; }
    haystack.windows(needle.len()).any(|w| w == needle)
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N',
    }).collect()
}

fn which_exists(cmd: &str) -> bool {
    Command::new("which").arg(cmd).output().map(|o| o.status.success()).unwrap_or(false)
}

fn load_chain_sites(path: &str) -> Vec<Vec<u8>> {
    // Parse chains.json to extract site sequences
    // Format: look for "seq": "ACATCACAAAG" patterns
    let content = std::fs::read_to_string(path).expect("Cannot read chains file");
    let mut sites = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if let Some(start) = line.find("\"seq\"") {
            // Extract the sequence value
            if let Some(colon) = line[start..].find(':') {
                let after_colon = &line[start + colon + 1..];
                if let Some(q1) = after_colon.find('"') {
                    if let Some(q2) = after_colon[q1 + 1..].find('"') {
                        let seq = &after_colon[q1 + 1..q1 + 1 + q2];
                        let seq_bytes = seq.as_bytes().to_vec();
                        if seq_bytes.len() >= 6 && !sites.contains(&seq_bytes) {
                            sites.push(seq_bytes);
                        }
                    }
                }
            }
        }
    }
    sites
}

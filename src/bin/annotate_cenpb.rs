use std::io::{BufRead, BufReader, Write};
use alphasplitter::monomer::revcomp;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: annotate_cenpb <monomers.tsv> <output.tsv>");
        std::process::exit(1);
    }

    // CENP-B box: nTTCGnnnnAnnCGGGn (17bp)
    let fixed: [(usize, u8); 9] = [
        (1, b'T'), (2, b'T'), (3, b'C'), (4, b'G'),
        (9, b'A'),
        (12, b'C'), (13, b'G'), (14, b'G'), (15, b'G'),
    ];
    let _var_pos: [usize; 8] = [0, 5, 6, 7, 8, 10, 11, 16];

    let file = std::fs::File::open(&args[1]).unwrap();
    let reader = BufReader::new(file);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&args[2]).unwrap());

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        if lineno == 0 {
            // Header
            writeln!(out, "{}\tcenpb\tcenpb_score\tcenpb_pos\tcenpb_strand\tcenpb_cigar9\tcenpb_cigar17\tcenpb_seq", line).unwrap();
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let seq_col = fields.len() - 1; // sequence is last column
        let seq = fields[seq_col].as_bytes();

        // Search both strands
        let rc = revcomp(seq);
        let mut best_score: u32 = 0;
        let mut best_pos: usize = 0;
        let mut best_strand = '+';
        let mut best_cigar9 = String::new();
        let mut best_cigar17 = String::new();
        let mut best_seq = String::new();

        for (strand_seq, strand) in [(&seq[..], '+'), (&rc[..], '-')] {
            if strand_seq.len() < 17 { continue; }
            for i in 0..strand_seq.len() - 16 {
                let mut score = 0u32;
                for &(pos, base) in &fixed {
                    if strand_seq[i + pos].to_ascii_uppercase() == base {
                        score += 1;
                    }
                }
                if score > best_score {
                    best_score = score;
                    best_pos = if strand == '+' { i } else { seq.len() - i - 17 };
                    best_strand = strand;

                    // Build CIGAR9 (fixed positions only)
                    let mut c9 = String::new();
                    for &(pos, base) in &fixed {
                        let actual = strand_seq[i + pos].to_ascii_uppercase();
                        if actual == base {
                            c9.push('=');
                        } else {
                            c9.push((actual as char).to_ascii_lowercase());
                        }
                    }
                    best_cigar9 = c9;

                    // Build CIGAR17 (all positions)
                    let mut c17 = String::new();
                    let box_seq = &strand_seq[i..i + 17];
                    for p in 0..17 {
                        let actual = box_seq[p].to_ascii_uppercase() as char;
                        let is_fixed = fixed.iter().any(|&(fp, _)| fp == p);
                        if is_fixed {
                            let expected = fixed.iter().find(|&&(fp, _)| fp == p).unwrap().1 as char;
                            if actual == expected {
                                c17.push('=');
                            } else {
                                c17.push(actual.to_ascii_lowercase());
                            }
                        } else {
                            c17.push(actual);
                        }
                    }
                    best_cigar17 = c17;
                    best_seq = String::from_utf8(box_seq.to_vec()).unwrap_or_default();
                }
            }
        }

        let cenpb_label = if best_score >= 7 { "B+" } else if best_score == 6 { "B?" } else { "B-" };
        let pos_str = if best_score >= 6 { format!("{}", best_pos) } else { ".".to_string() };
        let strand_str = if best_score >= 6 { format!("{}", best_strand) } else { ".".to_string() };
        let c9_str = if best_score >= 6 { &best_cigar9 } else { "." };
        let c17_str = if best_score >= 6 { &best_cigar17 } else { "." };
        let seq_str = if best_score >= 6 { &best_seq } else { "." };

        writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            line, cenpb_label, best_score, pos_str, strand_str, c9_str, c17_str, seq_str).unwrap();
    }
}


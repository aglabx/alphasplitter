# AlphaSplitter

## What is this

Ab initio satellite DNA monomer alphabet discovery tool. Rust, CLI.

## Build & Run

```bash
cargo build --release
./target/release/discover_chains input.fasta -o chains.json -t 96
./target/release/motif_cut input.fasta -m chains.json -o monomers.tsv -t 96
./target/release/annotate_cenpb monomers.tsv annotated.tsv
```

## Architecture

- `src/lib.rs` — shared modules (io, kmer, monomer)
- `src/bin/discover_chains.rs` — main chain discovery pipeline
- `src/bin/motif_cut.rs` — monomer cutting + classification
- `src/bin/annotate_cenpb.rs` — CENP-B box annotation
- `src/bin/find_box.rs` — generic box pattern search
- `src/bin/cenpb_spacing.rs` — spacing analysis
- `src/bin/find_periodic_boxes.rs` — de novo periodic box discovery

## Key Design Decisions

- **Indels, not substitutions**: Letters defined by site presence/absence + spacer distances
- **Chain-first**: Sites discovered by periodicity, then assembled into chain
- **Distance = presence**: If neighbor distance matches, intermediate site is present (no fuzzy matching)
- **Both strands**: CENP-B searched on both strands (critical fix)

## Server

```
Host: 77.234.216.99
Data: /mnt/data/claude/2026-03-30_alphasplitter/
Results: results/R1/{CHM13,zebrafinch,mouse}/
```

## Testing

End-to-end test: run full pipeline on a small test FASTA, verify output format and known letters.

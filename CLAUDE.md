# AlphaSplitter

## What is this

Ab initio satellite DNA monomer alphabet discovery tool. Rust, single CLI.

## Build & Run

```bash
cargo build --release

# End-to-end: discover → cut → annotate
./target/release/alphasplitter run input.fasta -o out/ -t 96

# Or stage-by-stage
./target/release/alphasplitter discover input.fasta -o chains.json -t 96
./target/release/alphasplitter cut      input.fasta -m chains.json -o monomers.tsv -t 96
./target/release/alphasplitter annotate monomers.tsv annotated.tsv
```

See `alphasplitter --help` for the full list of subcommands.

## Architecture

- `src/lib.rs` — shared modules (io, kmer, monomer, cmd)
- `src/main.rs` — CLI dispatcher (clap)
- `src/cmd/discover_chains.rs` — `discover`: chain discovery pipeline
- `src/cmd/motif_cut.rs` — `cut`: monomer cutting + classification
- `src/cmd/annotate_cenpb.rs` — `annotate`: CENP-B box annotation
- `src/cmd/find_box.rs` — `find-box`: generic box pattern search
- `src/cmd/cenpb_spacing.rs` — `spacing`: CENP-B spacing analysis
- `src/cmd/find_periodic_boxes.rs` — `find-periodic-boxes`: de novo box discovery
- `src/cmd/{build_hpc_hmms,scan_reads,classify_reads,reads_extract,reads_alphabet}.rs`
  — `reads {build-hmms,scan,classify,extract,alphabet}`: ONT HPC-HMM pipeline
- `src/cmd/{find_phase,motif_graph}.rs` — hidden `dev {find-phase,motif-graph}` (research)

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

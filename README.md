# AlphaSplitter

Ab initio discovery of satellite DNA monomer alphabets from tandem repeat arrays.

AlphaSplitter identifies periodic conserved anchor sites in satellite DNA, assembles them into chain templates, and classifies monomers by structural indels — without reference HMMs, curated databases, or prior knowledge of the satellite family.

## Quick Start

```bash
# Build
cargo build --release

# Full auto-satellome (discovers all periods and subfamilies)
alphasplitter discover-chains input.10kb.fasta -o chains.json -t 96

# Cut monomers and classify
alphasplitter motif-cut input.10kb.fasta -m chains.json -o monomers.tsv -t 96

# Annotate CENP-B boxes
alphasplitter annotate-cenpb monomers.tsv annotated.tsv
```

## What it does

**Input:** FASTA file of satellite DNA arrays (e.g., from [Satellome](https://github.com/ad3002/satellome) or TRASH)

**Output:**
- `chains.json` — discovered anchor sites and chain structure per satellite family
- `monomers.tsv` — per-monomer annotation: letter, subtype, site_order, site_structure, CENP-B box
- HOR strings showing higher-order repeat patterns

## Pipeline

```
satellome .10kb.fasta
  │
  ▼
discover_chains (per period, iterative SF discovery)
  ├── Scan 8-mers → filter by periodicity + support + complexity
  ├── Grow by Shannon entropy (cap 11bp)
  ├── Discover links (exact co-occurrence, stable spacing)
  ├── Assemble chain → enrichment (pairwise distances + occupancy ≥30%)
  └── → chains.json
  │
  ▼
motif_cut -m chains.json
  ├── Scan arrays for anchor sites (both strands)
  ├── Cut at first-site positions → monomers
  ├── Classify: letter (indel-defined) + subtype (site variants)
  ├── site_order: abcde, ab-de, abdebcde (dimers)
  ├── site_structure: 0:a18b24c43d18e:57
  └── → monomers.tsv + HOR strings
  │
  ▼
annotate_cenpb
  ├── CENP-B box scoring (nTTCGnnnnAnnCGGGn, both strands)
  ├── 9-char and 17-char CIGAR notation
  └── → annotated.tsv
```

## Core Model

A satellite monomer is a chain of conserved anchor sites separated by variable spacers:

```
[--site_A--]---spacer---[--site_B--]---spacer---[--site_C--]...
  conserved    LENGTH      conserved    LENGTH     conserved
  sequence     matters!    sequence     matters!   sequence
```

**Letter** = defined by indels only (which sites present/deleted + spacer length changes), NOT by substitutions.

**Site presence detection:** A site is present if the distance between its neighbors matches expected — no fuzzy sequence matching needed.

## Binaries

| Binary | Purpose |
|--------|---------|
| `discover_chains` | Chain-first anchor discovery + enrichment + multi-SF/period |
| `motif_cut` | Array cutting + alphabet + HOR + site_structure |
| `annotate_cenpb` | Add CENP-B box annotation to monomer TSV |
| `find_box` | Search for any degenerate box pattern (CENP-B, TIGD4, etc.) |
| `cenpb_spacing` | CENP-B box spacing analysis at single-base resolution |
| `find_periodic_boxes` | De novo periodic box discovery from paired k-mers |
| `discover_motifs` | v1 conservation-profile motif discovery |
| `motif_graph` | Motif transition graph (diagnostic) |
| `find_phase` | Phase optimization experiment |

## Validated on

| Species | Satellite | Period | Sites | Arrays | Monomers | Letters |
|---------|-----------|--------|-------|--------|----------|---------|
| Human (CHM13) | Alpha satellite | 171bp | 12 | 522/1223 | 411K | 585 |
| Human (CHM13) | Human sat I | 2420bp | 87 | 70 | — | — |
| Mouse (C57BL/6J) | Minor satellite | 120bp | 8 | 17/28 | 44K | 105 |
| Zebra finch | 191bp satellite | 191bp | 11 | 262 | 104K | 279 |
| Zebra finch | 717bp satellite | 717bp | 33 | 70 | 82K | 127 |
| Zebra finch | 3772bp satellite | 3772bp | 131 | 108 | — | — |
| 7 primate species | Alpha satellite | 171bp | 5 (v1) | — | 8M | 26 |

Alexandrov validation: 38/39 SF types recovered, 96.3% of monomers assigned.

## Key Findings

- **585 structural letter types** in human alpha satellite (15× more than HumAS-HMMER's 39)
- **CENP-B variable positions** are chromosome-specific suprachromosomal family markers
- **CENP-B spacing within CDR** = chromosome-specific nucleosome positioning signature
- **Mouse centromeric satellite saturated** with CENP-B boxes (89% of monomers)
- **Zebra finch has two alternative box systems** (BoxA-D at 355bp, TIGD4 at 715bp) on different chromosomes
- **Convergent architecture**: human (340bp), ZF (355bp), mouse (120bp) → 1 box ≈ 1 nucleosome repeat

## Output Format (annotated TSV)

```
#alphasplitter v1.0.2
#chain_sites: a=ACATCACAAAG b=AGAATGCTTCT c=GAAGATATTTC d=TCCACTTGCAG e=AAAGAGTGTTT
#CENP-B_box_ref: nTTCGnnnnAnnCGGGn (17bp), searched on BOTH strands
```

| Column | Description |
|--------|-------------|
| chr/array_id | Chromosome or array identifier |
| start, end | Genomic/array coordinates |
| letter | Structural type (indel-defined) |
| site_order | Chain sites: abcde, ab-de, abdebcde |
| site_structure | gap:sites_with_distances:gap |
| cenpb | B+ (≥7/9), B? (6/9), B- (≤5/9) |
| cenpb_cigar17 | Full 17-char CIGAR |
| sequence | Monomer DNA sequence |

## Requirements

- Rust 1.70+
- Input: satellite arrays in FASTA (recommend ≥10kb arrays)

## Citation

Komissarov A. (2026). AlphaSplitter: Ab initio discovery of satellite DNA monomer alphabets from tandem repeat arrays. *In preparation.*

## License

MIT

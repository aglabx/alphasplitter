# satkmer_compare — satellite k-mer comparison of two BAM files

A turnkey pipeline: give it **two BAM files** (e.g. a matched normal/tumour pair, or any
two samples of the same individual) and it produces a **coverage-robust comparison of their
satellite-DNA content**, including a per-chromosome "satellite karyotype" that picks up
whole-chromosome gains/losses from satellite k-mers alone.

It was built on a kidney-chromophobe-carcinoma (KICH) normal/tumour pair, where it
recovers the known KICH whole-chromosome-loss signature (chr 1, 2, 6, 10, 13, 17, 21) from
the satellite reads that are usually discarded.

---

## What it does

1. **`kmc -k23 -ci3 -fbam`** on each BAM — counts **canonical 23-mers directly from the BAM**
   (no FASTQ unpacking), keeping k-mers seen ≥ 3 times, no upper cap. → two KMC databases.
2. **Rank comparison of high-copy k-mers** — for the abundant 23-mers (count ≥ 50 000 by
   default), order them by count in each sample and compare **ranks**, not counts.
   `rank_delta = rank₂ − rank₁`: negative ⇒ relatively *up* in sample 2, positive ⇒ relatively
   *down*. Surfaces which repeat families (HSat, alpha, microsatellites, …) shift — and any
   technical artifacts (e.g. an Illumina adapter k-mer exploding in rank flags a library-prep
   difference).
3. **Per-chromosome "satellite karyotype"** *(needs the reference pack, `--refpack`)* — for
   each chromosome, sum the donor counts of its **array-specific alpha-satellite marker
   23-mers** (1.26 M markers, each unique to one CHM13 alpha array → one chromosome), express
   as a fraction of the total, and report the **ratio = frac₂ / frac₁** + the rank shift.
   `ratio < 1` ⇒ that chromosome's centromeric alpha is relatively *lost* in sample 2 (a
   candidate whole-chromosome loss); `ratio > 1` ⇒ relatively *gained* (or compensatory rise
   when other chromosomes are lost).

**Why ranks / fractions, not raw counts?** Two libraries almost always differ in depth and
prep, so raw k-mer counts are not comparable between samples. Rank order and normalised
fractions are invariant to a global coverage scale.

---

## Requirements

- [`kmc`](https://github.com/refresh-bio/KMC) ≥ 3.x (provides `kmc` and `kmc_tools`) — reads BAM natively (`-fbam`)
- standard Unix: `bash`, `awk`, `sort`, `comm`, `join`, `column`
- `samtools` is **not** required.
- `blastn` only if you also want to run the optional alpha-read extraction (not part of this script).

Resources: KMC on a ~85 GB BAM ≈ 7 min / ~13 GB temp; a ~180 GB BAM ≈ 18 min / ~200 GB temp,
on 60 threads. Put the temp on a fast disk via `TMPDIR=...`. The KMC databases (~19 GB each
for human WGS at high coverage) are kept and reused on re-runs.

---

## Quick start

```bash
# 1. get the reference pack (once)  -- ~24 MB compressed, ~165 MB unpacked
#    (T2T-CHM13v2.0-derived: alpha-satellite array-specific markers, canonical alpha monomer,
#     clean alpha 23-mer sets)
#    -> from your hosting / shared drive, or scp from the analysis server:
scp <server>:/mnt/data/claude/2026-03-30_alphasplitter/work/refpack.tar.gz .
tar xzf refpack.tar.gz          # -> ./refpack/

# 2. run on a pair
./satkmer_compare.sh \
    -o results_DO220554 \
    -t 60 \
    --names Normal,Tumour \
    --refpack ./refpack \
    /path/DO220554_Normal.bam \
    /path/DO220554_Primary.tumour.bam

# 3. read the summary
less results_DO220554/SUMMARY.md
```

Options: `-o DIR` `-t THREADS` `-k 23` `-c 3` (KMC -ci) `--hicopy 50000` `--topn 2000`
`--names A,B` `--refpack DIR`. Run with `-h` for the full list.

---

## Outputs (in the `-o` directory)

| file | what |
|---|---|
| `<name>_k23.kmc_pre` / `.kmc_suf` | KMC k-mer databases (kept; reusable, e.g. with `kmc_tools`) |
| `<name>_hicopy_ranked.tsv` | the sample's high-copy 23-mers, sorted by count, with ranks |
| `rank_compare.tsv` | union of top-N high-copy 23-mers: `kmer · count₁ · rank₁ · count₂ · rank₂ · rank_delta` |
| `per_chromosome_satellite_cn.tsv` | *(with `--refpack`)* `chr · n_markers · det₁ · det₂ · sum₁ · sum₂ · frac₁(%) · frac₂(%) · ratio(2/1) · rank₁ · rank₂ · rank_delta` |
| `marker_counts.tsv` | *(with `--refpack`)* per-marker: `kmer · chr · count₁ · count₂` |
| `SUMMARY.md` | human-readable summary: KMC DB sizes, largest rank shifts, the per-chromosome table, candidate chromosome losses |
| `run.log` | full log |

See `example_output/` for the SUMMARY and per-chromosome table from the KICH test pair.

---

## Reference pack contents (`refpack/`)

| file | description |
|---|---|
| `marker_chr_map.tsv` | `canonical_23mer ⟶ field_idx ⟶ CHM13_array_id ⟶ chromosome ⟶ count_in_CHM13` — 1 260 118 array-specific alpha-satellite marker 23-mers (DF = 1 over the 65 BLAST-verified CHM13 alpha-satellite arrays > 100 kb). The map column 4 (chromosome) is what drives the per-chromosome aggregation. |
| `markers_k23.kmc_*` | auto-built by the script on the first `--refpack` run (a KMC DB of the marker k-mers, for fast intersection with the sample DBs). |
| `canonical_alpha.fa` | one canonical alpha-satellite monomer (~170 bp) — for BLAST sanity checks. |
| `clean_alpha_k23_DF1.txt` | all clean alpha 23-mers (1 530 806 canonical → 3 061 612 with reverse complements, one per line). |
| `clean_alpha_k23_DF4.txt` | conserved subset (present in ≥ 4 of the 65 alpha arrays): 60 341 → 120 682. Good precision/recall balance for read extraction. |
| `clean_alpha_k23_DF8.txt` | conserved core (≥ 8 arrays): 12 990 → 25 980. Highest precision, lowest recall. |

The `clean_alpha_k23_*.txt` files are forward+reverse-complement so they can be used with a
read extractor that does not reverse-complement reads (e.g. `V5_subset.exe`). They were
derived by: BLAST-classifying the > 100 kb CHM13 satellome arrays against one canonical alpha
monomer (→ 65 ALPHA arrays), taking all their canonical 23-mers, then subtracting 23-mers
that occur in the non-alpha satellite arrays (HSat2/3, DYZ, …), low-complexity 23-mers
(homopolymer run > 6, single base > 16/23, < 5 distinct dinucleotides), and 23-mers that are
ultra-abundant (> 2 M) in a control WGS BAM. Read-extraction precision on the KICH pair was
**97–99 %** against a clean alpha reference.

---

## Caveats / limitations

- **Ranks and fractions are coverage-robust but relative** — you cannot turn them into absolute
  copy numbers without a normalising assumption.
- **Tumour purity < 100 %** (tumour + stroma) dilutes losses: a heterozygous whole-chromosome
  loss reads as `ratio ≈ 0.5–0.8`, not 0.5.
- **CHM13-derived markers approximate the patient.** "Array-specific in CHM13" ≈ "array-specific
  in the patient" only for the dominant HOR variants; very divergent / patient-private alpha is
  invisible to these markers, and some patient alpha may sit on different chromosomes than in CHM13.
- The signal is **centromeric/pericentromeric alpha-array copy number**, which tracks
  whole-chromosome loss in carcinomas like KICH, but a focal centromeric change without an
  arm-level change would look the same.
- **`-ci 3`** removes most sequencing-error k-mers but also drops genuinely very-rare satellite
  variants. **k = 23** is specific but substitution-intolerant — the most divergent alpha is missed.
- **Library/technical differences** between the two BAMs (adapter content, trimming, chemistry)
  can produce spurious rank shifts; inspect `rank_compare.tsv` for adapter / known-artifact
  k-mers (`GATCGGAAGAGC…`, `AGATCGGAAGAGC…`, …) before interpreting.
- Tested on n = 1 patient — treat the biology as hypothesis-generating until replicated.

---

## Citation / provenance

Reference pack derived from **T2T-CHM13v2.0** (the complete human reference genome) — alpha-satellite
arrays from a satellome (tandem-repeat) annotation, BLAST-verified against a canonical alpha monomer.
Tools used: KMC 3.2.1, NCBI blastn 2.12, samtools 1.13.

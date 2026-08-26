# orf_bls_pipeline

Compute Branch Length Scores (BLS) and ORF conservation ages across species
using CodAlignView multi-species alignments and PRANK ancestral sequence
reconstruction. Optionally search against protein databases with BLAST and
render per-ORF phylogenetic trees.

## Installation

```bash
mamba create -n bls2026 python=3.11
mamba activate bls2026
mamba install -c conda-forge ete4 biopython matplotlib
mamba install -c bioconda prank iqtree ucsc-liftover blast
```

`matplotlib` is only required when using `--plot-trees`. Everything else
works without it.

## Quick start

```bash
mamba activate bls2026

python orf_bls_pipeline.py \
    --gtf orfs.gtf \
    --alnset hg38_120mammals \
    --ref-species Human \
    --outdir results/ \
    --threads 8
```

## Input

**GTF file** with ORFs. Each ORF needs an `orf_id` attribute (falls back to
`gene_id`). Accepted feature types: `ORF`, `CDS`, `exon`. Multi-exon ORFs
are joined by `orf_id`.

**Alignment set** from CodAlignView. Available sets include:
`hg38_100_tetrapod`, `hg38_120mammals`, `hg38_243primates`,
`hg38_470mammals`, `mm10`, `mm10_60`, `mm39_35_placental`,
`mm39_35_vertebrate`, `rn6`, `rn6_20`, `danRer7`, `sacCer3`,
`dm6_124insects`.
Full list: https://data.broadinstitute.org/compbio1/cav.php?Alnsets

**Species tree** in Newick format (optional — auto-downloaded if omitted).

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--gtf` | *required* | GTF file with ORFs |
| `--alnset` | *required* | CodAlignView alignment set |
| `--ref-species` | *required* | Reference species (e.g. `Human`, `Mouse`) |
| `--outdir` | *required* | Output directory |
| `--tree` | auto | Species tree (Newick) |
| `--stop-cutoff` | 0.7 | Min fraction of ORF that must be stop-codon-free |
| `--invalid-cutoff` | 0.5 | Max fraction of Ns in ungapped sequence |
| `--start-mode` | `same` | Start codon criterion (see below) |
| `--initiation-offset` | 3 | Codons from 5' end to search for start codon |
| `--min-identity` | 0.0 | Min AA identity vs ref to count as conserved (0=off) |
| `--repeats` | off | RepeatMasker `.out` file for repeat overlap annotation |
| `--blast` | off | Protein FASTA file for BLASTP search |
| `--blast-evalue` | 1e-4 | E-value threshold for BLAST hits |
| `--plot-trees` | off | Render per-ORF phylogenetic trees (naive + ASR) as PNG |
| `--threads` | 1 | Parallel ORFs |
| `--force` | off | Re-download and re-run everything |
| `-v` | off | Verbose logging |

### Start mode options

| Mode | Description |
|---|---|
| `same` | Start codon must match the reference's first codon or be ATG (default) |
| `nearcognate` | Any of ATG/CTG/GTG/TTG within `--initiation-offset` codons |
| `atg` | Only ATG within `--initiation-offset` codons |
| `none` | Start codon not checked |

## Output

### `bls_results.tsv`

Main results — **every ORF appears here**, including failed ones.

| Column | Description |
|---|---|
| `orf_id` | ORF identifier from GTF |
| `chrom`, `strand`, `n_exons` | Genomic coordinates |
| `origin_manner` | `denovo`, `nondenovo`, `too_few_species`, or failure reason |
| `origin_age_local` | Origination node in the ASR (local) tree |
| `origin_age_global` | Origination node in the full species tree |
| `n_species_with_orf_local` | Number of species with conserved ORF (from local/ASR tree) |
| `species_with_orf_local` | Semicolon-separated species list (local tree) |
| `most_distant_orf_species_local` | Most distant conserved species by PRANK tree branch length |
| `identity_conserved` | Mean AA identity of conserved species vs reference |
| `identity_notconserved` | Mean AA identity of non-conserved species vs reference |
| `repeat_overlaps` | RepeatMasker class/family overlaps (`;`-separated) or `NA` |
| `blast_matches` | BLASTP hit subject IDs (`;`-separated) or `none`/`NA` |
| `origin_species` | All species under the origination node |
| `bl_all` .. `bls_local_naive` | BLS fields (two methods) |

### `bls_results_noreconstructed.tsv`

Failed ORFs with reason codes.

### `fasta/` directory

Per-ORF FASTA files:
- `{orf_id}.aln.fa` — nucleotide alignment, all species (`>Species|status=1/0/-1`)
- `{orf_id}.conserved.fa` — conserved species only + reference
- `{orf_id}.prot.fa` — translated proteins
- `{orf_id}.asr.fa` — PRANK ancestral reconstructions (nucleotide).
  Headers include descendant species: `>Node17|status=1|species=Bonobo;Chimp;Gorilla;Human`

### `trees/` directory (only with `--plot-trees`)

Two PNG files per ORF:
- `{orf_id}.naive.png` — full species tree. Each leaf gets a colored dot:
  blue = ORF conserved (species in `species_with_orf_local`), red = aligned
  but not conserved, grey = species not in the local alignment.
- `{orf_id}.asr.png` — PRANK (local) tree with both extant leaves and
  ancestral `Node*` labels. Dots are colored the same way, using the
  conservation status computed on each sequence (ancestral reconstructions
  included). Grey here means the sequence was invalid (too many Ns, wrong
  start/end context).

Each plot includes a legend. Figure height scales with the number of leaves
(~0.13 inches per leaf) so 120-species trees render legibly.

## Examples

```bash
# Full pipeline: BLS + BLAST + repeats + tree plots
python orf_bls_pipeline.py \
    --gtf orfs.gtf \
    --alnset hg38_120mammals \
    --ref-species Human \
    --outdir results/ \
    --blast swissprot.fa \
    --blast-evalue 1e-4 \
    --repeats repeats/hg38.sorted.fa.out \
    --plot-trees \
    --threads 8

# Mouse ORFs, BLS only
python orf_bls_pipeline.py \
    --gtf orfs.gtf \
    --alnset mm10_60 \
    --ref-species Mouse \
    --outdir results/ \
    --threads 4

# Strict: ATG-only, 30% identity
python orf_bls_pipeline.py \
    --gtf orfs.gtf \
    --alnset hg38_120mammals \
    --ref-species Human \
    --outdir results/ \
    --start-mode atg \
    --min-identity 0.3 \
    --threads 8
```

## How it works

1. **GTF parsing**: Extracts ORF coordinates, joins multi-exon ORFs.
2. **Alignment download**: Fetches codon-aligned sequences from CodAlignView. Species with `.` (frameshift) or `|` (splice marker) are excluded. All-gap columns removed.
3. **PRANK**: Ancestral reconstruction with `-DNA -once -showanc`.
4. **Conservation**: Checks start codon (configurable), premature stops (`--stop-cutoff`), sequence validity, and optional minimum identity (`--min-identity`).
5. **BLS**: Two methods — ASR-based origination node, and naive MRCA.
6. **BLAST** (optional): Searches each ORF's reference protein sequence against a user-provided database with BLASTP.
7. **Repeat annotation** (optional): Overlaps ORF exons with RepeatMasker elements.
8. **Tree plotting** (optional): Per-ORF naive and ASR trees as PNG with conservation-colored dots and legends.

## Citation

Based on methods from Sandmann et al. 2023 (Mol Cell), Vakirlis et al. 2022
(Cell Reports), and CodAlignView (Jungreis/Broad Institute).

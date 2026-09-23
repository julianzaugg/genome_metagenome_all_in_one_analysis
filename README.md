# gmaio

Downstream analysis of [genome_metagenome_all_in_one](https://github.com/julianzaugg/genome_metagenome_all_in_one) pipeline output for metagenomes and isolate genomes.

`gmaio` is an R package plus project templates.
The package holds all the logic: readers for every pipeline output, sample linking, a common feature-profile data model, palettes, tables, figures and statistics.
Each project gets a `config.yml`, a `code/main.R` entry point that does the bulk processing, and one script per analysis.

## Install

```r
install.packages("remotes")
remotes::install_github("julianzaugg/genome_metagenome_all_in_one_analysis")
```

Optional packages used by specific analyses:

```r
install.packages(c("Nonpareil", "MicrobiomeStat", "svglite", "ape"))
BiocManager::install(c("maaslin3", "mixOmics", "ggtree"))
```

## Start a project

```r
gmaio::create_project("~/projects/J12345", mode = "metagenome")   # or mode = "isolate"
```

This creates:

```
J12345/
  config.yml        pipeline results directory, metadata, groups, colours
  data/             put metadata.xlsx (or .csv/.tsv) here
  code/main.R       bulk processing, run first
  code/*.R          one script per analysis
```

Edit `config.yml`, add the metadata, then from the project directory:

```bash
Rscript code/main.R
```

```bash
Rscript code/barcharts.R
```

`main.R` writes the standard tables to `Result_tables/`, the colour assignments to `palettes.yml`, and the processed project to `Result_other/processed.rds`.
Analysis scripts read that file, so they can be run in any order and rerun individually.

## Metadata and sample linking

The metadata must have one row per sample, with a column matching the pipeline samplesheet `sample` column (`metadata: sample_id_column`).
Tools label samples differently (`S1_R1`, `S1.clean_1.fastq.gz`, `S1.metabat2.5`, `S1___NODE_1_2`); gmaio maps them to metadata IDs by exact match, then by the longest ID followed by `.`, `_` or `-`.
Any pipeline sample that cannot be matched is an error.
Samples can be excluded from analyses (but kept in tables) with `exclude_column` or `exclude_samples`.

## Pipeline outputs

Outputs are found in `pipeline_results` by directory name suffix (`*_sylph`, `*_checkm2`, ...), so Illumina, Nanopore and isolate numbering all work.
`gmaio::pipeline_registry()` lists every output gmaio reads.
Any entry can be overridden in `config.yml` under `files:` with a file, glob or directory, for example when outputs have been copied elsewhere.
Outputs that were not produced (skipped pipeline steps) are reported and skipped.

## Colours

Metadata variables (group variables, discrete covariates, `colour_variables`, samples) and taxa get deterministic colours.
Colours set in `config.yml` under `colours:` win, then colours saved in `palettes.yml` from earlier runs, then generated ones.
Every generated colour is written to `palettes.yml`, so colours stay the same across reruns and you can edit any of them.
Taxa below phylum get distinct colours ordered by abundance (`outputs: taxa_colour_scheme: distinct`), or shades of their phylum colour (`hierarchical`).
`Other` and `Unassigned` are always grey.

## Development

See [CLAUDE.md](CLAUDE.md) for conventions.

```bash
Rscript -e 'devtools::test()'
Rscript -e 'devtools::check()'
```

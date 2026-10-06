# gmaio

Downstream analysis of [genome_metagenome_all_in_one](https://github.com/julianzaugg/genome_metagenome_all_in_one) pipeline output for metagenomes and isolate genomes.

`gmaio` is an R package plus project templates.
The package holds all the logic: readers for every pipeline output, sample linking, a common feature-profile data model, palettes, tables, figures and statistics.
Each project gets a `config.yml`, a `code/main.R` entry point that does the bulk processing, and one script per analysis.

## Install

```r
install.packages("remotes")
remotes::install_github("julianzaugg/genome_metagenome_all_in_one_analysis", build_vignettes = TRUE)
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
  config.yml              pipeline results directory, metadata, groups, colours
  data/                   put metadata.xlsx (or .csv/.tsv) here
  code/fetch_results.R    copies the pipeline outputs gmaio reads from a server
  code/main.R             bulk processing, run first
  code/*.R                one script per analysis
```

Edit `config.yml`, add the metadata, get the pipeline outputs (below), then from the project directory:

```bash
Rscript code/main.R
```

```bash
Rscript code/barcharts.R
```

`main.R` writes the standard tables to `Result_tables/`, the colour assignments to `palettes.yml`, and the processed project to `Result_other/processed.rds`.
Analysis scripts read that file, so they can be run in any order and rerun individually.

## Getting the pipeline outputs

gmaio reads only small summary tables from the pipeline `--outdir`, not reads, BAMs, assemblies or bins.
If the pipeline ran on the same machine, set `pipeline_results` in `config.yml` to its results directory.

If it ran on a server, set both paths in `config.yml`:

```yaml
pipeline_results: pipeline_results                        # local copy, inside the project
pipeline_remote: user@server:/path/to/pipeline/results    # where the pipeline wrote its output
```

Then, from the project directory in a terminal:

```bash
Rscript code/fetch_results.R
```

This copies only the files gmaio reads, keeping the pipeline's directory layout, and reports which outputs were found.
It needs rsync on both machines and ssh access to the server; run it from a terminal rather than RStudio if ssh asks for a password.
A run of a few hundred GB usually transfers tens to a few hundred MB.

The pipeline does not have to be finished.
Outputs that do not exist yet, or steps that were skipped, are reported and skipped by `main.R`, and analysis scripts without their inputs stop with a message.
Rerun `fetch_results.R` and `main.R` as the run progresses; only new or changed files are copied.
Some outputs, such as `read_stat_report.tsv` and `software_versions.tsv`, are only written when the run finishes.

To copy with rsync directly instead:

```r
gmaio::write_rsync_filter("gmaio_filter.txt", mode = "metagenome")
```

```bash
rsync -av --prune-empty-dirs --include-from=gmaio_filter.txt user@server:/path/to/results/ pipeline_results/
```

`gmaio::check_inputs(gmaio::read_config())` reports what is available at any time; `main.R` runs it first.

## Analysis scripts

Metagenome projects (`mode: metagenome`):

| Script | What it does |
|---|---|
| `main.R` | Loads read statistics, sylph, SingleM, MAGs and reference genomes (CheckM, GTDB-Tk, CoverM, DRAM distillate), gene catalogue functions (DRAM + RPKM, base and expanded catalogue), Nonpareil, GenomeSPOT, geNomad/CheckV, the marker gene tree, strain comparison (inStrain, TRACS) and comparison samples; writes the standard tables and palettes |
| `barcharts.R` | Stacked barcharts of the top taxa per sample for each profile and rank |
| `heatmaps.R` | Abundance heatmaps (log10 relative abundance) with phylum annotation |
| `ordination.R` | rclr PCA and Jaccard PCoA of taxonomic, MAG and functional profiles, PERMANOVA, PERMDISP, loadings |
| `diversity.R` | Alpha diversity (richness, Shannon, Simpson) with group tests |
| `differential_abundance.R` | MaAsLin3, LinDA and sPLS-DA with a consensus table, effect plots and boxplots |
| `nonpareil.R` | Nonpareil coverage curves and metrics by group |
| `genomespot.R` | GenomeSPOT trait predictions for HQ MAGs and abundance-weighted community traits |
| `mag_summary.R` | MAG quality, yield per sample, representative taxonomy, how much of each sample its own MAGs and the pooled catalogue explain, DRAM module heatmaps, marker gene trees |
| `functional_profiles.R` | KEGG, CAZy and peptidase heatmaps, CAZy substrate barcharts, for the study's gene catalogue or the expanded one |
| `comparison.R` | The study's samples next to external comparison samples: combined profiles, ordinations with PERMANOVA and PERMDISP, distances from each study sample to each group, composition, alpha diversity and detection of the study's MAGs |
| `strains.R` | Strains shared between samples (inStrain popANI), sharing within versus between groups, TRACS SNP distances and per-sample strain diversity |
| `mobile_elements.R` | Viral and plasmid cluster prevalence, richness, richness against depth, Jaccard PCoA, virus taxonomy |

Isolate projects (`mode: isolate`):

| Script | What it does |
|---|---|
| `main.R` | Loads CheckM2, GTDB-Tk, Bakta, MLST, AMRFinderPlus, ISEScan, comparison groups (Panaroo, fastANI, chewBBACA, IQ-TREE), geNomad/CheckV; writes the genome summary and comparison tables |
| `genome_summary.R` | Assembly, quality and annotation summary figure |
| `ani.R` | ANI heatmaps and PCoA per comparison group |
| `pangenome.R` | Gene categories, accessory gene heatmap and PCoA, gene associations with the group variable |
| `amr_mobile_elements.R` | AMR gene presence, AMR class and IS family counts, virus taxonomy, AMR gene associations |
| `cgmlst.R` | cgMLST allele distance heatmaps |
| `phylogeny.R` | Trees with tips coloured by the group variable |

Each script has a short parameter block at the top (profiles, ranks, thresholds) that can be edited per project.
The scripts set only the plot options that differ from the defaults; every option is listed on the function's help page (`?gmaio::plot_heatmap`).

## Documentation

The vignettes are worked examples on synthetic example projects (`gmaio::example_project()`, and `mode = "isolate"`):

| Vignette | Covers |
|---|---|
| `vignette("gmaio", package = "gmaio")` | How a project fits together: the processed project, profiles, palettes, saving, which functions each script uses |
| `vignette("outputs", package = "gmaio")` | Every table and figure the scripts write, and the MAG genome sets (`derep_bins`, `hq_bins`, `hq_derep_bins`, `ws_*`, ...) and which to use |
| `vignette("heatmaps", package = "gmaio")` | `plot_heatmap()`: choosing the data, annotations and splits, transforms, colours and legend breaks, ordering, cell values, sizes and fonts |
| `vignette("ordination", package = "gmaio")` | `run_ordination()`, PERMANOVA (overall and pairwise) and PERMDISP, and `plot_ordination()`: shapes, ellipses, hulls, spiders, loadings, envfit, combining panels |
| `vignette("barcharts", package = "gmaio")` | `plot_stacked_barchart()`: which taxa are shown, annotation strips, legend and bar options; where colours come from and how to change them |
| `vignette("diversity", package = "gmaio")` | `alpha_diversity()` and `plot_alpha_diversity()`: measures, rarefying, group and Dunn pairwise tests with brackets, figure options |
| `vignette("differential_abundance", package = "gmaio")` | MaAsLin3, LinDA and sPLS-DA, the consensus, effect, heatmap and boxplot figures, covariates and random effects |
| `vignette("comparison", package = "gmaio")` | Comparison samples: their metadata, combined profiles and metadata, ordinations, distances to each group and MAG detection |
| `vignette("strains", package = "gmaio")` | Strain sharing (inStrain), within versus between group tests, TRACS distances and strain diversity |
| `vignette("isolates", package = "gmaio")` | Isolate mode: genome summary, AMR genes and insertion sequences, ANI, cgMLST, pangenome and trees |

Each workbook gmaio writes also has an `About` sheet describing its sheets and columns, and `Result_tables/README.md` and `Result_figures/README.md` list the files present.

`browseVignettes("gmaio")` lists them.
They are only installed when the package is installed with `build_vignettes = TRUE`; without it, `vignette("gmaio", package = "gmaio")` reports that the vignette is not found.
The differential abundance vignette needs maaslin3, MicrobiomeStat and mixOmics installed to show its results.

## Working with the package directly

```r
project.l <- gmaio::load_processed()
genus.p <- gmaio::aggregate_profile(gmaio::get_profile(project.l, "singlem_relative"), rank = "genus")
genus.p <- gmaio::filter_profile(genus.p, min_abundance = 0.5, min_prevalence = 3)
ordination <- gmaio::run_ordination(genus.p, "pca_rclr")
gmaio::run_permanova(ordination$distance, gmaio::analysis_metadata(project.l), "Treatment")
```

Every abundance or function table is a `gm_profile`: a features x samples matrix plus feature annotations, so the same functions work for sylph, SingleM, MAGs, DRAM functions, AMR genes and pangenomes.
`names(project.l$profiles)` lists what was loaded.

## MAG sets and read statistics

MAG abundance profiles are named `mags_<set>_<type>`.
The cross-sample sets (`derep_bins`, `hq_bins`, `hq_derep_bins`, `hq_ref_bins`) map every sample to one catalogue of genomes pooled across samples; `hq_derep_bins` is the one the analysis scripts use.
With the pipeline's `--within_sample_dereplication`, gmaio also loads `ws_derep_bins` and `ws_hq_bins`, where each sample is mapped only to its own bins.
A genome can then only be detected in the sample it came from, so these are for per-sample summaries (`within_sample_mag_summary()`, `plot_mag_mapping()`), not ordination or differential abundance; those functions warn if given one.
`<type>` is `relative_abundance` (each sample sums to 100), `coverm_relative_abundance` (CoverM's value: % of the reads after QC and host removal, not rescaled), `coverage`, `read_count` or `covered_fraction` (share of each genome covered by reads, for detection).

Genomes from the pipeline's `--reference_genomes` are kept apart from the MAGs.
`project.l$tables$bin_summary` (the `Bins` sheet of `MAG_summary.xlsx`) has the MAGs only, and `project.l$tables$reference_genomes` (the `Reference_genomes` sheet) the references, with their CheckM2 quality from the pipeline's reference report and their GTDB-Tk taxonomy.
A genome is a reference when it is in that report, or when GTDB-Tk or a dereplication cluster lists it but no MAG CheckM report does.
References are never assigned to a sample, get `REF` short IDs, and are marked `Reference` in the `Genome_type` column of the `hq_ref_bins` profiles; MAG figures and counts leave them out.

`Metadata_and_read_stats.xlsx` has the read funnel in `Read_stats`: raw reads, each QC step, and reads mapped to each genome set, all as % of raw reads.
For Nanopore runs it adds `Bases_mapped_<set>_percent`, the full length of mapped reads as % of raw sequenced bases (the pipeline's metric A).
`Read_stats_full` is the pipeline report unchanged (including metrics B and C), and `Column_descriptions` explains every column.
Mapping percentages against raw reads include host reads in the denominator, so they can be far lower than CoverM's percentages when host removal discards most reads.

## Metadata and sample linking

The metadata must have one row per sample, with a column matching the pipeline samplesheet `sample` column (`metadata: sample_id_column`).
Tools label samples differently (`S1_R1`, `S1.clean_1.fastq.gz`, `S1.metabat2.5`, `S1___NODE_1_2`); gmaio maps them to metadata IDs by exact match, then by the longest ID followed by `.`, `_` or `-`.
Any pipeline sample that cannot be matched is an error.
Samples can be excluded from analyses (but kept in tables) with `exclude_column` or `exclude_samples`.

Comparison samples (the pipeline's `--comparison_reads`, e.g. a published survey) are linked to their own metadata, set in the `comparison:` section of `config.yml` (file, sheet, ID and label columns, `exclude_samples`, a `group_variable` and the two dataset labels).
The metadata can list more samples than were sequenced; only those in the pipeline outputs are used, and a comparison sample missing from it is an error.
`combined_metadata()` and `combined_profile()` put both sets of samples side by side, with `Dataset` and `Comparison_group` (the study's primary group variable, then the comparison groups).
A comparison group or label that repeats a study one is suffixed, so groups from the two datasets are never pooled.

## Pipeline outputs

Outputs are found in `pipeline_results` by directory name suffix (`*_sylph`, `*_checkm2`, ...), so Illumina, Nanopore and isolate numbering all work.
`gmaio::pipeline_registry()` lists every output gmaio reads.
Any entry can be overridden in `config.yml` under `files:` with a file, glob or directory, for example when outputs have been copied elsewhere.
Outputs that were not produced (skipped pipeline steps or an unfinished run) are reported by `gmaio::check_inputs()` and skipped.
When several files compete for one output, `check_inputs()` lists them with their dates; files from an older run into the same `--outdir` (for example per-bin DRAM `distilled/` folders from an earlier pipeline version) are the usual cause, and should be removed.
gmaio only reads the DRAM distillate of all bins together (`<n>_dram_bins/all_bins/distilled/`).

The expanded gene catalogue (`13_dram_expanded`, `23_rpkm_expanded`) adds the genes of the reference genomes (with `--reference_genomes_in_catalogue`) and/or comparison assemblies to the study's; its function profiles are named `functions_expanded_*`.
Comparison reads are mapped to the expanded catalogue when there is one, so `combined_profile(project.l, "functions_ko")` joins them with the study's expanded profile.

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

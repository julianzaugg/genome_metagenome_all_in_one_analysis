#' Genome sets used for MAG abundance
#'
#' The pipeline maps reads to up to six sets of genomes, which differ in which bins they
#' contain and whether genomes are pooled across samples. Each set gives one
#' `MAG_<Set>_abundances.xlsx` and one `MAG_<Set>_taxonomy_relative_abundance.xlsx` file,
#' one `mags_<Set>_*` family of profiles, and one column pair in the read statistics.
#'
#' @return Data frame with `Set` (as in file and profile names), `Name` (as in figures),
#'   `Scope` (`"All samples"` or `"Within sample"`), `Genomes` (how the set is built),
#'   `Use_for`, `Pipeline_clusters` and `Pipeline_mapping` (pipeline output folders) and
#'   `Read_stats_label` (the label in `Reads_mapped_<label>_percent` columns).
#' @export
mag_set_catalogue <- function(){
  # nolint start: line_length_linter. Descriptions are prose, one sentence per string.
  hq.s <- "HQ means completeness - 3 x contamination >= 50 in CheckM2 or CheckM1 (the pipeline's --hq_quality_source)."
  sets.df <- data.frame(
    Set = c("derep_bins", "hq_bins", "hq_derep_bins", "hq_ref_bins", "ws_derep_bins", "ws_hq_bins"),
    Scope = c(rep("All samples", 4), rep("Within sample", 2)),
    Genomes = c(
      "All bins of any quality from all samples, dereplicated together at 95% ANI; one representative per cluster.",
      paste("The HQ representatives of the derep_bins catalogue: bins are dereplicated first, then only representatives",
            "that are HQ are kept. A cluster whose chosen representative is not HQ contributes no genome, even if it",
            "contains HQ bins.", hq.s),
      paste("HQ bins from all samples, dereplicated together at 95% ANI, so every cluster is represented by an HQ",
            "genome.", hq.s),
      "HQ bins from all samples dereplicated together with the external --reference_genomes at 95% ANI.",
      "Each sample's bins of any quality, dereplicated within that sample only (pipeline --within_sample_dereplication).",
      "Each sample's HQ bins, dereplicated within that sample only (pipeline --within_sample_dereplication)."),
    Use_for = c(
      "The fullest picture of what was assembled; low-quality genomes are included, so read their taxonomy and functions with care.",
      "Comparing with older results; prefer hq_derep_bins, which keeps every HQ cluster (the two are often nearly the same).",
      "Comparing samples with MAGs: ordination, differential abundance, barcharts and heatmaps (the template default).",
      "As hq_derep_bins, when reference genomes should be counted alongside the MAGs.",
      "How much of each sample its own genomes explain; each sample is mapped only to its own genomes, so do not compare samples with it.",
      "As ws_derep_bins, for each sample's own HQ MAGs; not for comparing samples."),
    Pipeline_clusters = c("11_dereplicated_bins", "11_dereplicated_bins/high_quality_representatives",
                          "11_dereplicated_hq_bins", "11_dereplicated_hq_ref_bins",
                          "11_within_sample_dereplicated_bins/<sample>", "11_within_sample_dereplicated_hq_bins/<sample>"),
    Pipeline_mapping = c("12_coverm_bins", "12_coverm_hq_bins", "12_coverm_hq_derep_bins", "12_coverm_hq_ref_bins",
                         "12_coverm_within_sample_derep_bins", "12_coverm_within_sample_hq_bins"),
    stringsAsFactors = FALSE)
  # nolint end
  sets.df$Name <- unname(mag_set_descriptions()[sets.df$Set])
  sets.df$Read_stats_label <- unname(mag_sets()[sets.df$Set])
  sets.df[, c("Set", "Name", "Scope", "Genomes", "Use_for", "Pipeline_clusters", "Pipeline_mapping", "Read_stats_label")]
}

# One catalogue entry: file name (display and pattern), mode, writer, description, sheets and key columns
output_entry <- function(file, pattern, mode, written_by, description, sheets, columns = NULL){
  sheets.df <- data.frame(Sheet = names(sheets), Description = unname(sheets), stringsAsFactors = FALSE)
  sheets.df$Pattern <- paste0("^", gsub("<[^>]+>", ".+", gsub("([.()])", "\\\\\\1", sheets.df$Sheet)), "$")
  columns.df <- if (is.null(columns)) data.frame(Sheet = character(), Column = character(), Description = character()) else
    do.call(rbind, lapply(names(columns), function(sheet.s){
      data.frame(Sheet = sheet.s, Column = names(columns[[sheet.s]]), Description = unname(columns[[sheet.s]]),
                 stringsAsFactors = FALSE)
    }))
  list(File = file, Pattern = pattern, Mode = mode, Written_by = written_by, Description = description,
       Sheets = sheets.df, Columns = columns.df)
}

output_entries <- function(){
  # nolint start: line_length_linter. Descriptions are prose, one sentence per string.
  sets.s <- paste(names(mag_sets()), collapse = "|")
  samples.s <- "One column per sample, named by Sample_ID (the Metadata sheet of Metadata_and_read_stats.xlsx has the labels)."
  lineage.s <- "The full GTDB lineage (d__ to the rank), used as the feature ID."
  ranks.l <- c(phylum = "Summed to phylum.", class = "Summed to class.", order = "Summed to order.",
               family = "Summed to family.", genus = "Summed to genus.", species = "Summed to species.")
  rank_columns.l <- c(`<rank>` = lineage.s, Label = "Short display name: the lowest named rank with its phylum.",
                      `<Sample_ID>` = samples.s)
  bin_columns.v <- c(
    Bin_ID = "Bin name from the pipeline (sample, binner and bin number).",
    Sample_ID = "Sample the bin was assembled from.",
    `Completeness_<CheckM>` = "Estimated completeness (%) from CheckM2 or CheckM1.",
    `Contamination_<CheckM>` = "Estimated contamination (%) from CheckM2 or CheckM1.",
    `Quality_score_<CheckM>` = "Completeness - quality_weight x contamination (config mags:).",
    `Passes_HQ_<CheckM>` = "Whether the quality score reaches quality_threshold (config mags:).",
    High_quality = "HQ by the rule above, in the CheckM report(s) chosen by mags: quality_source.",
    MIMAG_tier = paste("High (> 90% complete, < 5% contaminated), Medium (>= 50%, < 10%) or Low, from CheckM2 (else CheckM1).",
                       "rRNA and tRNA genes are not checked, so this approximates the MIMAG standard; empty when CheckM could not assess the bin."),
    Classification = "GTDB-Tk classification; Domain to Species are its ranks.",
    closest_genome_reference = "GTDB-Tk closest reference genome; closest_genome_ani is its ANI (%).",
    `Representative_<Set>` = "The representative of the bin's cluster in that genome set (see MAG sets).",
    `Is_representative_<Set>` = "TRUE when the bin is itself the representative, i.e. it is in that set.",
    Short_ID = "Short name (MAG001, ...) used in figures.",
    Label = "Figure label: lowest named taxon and Short_ID.")
  test_columns.v <- c(Test = "Wilcoxon rank-sum (two groups) or Kruskal-Wallis (more).", Statistic = "Test statistic.",
                      P_value = "P value.", P_adjusted = "BH-adjusted across the measures in the sheet.",
                      N = "Samples tested.")
  pairwise_columns.v <- c(`Group_1, Group_2` = "The pair of groups compared.", `N_1, N_2` = "Samples per group.",
                          `Mean_1, Mean_2, Median_1, Median_2` = "Group means and medians.",
                          Test = "Dunn (after Kruskal-Wallis) or Wilcoxon.",
                          Statistic = "z for Dunn (positive when Group_2 ranks higher), W for Wilcoxon.",
                          P_adjusted = "BH-adjusted within each measure.", Significance = "*** < 0.001, ** < 0.01, * < 0.05, ns.")
  da_columns.v <- c(Method = "MaAsLin3, LinDA or sPLS-DA.",
                    Model = "MaAsLin3 abundance or prevalence model; sPLS-DA component.",
                    `Contrast, Reference` = "Group compared, and the reference group it is compared with.",
                    `Effect, Effect_type, SE` = "Effect size and its kind: log2 coefficient (MaAsLin3), log2 fold change (LinDA), loading (sPLS-DA).",
                    `P_value, Q_value` = "Raw and FDR-adjusted p values.",
                    Stability = "How often sPLS-DA selected the feature in cross-validation.",
                    Enriched_in = "Group in which the feature is higher.")

  list(
    output_entry("Metadata_and_read_stats.xlsx", "^Metadata_and_read_stats\\.xlsx$", "metagenome, isolate", "code/main.R",
      "The metadata as linked to the pipeline samples, read counts through QC and host removal, how many reads map to each genome set, and the pipeline run parameters.",
      c(Metadata = "Metadata with the Sample_ID, Sample_label and Excluded columns gmaio adds.",
        Read_stats = "Reads (and, for long reads, bases) through each QC step, and the share mapped to each genome set.",
        Read_stats_full = "The pipeline read statistics report, unchanged.",
        Column_descriptions = "What each Read_stats and Read_stats_full column means.",
        Run_parameters = "Pipeline parameters of the run.")),
    output_entry("Sylph_taxonomic_abundance.xlsx", "^Sylph_taxonomic_abundance\\.xlsx$", "metagenome", "code/main.R",
      "sylph taxonomic abundance (%): the share of genomes (cells) of each taxon, from reads matched to GTDB genomes. Each sample sums to 100.",
      c(genome = "One row per sylph genome (GTDB species representative).", ranks.l), list(`<rank>` = rank_columns.l)),
    output_entry("Sylph_sequence_abundance.xlsx", "^Sylph_sequence_abundance\\.xlsx$", "metagenome", "code/main.R",
      "sylph sequence abundance (%): the share of reads of each taxon; larger genomes get more reads per cell than in the taxonomic abundance.",
      c(genome = "One row per sylph genome (GTDB species representative).", ranks.l), list(`<rank>` = rank_columns.l)),
    output_entry("SingleM_coverage.xlsx", "^SingleM_coverage\\.xlsx$", "metagenome", "code/main.R",
      "SingleM coverage: estimated read depth (x) of each taxon from single-copy marker genes; not limited to taxa with reference genomes.",
      c(lineage = "One row per SingleM lineage.", ranks.l), list(`<rank>` = rank_columns.l)),
    output_entry("SingleM_relative_abundance.xlsx", "^SingleM_relative_abundance\\.xlsx$", "metagenome", "code/main.R",
      "SingleM coverage as a share (%) of the prokaryotic community; each sample sums to 100.",
      c(lineage = "One row per SingleM lineage.", ranks.l), list(`<rank>` = rank_columns.l)),
    output_entry("SingleM_scaled_abundance.xlsx", "^SingleM_scaled_abundance\\.xlsx$", "metagenome", "code/main.R",
      "SingleM relative abundance scaled by the prokaryotic fraction of reads: % of all reads, so samples with more host or eukaryotic DNA sum to less than 100.",
      c(lineage = "One row per SingleM lineage.", ranks.l,
        Prokaryotic_fraction = "SingleM prokaryotic fraction: read_fraction is the % of reads from bacteria and archaea in each sample."),
      list(`<rank>` = rank_columns.l)),
    output_entry("MAG_summary.xlsx", "^MAG_summary\\.xlsx$", "metagenome", "code/main.R",
      "Every bin from every sample with its quality, taxonomy and membership of each genome set. Start here to look up a MAG.",
      c(Bins = "One row per bin.", `Clusters_<Set>` = "Dereplication clusters of a genome set: each bin and its cluster's representative."),
      list(Bins = bin_columns.v,
           `Clusters_<Set>` = c(Representative = "Bin chosen to represent the cluster.", Bin_ID = "A bin in the cluster."))),
    output_entry("MAG_<Set>_abundances.xlsx", paste0("^MAG_(", sets.s, ")_abundances\\.xlsx$"), "metagenome", "code/main.R",
      "Abundance of each genome of one genome set (see MAG sets) in each sample, from the pipeline's CoverM mapping of all reads to that set.",
      c(relative_abundance = "Mean coverage as a share (%) of all genomes in the set; each sample sums to 100. Use this to compare composition.",
        coverm_relative_abundance = "CoverM relative abundance: % of the sample's reads (after QC and host removal) mapped to the genome. The rest, unmapped reads, is not shown, so samples sum to less than 100.",
        coverage = "Mean read depth over the genome (x).",
        read_count = "Reads mapped to the genome."),
      list(`All sheets` = c(genome = "Bin_ID of the genome (see MAG_summary.xlsx).", Short_ID = "Short name used in figures.",
                         Label = "Figure label.", `<Sample_ID>` = samples.s))),
    output_entry("MAG_<Set>_taxonomy_relative_abundance.xlsx", paste0("^MAG_(", sets.s, ")_taxonomy_relative_abundance\\.xlsx$"),
      "metagenome", "code/main.R",
      "The relative_abundance sheet of MAG_<Set>_abundances.xlsx summed to each GTDB rank; each sample sums to 100.",
      ranks.l, list(`<rank>` = rank_columns.l)),
    output_entry("MAG_per_sample_recovery.xlsx", "^MAG_per_sample_recovery\\.xlsx$", "metagenome", "code/mag_summary.R",
      "How many MAGs each sample yielded and how much of each sample's reads each genome set explains.",
      c(Samples = "One row per sample: bins and HQ bins assembled from it, and its own within-sample representatives.",
        Mapping = "One row per sample and genome set: the share of the sample's reads (or bases) mapped to the set."),
      list(Samples = c(Bins = "Bins assembled from the sample.", HQ_bins = "Of these, HQ.",
                       Own_representatives = "Its within-sample dereplicated genomes (ws_derep_bins).",
                       Own_HQ_MAGs = "Its within-sample HQ genomes (ws_hq_bins)."),
           Mapping = c(Set = "Genome set (see MAG sets).", Genome_set = "Its name in figures.",
                       Scope = "All samples (pooled catalogue) or Within sample (own genomes only).",
                       Input_reads_mapped_percent = "% of reads after QC and host removal mapped to the set (CoverM).",
                       Reads_mapped_percent = "% of raw reads mapped (pipeline read statistics).",
                       Bases_mapped_percent = "% of raw bases mapped (long reads; metric A of the read statistics)."))),
    output_entry("Gene_catalogue_functions_normalised_rpkm.xlsx", "^Gene_catalogue_functions_normalised_rpkm\\.xlsx$", "metagenome",
      "code/main.R",
      "Functions of the gene catalogue (all assembled genes, clustered across samples) as normalised RPKM, summed per function.",
      c(KEGG_KO = "KEGG orthologues.", CAZy_family = "CAZyme families.", CAZy_substrate = "CAZyme substrates.",
        Peptidase = "MEROPS peptidase families.", Annotation_coverage = "Share of genes and RPKM with each annotation type.",
        CAZy_substrate_map = "Which CAZy families count towards each substrate.")),
    output_entry("Nonpareil_summary.xlsx", "^Nonpareil_summary\\.xlsx$", "metagenome", "code/main.R",
      "Nonpareil estimates of sequencing coverage and diversity per sample.",
      c(Summary = "Coverage at the sequenced depth (C), Nonpareil diversity (diversity) and the sequencing effort for 95% coverage (LRstar).")),
    output_entry("GenomeSPOT_predictions.xlsx", "^GenomeSPOT_predictions\\.xlsx$", "metagenome, isolate", "code/main.R",
      "GenomeSPOT predictions of oxygen, temperature, salinity and pH preferences for each genome.",
      c(Predictions = "One row per genome, one column per trait.", Long = "The same, one row per genome and trait.")),
    output_entry("Mobile_elements.xlsx", "^Mobile_elements\\.xlsx$", "metagenome, isolate", "code/main.R",
      "Viruses and plasmids found by geNomad, their CheckV quality and their clusters across samples.",
      c(Virus_summary = "geNomad virus predictions.", Plasmid_summary = "geNomad plasmid predictions.",
        CheckV_quality = "CheckV completeness and quality of the viruses.",
        Virus_clusters = "Virus clusters across samples.", Plasmid_clusters = "Plasmid clusters across samples.")),
    output_entry("Alpha_diversity.xlsx", "^Alpha_diversity\\.xlsx$", "metagenome", "code/diversity.R",
      "Alpha diversity of each dataset (profile and rank) and tests between groups of the primary group variable.",
      c(`<dataset>_values` = "Richness, Shannon, Simpson (Gini-Simpson) and Pielou evenness per sample.",
        `<dataset>_tests` = "Group test per measure.",
        `<dataset>_pairwise` = "Pairwise tests between groups, when there are more than two."),
      list(`<dataset>_tests` = test_columns.v, `<dataset>_pairwise` = pairwise_columns.v)),
    output_entry("Ordination_statistics.xlsx", "^Ordination_statistics\\.xlsx$", "metagenome", "code/ordination.R",
      "PERMANOVA and PERMDISP tests of the primary group variable on each ordination (dataset and method, in the Label column).",
      c(PERMANOVA = "Overall test: how much of the variation between samples the groups explain (R2).",
        PERMANOVA_pairwise = "PERMANOVA for each pair of groups on just those samples, BH adjusted (more than two groups).",
        PERMDISP = "Whether groups differ in spread; a difference can by itself make PERMANOVA significant.",
        PERMDISP_pairwise = "PERMDISP for each pair of groups."),
      list(PERMANOVA = c(Term = "Model term.", R2 = "Share of variation explained.", F = "Pseudo-F.",
                         P_value = "Permutation p value.", N = "Samples.", Permutations = "Permutations run."),
           PERMANOVA_pairwise = c(`Group_1, Group_2` = "The pair.", `N_1, N_2` = "Samples per group.",
                                  P_adjusted = "BH-adjusted across the pairs."))),
    output_entry("Differential_abundance_all_results.xlsx", "^Differential_abundance_all_results\\.xlsx$", "metagenome",
      "code/differential_abundance.R",
      "Every feature from every method (MaAsLin3, LinDA, sPLS-DA), significant or not, one sheet per dataset.",
      c(`<dataset>` = "All results for a profile and rank."), list(`<dataset>` = da_columns.v)),
    output_entry("Differential_abundance_consensus.xlsx", "^Differential_abundance_consensus\\.xlsx$", "metagenome",
      "code/differential_abundance.R",
      "Features called significant by at least one method, with the methods that agree; report those found by two or more.",
      c(`<dataset>` = "Consensus for a profile and rank."),
      list(`<dataset>` = c(Enriched_in = "Group in which the feature is higher.",
                           `<Method>_effect` = "Effect from each method (empty when that method did not call it).",
                           `MaAsLin3_q, LinDA_q, sPLSDA_stability` = "Evidence from each method.",
                           N_methods = "Methods calling the feature higher in this group.",
                           Conflicting_direction = "Methods call the feature higher in different groups."))),
    output_entry("Nonpareil_statistics.xlsx", "^Nonpareil_statistics\\.xlsx$", "metagenome", "code/nonpareil.R",
      "Nonpareil estimates per sample and tests between groups.",
      c(Summary = "Nonpareil estimates per sample.", Group_tests = "Group test per metric."),
      list(Group_tests = test_columns.v)),
    output_entry("Mobile_elements_analysis.xlsx", "^Mobile_elements_analysis\\.xlsx$", "metagenome", "code/mobile_elements.R",
      "Virus and plasmid cluster richness and prevalence per sample, with tests between groups.",
      c(`<type>_richness` = "Clusters per sample.", `<type>_prevalence` = "Samples each cluster is found in.",
        `<type>_depth_correlation` = "Correlation of richness with sequencing depth.",
        `<type>_richness_tests` = "Group test of richness.", `<type>_permanova` = "PERMANOVA on cluster presence (Jaccard).",
        `<type>_permanova_pairwise` = "Pairwise PERMANOVA (more than two groups).")),
    output_entry("GenomeSPOT_analysis.xlsx", "^GenomeSPOT_analysis\\.xlsx$", "metagenome", "code/genomespot.R",
      "GenomeSPOT traits of the HQ MAGs and community traits weighted by MAG abundance.",
      c(HQ_MAG_predictions = "Predicted traits of each HQ MAG.",
        Community_weighted = "Abundance-weighted mean of each trait per sample.", Community_tests = "Group test per trait.")),
    output_entry("Genome_summary.xlsx", "^Genome_summary\\.xlsx$", "isolate", "code/main.R",
      "Quality, taxonomy, typing and annotation of each isolate.",
      c(Summary = "One row per isolate: assembly and quality metrics, species, MLST type, AMR and IS counts.",
        CheckM2 = "CheckM2 report.", GTDBTk = "GTDB-Tk classification.", Bakta = "Bakta annotation statistics.",
        MLST = "MLST scheme, sequence type and alleles.", AMRFinderPlus = "AMR, stress and virulence genes found.",
        ISEScan = "Insertion sequences found.")),
    output_entry("AMR_and_IS_profiles.xlsx", "^AMR_and_IS_profiles\\.xlsx$", "isolate", "code/main.R",
      "Isolate x feature tables of AMR genes and insertion sequences.",
      c(amr_genes = "Presence (1) or absence (0) of each AMR gene.", amr_classes = "AMR genes per drug class.",
        is_families = "Insertion sequences per family.")),
    output_entry("Comparison_<group>.xlsx", "^Comparison_.+\\.xlsx$", "isolate", "code/main.R",
      "Pangenome and ANI of one comparison group (isolates plus reference genomes).",
      c(Entries = "Genomes in the group: samples and references.", Gene_categories = "Gene clusters per pangenome category.",
        Presence_absence = "Panaroo gene clusters: presence per genome, prevalence and category (core, soft core, shell, cloud).",
        fastANI = "Pairwise ANI.")),
    output_entry("AMR_associations.xlsx", "^AMR_associations\\.xlsx$", "isolate", "code/amr_mobile_elements.R",
      "Fisher's exact tests of each AMR gene against the primary group variable.",
      c(AMR_gene_associations = "Presence per group, odds ratio (two groups), p and BH q values.")),
    output_entry("Pangenome_associations.xlsx", "^Pangenome_associations\\.xlsx$", "isolate", "code/pangenome.R",
      "Fisher's exact tests of each accessory gene against the primary group variable.",
      c(`<group>_associations` = "One sheet per comparison group.")),
    output_entry("cgMLST_distances.xlsx", "^cgMLST_distances\\.xlsx$", "isolate", "code/cgmlst.R",
      "Allele differences between genomes from chewBBACA cgMLST, ignoring loci missing in either genome.",
      c(`<group>_<schema>` = "Distance matrix for a comparison group and cgMLST threshold (e.g. cgMLST95: loci in 95% of genomes)."))
  )
}

  # nolint end
#' Table files written by gmaio
#'
#' Lists every table file `main.R` and the analysis scripts write, with what it contains.
#' The same descriptions are written to an `About` sheet in each file and to `README.md` in
#' the table directory; `vignette("outputs", package = "gmaio")` explains them with examples.
#'
#' @param mode `"all"`, `"metagenome"` or `"isolate"`.
#' @return Data frame with `File`, `Mode`, `Written_by` and `Description`.
#' @export
output_catalogue <- function(mode = c("all", "metagenome", "isolate")){
  mode <- match.arg(mode)
  entries.l <- output_entries()
  catalogue.df <- data.frame(File = vapply(entries.l, `[[`, character(1), "File"),
                             Mode = vapply(entries.l, `[[`, character(1), "Mode"),
                             Written_by = vapply(entries.l, `[[`, character(1), "Written_by"),
                             Description = vapply(entries.l, `[[`, character(1), "Description"), stringsAsFactors = FALSE)
  if (mode != "all") catalogue.df <- catalogue.df[grepl(mode, catalogue.df$Mode), , drop = FALSE]
  rownames(catalogue.df) <- NULL
  catalogue.df
}

#' Describe a table file written by gmaio
#'
#' @param file File name or path of a gmaio table (e.g. `"MAG_hq_derep_bins_abundances.xlsx"`).
#' @return `NULL` for files gmaio does not write; otherwise a list with `File`, `Description`,
#'   `Written_by`, `Sheets` and `Columns` (data frames) and, for MAG files, `Set` (a row of
#'   [mag_set_catalogue()]).
#' @export
describe_output <- function(file){
  file.s <- basename(file)
  entries.l <- output_entries()
  hit.v <- vapply(entries.l, function(e) grepl(e$Pattern, file.s), logical(1))
  if (!any(hit.v)) return(NULL)
  entry.l <- entries.l[[which(hit.v)[1]]]
  entry.l$Set <- NULL
  set.s <- regmatches(file.s, regexec(paste0("^MAG_(", paste(names(mag_sets()), collapse = "|"), ")_"), file.s))[[1]][2]
  if (!is.na(set.s)){
    sets.df <- mag_set_catalogue()
    entry.l$Set <- sets.df[sets.df$Set == set.s, , drop = FALSE]
    entry.l$File <- file.s
  }
  entry.l
}

# The About sheet: what the file is, each sheet present, key columns and, for MAG files, the genome set
about_rows <- function(entry.l, sheet_names.v){
  sheets.df <- entry.l$Sheets
  present.v <- vapply(sheet_names.v, function(sheet.s){
    hit.v <- which(vapply(sheets.df$Pattern, function(p) grepl(p, sheet.s), logical(1)))
    if (length(hit.v) == 0) "" else sheets.df$Description[hit.v[1]]
  }, character(1))
  sections.l <- list(
    list(title = entry.l$File, rows = data.frame(A = c(entry.l$Description, paste("Written by", entry.l$Written_by)), B = "")),
    list(title = "Sheets", rows = data.frame(A = sheet_names.v, B = present.v)))
  if (!is.null(entry.l$Set)){
    set.df <- entry.l$Set
    sections.l <- c(sections.l, list(list(title = paste("Genome set:", set.df$Set), rows = data.frame(
      A = c("Name in figures", "Genomes", "Use for", "Pipeline folders"),
      B = c(set.df$Name, set.df$Genomes, set.df$Use_for, paste(set.df$Pipeline_clusters, "and", set.df$Pipeline_mapping))))))
  }
  if (nrow(entry.l$Columns) > 0){
    columns.df <- entry.l$Columns
    columns_rows.df <- data.frame(A = paste0(columns.df$Sheet, ": ", columns.df$Column), B = columns.df$Description)
    sections.l <- c(sections.l, list(list(title = "Columns", rows = columns_rows.df)))
  }
  sections.l <- c(sections.l, list(list(title = "More", rows = data.frame(
    A = "README.md in this folder lists every table; in R, vignette(\"outputs\", package = \"gmaio\") explains them.", B = ""))))
  sections.l
}

add_about_sheet <- function(workbook, entry.l, sheet_names.v){
  openxlsx::addWorksheet(workbook, "About")
  title_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12)
  wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  row.n <- 1
  for (section.l in about_rows(entry.l, sheet_names.v)){
    openxlsx::writeData(workbook, "About", section.l$title, startRow = row.n)
    openxlsx::addStyle(workbook, "About", title_style, rows = row.n, cols = 1)
    rows.df <- section.l$rows
    openxlsx::writeData(workbook, "About", rows.df, startRow = row.n + 1, colNames = FALSE)
    openxlsx::addStyle(workbook, "About", wrap_style, rows = row.n + seq_len(nrow(rows.df)), cols = 1:2, gridExpand = TRUE)
    row.n <- row.n + nrow(rows.df) + 2
  }
  openxlsx::setColWidths(workbook, "About", cols = 1:2, widths = c(45, 100))
  openxlsx::worksheetOrder(workbook) <- c(length(openxlsx::sheets(workbook)), seq_len(length(openxlsx::sheets(workbook)) - 1))
  invisible(workbook)
}

#' Write the README.md index of a table directory
#'
#' Lists every table file in the directory with what it contains and, when MAG files are
#' present, which genome set to use for what. Written automatically whenever gmaio writes
#' one of its tables, so it matches the files present.
#'
#' @param dir Table directory (e.g. `Result_tables`).
#' @return The README path, invisibly.
#' @export
write_output_index <- function(dir){
  # nolint start: line_length_linter. README text is prose, one sentence per string.
  files.v <- sort(list.files(dir, pattern = "\\.xlsx$"))
  files.v <- files.v[!startsWith(files.v, "~$")]
  described.l <- lapply(stats::setNames(files.v, files.v), describe_output)
  known.v <- files.v[!vapply(described.l, is.null, logical(1))]
  lines.v <- c("# Result tables", "",
               "Written by gmaio; this file is regenerated whenever a gmaio table is written.",
               "Each workbook also has an `About` sheet describing its sheets and columns.",
               "In R, `vignette(\"outputs\", package = \"gmaio\")` explains the files with examples.", "")
  if (any(grepl("^MAG_", known.v))){
    sets.df <- mag_set_catalogue()
    present.v <- sets.df$Set[vapply(sets.df$Set, function(set.s) any(startsWith(files.v, paste0("MAG_", set.s, "_"))),
                                    logical(1))]
    if (length(present.v) > 0) sets.df <- sets.df[sets.df$Set %in% present.v, , drop = FALSE]
    lines.v <- c(lines.v, "## Which MAG file to use", "",
                 "`MAG_<set>_abundances.xlsx` and `MAG_<set>_taxonomy_relative_abundance.xlsx` exist for each genome set the pipeline mapped reads to.",
                 "The sets differ in which bins they contain and whether genomes are pooled across samples:", "",
                 "| Set | In figures | Genomes | Use for |", "|---|---|---|---|",
                 sprintf("| `%s` | %s | %s | %s |", sets.df$Set, sets.df$Name, sets.df$Genomes, sets.df$Use_for), "",
                 "Within each abundance file, `relative_abundance` (genome copies, each sample sums to 100) is the one to compare composition with.", "")
  }
  lines.v <- c(lines.v, "## Files", "")
  for (file.s in known.v){
    entry.l <- described.l[[file.s]]
    lines.v <- c(lines.v, paste0("### ", file.s), "", entry.l$Description,
                 if (!is.null(entry.l$Set)) paste0("Genome set `", entry.l$Set$Set, "`: ", entry.l$Set$Genomes),
                 paste0("Written by `", entry.l$Written_by, "`."), "",
                 paste0("- `", entry.l$Sheets$Sheet, "`: ", entry.l$Sheets$Description), "")
  }
  other.v <- setdiff(files.v, known.v)
  if (length(other.v) > 0) lines.v <- c(lines.v, "## Other files", "", "Not written by gmaio:", "", paste0("- `", other.v, "`"), "")
  # nolint end
  path.s <- file.path(dir, "README.md")
  writeLines(lines.v, path.s)
  invisible(path.s)
}

#' Figure folders written by gmaio
#'
#' Each analysis script saves its figures to one folder of the project figure directory.
#' `save_plot()` keeps a `README.md` index of them in the figure directory.
#'
#' @param mode `"all"`, `"metagenome"` or `"isolate"`.
#' @return Data frame with `Folder`, `Mode`, `Written_by`, `Files` (file names, with
#'   `<dataset>` for a profile and rank such as `sylph_taxonomic__genus`) and `Description`.
#' @export
figure_catalogue <- function(mode = c("all", "metagenome", "isolate")){
  # nolint start: line_length_linter. Descriptions are prose, one sentence per string.
  mode <- match.arg(mode)
  row.f <- function(folder, mode, script, files, description){
    data.frame(Folder = folder, Mode = mode, Written_by = paste0("code/", script), Files = files, Description = description,
               stringsAsFactors = FALSE)
  }
  figures.df <- rbind(
    row.f("barcharts", "metagenome", "barcharts.R", "<dataset>",
          "Stacked barcharts of the most abundant taxa per sample, split by the primary group."),
    row.f("heatmaps", "metagenome", "heatmaps.R", "<dataset>", "Heatmaps of the most abundant taxa across samples."),
    row.f("ordination", "metagenome", "ordination.R", "<dataset>__pca_rclr, <dataset>__jaccard, <dataset>__pca_rclr_PC1_loadings, <dataset>__pca_rclr_PC2_loadings",
          "Ordinations (pca_rclr: PCA of robust centred log-ratios; jaccard: PCoA of presence/absence) with the PERMANOVA and PERMDISP result in the caption, and the taxa loading most on each axis."),
    row.f("diversity", "metagenome", "diversity.R", "<dataset>",
          "Alpha diversity by group with the group test and, for more than two groups, brackets for the pairs that differ."),
    row.f("differential_abundance", "metagenome", "differential_abundance.R",
          "<dataset>__consensus, <dataset>__heatmap, <dataset>__boxplots (or _page1, _page2, ...)",
          "Effects of features called by at least two methods, a heatmap of every feature called with a track per method, and boxplots of the consensus features."),
    row.f("mags", "metagenome", "mag_summary.R",
          "bin_quality_checkm2, bin_quality_checkm1, bins_per_sample, representative_taxonomy_phylum, representative_taxonomy_family, reads_mapped_to_mags, own_hq_mag_composition_phylum, own_hq_mag_composition_genus, dram_module_completeness, dram_function_presence",
          "MAG quality, bins per sample by quality tier, taxonomy of the representatives, how much of each sample its own MAGs and the pooled catalogue explain, and DRAM heatmaps."),
    row.f("functions", "metagenome", "functional_profiles.R", "<profile>__top_heatmap, cazy_substrates",
          "Heatmaps of the most abundant gene catalogue functions and CAZyme substrates per sample."),
    row.f("mobile_elements", "metagenome", "mobile_elements.R",
          "virus_prevalence, virus_richness, virus_richness_vs_depth, virus_jaccard_pcoa (and plasmid_ ...), virus_taxonomy",
          "Virus and plasmid cluster prevalence, richness by group and against sequencing depth, ordination and virus taxonomy."),
    row.f("nonpareil", "metagenome", "nonpareil.R", "nonpareil_curves_by_sample, nonpareil_curves_by_group, nonpareil_metrics",
          "Nonpareil coverage curves and coverage, diversity and required effort by group."),
    row.f("genomespot", "metagenome", "genomespot.R", "hq_mag_traits_by_phylum, hq_mag_oxygen_by_phylum, community_weighted_traits",
          "GenomeSPOT traits of the HQ MAGs by phylum, and abundance-weighted community traits by group."),
    row.f("genome_summary", "isolate", "genome_summary.R", "genome_summary", "Assembly and quality metrics of each isolate."),
    row.f("amr_mobile_elements", "isolate", "amr_mobile_elements.R", "amr_gene_presence, amr_classes_counts, is_families_counts, virus_taxonomy",
          "AMR gene presence by drug class, AMR and IS counts per isolate, and prophage taxonomy."),
    row.f("ani", "isolate", "ani.R", "<group>__ani_heatmap, <group>__ani_pcoa", "ANI between the genomes of each comparison group."),
    row.f("cgmlst", "isolate", "cgmlst.R", "<group>__<schema>", "cgMLST allele differences between genomes."),
    row.f("pangenome", "isolate", "pangenome.R", "<group>__gene_categories, <group>__accessory_heatmap, <group>__accessory_pcoa",
          "Pangenome categories, accessory gene presence and ordination."),
    row.f("phylogeny", "isolate", "phylogeny.R", "<group>__tree", "Phylogenetic tree of each comparison group."))
  # nolint end
  if (mode != "all") figures.df <- figures.df[figures.df$Mode == mode, , drop = FALSE]
  rownames(figures.df) <- NULL
  figures.df
}

#' Write the README.md index of a figure directory
#'
#' @param dir Figure directory (e.g. `Result_figures`).
#' @return The README path, invisibly.
#' @export
write_figure_index <- function(dir){
  # nolint start: line_length_linter. Descriptions are prose, one sentence per string.
  folders.v <- sort(list.dirs(dir, full.names = FALSE, recursive = FALSE))
  figures.df <- figure_catalogue()
  figures.df <- figures.df[!duplicated(figures.df$Folder), , drop = FALSE]
  known.v <- intersect(folders.v, figures.df$Folder)
  lines.v <- c("# Result figures", "",
               "Written by gmaio; this file is regenerated whenever a gmaio figure is saved.",
               "`<dataset>` is a profile and rank, e.g. `sylph_taxonomic__genus`, or a MAG profile such as `mags_hq_derep_bins_relative_abundance`.",
               "In R, `vignette(\"outputs\", package = \"gmaio\")` explains the figures and tables.", "")
  # nolint end
  for (folder.s in known.v){
    row.df <- figures.df[figures.df$Folder == folder.s, , drop = FALSE]
    lines.v <- c(lines.v, paste0("## ", folder.s, "/"), "", row.df$Description,
                 paste0("Written by `", row.df$Written_by, "`: ", row.df$Files, "."), "")
  }
  other.v <- setdiff(folders.v, known.v)
  if (length(other.v) > 0) lines.v <- c(lines.v, "## Other folders", "", paste0("- `", other.v, "/`"), "")
  path.s <- file.path(dir, "README.md")
  writeLines(lines.v, path.s)
  invisible(path.s)
}

# Builds the synthetic pipeline output used by the tests and examples:
# inst/extdata/metagenome_results (plus a metadata table). Headers follow the real
# pipeline outputs; values are made up. Run from the package root:
#   Rscript data-raw/make_fixtures.R

set.seed(42)
root.s <- "inst/extdata/metagenome"
unlink(root.s, recursive = TRUE)
results.s <- file.path(root.s, "results")
write_tsv <- function(x.df, ...){
  path.s <- file.path(results.s, ...)
  dir.create(dirname(path.s), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x.df, path.s, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

samples.v <- c("S1", "S2", "S3", "S10", "S11", "S12")
groups.v <- c("Control", "Control", "Control", "Treated", "Treated", "Treated")
dir.create(root.s, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(Sample_ID = samples.v, Name = paste0("Mouse_", seq_along(samples.v)), Treatment = groups.v,
                            Exclude = c("no", "no", "no", "no", "no", "yes"), Age_weeks = c(8, 9, 8, 10, 9, 8)),
                 file.path(root.s, "metadata.csv"), row.names = FALSE)

dir.create(file.path(results.s, "pipeline_info"), recursive = TRUE)
writeLines('{"mode": "illumina_metagenome", "outdir": "results", "hq_quality_source": "both"}',
           file.path(results.s, "pipeline_info", "run_params.json"))

write_tsv(data.frame(Sample_ID = samples.v, GBbp = round(runif(6, 1, 3), 3), Raw_count = 1e6 + 1:6,
                     Fastp_count = 9e5 + 1:6, Fastp_percent = 90, Reads_mapped_HQ_Derep_MAGs_percent = round(runif(6, 20, 60), 2)),
          "00_read_stats", "read_stat_report.tsv")

lineages.v <- c(
  "d__Bacteria|p__Bacillota_A|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia obeum|t__GCF_000001.1",
  "d__Bacteria|p__Bacillota_A|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia wexlerae|t__GCF_000002.1",
  "d__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides fragilis|t__GCF_000003.1",
  "d__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola vulgatus|t__GCF_000004.1",
  "d__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia coli|t__GCF_000005.1",
  "d__Archaea|p__Methanobacteriota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter smithii|t__GCF_000006.1"
)
leaf.m <- matrix(rexp(6 * 6), nrow = 6)
leaf.m[5, 1:3] <- 0
leaf.m[6, 4:6] <- 0
leaf.m <- sweep(leaf.m, 2, colSums(leaf.m), "/") * 100
sylph_rows.l <- list()
for (depth.n in 1:8){
  prefix.v <- vapply(strsplit(lineages.v, "|", fixed = TRUE), function(x) paste(x[seq_len(depth.n)], collapse = "|"), character(1))
  summed.m <- rowsum(leaf.m, prefix.v, reorder = FALSE)
  sylph_rows.l[[depth.n]] <- data.frame(clade_name = rownames(summed.m), round(summed.m, 4), check.names = FALSE)
}
sylph.df <- do.call(rbind, sylph_rows.l)
names(sylph.df)[-1] <- paste0("/work/reads/", samples.v, ".clean_1.fastq.gz")
write_tsv(sylph.df, "02_sylph", "merged_relative_abundance.tsv")
write_tsv(sylph.df, "02_sylph", "merged_sequence_abundance.tsv")

singlem_lineages.v <- c("Root; d__Bacteria", "Root; d__Bacteria; p__Bacillota_A; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae",
                        "Root; d__Bacteria; p__Bacillota_A; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae; g__Blautia",
                        "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides",
                        "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides fragilis",
                        "Root; d__Archaea; p__Methanobacteriota")
singlem.df <- expand.grid(taxonomy = singlem_lineages.v, sample = paste0(samples.v, "_R1"), stringsAsFactors = FALSE)
singlem.df$coverage <- round(rexp(nrow(singlem.df)) * 5, 2)
write_tsv(singlem.df[, c("sample", "coverage", "taxonomy")], "03_singlem", "metagenome.condensed.tsv")
write_tsv(data.frame(sample = paste0(samples.v, "_R1"), bacterial_archaeal_bases = 1e8, metagenome_size = 2e8,
                     read_fraction = round(runif(6, 40, 90), 2), average_bacterial_archaeal_genome_size = 3e6, warning = "",
                     domain_relative_abundance = 1),
          "03_singlem", "metagenome.prokaryotic_fraction.tsv")

bins.v <- c("S1.metabat2.1", "S1.semibin.2", "S2.metabat2.1", "S10.vamb.3", "S11.metabat2.4", "S12.rosella.5")
write_tsv(data.frame(Name = bins.v, Completeness = c(98, 72, 95, 60, 91, 88), Contamination = c(1, 8, 2, 12, 0.5, 3),
                     Completeness_Model_Used = "Gradient Boost (General Model)", Translation_Table_Used = 11,
                     Coding_Density = 0.9, Contig_N50 = 25000, Average_Gene_Length = 300, Genome_Size = 3e6,
                     GC_Content = 0.45, Total_Coding_Sequences = 3000, Total_Contigs = 120, Max_Contig_Length = 1e5,
                     Additional_Notes = "None"),
          "07_checkm2", "all_bins_checkm2_report.tsv")
checkm1.df <- data.frame(bins.v, "o__Clostridiales (UID1)", 100, 500, 300, c(96, 70, 90, 55, 92, 85), c(1.5, 9, 2, 15, 1, 4), 0)
names(checkm1.df) <- c("Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "Completeness",
                       "Contamination", "Strain heterogeneity")
write_tsv(checkm1.df, "16_checkm1", "checkm_qa_results.tsv")
write_tsv(data.frame(user_genome = bins.v,
                     classification = c(gsub("\\|", ";", sub("\\|t__.*", "", lineages.v[c(1, 1, 3, 4, 5)])),
                                        "d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__;s__"),
                     classification_method = "ANI", note = "", msa_percent = 90, red_value = "N/A", warnings = "N/A"),
          "15_gtdbtk", "all_genomes", "classify", "all_genomes.bac120.summary.tsv")
write_tsv(data.frame(V1 = paste0("bins/", bins.v[c(1, 1, 3, 4, 5, 6)], ".fasta"),
                     V2 = paste0("bins/", bins.v[c(1, 2, 3, 4, 5, 6)], ".fasta")),
          "08_dereplicated_bins", "cluster_definition.tsv")
file_lines.v <- readLines(file.path(results.s, "08_dereplicated_bins", "cluster_definition.tsv"))[-1]
writeLines(file_lines.v, file.path(results.s, "08_dereplicated_bins", "cluster_definition.tsv"))
hq.v <- bins.v[c(1, 3, 5, 6)]
dir.create(file.path(results.s, "08_dereplicated_hq_bins"), showWarnings = FALSE)
writeLines(paste0("bins/", hq.v, ".fasta\tbins/", hq.v, ".fasta"), file.path(results.s, "08_dereplicated_hq_bins", "cluster_definition.tsv"))

for (i in seq_along(samples.v)){
  label.s <- paste0(samples.v[i], ".clean_1.fastq.gz")
  coverage.v <- round(rexp(length(hq.v)) * 5, 4)
  relative.v <- round(coverage.v / sum(coverage.v) * 60, 4)
  coverm.df <- data.frame(Genome = c("unmapped", hq.v), c(40, relative.v), c(NA, rep(0.9, length(hq.v))),
                          c(NA, coverage.v), c(NA, round(coverage.v * 1000)), check.names = FALSE)
  names(coverm.df)[-1] <- paste(label.s, c("Relative Abundance (%)", "Covered Fraction", "Mean", "Read Count"))
  write_tsv(coverm.df, "09_coverm_hq_derep_bins", paste0(samples.v[i], "_abundances.tsv"))
}

spot_targets.v <- c("ph_optimum", "ph_min", "ph_max", "salinity_optimum", "temperature_optimum", "oxygen")
genomespot.df <- expand.grid(target = spot_targets.v, bin_id = bins.v, stringsAsFactors = FALSE)
genomespot.df$value <- ifelse(genomespot.df$target == "oxygen", sample(c("tolerant", "not tolerant"), nrow(genomespot.df), TRUE),
                              round(runif(nrow(genomespot.df), 5, 40), 3))
genomespot.df$value[genomespot.df$bin_id == bins.v[4] & genomespot.df$target != "oxygen"] <- "None"
genomespot.df$error <- round(runif(nrow(genomespot.df)), 3)
genomespot.df$units <- "unit"
genomespot.df$is_novel <- "False"
genomespot.df$warning <- "None"
write_tsv(genomespot.df[, c("bin_id", "target", "value", "error", "units", "is_novel", "warning")],
          "18_genomespot", "genomespot_combined.tsv")

genes.v <- paste0(rep(samples.v, each = 5), "___NODE_", 1:30, "_length_1000_cov_2.0_1")
dram.df <- data.frame(V1 = paste0("chunk.001_", genes.v), fasta = "chunk.001", rank = "C",
                      ko_id = sample(c("K00001", "K00002", "K00003", NA), 30, TRUE),
                      kegg_hit = "alcohol dehydrogenase [EC:1.1.1.1]",
                      cazy_ids = sample(c("GH13; CBM48", "GT2", NA, "GH23"), 30, TRUE), cazy_best_hit = "",
                      peptidase_family = sample(c("M23", NA), 30, TRUE))
names(dram.df)[1] <- ""
write_tsv(dram.df, "13_dram", "dram_annotations", "annotations.tsv")
rpkm.df <- data.frame(Gene_ID = genes.v, matrix(round(rexp(30 * 6), 4), nrow = 30, dimnames = list(NULL, samples.v)), check.names = FALSE)
write_tsv(rpkm.df, "23_rpkm", "gene_catalogue_rpkm_per_gene_normalised.tsv")

virus.v <- paste0(samples.v, ".scaffolds__genomad__NODE_", 1:6, "_length_5000_cov_3.0")
write_tsv(data.frame(seq_name = virus.v, length = 5000, topology = "No terminal repeats", coordinates = NA, n_genes = 6,
                     genetic_code = 11, virus_score = 0.9, fdr = NA, n_hallmarks = 2, marker_enrichment = 5,
                     taxonomy = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"),
          "20_genomad", "genomad_virus_summary.tsv")
dir.create(file.path(results.s, "22_checkv_clustering"), showWarnings = FALSE)
writeLines(c(paste(virus.v[1], paste(virus.v[c(1, 4)], collapse = ","), sep = "\t"),
             paste(virus.v[2], paste(virus.v[c(2, 3, 5)], collapse = ","), sep = "\t"),
             paste(virus.v[6], virus.v[6], sep = "\t")),
           file.path(results.s, "22_checkv_clustering", "virus_clusters.tsv"))

npo_header.v <- c("# @impl: Nonpareil", "# @ksize: 24", "# @version: 3.5.5", "# @L: 150.000", "# @AL: 120.000",
                  "# @R: 2000000", "# @overlap: 50.00", "# @divide: 0.70")
for (i in seq_along(samples.v)){
  effort.v <- unique(round(0.7^(40:0) * 2e6))
  coverage.v <- 1 - exp(-effort.v / (2e5 * i))
  npo.df <- data.frame(effort.v, round(coverage.v * 0.9, 5), 0.02, round(coverage.v * 0.85, 5), round(coverage.v * 0.9, 5),
                       round(pmin(coverage.v * 0.95, 1), 5))
  path.s <- file.path(results.s, "17_nonpareil", paste0(samples.v[i], ".npo"))
  dir.create(dirname(path.s), showWarnings = FALSE)
  writeLines(npo_header.v, path.s)
  utils::write.table(npo.df, path.s, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}
cat("Fixtures written to", root.s, "\n")

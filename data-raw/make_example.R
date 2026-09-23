# Builds the synthetic example project used by the vignettes and example_project():
# inst/extdata/example (pipeline results plus a metadata table). Larger than the test
# fixtures so figures look like a real study: 16 mouse samples (Control vs Treated, two
# time points), about 40 sylph genomes across six phyla with treatment and time effects,
# and SingleM. Only sylph and SingleM are included, like an early partial pipeline run.
# Values are made up. Run from the package root:
#   Rscript data-raw/make_example.R

set.seed(2026)
root.s <- "inst/extdata/example"
unlink(root.s, recursive = TRUE)
results.s <- file.path(root.s, "results")
write_tsv <- function(x.df, ...){
  path.s <- file.path(results.s, ...)
  dir.create(dirname(path.s), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x.df, path.s, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

metadata.df <- expand.grid(Mouse = 1:4, Treatment = c("Control", "Treated"), Time = c("Week_0", "Week_4"),
                           stringsAsFactors = FALSE)
metadata.df$Sample_ID <- sprintf("S%02d", seq_len(nrow(metadata.df)))
metadata.df$Name <- sprintf("%s_M%d_%s", substr(metadata.df$Treatment, 1, 3), metadata.df$Mouse + 4 * (metadata.df$Treatment == "Treated"),
                            sub("Week_", "W", metadata.df$Time))
metadata.df$Age_weeks <- 8 + 4 * (metadata.df$Time == "Week_4") + sample(0:1, nrow(metadata.df), replace = TRUE)
metadata.df$Cage <- paste0("Cage_", 1 + (metadata.df$Mouse > 2) + 2 * (metadata.df$Treatment == "Treated"))
metadata.df$Exclude <- "no"
samples.v <- metadata.df$Sample_ID
dir.create(root.s, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(metadata.df[c("Sample_ID", "Name", "Treatment", "Time", "Age_weeks", "Cage", "Exclude")],
                 file.path(root.s, "metadata.csv"), row.names = FALSE)

dir.create(file.path(results.s, "pipeline_info"), recursive = TRUE)
writeLines('{"mode": "illumina_metagenome", "outdir": "results", "hq_quality_source": "both"}',
           file.path(results.s, "pipeline_info", "run_params.json"))

# Genus, family, order, class, phylum, domain; one or two genomes per genus
taxa.df <- read.table(text = "
Muribaculum Muribaculaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Duncaniella Muribaculaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Paramuribaculum Muribaculaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
CAG-485 Muribaculaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Bacteroides Bacteroidaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Phocaeicola Bacteroidaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Alistipes Rikenellaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Parabacteroides Tannerellaceae Bacteroidales Bacteroidia Bacteroidota Bacteria
Lactobacillus Lactobacillaceae Lactobacillales Bacilli Bacillota Bacteria
Ligilactobacillus Lactobacillaceae Lactobacillales Bacilli Bacillota Bacteria
Limosilactobacillus Lactobacillaceae Lactobacillales Bacilli Bacillota Bacteria
Turicibacter Turicibacteraceae Erysipelotrichales Bacilli Bacillota Bacteria
Dubosiella Erysipelotrichaceae Erysipelotrichales Bacilli Bacillota Bacteria
Faecalibaculum Erysipelotrichaceae Erysipelotrichales Bacilli Bacillota Bacteria
Enterococcus Enterococcaceae Lactobacillales Bacilli Bacillota Bacteria
Lachnospiraceae_NK4A136 Lachnospiraceae Lachnospirales Clostridia Bacillota_A Bacteria
Roseburia Lachnospiraceae Lachnospirales Clostridia Bacillota_A Bacteria
Blautia Lachnospiraceae Lachnospirales Clostridia Bacillota_A Bacteria
Acetatifactor Lachnospiraceae Lachnospirales Clostridia Bacillota_A Bacteria
Kineothrix Lachnospiraceae Lachnospirales Clostridia Bacillota_A Bacteria
Schaedlerella Lachnospiraceae Lachnospirales Clostridia Bacillota_A Bacteria
Oscillibacter Oscillospiraceae Oscillospirales Clostridia Bacillota_A Bacteria
Flavonifractor Oscillospiraceae Oscillospirales Clostridia Bacillota_A Bacteria
Ruminococcus Ruminococcaceae Oscillospirales Clostridia Bacillota_A Bacteria
Anaerotruncus Ruminococcaceae Oscillospirales Clostridia Bacillota_A Bacteria
Clostridium Clostridiaceae Clostridiales Clostridia Bacillota_A Bacteria
Eubacterium Eubacteriaceae Eubacteriales Clostridia Bacillota_A Bacteria
Bifidobacterium Bifidobacteriaceae Bifidobacteriales Actinomycetes Actinomycetota Bacteria
Adlercreutzia Eggerthellaceae Coriobacteriales Coriobacteriia Actinomycetota Bacteria
Enterorhabdus Eggerthellaceae Coriobacteriales Coriobacteriia Actinomycetota Bacteria
Akkermansia Akkermansiaceae Verrucomicrobiales Verrucomicrobiae Verrucomicrobiota Bacteria
Escherichia Enterobacteriaceae Enterobacterales Gammaproteobacteria Pseudomonadota Bacteria
Parasutterella Burkholderiaceae Burkholderiales Gammaproteobacteria Pseudomonadota Bacteria
Desulfovibrio Desulfovibrionaceae Desulfovibrionales Desulfovibrionia Desulfobacterota Bacteria
Mucispirillum Mucispirillaceae Deferribacterales Deferribacteres Deferribacterota Bacteria
Methanobrevibacter Methanobacteriaceae Methanobacteriales Methanobacteria Methanobacteriota Archaea
", col.names = c("Genus", "Family", "Order", "Class", "Phylum", "Domain"), stringsAsFactors = FALSE)
taxa.df <- taxa.df[rep(seq_len(nrow(taxa.df)), ifelse(taxa.df$Genus %in% c("Muribaculum", "Duncaniella", "Lactobacillus",
                                                                           "Bacteroides", "Blautia"), 2, 1)), ]
taxa.df$Species <- paste(taxa.df$Genus, sprintf("sp%06d", sample(1e5:9e5, nrow(taxa.df))))
taxa.df$Accession <- sprintf("GCF_%09d.1", sample(1e8:9e8, nrow(taxa.df)))
lineages.v <- sprintf("d__%s|p__%s|c__%s|o__%s|f__%s|g__%s|s__%s|t__%s", taxa.df$Domain, taxa.df$Phylum, taxa.df$Class,
                      taxa.df$Order, taxa.df$Family, taxa.df$Genus, taxa.df$Species, taxa.df$Accession)
n_genomes.n <- length(lineages.v)

# Log-normal baseline abundances (Muribaculaceae and Lactobacillus dominate the mouse gut),
# a treatment effect on some genera, a time effect on others, then per-sample noise and dropout
baseline.v <- stats::rnorm(n_genomes.n, mean = 0, sd = 1.4)
baseline.v[taxa.df$Family %in% c("Muribaculaceae", "Lactobacillaceae")] <- baseline.v[taxa.df$Family %in% c("Muribaculaceae", "Lactobacillaceae")] + 2.5
treated.v <- setNames(rep(0, n_genomes.n), taxa.df$Genus)
treated.v[taxa.df$Genus %in% c("Akkermansia", "Bacteroides", "Escherichia", "Enterococcus", "Parasutterella")] <- 2.2
treated.v[taxa.df$Genus %in% c("Lachnospiraceae_NK4A136", "Roseburia", "Dubosiella", "Turicibacter", "Kineothrix")] <- -2.4
time.v <- setNames(rep(0, n_genomes.n), taxa.df$Genus)
time.v[taxa.df$Genus %in% c("Bifidobacterium", "Ligilactobacillus", "Faecalibaculum")] <- -1.5
time.v[taxa.df$Genus %in% c("Oscillibacter", "Alistipes", "Desulfovibrio")] <- 1.3
treated_after.v <- (metadata.df$Treatment == "Treated") * (metadata.df$Time == "Week_4")
log.m <- outer(baseline.v, rep(1, length(samples.v))) + outer(treated.v, treated_after.v) +
  outer(time.v, metadata.df$Time == "Week_4") + matrix(stats::rnorm(n_genomes.n * length(samples.v), sd = 0.7), n_genomes.n)
abundance.m <- exp(log.m)
# Rare genomes drop below detection more often
detect.m <- matrix(stats::runif(length(abundance.m)), n_genomes.n) < stats::plogis(log.m + 1)
abundance.m <- abundance.m * detect.m
leaf.m <- sweep(abundance.m, 2, colSums(abundance.m), "/") * 100

sylph_rows.l <- list()
for (depth.n in 1:8){
  prefix.v <- vapply(strsplit(lineages.v, "|", fixed = TRUE), function(x) paste(x[seq_len(depth.n)], collapse = "|"), character(1))
  summed.m <- rowsum(leaf.m, prefix.v, reorder = FALSE)
  sylph_rows.l[[depth.n]] <- data.frame(clade_name = rownames(summed.m), round(summed.m, 5), check.names = FALSE)
}
sylph.df <- do.call(rbind, sylph_rows.l)
sylph.df <- sylph.df[rowSums(sylph.df[-1]) > 0, ]
names(sylph.df)[-1] <- paste0("/work/reads/", samples.v, ".clean_1.fastq.gz")
write_tsv(sylph.df, "04_sylph", "merged_relative_abundance.tsv")
write_tsv(sylph.df, "04_sylph", "merged_sequence_abundance.tsv")

# SingleM: coverage per genus (and family for some), from the same community with extra noise
genus_lineage.v <- sprintf("Root; d__%s; p__%s; c__%s; o__%s; f__%s; g__%s", taxa.df$Domain, taxa.df$Phylum, taxa.df$Class,
                           taxa.df$Order, taxa.df$Family, taxa.df$Genus)
coverage.m <- rowsum(leaf.m, genus_lineage.v) * stats::runif(length(samples.v), 0.2, 0.6)
coverage.m <- coverage.m * matrix(exp(stats::rnorm(length(coverage.m), sd = 0.3)), nrow(coverage.m))
singlem.df <- data.frame(sample = rep(paste0(samples.v, "_R1"), each = nrow(coverage.m)),
                         coverage = round(as.vector(coverage.m), 2), taxonomy = rep(rownames(coverage.m), length(samples.v)))
write_tsv(singlem.df[singlem.df$coverage > 0, ], "05_singlem", "metagenome.condensed.tsv")
write_tsv(data.frame(sample = paste0(samples.v, "_R1"), bacterial_archaeal_bases = 1e9, metagenome_size = 1.2e9,
                     read_fraction = round(stats::runif(length(samples.v), 75, 95), 2),
                     average_bacterial_archaeal_genome_size = 3.2e6, warning = "", domain_relative_abundance = 1),
          "05_singlem", "metagenome.prokaryotic_fraction.tsv")

# Builds the synthetic example project used by the vignettes and example_project():
# inst/extdata/example (pipeline results plus a metadata table). Larger than the test
# fixtures so figures look like a real study: 24 samples from 12 mice (Control vs Treated, two
# time points), about 40 sylph genomes across six phyla with treatment and time effects
# (treated mice at week 4 also lose rare genomes, so their richness is lower),
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

metadata.df <- expand.grid(Mouse = 1:6, Treatment = c("Control", "Treated"), Time = c("Week_0", "Week_4"),
                           stringsAsFactors = FALSE)
metadata.df$Sample_ID <- sprintf("S%02d", seq_len(nrow(metadata.df)))
metadata.df$Name <- sprintf("%s_M%d_%s", substr(metadata.df$Treatment, 1, 3), metadata.df$Mouse + 6 * (metadata.df$Treatment == "Treated"),
                            sub("Week_", "W", metadata.df$Time))
metadata.df$Age_weeks <- 8 + 4 * (metadata.df$Time == "Week_4") + sample(0:1, nrow(metadata.df), replace = TRUE)
metadata.df$Cage <- paste0("Cage_", 1 + (metadata.df$Mouse > 3) + 2 * (metadata.df$Treatment == "Treated"))
metadata.df$Mouse <- sprintf("M%d", metadata.df$Mouse + 6 * (metadata.df$Treatment == "Treated"))
metadata.df$Exclude <- "no"
samples.v <- metadata.df$Sample_ID
dir.create(root.s, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(metadata.df[c("Sample_ID", "Name", "Treatment", "Time", "Mouse", "Age_weeks", "Cage", "Exclude")],
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
treated.v[taxa.df$Genus %in% c("Akkermansia", "Bacteroides", "Escherichia", "Enterococcus", "Parasutterella")] <- 2.5
treated.v[taxa.df$Genus %in% c("Lachnospiraceae_NK4A136", "Roseburia", "Dubosiella", "Turicibacter", "Kineothrix")] <- -2.5
time.v <- setNames(rep(0, n_genomes.n), taxa.df$Genus)
time.v[taxa.df$Genus %in% c("Bifidobacterium", "Ligilactobacillus", "Faecalibaculum")] <- -1.5
time.v[taxa.df$Genus %in% c("Oscillibacter", "Alistipes", "Desulfovibrio")] <- 1.3
# Responding genera start at moderate abundance, so their changes are measurable
responders.v <- treated.v != 0 | time.v != 0
baseline.v[responders.v] <- stats::rnorm(sum(responders.v), mean = 1.2, sd = 0.4)
treated_after.v <- (metadata.df$Treatment == "Treated") * (metadata.df$Time == "Week_4")
log.m <- outer(baseline.v, rep(1, length(samples.v))) + outer(treated.v, treated_after.v) +
  outer(time.v, metadata.df$Time == "Week_4") + matrix(stats::rnorm(n_genomes.n * length(samples.v), sd = 0.7), n_genomes.n)
abundance.m <- exp(log.m)
# Rare genomes drop below detection more often, and in treated mice at week 4 so do genomes the
# treatment does not favour (lower richness)
detect.m <- matrix(stats::runif(length(abundance.m)), n_genomes.n) <
  stats::plogis(log.m + 1 - outer(2.5 * (treated.v <= 0), treated_after.v))
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

# ---- Isolate example ---------------------------------------------------------------
# inst/extdata/example_isolate: 20 Escherichia coli isolates from four sequence types (clades)
# and three sources, plus two reference genomes, with clade-linked AMR genes, insertion
# sequences, accessory genes and cgMLST alleles, and a tree built from the allele distances.
set.seed(2027)
root.s <- "inst/extdata/example_isolate"
unlink(root.s, recursive = TRUE)
# Short directory name keeps package paths within the 100-byte tar limit
results.s <- file.path(root.s, "res")
dir.create(file.path(results.s, "pipeline_info"), recursive = TRUE)
writeLines('{"mode": "illumina_isolate", "outdir": "results"}', file.path(results.s, "pipeline_info", "run_params.json"))

isolates.v <- sprintf("ISO%02d", 1:20)
clade.v <- rep(c("ST131", "ST73", "ST10", "ST69"), c(8, 5, 5, 2))
source.v <- c("Clinical", "Clinical", "Clinical", "Clinical", "Clinical", "Clinical", "Animal", "Environmental",
              "Clinical", "Clinical", "Animal", "Environmental", "Clinical",
              "Animal", "Animal", "Environmental", "Environmental", "Animal",
              "Clinical", "Environmental")
references.v <- c("GCF_000005845.2", "GCF_000010485.1")
genomes.v <- c(isolates.v, references.v)
genome_clade.v <- c(clade.v, "ST10", "ST131")
dir.create(root.s, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(Sample_ID = isolates.v, Name = sprintf("EC%02d", 1:20), Source = source.v,
                            Year = sample(2019:2023, 20, replace = TRUE), Exclude = "no"),
                 file.path(root.s, "metadata.csv"), row.names = FALSE)

write_tsv(data.frame(Sample_ID = isolates.v, GBbp = round(stats::runif(20, 0.3, 0.9), 3),
                     Raw_count = round(stats::runif(20, 1.5e6, 4e6)), Covered_fraction = round(stats::runif(20, 0.985, 0.999), 3),
                     Mean_coverage = round(stats::runif(20, 45, 160), 1), Read_count = round(stats::runif(20, 1.4e6, 3.8e6)),
                     Read_count_percent = round(stats::runif(20, 94, 99), 1)),
          "00_read_stats", "read_stat_report.tsv")
genome_size.v <- round(c(ST131 = 5.15e6, ST73 = 5.05e6, ST10 = 4.75e6, ST69 = 5.0e6)[clade.v] + stats::rnorm(20, 0, 6e4))
for (i in seq_along(isolates.v)){
  write_tsv(data.frame(Name = paste0(isolates.v[i], ".scaffolds"), Completeness = round(stats::runif(1, 98.8, 100), 2),
                       Contamination = round(stats::runif(1, 0.05, 1.8), 2),
                       Completeness_Model_Used = "Neural Network (Specific Model)", Translation_Table_Used = 11,
                       Coding_Density = 0.87, Contig_N50 = round(stats::runif(1, 9e4, 3.2e5)), Average_Gene_Length = 310,
                       Genome_Size = genome_size.v[i], GC_Content = round(stats::rnorm(1, 0.506, 0.002), 3),
                       Total_Coding_Sequences = round(genome_size.v[i] / 1070), Total_Contigs = round(stats::runif(1, 40, 180)),
                       Max_Contig_Length = 5e5, Additional_Notes = "None"),
            "07_checkm2", paste0(isolates.v[i], "_checkm2_report.tsv"))
}
write_tsv(data.frame(user_genome = paste0(isolates.v, ".scaffolds"),
                     classification = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
                     classification_method = "ANI"),
          "15_gtdbtk", "all_genomes", "classify", "all_genomes.bac120.summary.tsv")
write_tsv(data.frame(Genome = isolates.v, Taxon = "Escherichia coli", Complete = "False", `Translation table` = 11,
                     `# Sequences` = 100, Size = genome_size.v, GC = 0.506, `N ratio` = 0, N50 = 150000, `Coding ratio` = 0.87,
                     tRNA = sample(78:88, 20, replace = TRUE), tmRNA = 1, rRNA = sample(5:8, 20, replace = TRUE), ncRNA = 70,
                     `ncRNA region` = 30, CRISPR = sample(0:2, 20, replace = TRUE), CDS = round(genome_size.v / 1080),
                     `CDS hypothetical` = 110, `CDS pseudogene` = 25, sORF = 5, GAP = 0, oriC = 1, oriV = 0, oriT = 0,
                     check.names = FALSE),
          "10_bakta", "bakta_annotation_stats.tsv")
st_alleles.l <- list(ST131 = "adk(53),fumC(40),gyrB(47),icd(13),mdh(36),purA(28),recA(29)",
                     ST73 = "adk(36),fumC(24),gyrB(9),icd(13),mdh(17),purA(11),recA(25)",
                     ST10 = "adk(10),fumC(11),gyrB(4),icd(8),mdh(8),purA(8),recA(2)",
                     ST69 = "adk(21),fumC(35),gyrB(27),icd(6),mdh(5),purA(5),recA(4)")
dir.create(file.path(results.s, "11_mlst"))
writeLines(paste0("/work/", isolates.v, ".scaffolds.fasta,ecoli_achtman_4,", sub("ST", "", clade.v), ",",
                  unlist(st_alleles.l[clade.v])),
           file.path(results.s, "11_mlst", "mlst_summary.csv"))

# AMR genes: probability of carriage per clade (ST131 multidrug resistant, ST10 from animals with tet/sul)
amr_genes.df <- read.table(text = "
blaCTX-M-15 BETA-LACTAM 0.85 0.1 0.1 0
blaTEM-1 BETA-LACTAM 0.4 0.4 0.5 0.5
aac(6')-Ib-cr AMINOGLYCOSIDE/QUINOLONE 0.7 0 0.05 0
aadA5 AMINOGLYCOSIDE 0.6 0.2 0.1 0
dfrA17 TRIMETHOPRIM 0.7 0.2 0.1 0
sul1 SULFONAMIDE 0.6 0.2 0.5 0
sul2 SULFONAMIDE 0.3 0.4 0.7 0.5
tet(A) TETRACYCLINE 0.3 0.2 0.8 0
tet(B) TETRACYCLINE 0 0.1 0.4 0.5
mph(A) MACROLIDE 0.6 0.1 0 0
qnrS1 QUINOLONE 0.1 0 0.3 0
floR PHENICOL 0 0 0.4 0
mdtM EFFLUX 1 1 1 1
", col.names = c("Symbol", "Class", "ST131", "ST73", "ST10", "ST69"), stringsAsFactors = FALSE)
amr_rows.l <- list()
for (i in seq_along(isolates.v)){
  carried.v <- amr_genes.df$Symbol[stats::runif(nrow(amr_genes.df)) < amr_genes.df[[clade.v[i]]]]
  if (length(carried.v) == 0) next
  amr_rows.l[[i]] <- data.frame(Sample = isolates.v[i], `Protein id` = paste0(isolates.v[i], "_", seq_along(carried.v)),
                                `Contig id` = "contig_1", Start = 100, Stop = 900, Strand = "+", `Element symbol` = carried.v,
                                `Element name` = "resistance gene", Scope = ifelse(carried.v == "mdtM", "plus", "core"),
                                Type = "AMR", Subtype = "AMR", Class = amr_genes.df$Class[match(carried.v, amr_genes.df$Symbol)],
                                Subclass = "", Method = "EXACTX", check.names = FALSE)
}
write_tsv(do.call(rbind, amr_rows.l), "12_amrfinder", "amrfinder_all.tsv")

# Insertion sequences: copies per family, with more IS in ST131
is_rate.m <- rbind(IS1 = c(4, 2, 3, 2), IS3 = c(6, 3, 2, 3), IS30 = c(2, 1, 1, 0), IS66 = c(3, 0, 1, 1), IS110 = c(1, 1, 0, 1))
colnames(is_rate.m) <- c("ST131", "ST73", "ST10", "ST69")
is_rows.l <- list()
for (i in seq_along(isolates.v)){
  copies.v <- stats::rpois(nrow(is_rate.m), is_rate.m[, clade.v[i]])
  families.v <- rep(rownames(is_rate.m), copies.v)
  if (length(families.v) == 0) next
  is_rows.l[[i]] <- data.frame(Genome_ID = isolates.v[i], seqID = "contig_1", family = families.v, cluster = "c1",
                               isBegin = 100, isEnd = 1400, isLen = 1300, ncopy4is = 1, type = "c")
}
write_tsv(do.call(rbind, is_rows.l), "13_isescan", "ISEScan_summary.tsv")

group_dir.s <- file.path("14_comparison_groups", "comparison_groups", "all")
write_tsv(data.frame(comparison_group = "all", entry_type = c(rep("sample", 20), "reference", "reference"), id = genomes.v,
                     fasta = paste0("/work/comparison_groups/all/fastas/", genomes.v, ".fasta"), gff = "", faa = "",
                     parsnp_reference = c(rep("false", 20), "true", "false")),
          group_dir.s, "entries.tsv")

# Pangenome: core genes, clade-specific shell genes, scattered cloud genes and some mobile genes
core.v <- sprintf("core_%03d", 1:150)
clade_genes.l <- lapply(stats::setNames(c("ST131", "ST73", "ST10", "ST69"), c("ST131", "ST73", "ST10", "ST69")),
                        function(clade.s) sprintf("%s_group_%02d", clade.s, 1:15))
cloud.v <- sprintf("group_%04d", sample(1000:9999, 40))
mobile.v <- c("intI1", "traA", "traD", "repFIB", "repFII", "phage_int", "papC", "papG", "hlyA", "cnf1", "iutA", "fyuA")
genes.v <- c(core.v, unlist(clade_genes.l), cloud.v, mobile.v)
presence.m <- matrix(0, nrow = length(genes.v), ncol = length(genomes.v), dimnames = list(genes.v, genomes.v))
presence.m[core.v, ] <- 1
presence.m[core.v, ][stats::runif(length(core.v) * length(genomes.v)) < 0.01] <- 0
for (j in seq_along(genomes.v)){
  own.v <- clade_genes.l[[genome_clade.v[j]]]
  presence.m[own.v, j] <- as.numeric(stats::runif(length(own.v)) < 0.93)
  presence.m[cloud.v, j] <- as.numeric(stats::runif(length(cloud.v)) < 0.08)
  mobile_rate.v <- if (genome_clade.v[j] %in% c("ST131", "ST73")) 0.6 else 0.25
  presence.m[mobile.v, j] <- as.numeric(stats::runif(length(mobile.v)) < mobile_rate.v)
}
presence.m <- presence.m[rowSums(presence.m) > 0, ]
annotation.v <- ifelse(grepl("^core_", rownames(presence.m)), "conserved protein",
                       ifelse(rownames(presence.m) %in% mobile.v, "mobile element or virulence protein", "hypothetical protein"))
panaroo.df <- data.frame(Gene = rownames(presence.m), `Non-unique Gene name` = "", Annotation = annotation.v, check.names = FALSE)
for (j in seq_along(genomes.v)){
  panaroo.df[[genomes.v[j]]] <- ifelse(presence.m[, j] == 1, paste0(genomes.v[j], "_", sprintf("%04d", seq_len(nrow(presence.m)))), "")
}
path.s <- file.path(results.s, "15_panaroo", "all_panaroo", "gene_presence_absence.csv")
dir.create(dirname(path.s), recursive = TRUE)
utils::write.csv(panaroo.df, path.s, row.names = FALSE, quote = FALSE)

# cgMLST: 300 loci; clades differ at most loci, isolates within a clade at a few
n_loci.n <- 300
clade_profiles.l <- lapply(stats::setNames(unique(genome_clade.v), unique(genome_clade.v)),
                           function(x) sample(1:40, n_loci.n, replace = TRUE))
alleles.m <- t(vapply(genome_clade.v, function(clade.s){
  alleles.v <- clade_profiles.l[[clade.s]]
  changed.v <- stats::runif(n_loci.n) < 0.03
  alleles.v[changed.v] <- sample(41:99, sum(changed.v), replace = TRUE)
  as.character(alleles.v)
}, character(n_loci.n)))
alleles.m[sample(length(alleles.m), 60)] <- "LNF"
alleles.m[sample(length(alleles.m), 20)] <- paste0("INF-", sample(100:200, 20, replace = TRUE))
rownames(alleles.m) <- genomes.v
cgmlst.df <- data.frame(FILE = paste0(genomes.v, ".fasta"), alleles.m, check.names = FALSE)
names(cgmlst.df)[-1] <- paste0("EC_locus_", sprintf("%04d", seq_len(n_loci.n)), ".fasta")
write_tsv(cgmlst.df, "19_chewbacca", "all_chewbbaca", "cgMLST", "cgMLST95.tsv")
write_tsv(data.frame(original_id = genomes.v, chewbbaca_id = genomes.v), "19_chewbacca", "all_chewbbaca", "genome_hash_map.tsv")

# fastANI from the allele distances: 99.9% within a sequence type, about 97.5-98.8% between
calls.m <- matrix(suppressWarnings(as.numeric(sub("^INF-", "", alleles.m))), nrow(alleles.m), dimnames = dimnames(alleles.m))
difference.m <- matrix(0, length(genomes.v), length(genomes.v), dimnames = list(genomes.v, genomes.v))
for (i in seq_along(genomes.v)) for (j in seq_along(genomes.v)){
  both.v <- !is.na(calls.m[i, ]) & !is.na(calls.m[j, ])
  difference.m[i, j] <- mean(calls.m[i, both.v] != calls.m[j, both.v])
}
ani.m <- 100 - 2.4 * difference.m + matrix(stats::runif(length(difference.m), -0.03, 0.03), nrow(difference.m))
diag(ani.m) <- 100
pairs.df <- expand.grid(q = seq_along(genomes.v), r = seq_along(genomes.v))
dir.create(file.path(results.s, "18_fastani"))
utils::write.table(data.frame(paste0("/work/fastas/", genomes.v[pairs.df$q], ".fasta"),
                              paste0("/work/fastas/", genomes.v[pairs.df$r], ".fasta"),
                              round(ani.m[cbind(pairs.df$q, pairs.df$r)], 4), 1650, 1700),
                   file.path(results.s, "18_fastani", "all.fastani.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
                   col.names = FALSE)

# Tree: neighbour joining on the allele differences, tips named like IQ-TREE output
tree <- ape::nj(stats::as.dist(difference.m))
tree$tip.label <- ifelse(tree$tip.label %in% references.v, paste0(tree$tip.label, ".fasta.ref"), paste0(tree$tip.label, ".fasta"))
tree$edge.length <- pmax(tree$edge.length, 1e-5)
dir.create(file.path(results.s, "20_tree"))
ape::write.tree(tree, file.path(results.s, "20_tree", "all.treefile"))

prophage.v <- isolates.v[stats::runif(20) < 0.6]
write_tsv(data.frame(seq_name = paste0(prophage.v, ".scaffolds__genomad__contig_", seq_along(prophage.v)),
                     length = round(stats::runif(length(prophage.v), 2e4, 5e4)), topology = "Provirus", coordinates = "1-40000",
                     n_genes = 50, genetic_code = 11, virus_score = 0.95, fdr = NA, n_hallmarks = 5, marker_enrichment = 20,
                     taxonomy = sample(c("Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;",
                                         "Viruses;Monodnaviria;Loebvirae;Hofneiviricota;Faserviricetes;Tubulavirales;Inoviridae"),
                                       length(prophage.v), replace = TRUE, prob = c(0.8, 0.2))),
          "20_genomad", "genomad_virus_summary.tsv")
cat("Example projects written to inst/extdata/example and", root.s, "\n")

# ---- Metagenome example: MAGs, strains and comparison samples -----------------------
# Ten MAGs with CheckM2 and GTDB-Tk results; inStrain and TRACS comparisons in which mice sharing a cage
# share strains; and 16 mice from two other facilities as comparison samples (pipeline --comparison_reads),
# profiled with sylph and SingleM, with their own metadata. A separate seed keeps the data above unchanged.
set.seed(2028)
root.s <- "inst/extdata/example"
results.s <- file.path(root.s, "results")
metadata.df <- utils::read.csv(file.path(root.s, "metadata.csv"), stringsAsFactors = FALSE)
samples.v <- metadata.df$Sample_ID
sylph.df <- utils::read.delim(file.path(results.s, "04_sylph", "merged_relative_abundance.tsv"), check.names = FALSE,
                              stringsAsFactors = FALSE)
leaf.df <- sylph.df[grepl("\\|t__", sylph.df$clade_name), , drop = FALSE]
lineages.v <- leaf.df$clade_name

# MAGs: one per abundant genome, from the sample where it is most abundant
mag_rows.v <- utils::head(order(-rowMeans(leaf.df[-1])), 10)
mag_lineages.v <- lineages.v[mag_rows.v]
mags.v <- sprintf("%s.metabat2.%d", samples.v[apply(leaf.df[mag_rows.v, -1], 1, which.max)], seq_along(mag_rows.v))
write_tsv(data.frame(Name = mags.v, Completeness = round(stats::runif(10, 72, 99.5), 2),
                     Contamination = round(stats::runif(10, 0.2, 6), 2), Completeness_Model_Used = "Gradient Boost (General Model)",
                     Translation_Table_Used = 11, Coding_Density = 0.89, Contig_N50 = round(stats::runif(10, 1e4, 9e4)),
                     Average_Gene_Length = 310, Genome_Size = round(stats::runif(10, 2e6, 4e6)), GC_Content = 0.46,
                     Total_Coding_Sequences = 2800, Total_Contigs = round(stats::runif(10, 40, 300)), Max_Contig_Length = 3e5,
                     Additional_Notes = "None"),
          "07_checkm2", "all_bins_checkm2_report.tsv")
write_tsv(data.frame(user_genome = mags.v, classification = gsub("\\|", ";", sub("\\|t__.*", "", mag_lineages.v)),
                     classification_method = "ANI", msa_percent = 95, red_value = "N/A", warnings = "N/A"),
          "15_gtdbtk", "all_genomes", "classify", "all_genomes.bac120.summary.tsv")

# inStrain and TRACS: mice in the same cage usually carry the same strain, mice in different cages rarely
cage.v <- stats::setNames(metadata.df$Cage, samples.v)
pairs.df <- data.frame(t(utils::combn(samples.v, 2)), stringsAsFactors = FALSE)
names(pairs.df) <- c("sample_a", "sample_b")
strain_mags.v <- mags.v[1:6]
sharing.df <- do.call(rbind, lapply(strain_mags.v, function(mag.s){
  compared.df <- pairs.df[stats::runif(nrow(pairs.df)) < 0.6, , drop = FALSE]
  same_cage.v <- cage.v[compared.df$sample_a] == cage.v[compared.df$sample_b]
  same_strain.v <- stats::runif(nrow(compared.df)) < ifelse(same_cage.v, 0.7, 0.05)
  popani.v <- ifelse(same_strain.v, stats::runif(nrow(compared.df), 0.999991, 1), stats::runif(nrow(compared.df), 0.9985, 0.99998))
  data.frame(genome = mag.s, sample_a = paste0(compared.df$sample_a, ".bam"), sample_b = paste0(compared.df$sample_b, ".bam"),
             popANI = round(popani.v, 7), conANI = round(popani.v - 1e-5, 7), percent_genome_compared = "",
             coverage_overlap = round(stats::runif(nrow(compared.df), 0.5, 1), 3), cluster_a = "1_1", cluster_b = "1_2",
             same_strain = ifelse(popani.v >= 0.99999, "TRUE", "FALSE"))
}))
write_tsv(sharing.df, "27_instrain", "summary", "strain_sharing_summary.tsv")
counts.df <- stats::aggregate(sharing.df$same_strain == "TRUE", list(sample_a = sharing.df$sample_a, sample_b = sharing.df$sample_b),
                              function(x) c(length(x), sum(x)))
write_tsv(data.frame(sample_a = counts.df$sample_a, sample_b = counts.df$sample_b, n_genomes_compared = counts.df$x[, 1],
                     n_same_strain = counts.df$x[, 2]),
          "27_instrain", "summary", "strain_sharing_counts.tsv")
treated_week4.v <- stats::setNames(metadata.df$Treatment == "Treated" & metadata.df$Time == "Week_4", samples.v)
for (sample.s in samples.v){
  # Treatment at week 4 narrows the strain populations (lower nucleotide diversity)
  diversity.n <- if (treated_week4.v[[sample.s]]) 0.0008 else 0.0016
  write_tsv(data.frame(genome = mags.v, coverage = round(stats::runif(10, 3, 60), 2), breadth = round(stats::runif(10, 0.7, 1), 3),
                       nucl_diversity = round(stats::rlnorm(10, log(diversity.n), 0.3), 6), length = 3e6,
                       breadth_minCov = round(stats::runif(10, 0.2, 0.98), 3),
                       popANI_reference = round(stats::runif(10, 0.997, 0.9999), 6),
                       SNV_count = stats::rpois(10, diversity.n * 3e6)),
            "27_instrain", "profiles", paste0(sample.s, ".IS"), "output", paste0(sample.s, ".IS_genome_info.tsv"))
}
tracs_pairs.df <- sharing.df[, c("genome", "sample_a", "sample_b", "same_strain")]
tracs.df <- data.frame(sampleA = sub("\\.bam$", "", tracs_pairs.df$sample_a), sampleB = sub("\\.bam$", "", tracs_pairs.df$sample_b),
                       `date difference` = NA, check.names = FALSE)
tracs.df$`SNP distance` <- ifelse(tracs_pairs.df$same_strain == "TRUE", stats::rpois(nrow(tracs.df), 3), stats::rpois(nrow(tracs.df), 400))
tracs.df$`transmission distance` <- NA
tracs.df$`expected K` <- NA
tracs.df$`filtered SNP distance` <- round(tracs.df$`SNP distance` * stats::runif(nrow(tracs.df), 0.6, 0.9))
tracs.df$`sites considered` <- round(stats::runif(nrow(tracs.df), 1e6, 2.5e6))
tracs.df$`MSA file` <- tracs_pairs.df$genome
dir.create(file.path(results.s, "28_tracs"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tracs.df, file.path(results.s, "28_tracs", "transmission_distances.csv"), row.names = FALSE, na = "NA")

# Comparison samples: 16 mice from two facilities, with the same taxa at facility-specific abundances,
# two genomes the study mice lack and some study genomes missing
comparison.df <- data.frame(Sample_ID = sprintf("F%02d", 1:16), Facility = rep(c("Facility_A", "Facility_B"), each = 8),
                            stringsAsFactors = FALSE)
comparison.df$Name <- sprintf("%s_%d", sub("Facility_", "Fac", comparison.df$Facility), rep(1:8, 2))
utils::write.csv(comparison.df[c("Sample_ID", "Name", "Facility")], file.path(root.s, "comparison_metadata.csv"), row.names = FALSE)
extra_lineages.v <- c(
  "d__Bacteria|p__Bacillota_A|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Dorea|s__Dorea sp100001|t__GCF_900000001.1",
  "d__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Muribaculaceae|g__Sodaliphilus|s__Sodaliphilus sp100002|t__GCF_900000002.1")
comparison_lineages.v <- c(lineages.v, extra_lineages.v)
baseline.v <- log(rowMeans(as.matrix(leaf.df[-1])) + 0.05)
facility.m <- cbind(stats::rnorm(length(lineages.v), sd = 0.8), stats::rnorm(length(lineages.v), sd = 0.8))
log.m <- sapply(seq_len(nrow(comparison.df)), function(j){
  c(baseline.v + facility.m[, 1 + (comparison.df$Facility[j] == "Facility_B")], c(1.5, 1.2)) + stats::rnorm(length(comparison_lineages.v), sd = 0.6)
})
abundance.m <- exp(log.m) * (matrix(stats::runif(length(log.m)), nrow(log.m)) < stats::plogis(log.m + 2))
abundance.m <- sweep(abundance.m, 2, colSums(abundance.m), "/") * 100
rows.l <- list()
for (depth.n in 1:8){
  prefix.v <- vapply(strsplit(comparison_lineages.v, "|", fixed = TRUE), function(x) paste(x[seq_len(depth.n)], collapse = "|"),
                     character(1))
  summed.m <- rowsum(abundance.m, prefix.v, reorder = FALSE)
  rows.l[[depth.n]] <- data.frame(clade_name = rownames(summed.m), round(summed.m, 5), check.names = FALSE)
}
comparison_sylph.df <- do.call(rbind, rows.l)
comparison_sylph.df <- comparison_sylph.df[rowSums(comparison_sylph.df[-1]) > 0, ]
names(comparison_sylph.df)[-1] <- paste0(comparison.df$Sample_ID, ".clean_1.fastq.gz")
write_tsv(comparison_sylph.df, "29_comparison_reads", "sylph", "merged_relative_abundance.tsv")
write_tsv(comparison_sylph.df, "29_comparison_reads", "sylph", "merged_sequence_abundance.tsv")
singlem_lineage.v <- vapply(strsplit(comparison_lineages.v, "|", fixed = TRUE),
                            function(x) paste(c("Root", x[1:6]), collapse = "; "), character(1))
coverage.m <- rowsum(abundance.m, singlem_lineage.v) * stats::runif(nrow(comparison.df), 0.2, 0.6)
write_tsv(data.frame(sample = rep(paste0(comparison.df$Sample_ID, "_R1"), each = nrow(coverage.m)),
                     coverage = round(as.vector(coverage.m), 2), taxonomy = rep(rownames(coverage.m), nrow(comparison.df)))[
            as.vector(coverage.m) > 0, ],
          "29_comparison_reads", "singlem", "metagenome.condensed.tsv")
cat("MAG, strain and comparison examples written to", root.s, "\n")

# ---- Metagenome example: comparison samples in the strain comparison ----------------------------------------
# As with --strain_include_comparison_reads: the 16 comparison mice were profiled by inStrain and TRACS too. Mice of
# one facility often share strains, mice of different facilities rarely, and a few study mice carry a Facility_A
# strain. A separate seed keeps the data above unchanged.
set.seed(2029)
read_example <- function(...) utils::read.delim(file.path(results.s, ...), check.names = FALSE, stringsAsFactors = FALSE,
                                                colClasses = "character")
facility.v <- stats::setNames(comparison.df$Facility, comparison.df$Sample_ID)
all_samples.v <- c(samples.v, comparison.df$Sample_ID)
write_tsv(data.frame(sample = all_samples.v, cohort = rep(c("main", "comparison"), c(length(samples.v), nrow(comparison.df)))),
          "26_strain_reference", "strain_samples.tsv")
pairs.df <- data.frame(t(utils::combn(all_samples.v, 2)), stringsAsFactors = FALSE)
names(pairs.df) <- c("sample_a", "sample_b")
pairs.df <- pairs.df[pairs.df$sample_b %in% comparison.df$Sample_ID, , drop = FALSE]
new_sharing.df <- do.call(rbind, lapply(strain_mags.v, function(mag.s){
  compared.df <- pairs.df[stats::runif(nrow(pairs.df)) < 0.5, , drop = FALSE]
  facility_a.v <- facility.v[compared.df$sample_a]
  same_facility.v <- !is.na(facility_a.v) & facility_a.v == facility.v[compared.df$sample_b]
  study_to_a.v <- is.na(facility_a.v) & facility.v[compared.df$sample_b] == "Facility_A" & mag.s == strain_mags.v[1]
  same_strain.v <- stats::runif(nrow(compared.df)) < ifelse(same_facility.v, 0.4, ifelse(study_to_a.v, 0.15, 0.01))
  popani.v <- ifelse(same_strain.v, stats::runif(nrow(compared.df), 0.999991, 1), stats::runif(nrow(compared.df), 0.9985, 0.99998))
  data.frame(genome = mag.s, sample_a = paste0(compared.df$sample_a, ".bam"), sample_b = paste0(compared.df$sample_b, ".bam"),
             popANI = round(popani.v, 7), conANI = round(popani.v - 1e-5, 7), percent_genome_compared = "",
             coverage_overlap = round(stats::runif(nrow(compared.df), 0.5, 1), 3), cluster_a = "1_1", cluster_b = "1_2",
             same_strain = ifelse(popani.v >= 0.99999, "TRUE", "FALSE"))
}))
sharing.df <- rbind(read_example("27_instrain", "summary", "strain_sharing_summary.tsv"), new_sharing.df)
write_tsv(sharing.df, "27_instrain", "summary", "strain_sharing_summary.tsv")
counts.df <- stats::aggregate(sharing.df$same_strain == "TRUE", list(sample_a = sharing.df$sample_a, sample_b = sharing.df$sample_b),
                              function(x) c(length(x), sum(x)))
write_tsv(data.frame(sample_a = counts.df$sample_a, sample_b = counts.df$sample_b, n_genomes_compared = counts.df$x[, 1],
                     n_same_strain = counts.df$x[, 2]),
          "27_instrain", "summary", "strain_sharing_counts.tsv")
for (sample.s in comparison.df$Sample_ID){
  write_tsv(data.frame(genome = mags.v, coverage = round(stats::runif(10, 3, 60), 2), breadth = round(stats::runif(10, 0.7, 1), 3),
                       nucl_diversity = round(stats::rlnorm(10, log(0.0014), 0.3), 6), length = 3e6,
                       breadth_minCov = round(stats::runif(10, 0.2, 0.98), 3),
                       popANI_reference = round(stats::runif(10, 0.997, 0.9999), 6), SNV_count = stats::rpois(10, 0.0014 * 3e6)),
            "27_instrain", "profiles", paste0(sample.s, ".IS"), "output", paste0(sample.s, ".IS_genome_info.tsv"))
}
new_tracs.df <- data.frame(sampleA = sub("\\.bam$", "", new_sharing.df$sample_a), sampleB = sub("\\.bam$", "", new_sharing.df$sample_b),
                           `date difference` = NA, check.names = FALSE)
new_tracs.df$`SNP distance` <- ifelse(new_sharing.df$same_strain == "TRUE", stats::rpois(nrow(new_tracs.df), 3),
                                      stats::rpois(nrow(new_tracs.df), 400))
new_tracs.df$`transmission distance` <- NA
new_tracs.df$`expected K` <- NA
new_tracs.df$`filtered SNP distance` <- round(new_tracs.df$`SNP distance` * stats::runif(nrow(new_tracs.df), 0.6, 0.9))
new_tracs.df$`sites considered` <- round(stats::runif(nrow(new_tracs.df), 1e6, 2.5e6))
new_tracs.df$`MSA file` <- new_sharing.df$genome
tracs_path.s <- file.path(results.s, "28_tracs", "transmission_distances.csv")
tracs.df <- utils::read.csv(tracs_path.s, check.names = FALSE, stringsAsFactors = FALSE)
utils::write.csv(rbind(tracs.df, new_tracs.df), tracs_path.s, row.names = FALSE, na = "NA")
cat("Comparison samples added to the strain examples in", root.s, "\n")

# One realistic path per registry entry plus decoys (large files gmaio does not read)
example_paths <- function(mode){
  shared.v <- c("pipeline_info/run_params.json", "pipeline_info/software_versions.tsv", "pipeline_info/report.html",
                "00_read_stats/read_stat_report.tsv", "00_read_stats/mapping_assessment.tsv",
                "00_read_stats/S1.raw_long.seqkit_stats.tsv",
                "18_gtdbtk/all_genomes/all_genomes.bac120.summary.tsv",
                "18_gtdbtk/all_genomes/classify/all_genomes.bac120.summary.tsv",
                "18_gtdbtk/all_genomes/classify/all_genomes.ar53.summary.tsv", "18_gtdbtk/all_genomes/align/x.fasta.gz",
                "20_genomespot/genomespot_combined.tsv", "20_genomespot/S1.predictions.tsv",
                "24_genomad/genomad_virus_summary.tsv", "24_genomad/genomad_plasmid_summary.tsv",
                "24_genomad/S1/S1_summary/S1_virus.fna", "25_checkv/S1/quality_summary.tsv", "25_checkv/S1/viruses.fna",
                "26_checkv_clustering/virus_clusters.tsv", "26_checkv_clustering/plasmid_clusters.tsv",
                "26_checkv_clustering/virus_ani.tsv", "cleanifier/S1.clean.fastq.gz")
  if (mode == "isolate"){
    return(c(shared.v, "07_checkm2/ISO1_checkm2_report.tsv", "07_checkm2/ISO1/protein_files/ISO1.faa",
             "10_bakta/bakta_annotation_stats.tsv", "10_bakta/ISO1/ISO1.gff3", "11_mlst/mlst_summary.csv",
             "12_amrfinder/amrfinder_all.tsv", "12_amrfinder/ISO1.tsv", "13_isescan/ISEScan_summary.tsv",
             "13_isescan/ISO1/ISO1.fasta.is.fna", "14_comparison_groups/all/entries.tsv",
             "15_panaroo/all_panaroo/gene_presence_absence.csv", "15_panaroo/all_panaroo/pan_genome_reference.fa",
             "18_fastani/all.fastani.tsv", "19_chewbacca/all_chewbbaca/cgMLST/cgMLST95.tsv",
             "19_chewbacca/all_chewbbaca/genome_hash_map.tsv", "19_chewbacca/all_chewbbaca/results_alleles.tsv",
             "20_tree/all.treefile", "20_tree/all.iqtree"))
  }
  c(shared.v, "04_sylph/sylph_profile.tsv", "04_sylph/merged_relative_abundance.tsv",
    "04_sylph/merged_sequence_abundance.tsv", "04_sylph/S1.sylsp", "05_singlem/metagenome.condensed.tsv",
    "05_singlem/metagenome.prokaryotic_fraction.tsv", "05_singlem/metagenome.otu_table.tsv", "05_singlem/S1.json",
    "06_nonpareil/S1.npo", "06_nonpareil/S1.npa", "07_myloasm/S1/assembly.fasta", "09_aviary/S1/bins/final_bins/b1.fna",
    "10_checkm2/all_bins_checkm2_report.tsv", "10_checkm2/all_bins/protein_files/S1.b1.faa",
    "19_checkm1/checkm_qa_results.tsv", "19_checkm1/checkm1/bins/S1.b1/genes.faa",
    "11_dereplicated_bins/cluster_definition.tsv", "11_dereplicated_bins/high_quality_representatives/cluster_definition.tsv",
    "11_dereplicated_bins/representatives/S1.b1.fna", "11_dereplicated_hq_bins/cluster_definition.tsv",
    "11_dereplicated_hq_ref_bins/cluster_definition.tsv", "12_coverm_bins/S1_abundances.tsv", "12_coverm_bins/S1.bam",
    "12_coverm_hq_bins/S1_abundances.tsv", "12_coverm_hq_derep_bins/S1_abundances.tsv",
    "12_coverm_hq_ref_bins/S1_abundances.tsv", "13_coverm_scaffolds_nanopore/S1_counts.tsv",
    "13_coverm_scaffolds_nanopore/S1.bam", "17_dram_bins/all_bins/distilled/product.tsv",
    "17_dram_bins/all_bins/distilled/product.html", "17_dram_bins/all_bins/distilled/metabolism_summary.xlsx",
    "17_dram_bins/all_bins/distilled/genome_stats.tsv", "17_dram_bins/all_bins/combined_annotations.tsv",
    "17_dram_bins/S1.b1/dram_annotations/annotations.tsv", "21_gene_catalogue/gene_catalogue_membership.tsv",
    "21_gene_catalogue/gene_catalogue.fna", "22_dram/dram_annotations/annotations.tsv", "22_dram/dram_annotations/genes.faa",
    "23_rpkm/gene_catalogue_rpkm_per_gene_normalised.tsv", "23_rpkm/gene_catalogue_mapped_reads_per_gene.tsv",
    "23_rpkm/singlem_rpkm_means.tsv")
}

make_tree <- function(root.s, paths.v){
  for (path.s in file.path(root.s, paths.v)){
    dir.create(dirname(path.s), recursive = TRUE, showWarnings = FALSE)
    writeLines("x", path.s)
  }
  root.s
}

relative_files <- function(root.s) sort(list.files(root.s, recursive = TRUE))

run_rsync <- function(filter.v, source.s, dest.s){
  filter.s <- withr::local_tempfile(lines = filter.v)
  system2(Sys.which("rsync"), shQuote(c("-a", "--prune-empty-dirs", paste0("--include-from=", filter.s), source.s, dest.s)))
}

test_that("the rsync filter copies exactly the registry outputs, with or without a trailing slash", {
  skip_if(!nzchar(Sys.which("rsync")), "rsync not available")
  for (mode.s in c("metagenome", "isolate")){
    results.s <- make_tree(file.path(withr::local_tempdir(), "results"), example_paths(mode.s))
    registry.df <- registry_for_mode(mode.s)
    matched.l <- lapply(seq_len(nrow(registry.df)), function(i){
      sub(paste0("^", results.s, "/"), "", search_results_dir(results.s, registry.df$dir[i], registry.df$file[i]))
    })
    covered.v <- lengths(matched.l) > 0
    expect_true(all(covered.v), info = paste(mode.s, "entries without an example:", toString(registry.df$key[!covered.v])))
    expected.v <- sort(unique(unlist(matched.l[!is.na(registry.df$transfer)])))

    slash.s <- withr::local_tempdir()
    run_rsync(rsync_filter(mode.s), paste0(results.s, "/"), slash.s)
    expect_equal(relative_files(slash.s), expected.v, info = mode.s)

    no_slash.s <- withr::local_tempdir()
    run_rsync(rsync_filter(mode.s), results.s, no_slash.s)
    expect_equal(relative_files(file.path(no_slash.s, "results")), expected.v, info = mode.s)

    all.s <- withr::local_tempdir()
    run_rsync(rsync_filter("all"), paste0(results.s, "/"), all.s)
    expect_true(all(expected.v %in% relative_files(all.s)), info = mode.s)
  }
})

test_that("write_rsync_filter writes the filter and explains the copy command", {
  path.s <- withr::local_tempfile(fileext = ".txt")
  expect_message(write_rsync_filter(path.s, mode = "metagenome"), "rsync -av --prune-empty-dirs")
  lines.v <- readLines(path.s)
  expect_true("+ [0-9]*_sylph/merged_relative_abundance.tsv" %in% lines.v)
  expect_false(any(grepl("gene_catalogue_membership|bakta", lines.v)))
  expect_equal(utils::tail(lines.v, 2), c("+ */", "- *"))
})

test_that("check_inputs reports complete, partial and missing steps", {
  inputs.df <- check_inputs(make_test_project(), verbose = FALSE)
  expect_equal(inputs.df$status[inputs.df$key == "sylph_relative"], "found")
  expect_equal(inputs.df$n_files[inputs.df$key == "coverm_hq_derep_bins"], 6)
  expect_false(any(c("gene_catalogue_membership", "rpkm_mapped_reads") %in% inputs.df$key))

  config.l <- make_test_project(pipeline_results = partial_results())
  expect_message(inputs.df <- check_inputs(config.l), "MAG abundance: not found")
  expect_equal(inputs.df$status[inputs.df$key == "sylph_relative"], "found")
  expect_equal(inputs.df$status[inputs.df$key == "read_stats"], "missing")
  expect_message(check_inputs(config.l), "Missing outputs are skipped")
})

test_that("a partial run processes what is there and skips the rest", {
  config.l <- make_test_project(pipeline_results = partial_results())
  project.l <- suppressMessages(start_project(config.l))
  expect_message(project.l <- add_sylph(project.l), NA)
  expect_message(project.l <- add_singlem(project.l), "Skipping SingleM")
  project.l <- suppressMessages(add_mags(project.l))
  project.l <- suppressMessages(add_gene_catalogue_functions(project.l))
  project.l <- suppressMessages(add_nonpareil(project.l))
  project.l <- suppressMessages(add_genomespot(project.l))
  project.l <- suppressMessages(add_mobile_elements(project.l))
  project.l <- suppressMessages(add_palettes(project.l))
  expect_setequal(names(project.l$profiles), c("sylph_taxonomic", "sylph_sequence"))
  expect_gt(nrow(project.l$tables$bin_summary), 0)
  expect_null(project.l$read_stats)
  written.v <- suppressMessages(write_project_tables(project.l))
  expect_true(all(file.exists(written.v)))
})

test_that("check_inputs points to outputs copied one level too deep", {
  config.l <- make_test_project(pipeline_results = fixture_dir())
  expect_message(check_inputs(config.l), "look like they are in")
})

test_that("fetch_pipeline_results copies the outputs into pipeline_results and can be rerun", {
  skip_if(!nzchar(Sys.which("rsync")), "rsync not available")
  source.s <- make_tree(file.path(withr::local_tempdir(), "results"), example_paths("metagenome"))
  project_dir.s <- withr::local_tempdir()
  config_path.s <- file.path(project_dir.s, "config.yml")
  yaml::write_yaml(list(mode = "metagenome", pipeline_results = "pipeline_results", pipeline_remote = source.s),
                   config_path.s)
  expect_error(read_config(config_path.s), "code/fetch_results.R")

  dest.s <- file.path(project_dir.s, "pipeline_results")
  suppressMessages(fetch_pipeline_results(config_path.s, dry_run = TRUE, rsync_args = "--quiet"))
  expect_length(list.files(dest.s, recursive = TRUE), 0)

  inputs.df <- suppressMessages(fetch_pipeline_results(config_path.s, rsync_args = "--quiet"))
  expect_true(all(inputs.df$status == "found"))
  expect_true(file.exists(file.path(dest.s, "04_sylph", "merged_relative_abundance.tsv")))
  expect_false(file.exists(file.path(dest.s, "12_coverm_bins", "S1.bam")))

  writeLines("x", file.path(source.s, "06_nonpareil", "S2.npo"))
  suppressMessages(fetch_pipeline_results(config_path.s, rsync_args = "--quiet"))
  expect_true(file.exists(file.path(dest.s, "06_nonpareil", "S2.npo")))
})

test_that("fetch_pipeline_results needs a source", {
  project_dir.s <- withr::local_tempdir()
  yaml::write_yaml(list(mode = "metagenome", pipeline_results = "pipeline_results"), file.path(project_dir.s, "config.yml"))
  expect_error(fetch_pipeline_results(file.path(project_dir.s, "config.yml")), "pipeline_remote")
})

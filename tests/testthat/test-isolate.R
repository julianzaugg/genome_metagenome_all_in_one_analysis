isolate_fixture_dir <- function() system.file("extdata", "isolate", package = "gmaio", mustWork = TRUE)

processed_isolate_project <- function(env = parent.frame()){
  project_dir.s <- withr::local_tempdir(.local_envir = env)
  yaml::write_yaml(list(mode = "isolate", pipeline_results = file.path(isolate_fixture_dir(), "results"),
                        metadata = list(file = file.path(isolate_fixture_dir(), "metadata.csv"), label_column = "Isolate"),
                        analysis = list(group_variables = list("Source"), reference_levels = list(Source = "Clinical"))),
                   file.path(project_dir.s, "config.yml"))
  project.l <- suppressMessages(start_project(read_config(file.path(project_dir.s, "config.yml"))))
  project.l <- suppressMessages(add_isolate_genomes(project.l))
  project.l <- suppressMessages(add_comparisons(project.l))
  project.l <- suppressMessages(add_mobile_elements(project.l))
  suppressMessages(add_palettes(project.l))
}

test_that("isolate readers link genomes and build the genome summary", {
  project.l <- processed_isolate_project()
  summary.df <- project.l$tables$genome_summary
  expect_equal(summary.df$Sample_ID, paste0("ISO", 1:6))
  expect_equal(summary.df$ST, c("131", "131", "73", "10", "10", "-"))
  expect_equal(summary.df$AMR_elements, c(3, 3, 1, 1, 1, 1))
  expect_equal(summary.df$IS_elements, c(2, 1, 1, 1, 2, 1))
  expect_true(all(c("Completeness", "Genome_size", "CDS", "Species", "Mean_coverage") %in% names(summary.df)))
  expect_equal(sort(rownames(project.l$profiles$amr_genes$values)), sort(c("blaCTX-M-15", "mdtM", "qnrS1", "sul1", "tet(A)")))
  expect_true("Class" %in% names(project.l$profiles$amr_genes$features))
})

test_that("comparison groups load pangenome, ANI, cgMLST and tree", {
  project.l <- processed_isolate_project()
  group.l <- project.l$tables$comparisons$all
  expect_equal(nrow(group.l$entries), 7)
  expect_equal(sum(group.l$pangenome$features$Category == "Core"), 8)
  ani.m <- ani_matrix(group.l$ani)
  expect_equal(unname(diag(ani.m)), rep(100, 7))
  expect_true(isSymmetric(unname(ani.m)))
  distance.m <- cgmlst_distances(group.l$cgmlst$cgMLST95)
  expect_equal(unname(distance.m["ISO4", "ISO5"]) <= 10, TRUE)
  expect_true(isSymmetric(distance.m))
  skip_if_not_installed("ape")
  expect_setequal(group.l$tree$tip.label, group.l$entries$Genome)
})

test_that("genome metadata marks references and gives them a Reference group", {
  project.l <- processed_isolate_project()
  genomes.df <- genome_metadata(project.l, "all")
  expect_equal(genomes.df$Entry_type, c(rep("Sample", 6), "Reference"))
  expect_equal(as.character(genomes.df$Source[7]), "Reference")
  expect_equal(genomes.df$Sample_label[1], "Strain_1")
  expect_equal(project.l$palettes$Source[["Reference"]], special_colours()[["Reference"]])
})

test_that("gene associations test presence against groups", {
  project.l <- processed_isolate_project()
  association.df <- gene_association(project.l$profiles$amr_genes, analysis_metadata(project.l), "Source", min_present = 1)
  expect_true(all(c("Log2_odds_ratio", "P_value", "Q_value", "Enriched_in") %in% names(association.df)))
  expect_equal(association.df$Enriched_in[association.df$Feature_ID == "blaCTX-M-15"], "Clinical")
})

test_that("isolate figures and tables are produced", {
  project.l <- processed_isolate_project()
  out.s <- withr::local_tempdir()
  group.l <- project.l$tables$comparisons$all
  genomes.df <- genome_metadata(project.l, "all")
  save_plot(plot_genome_matrix(ani_matrix(group.l$ani), genomes.df, project.l$palettes, c("Entry_type", "Source")),
            file.path(out.s, "ani.pdf"))
  save_plot(plot_genome_matrix(cgmlst_distances(group.l$cgmlst$cgMLST95), genomes.df, project.l$palettes, "Source",
                               type = "distance"), file.path(out.s, "cgmlst.pdf"))
  save_plot(plot_presence_heatmap(group.l$pangenome, genomes.df, project.l$palettes, "Source"), file.path(out.s, "pangenome.pdf"))
  save_plot(plot_pangenome_categories(group.l$pangenome), file.path(out.s, "categories.pdf"))
  save_plot(plot_genome_summary(project.l$tables$genome_summary, analysis_metadata(project.l), "Source",
                                project.l$palettes$Source), file.path(out.s, "summary.pdf"))
  ordination <- ordinate_distance(100 - ani_matrix(group.l$ani), "100 - ANI")
  save_plot(plot_ordination(ordination, genomes.df, "Source", project.l$palettes$Source), file.path(out.s, "pcoa.pdf"))
  if (requireNamespace("ggtree", quietly = TRUE)){
    save_plot(plot_tree(group.l$tree, genomes.df, "Source", project.l$palettes$Source), file.path(out.s, "tree.pdf"))
  }
  expect_true(all(file.size(list.files(out.s, full.names = TRUE)) > 1000))
  written.v <- suppressMessages(write_isolate_tables(project.l))
  expect_true(all(file.exists(written.v)))
  expect_true(any(grepl("Comparison_all.xlsx$", written.v)))
})

test_that("the isolate example links typing to genomes and drives the isolate figures", {
  project.l <- example_project(mode = "isolate")
  expect_equal(nrow(project.l$metadata), 20)
  group.l <- project.l$tables$comparisons$all
  ani.m <- ani_matrix(group.l$ani)
  genomes.df <- genome_metadata(project.l, "all", rownames(ani.m))
  expect_setequal(unique(genomes.df$ST), c("ST131", "ST73", "ST10", "ST69", "Reference"))
  expect_equal(genomes.df$ST[genomes.df$Entry_type == "Reference"], c("Reference", "Reference"))

  split.ht <- plot_genome_matrix(ani.m, genomes.df, project.l$palettes, "Source", split_by = "ST")
  expect_equal(levels(split.ht@matrix_param$row_split[[1]]), c("ST10", "ST69", "ST73", "ST131", "Reference"))
  expect_error(plot_genome_matrix(ani.m, genomes.df, split_by = "Nope"), "Nope")
  expect_error(plot_genome_matrix(ani.m, genomes.df, breaks = c(90, 100), colours = "red"), "one colour per break")
  presence.ht <- plot_presence_heatmap(get_profile(project.l, "amr_genes"), analysis_metadata(project.l), project.l$palettes,
                                       "Source", row_split = "Class", column_split = "Source")
  expect_s4_class(presence.ht, "Heatmap")

  skip_if_not_installed("ggtree")
  tree_genomes.df <- genome_metadata(project.l, "all", group.l$tree$tip.label)
  tree.gg <- plot_tree(group.l$tree, tree_genomes.df, "Source", project.l$palettes$Source, shape_by = "Entry_type",
                       align_labels = TRUE)
  expect_silent(built.l <- ggplot2::ggplot_build(tree.gg))
  shape_scale <- tree.gg$scales$get_scales("shape")
  expect_equal(shape_scale$palette(2)[c("Sample", "Reference")], c(Sample = 21, Reference = 23))
  expect_s3_class(plot_tree(group.l$tree, tree_genomes.df, layout = "circular", label_by = NULL), "ggplot")
})

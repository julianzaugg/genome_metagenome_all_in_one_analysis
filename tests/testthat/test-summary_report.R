test_that("the summary report is off by default and records nothing then", {
  config.l <- make_test_project()
  expect_false(report_enabled(config.l))
  expect_null(record_summary(config.l, "diversity", "Alpha diversity", tables = list(Tests = data.frame(x = 1))))
  expect_false(dir.exists(project_path(config.l, "Result_other", "summary")))
  expect_error(make_test_project(outputs = list(report = "yes")), "outputs\\$report")
})

test_that("recorded sections and the processed project make one self-contained report", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available("2.0"))
  project.l <- processed_test_project()
  project.l$config$outputs$report <- TRUE
  config.l <- project.l$config
  suppressMessages(save_processed(project.l))

  plot.gg <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(.data$x, .data$y)) + ggplot2::geom_point()
  figure.s <- save_plot(plot.gg, figure_path(config.l, "diversity", "sylph_taxonomic__species"), 8, 6)
  table.s <- output_path(config.l, "tables", "Alpha_diversity.xlsx")
  suppressMessages(write_xlsx_tables(list(Values = data.frame(x = 1)), table.s))
  # "-" (e.g. an unknown MLST type) must stay a cell, not become a markdown list
  tests.df <- data.frame(Measure = c("Richness", "Shannon"), Type = c("-", "1"), P_value = c(0.01234, 0.0002), Test_label = "p")
  record.s <- suppressMessages(record_summary(config.l, "diversity", "Alpha diversity",
                                              tables = list(`Group tests` = tests.df, Empty = data.frame()),
                                              figures = list(summary_figure(plot.gg, figure.s, "Diversity figure")),
                                              files = table.s))
  expect_true(file.exists(record.s))
  expect_true(file.exists(file.path(dirname(record.s), "figures", "diversity__01.png")))
  # A record older than processed.rds is flagged as possibly out of date
  stale.l <- list(section = "strains", title = "Strains", tables = list(Test = data.frame(P_value = 0.5)), script = "strains.R",
                  created = Sys.time() - 3600)
  saveRDS(stale.l, file.path(dirname(record.s), "strains.rds"))

  report.s <- suppressMessages(write_summary_report(project.l))
  expect_equal(basename(report.s), "Summary_report.html")
  html.s <- paste(readLines(report.s, warn = FALSE), collapse = "\n")
  # expect_match() would print the whole page, scripts included, on failure
  titles.v <- c("Overview", "Samples and reads", "Bins and MAGs", "Alpha diversity", "Strains", "Outputs")
  expected.v <- c(paste0(">", titles.v, "</h2>"),
                  "data:image/png;base64,", "href=\"../Result_figures/diversity/sylph_taxonomic__species.pdf\"",
                  "href=\"../Result_tables/Alpha_diversity.xlsx\"", "> 0.012 <", "> - <", "&lt; 0.001",
                  "Recorded before the latest main.R run; rerun strains.R",
                  "No results recorded from mag_summary.R, ordination.R")
  for (text.s in expected.v) expect_true(grepl(text.s, html.s, fixed = TRUE), label = text.s)
  body.s <- substring(html.s, regexpr("<div id=\"overview\"", html.s, fixed = TRUE))
  for (text.s in c("Test label", ">Empty<", "<li>")) expect_false(grepl(text.s, body.s, fixed = TRUE), label = text.s)
})

test_that("the read table has one median per genome set", {
  project.l <- list(config = list(analysis = list(group_variables = character()), project_dir = tempdir(),
                                  outputs = list(tables = "Result_tables")),
                    metadata = data.frame(Sample_ID = c("a", "b"), Excluded = FALSE, row.names = c("a", "b")),
                    read_stats = data.frame(Sample_ID = c("a", "b"), Raw_count = c(100, 200), Clean_count = c(90, 180),
                                            Reads_mapped_Set_A_percent = c(10, 20), Reads_mapped_Set_B_percent = c(50, 70)))
  reads.df <- report_reads(project.l)$tables$Reads
  mapped.v <- reads.df$All[match(c("Mapped to Set A (median %)", "Mapped to Set B (median %)"), reads.df$Measure)]
  expect_equal(mapped.v, c("15", "60"))
  expect_equal(reads.df$All[reads.df$Measure == "Kept after QC (median %)"], "90")
})

test_that("summary tables combine statistics compactly", {
  stacked.df <- stack_tables(list(a_tests = data.frame(P_value = 0.1), b_tests = data.frame(P_value = 0.2), c_tests = NULL),
                             strip = "_tests$")
  expect_equal(stacked.df$Dataset, c("a", "b"))
  expect_null(stack_tables(list()))

  permanova.df <- data.frame(Label = rep(c("x", "y"), each = 3), Term = rep(c("Group", "Residual", "Total"), 2), N = 10,
                             R2 = c(0.2, 0.8, 1, 0.3, 0.7, 1), F = c(2, NA, NA, 3, NA, NA),
                             P_value = c(0.01, NA, NA, 0.2, NA, NA))
  permdisp.df <- data.frame(Label = c("x", "x"), Term = c("Groups", "Residuals"), P_value = c(0.5, NA))
  summary.df <- summary_permanova(permanova.df, permdisp.df)
  expect_equal(summary.df$Dataset, c("x", "y"))
  expect_equal(summary.df$PERMDISP_P_value, c(0.5, NA))
  expect_true(is.na(summary_permanova(permanova.df, permdisp.df, permdisp_term = "Other")$PERMDISP_P_value[1]))

  consensus.l <- list(called = data.frame(Feature_ID = c("f1", "f2", "f3"), Label = c("A", "B", "C"),
                                          Enriched_in = "T", N_methods = c(1, 3, 2),
                                          Conflicting_direction = c(FALSE, FALSE, TRUE)),
                      none = data.frame())
  da.l <- summary_da(consensus.l, top_n = 2)
  expect_equal(da.l$Calls$Features_called, c(3, 0))
  expect_equal(da.l$Calls$By_2_or_more_methods, c(2, 0))
  expect_equal(da.l$Calls$Conflicting, c(1, 0))
  expect_equal(da.l$`Top features`$Label, c("B", "C"))
  expect_equal(summary_da(list()), list())
})

test_that("report values are formatted for reading", {
  formatted.df <- format_report_table(data.frame(P_adjusted = c(0.0004, 0.04567, NA), R2 = c(0.12345, 2, 1234567),
                                                 Flag = c(TRUE, FALSE, NA), Test_label = "x"))
  expect_equal(names(formatted.df), c("P adjusted", "R2", "Flag"))
  expect_equal(formatted.df$`P adjusted`, c("< 0.001", "0.046", ""))
  expect_equal(formatted.df$R2, c("0.123", "2", "1,234,567"))
  expect_equal(formatted.df$Flag, c("yes", "no", ""))
  expect_equal(format_count(5386892), "5.4M")
})

test_that("png sizes are read from the header", {
  path.s <- withr::local_tempfile(fileext = ".png")
  grDevices::png(path.s, width = 300, height = 100)
  grid::grid.newpage()
  grDevices::dev.off()
  expect_equal(png_aspect_ratio(readBin(path.s, "raw", file.info(path.s)$size)), 3)
})

test_that("links are relative to the report folder", {
  root.s <- withr::local_tempdir()
  dir.create(file.path(root.s, "Result_other"))
  dir.create(file.path(root.s, "Result_figures", "mags"), recursive = TRUE)
  expect_equal(relative_path(file.path(root.s, "Result_figures", "mags", "a.pdf"), file.path(root.s, "Result_other")),
               "../Result_figures/mags/a.pdf")
  expect_equal(relative_path(file.path(root.s, "Result_other", "b.png"), file.path(root.s, "Result_other")), "b.png")
})

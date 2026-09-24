test_that("every table and figure folder the templates write is in the catalogue", {
  scripts.v <- list.files(system.file("templates", package = "gmaio"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  code.s <- paste(unlist(lapply(scripts.v, readLines)), collapse = "\n")
  tables.v <- unique(regmatches(code.s, gregexpr('output_path\\(config.l, "tables", "[^"]+\\.xlsx"', code.s))[[1]])
  tables.v <- sub('.*"tables", "', "", sub('"$', "", tables.v))
  expect_gt(length(tables.v), 5)
  for (file.s in tables.v) expect_false(is.null(describe_output(file.s)), label = file.s)
  figure_calls.v <- regmatches(code.s, gregexpr('figure_path\\(config.l, "[a-z_]+"', code.s))[[1]]
  folders.v <- unique(sub('.*"([a-z_]+)"$', "\\1", figure_calls.v))
  expect_setequal(folders.v, figure_catalogue()$Folder)
  expect_equal(mag_set_catalogue()$Set, names(mag_sets()))
})

test_that("gmaio tables get an About sheet, described sheets and a README index", {
  project.l <- processed_test_project()
  written.v <- suppressMessages(write_project_tables(project.l))
  for (path.s in written.v){
    expect_equal(readxl::excel_sheets(path.s)[1], "About", label = basename(path.s))
    entry.l <- describe_output(path.s)
    data_sheets.v <- readxl::excel_sheets(path.s)[-1]
    described.v <- vapply(data_sheets.v, function(sheet.s) any(vapply(entry.l$Sheets$Pattern, grepl, logical(1), sheet.s)),
                          logical(1))
    expect_true(all(described.v), label = paste(basename(path.s), toString(data_sheets.v[!described.v])))
  }
  about.df <- readxl::read_excel(grep("MAG_hq_derep_bins_abundances", written.v, value = TRUE), "About", col_names = FALSE)
  expect_true("Genome set: hq_derep_bins" %in% about.df[[1]])
  readme.v <- readLines(file.path(dirname(written.v[1]), "README.md"))
  expect_true(all(paste0("### ", basename(written.v)) %in% readme.v))
  expect_true("## Which MAG file to use" %in% readme.v)
  expect_false(any(grepl("ws_hq_bins", readme.v[seq_len(which(readme.v == "## Files"))]) &
                 !any(grepl("^MAG_ws_hq_bins_", basename(written.v)))))
})

test_that("isolate tables are described too", {
  written.v <- suppressMessages(write_isolate_tables(processed_isolate_project()))
  expect_true(all(vapply(written.v, function(p) readxl::excel_sheets(p)[1] == "About", logical(1))))
})

test_that("other workbooks are left alone", {
  dir.s <- withr::local_tempdir()
  write_xlsx_tables(list(Data = data.frame(a = 1)), file.path(dir.s, "my_table.xlsx"))
  expect_equal(readxl::excel_sheets(file.path(dir.s, "my_table.xlsx")), "Data")
  expect_false(file.exists(file.path(dir.s, "README.md")))
  write_xlsx_tables(list(Data = data.frame(a = 1)), file.path(dir.s, "Alpha_diversity.xlsx"), about = FALSE)
  expect_false(file.exists(file.path(dir.s, "README.md")))
})

test_that("saving a figure to a gmaio folder keeps the figure README current", {
  dir.s <- withr::local_tempdir()
  figure.gg <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(.data$x, .data$y)) + ggplot2::geom_point()
  save_plot(figure.gg, file.path(dir.s, "diversity", "test.png"), width = 5, height = 5)
  readme.v <- readLines(file.path(dir.s, "README.md"))
  expect_true("## diversity/" %in% readme.v)
  # A folder of your own is listed under other folders the next time a gmaio figure is saved
  save_plot(figure.gg, file.path(dir.s, "my_folder", "test.png"), width = 5, height = 5)
  save_plot(figure.gg, file.path(dir.s, "diversity", "test.png"), width = 5, height = 5)
  expect_true("- `my_folder/`" %in% readLines(file.path(dir.s, "README.md")))
  other.s <- withr::local_tempdir()
  save_plot(figure.gg, file.path(other.s, "plot.png"), width = 5, height = 5)
  expect_false(file.exists(file.path(other.s, "README.md")))
})

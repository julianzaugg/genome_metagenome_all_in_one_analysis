test_that("write_xlsx_tables writes styled sheets including empty ones", {
  path.s <- withr::local_tempfile(fileext = ".xlsx")
  write_xlsx_tables(list(Data = data.frame(a = 1:3, b = c("x", "y", "z")), Empty = data.frame(a = numeric())), path.s)
  expect_equal(readxl::excel_sheets(path.s), c("Data", "Empty"))
  expect_equal(nrow(readxl::read_excel(path.s, "Data")), 3)
  expect_equal(names(readxl::read_excel(path.s, "Empty")), "a")
})

test_that("write_xlsx_tables rejects unnamed lists and truncation clashes", {
  path.s <- withr::local_tempfile(fileext = ".xlsx")
  expect_error(write_xlsx_tables(list(data.frame(a = 1)), path.s), "named")
  long.s <- strrep("a", 35)
  tables.l <- stats::setNames(list(data.frame(a = 1), data.frame(a = 2)), c(paste0(long.s, "1"), paste0(long.s, "2")))
  expect_error(write_xlsx_tables(tables.l, path.s), "not unique")
})

test_that("an empty table list writes nothing", {
  path.s <- file.path(withr::local_tempdir(), "empty.xlsx")
  expect_message(result <- write_xlsx_tables(list(), path.s), "No tables")
  expect_null(result)
  expect_false(file.exists(path.s))
})

test_that("read statistics are written as a readable funnel with the full report and descriptions", {
  report.df <- data.frame(Sample_ID = c("S1", "S2"), GBbp = c(9.6, 15.6), Raw_count = c(3e6, 6e6),
                          Porechop_count = c(2.9e6, 5.9e6), Porechop_percent = c(99.7, 99.7),
                          Reads_mapped_HQ_Derep_MAGs_count = c(19428, 770009),
                          Reads_mapped_HQ_Derep_MAGs_percent = c(0.64, 12.94),
                          Bases_mapped_A_HQ_Derep_MAGs_count = c(6e7, 2.3e9),
                          Bases_mapped_A_HQ_Derep_MAGs_percent = c(0.63, 14.81),
                          Bases_mapped_B_HQ_Derep_MAGs_percent = c(0.63, 14.77),
                          Bases_mapped_C_HQ_Derep_MAGs_percent = c(0.63, 14.79))
  metadata.df <- data.frame(Sample_ID = c("S1", "S2"), Sample_label = c("Mouse 1", "Mouse 2"))
  tables.l <- read_stats_tables(report.df, metadata.df)
  expect_equal(names(tables.l$Read_stats),
               c("Sample_ID", "Sample_label", "GBbp", "Raw_count", "Porechop_count", "Porechop_percent",
                 "Reads_mapped_HQ_Derep_MAGs_count", "Reads_mapped_HQ_Derep_MAGs_percent", "Bases_mapped_HQ_Derep_MAGs_percent"))
  expect_equal(tables.l$Read_stats$Bases_mapped_HQ_Derep_MAGs_percent, c(0.63, 14.81))
  expect_equal(tables.l$Read_stats_full, report.df)
  descriptions.df <- tables.l$Column_descriptions
  expect_setequal(descriptions.df$Column[descriptions.df$Sheet == "Read_stats_full"],
                  grep("^Bases_mapped_[ABC]_", names(report.df), value = TRUE))
  expect_false(any(descriptions.df$Description == "From the pipeline read statistics report"))
  expect_match(descriptions.df$Description[descriptions.df$Column == "Bases_mapped_HQ_Derep_MAGs_percent"],
               "full length.*% of raw sequenced bases")
  expect_match(descriptions.df$Description[descriptions.df$Column == "Porechop_percent"], "adapter trimming.*% of raw reads")

  short.l <- read_stats_tables(report.df, bases = FALSE)
  expect_false(any(grepl("^Bases_mapped", names(short.l$Read_stats))))
  plain.l <- read_stats_tables(report.df[, 1:5])
  expect_null(plain.l$Read_stats_full)
  expect_null(read_stats_tables(NULL))
})

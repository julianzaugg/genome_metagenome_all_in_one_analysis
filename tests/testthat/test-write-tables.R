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

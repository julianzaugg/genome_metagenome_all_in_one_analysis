local_build <- function(loaded, on_disk, env = parent.frame()){
  old.l <- as.list(gm_build)
  gm_build$loaded <- loaded
  withr::defer(list2env(old.l, envir = gm_build), envir = env)
  testthat::local_mocked_bindings(is_dev_session = function() FALSE, installed_build = function(...) on_disk, .env = env)
}

test_that("load checks are quiet when session, install and processed.rds agree", {
  local_build("0.1.0 aaaaaaa", "0.1.0 aaaaaaa")
  expect_no_warning(check_build(list(gmaio_build = "0.1.0 aaaaaaa")))
})

test_that("a reinstall after the session loaded gmaio asks for an R restart", {
  local_build("0.1.0 aaaaaaa", "0.1.0 bbbbbbb")
  expect_warning(check_build(list(gmaio_build = "0.1.0 aaaaaaa")), "Restart R")
})

test_that("processed.rds from another build asks for main.R to be rerun", {
  local_build("0.1.0 bbbbbbb", "0.1.0 bbbbbbb")
  expect_warning(check_build(list(gmaio_build = "0.1.0 aaaaaaa")), "main.R")
  expect_warning(check_build(list()), "an older gmaio")
})

test_that("load checks are skipped under devtools::load_all()", {
  local_build("0.1.0 aaaaaaa", "0.1.0 bbbbbbb")
  testthat::local_mocked_bindings(is_dev_session = function() TRUE)
  expect_no_warning(check_build(list()))
})

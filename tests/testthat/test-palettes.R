test_that("assign_colours is deterministic and respects precedence", {
  colours.v <- assign_colours(c("b", "a", "c"))
  expect_equal(names(colours.v), c("a", "b", "c"))
  expect_identical(colours.v, assign_colours(c("c", "a", "b")))
  fixed.v <- assign_colours(c("a", "b"), existing.v = c(a = "#000001", b = "#000002"), fixed.v = c(a = "#FFFFFF"))
  expect_equal(unname(fixed.v), c("#FFFFFF", "#000002"))
})

test_that("new values never reuse assigned colours", {
  colours.v <- assign_colours(c("a", "b", "c"), existing.v = c(a = "#E69F00"))
  expect_equal(anyDuplicated(toupper(colours.v)), 0L)
})

test_that("qualitative palettes contain no greys and no duplicates", {
  for (n in c(5, 10, 20, 30, 45)){
    palette.v <- qualitative_palette(n)
    expect_length(palette.v, n)
    expect_equal(anyDuplicated(palette.v), 0L)
    expect_false(any(toupper(palette.v) %in% toupper(special_colours())))
  }
})

test_that("taxa colours: Other and Unassigned are grey, lower ranks get colours", {
  taxa.df <- data.frame(Label = c("p__A", "p__B", "(p__A) g__X", "(p__A) g__Y", "Unclassified d__Bacteria"),
                        Phylum = c("p__A", "p__B", "p__A", "p__A", NA), Abundance = c(50, 20, 30, 10, 5))
  for (scheme.s in c("distinct", "hierarchical")){
    colours.v <- assign_taxa_colours(taxa.df, scheme = scheme.s)
    expect_equal(colours.v[["Other"]], special_colours()[["Other"]])
    expect_equal(colours.v[["Unclassified d__Bacteria"]], special_colours()[["Unassigned"]])
    expect_true(all(c("(p__A) g__X", "(p__A) g__Y") %in% names(colours.v)))
    expect_false(colours.v[["(p__A) g__X"]] == colours.v[["(p__A) g__Y"]])
  }
})

test_that("palettes round-trip through palettes.yml", {
  path.s <- withr::local_tempfile(fileext = ".yml")
  palettes.l <- list(Treatment = c(Control = "#123456", Treated = "#654321"), taxa = c(`p__A` = "#ABCDEF"))
  write_palettes(palettes.l, path.s)
  expect_equal(read_palettes(path.s), palettes.l)
})

test_that("metadata variables avoid each other's colours", {
  metadata.df <- data.frame(Sample_ID = paste0("S", 1:4), Sample_label = paste0("S", 1:4),
                            Treatment = c("Control", "Control", "Treated", "Treated"), Time = c("W0", "W4", "W0", "W4"))
  palettes.l <- build_palettes(metadata.df, c("Treatment", "Time"))
  expect_length(intersect(palettes.l$Treatment, palettes.l$Time), 0)
  saved.l <- build_palettes(metadata.df, c("Treatment", "Time"), existing.l = list(Time = c(W0 = "#E69F00", W4 = "#56B4E9")))
  expect_equal(unname(saved.l$Time), c("#E69F00", "#56B4E9"))
})

test_that("split_taxonomy handles GTDB, SingleM and sylph styles", {
  taxonomy.df <- split_taxonomy(c("d__Bacteria;p__Bacillota;g__",
                                  "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia",
                                  "d__Archaea|p__Methanobacteriota|c__M|o__M|f__M|g__Methanobrevibacter|s__M smithii|t__GCF_1"))
  expect_named(taxonomy.df, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  expect_equal(taxonomy.df$Phylum, c("p__Bacillota", "p__Bacteroidota", "p__Methanobacteriota"))
  expect_true(is.na(taxonomy.df$Genus[1]))
  expect_equal(taxonomy.df$Species[3], "s__M smithii")
})

test_that("unprefixed tokens take the next rank", {
  expect_equal(split_taxonomy("Root; d__Eukaryota; Fungi")$Phylum, "Fungi")
})

test_that("taxonomy_string keeps unresolved lineages under their parent", {
  taxonomy.df <- split_taxonomy(c("d__Bacteria;p__Bacillota;c__Clostridia", "d__Bacteria"))
  expect_equal(taxonomy_string(taxonomy.df, "class"),
               c("d__Bacteria;p__Bacillota;c__Clostridia", "d__Bacteria;Unassigned;Unassigned"))
})

test_that("taxon_label gives resolved, unclassified and unassigned labels", {
  labels.v <- taxon_label(c("d__Bacteria;p__Bacillota_A;c__C;o__O;f__Lachnospiraceae;g__Blautia",
                            "d__Bacteria;p__Bacillota_A;c__C;o__O;f__Lachnospiraceae;Unassigned",
                            "d__Bacteria;p__Bacillota_A",
                            "d__Bacteria;Unassigned",
                            "Other"))
  expect_equal(labels.v, c("(p__Bacillota_A) g__Blautia", "Unclassified (p__Bacillota_A) f__Lachnospiraceae",
                           "p__Bacillota_A", "Unclassified d__Bacteria", "Other"))
})

test_that("clashing labels are made unique", {
  labels.v <- taxon_label(c("d__A;p__X;c__C1;o__O;f__F;g__G", "d__A;p__X;c__C2;o__O;f__F;g__G"))
  expect_equal(anyDuplicated(labels.v), 0L)
})

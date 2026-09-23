# gmaio

R package for downstream analysis of [genome_metagenome_all_in_one](https://github.com/julianzaugg/genome_metagenome_all_in_one) pipeline output, plus project templates.

## Layout

- `R/` package functions; `inst/templates/{metagenome,isolate}/` project templates (`config.yml`, `code/main.R`, one script per analysis).
- `R/pipeline_registry.R` is the single map of pipeline outputs (directory and file patterns, matched by name suffix so mode-specific numbering does not matter).
- `inst/extdata/metagenome` (small, for tests) and `inst/extdata/example` (larger, for vignettes and `example_project()`) are synthetic pipeline output built by `data-raw/make_fixtures.R` and `data-raw/make_example.R`. Never commit client data.
- `vignettes/` worked examples on `example_project()`; update them when plot function options change.
- `code_examples_dump/` holds old ad-hoc scripts for reference; it is gitignored.

## Design rules

- All logic lives in the package; templates only call `gmaio::` functions with a short parameter block.
- `main.R` does bulk processing and saves `Result_other/processed.rds`; analysis scripts start with `gmaio::load_processed()` and are independent of each other.
- Every abundance or function table is a `gm_profile` (features x samples matrix plus a `features` data frame with `Feature_ID` and `Label`).
- Samples are linked to metadata only through `resolve_sample_ids()` / `link_samples()`; never add project-specific regexes.
- Colours come from the project palettes (`get_palette()`, `project.l$palettes`); precedence is config `colours:` > `palettes.yml` > generated.
- Tables are written with `write_xlsx_tables()`; figures with `save_plot()` (sizes in cm).
- Heavy optional packages (maaslin3, MicrobiomeStat, mixOmics, Nonpareil, ggtree, ape) are in Suggests and checked with `require_pkg()`.

## Style

- Object suffixes `.df`, `.m`, `.v`, `.l`, `.s`, `.n`, `.b`, `.p` (profile), `.gg`, `.ht`; snake_case functions; `function(x){` brace style.
- Native pipe `|>`; explicit `pkg::` calls; no `library()` in package code.
- Informative but brief comments; roxygen docs for every exported function.

## Commands

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::check()'
Rscript data-raw/make_fixtures.R
Rscript data-raw/make_example.R
```

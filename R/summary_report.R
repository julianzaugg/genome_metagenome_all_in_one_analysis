#' Whether the summary report is switched on
#'
#' The summary report is optional and off by default; set `outputs: report: true` in
#' `config.yml` to have the analysis scripts record their headline results for it.
#'
#' @param config.l A `gm_config`.
#' @return `TRUE` or `FALSE`.
#' @export
report_enabled <- function(config.l){
  isTRUE(config.l$outputs$report)
}

#' A figure to show in the summary report
#'
#' Holds a figure and the path it was saved to. [record_summary()] draws a small png of it
#' for the report, which links to the saved figure.
#'
#' @param plot A figure accepted by [save_plot()].
#' @param path Path the full figure was saved to (e.g. from [figure_path()]).
#' @param caption Caption; defaults to the file name.
#' @param width,height Size of the thumbnail in centimetres; `NULL` uses the figure's own size.
#' @return A `gm_summary_figure`.
#' @export
summary_figure <- function(plot, path, caption = NULL, width = NULL, height = NULL){
  structure(list(plot = plot, path = path, caption = caption %||% tools::file_path_sans_ext(basename(path)),
                 width = width, height = height), class = "gm_summary_figure")
}

#' Record the headline results of an analysis for the summary report
#'
#' Saves a small record to `Result_other/summary/<section>.rds`, replacing the one from an
#' earlier run, and a png thumbnail of each figure. [write_summary_report()] collects the records.
#' Does nothing unless the report is switched on ([report_enabled()]), so analysis scripts can
#' always call it.
#'
#' @param config.l A `gm_config`.
#' @param section Short name of the section, e.g. `"diversity"`; letters, digits and underscores.
#' @param title Section title.
#' @param tables Named list of data frames (names are the table captions).
#' @param values Named list of single values shown as a short list (e.g. counts).
#' @param figures List of [summary_figure()]s.
#' @param files Paths of the full output files, linked from the section.
#' @param notes Character vector of short notes shown under the section.
#' @return The record path invisibly, or `NULL` when the report is off.
#' @export
record_summary <- function(config.l, section, title, tables = list(), values = list(), figures = list(),
                           files = character(), notes = character()){
  if (!report_enabled(config.l)) return(invisible(NULL))
  if (!grepl("^[A-Za-z0-9_]+$", section)) cli::cli_abort("{.arg section} must be letters, digits and underscores")
  tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, tables)
  values <- Filter(function(x) length(x) == 1 && !is.na(x), values)
  record_path.s <- output_path(config.l, "other", "summary", paste0(section, ".rds"))
  figure_dir.s <- file.path(dirname(record_path.s), "figures")
  dir.create(figure_dir.s, showWarnings = FALSE)
  unlink(list.files(figure_dir.s, pattern = paste0("^", section, "__.*\\.png$"), full.names = TRUE))
  figures.l <- list()
  for (figure.l in Filter(function(x) inherits(x, "gm_summary_figure"), figures)){
    thumbnail.s <- file.path(figure_dir.s, sprintf("%s__%02d.png", section, length(figures.l) + 1))
    # A thumbnail that cannot be drawn should not stop the analysis that recorded it
    saved.s <- tryCatch(save_plot(figure.l$plot, thumbnail.s, figure.l$width, figure.l$height, dpi = 110),
                        error = function(e){
                          cli::cli_warn("Summary report thumbnail for {.path {figure.l$path}} failed: {conditionMessage(e)}")
                          NULL
                        })
    if (is.null(saved.s)) next
    figures.l[[length(figures.l) + 1]] <- list(caption = figure.l$caption, thumbnail = thumbnail.s,
                                               path = normalizePath(figure.l$path, mustWork = FALSE))
  }
  record.l <- list(section = section, title = title, tables = tables, values = values, figures = figures.l,
                   files = normalizePath(files, mustWork = FALSE), notes = notes, script = current_script(),
                   created = Sys.time(), gmaio_version = as.character(utils::packageVersion("gmaio")))
  saveRDS(record.l, record_path.s)
  cli::cli_inform(c("v" = "Recorded {.val {section}} for the summary report"))
  invisible(record_path.s)
}

# The script run with Rscript, for the report's "recorded by" line
current_script <- function(){
  file.v <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file.v) == 0) NA_character_ else basename(file.v[1])
}

#' Stack named tables into one, with the names as a column
#'
#' @param tables.l Named list of data frames with the same columns.
#' @param id Name of the column holding the list names.
#' @param strip Regular expression removed from the names (e.g. `"_tests$"`).
#' @return Data frame, or `NULL` when `tables.l` has no rows.
#' @export
stack_tables <- function(tables.l, id = "Dataset", strip = NULL){
  tables.l <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, tables.l)
  if (length(tables.l) == 0) return(NULL)
  names.v <- if (is.null(strip)) names(tables.l) else sub(strip, "", names(tables.l))
  stacked.df <- do.call(rbind, Map(function(x.df, name.s){
    cbind(stats::setNames(data.frame(name.s, stringsAsFactors = FALSE), id), x.df)
  }, tables.l, names.v))
  rownames(stacked.df) <- NULL
  stacked.df
}

#' PERMANOVA and PERMDISP results in one table for the summary report
#'
#' @param permanova.df Stacked [run_permanova()] results.
#' @param permdisp.df Stacked `overall` tables from [run_permdisp()].
#' @param permdisp_term The PERMANOVA term the PERMDISP tests belong to; `NULL` for every term.
#' @return Data frame with one row per dataset and term.
#' @export
summary_permanova <- function(permanova.df, permdisp.df = NULL, permdisp_term = NULL){
  if (is.null(permanova.df) || nrow(permanova.df) == 0) return(NULL)
  tests.df <- permanova.df[!permanova.df$Term %in% c("Residual", "Total"), , drop = FALSE]
  result.df <- data.frame(Dataset = tests.df$Label, Term = tests.df$Term, N = tests.df$N, R2 = tests.df$R2,
                          F = tests.df$F, P_value = tests.df$P_value, stringsAsFactors = FALSE)
  if (!is.null(permdisp.df) && nrow(permdisp.df) > 0){
    groups.df <- permdisp.df[permdisp.df$Term == "Groups", , drop = FALSE]
    result.df$PERMDISP_P_value <- groups.df$P_value[match(result.df$Dataset, groups.df$Label)]
    if (!is.null(permdisp_term)) result.df$PERMDISP_P_value[result.df$Term != permdisp_term] <- NA
  }
  result.df
}

#' Differential abundance consensus in two tables for the summary report
#'
#' @param consensus.l Named list of [da_consensus()] results, one per dataset.
#' @param top_n Features per dataset in the top features table (most methods first).
#' @return List with `Calls` (features called per dataset, by two or more methods and with
#'   conflicting directions) and `Top features`.
#' @export
summary_da <- function(consensus.l, top_n = 3){
  consensus.l <- Filter(is.data.frame, consensus.l)
  if (length(consensus.l) == 0) return(list())
  # Datasets with no significant features stay in the counts, with zeros
  calls.df <- do.call(rbind, lapply(names(consensus.l), function(name.s){
    x.df <- consensus.l[[name.s]]
    data.frame(Dataset = name.s, Features_called = length(unique(x.df$Feature_ID)),
               By_2_or_more_methods = length(unique(x.df$Feature_ID[x.df$N_methods >= 2])),
               Conflicting = length(unique(x.df$Feature_ID[x.df$Conflicting_direction %in% TRUE])))
  }))
  called.l <- Filter(function(x) nrow(x) > 0, consensus.l)
  top.df <- stack_tables(lapply(called.l, function(x.df){
    x.df <- x.df[order(-x.df$N_methods), c("Label", "Enriched_in", "N_methods", "Conflicting_direction"), drop = FALSE]
    utils::head(x.df, top_n)
  }))
  list(Calls = calls.df, `Top features` = top.df)
}

# Sections in report order, with the script that records each; built-in sections come from the processed project
report_sections <- function(mode = c("metagenome", "isolate")){
  mode <- match.arg(mode)
  if (mode == "metagenome"){
    data.frame(Section = c("overview", "reads", "mags", "diversity", "ordination", "differential_abundance", "strains",
                           "comparison"),
               Script = c(NA, NA, "mag_summary.R", "diversity.R", "ordination.R", "differential_abundance.R",
                          "strains.R", "comparison.R"), stringsAsFactors = FALSE)
  } else {
    data.frame(Section = c("overview", "reads", "genomes"), Script = c(NA, NA, "genome_summary.R"),
               stringsAsFactors = FALSE)
  }
}

#' Write the summary report
#'
#' A single self-contained HTML page with the headline numbers of the project (samples, reads,
#' bins or genomes, from the processed project) and the results recorded by the analysis scripts
#' with [record_summary()], with links to the full tables and figures.
#' Run it after the analysis scripts (`code/summary_report.R`); it can be rerun at any time.
#' Needs the rmarkdown package and pandoc (which comes with RStudio).
#'
#' @param project.l A `gm_project`.
#' @param file Output path; defaults to `Result_other/Summary_report.html`.
#' @param max_rows Rows shown per table; longer tables link to the full file.
#' @return The report path, invisibly.
#' @export
write_summary_report <- function(project.l, file = NULL, max_rows = 30){
  require_pkg("rmarkdown", "to write the summary report")
  if (!rmarkdown::pandoc_available("2.0")){
    cli::cli_abort(c("The summary report needs pandoc 2.0 or later",
                     "i" = "Install pandoc (it comes with RStudio) or set {.envvar RSTUDIO_PANDOC}"))
  }
  config.l <- project.l$config
  file <- file %||% output_path(config.l, "other", "Summary_report.html")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  file <- file.path(normalizePath(dirname(file)), basename(file))
  sections.l <- collect_report_sections(project.l)
  body.s <- report_body_html(sections.l, project.l, dirname(file), max_rows)

  work_dir.s <- tempfile("gmaio_report_")
  dir.create(work_dir.s)
  on.exit(unlink(work_dir.s, recursive = TRUE), add = TRUE)
  template.s <- file.path(work_dir.s, "summary_report.Rmd")
  file.copy(system.file("report", "summary_report.Rmd", package = "gmaio", mustWork = TRUE), template.s)
  title.s <- paste(config.l$project_name %||% "gmaio project", "summary report")
  rmarkdown::render(template.s, output_file = basename(file), output_dir = dirname(file),
                    intermediates_dir = work_dir.s, knit_root_dir = work_dir.s,
                    params = list(title = title.s, body = body.s), envir = new.env(parent = globalenv()), quiet = TRUE)
  cli::cli_inform(c("v" = "Wrote the summary report to {.path {file}}"))
  invisible(file)
}

# Built-in sections from the processed project, merged with the recorded ones, in report order
collect_report_sections <- function(project.l){
  config.l <- project.l$config
  builtin.l <- Filter(Negate(is.null), list(overview = report_overview(project.l), reads = report_reads(project.l),
                                            mags = report_mags(project.l), genomes = report_genomes(project.l)))
  summary_dir.s <- project_path(config.l, config.l$outputs$other, "summary")
  recorded.l <- lapply(list.files(summary_dir.s, pattern = "\\.rds$", full.names = TRUE), readRDS)
  names(recorded.l) <- vapply(recorded.l, `[[`, character(1), "section")
  processed.s <- project_path(config.l, config.l$outputs$other, "processed.rds")
  processed_time.x <- if (file.exists(processed.s)) file.mtime(processed.s)
  for (name.s in names(recorded.l)){
    record.l <- recorded.l[[name.s]]
    # Results recorded before the latest main.R run may not match the processed project
    if (!is.null(processed_time.x) && record.l$created < processed_time.x){
      record.l$notes <- c(record.l$notes, paste("Recorded before the latest main.R run; rerun",
                                                if_missing(record.l$script, "its script"), "to update."))
    }
    if (name.s %in% names(builtin.l)){
      merged.l <- builtin.l[[name.s]]
      for (field.s in c("tables", "values", "figures", "files", "notes")){
        merged.l[[field.s]] <- c(merged.l[[field.s]], record.l[[field.s]])
      }
      merged.l[c("script", "created")] <- record.l[c("script", "created")]
      builtin.l[[name.s]] <- merged.l
    } else {
      builtin.l[[name.s]] <- record.l
    }
  }
  order.v <- report_sections(config.l$mode)$Section
  builtin.l[c(intersect(order.v, names(builtin.l)), sort(setdiff(names(builtin.l), order.v)))]
}

if_missing <- function(x, y) if (is.null(x) || is.na(x)) y else x

report_overview <- function(project.l){
  metadata.df <- project.l$metadata
  analysed.df <- analysis_metadata(project.l)
  excluded.v <- metadata.df$Sample_ID[metadata.df$Excluded %in% TRUE]
  values.l <- list(Samples = nrow(analysed.df))
  if (length(excluded.v) > 0) values.l$Excluded <- length(excluded.v)
  bins.df <- project.l$tables$bin_summary
  if (!is.null(bins.df)){
    values.l$Bins <- nrow(bins.df)
    values.l$`High-quality bins` <- sum(bins.df$High_quality %in% TRUE)
  }
  if (!is.null(project.l$tables$reference_genomes)) values.l$`Reference genomes` <- nrow(project.l$tables$reference_genomes)
  if (!is.null(project.l$comparison)) values.l$`Comparison samples` <- nrow(project.l$comparison$metadata)
  tables.l <- list()
  for (group.s in project.l$config$analysis$group_variables){
    if (!group.s %in% names(analysed.df)) next
    counts.v <- table(analysed.df[[group.s]], useNA = "ifany")
    tables.l[[paste("Samples by", variable_label(group.s))]] <-
      stats::setNames(data.frame(names(counts.v), as.integer(counts.v), stringsAsFactors = FALSE),
                      c(variable_label(group.s), "Samples"))
  }
  notes.v <- if (length(excluded.v) > 0) paste("Excluded from analyses:", paste(excluded.v, collapse = ", "))
  list(section = "overview", title = "Overview", values = values.l, tables = tables.l, notes = notes.v,
       files = project_path(project.l$config, project.l$config$outputs$tables, "Metadata_and_read_stats.xlsx"))
}

report_reads <- function(project.l){
  stats.df <- project.l$read_stats
  if (is.null(stats.df) || !"Raw_count" %in% names(stats.df)) return(NULL)
  analysed.df <- analysis_metadata(project.l)
  stats.df <- stats.df[stats.df$Sample_ID %in% analysed.df$Sample_ID, , drop = FALSE]
  if (nrow(stats.df) == 0) return(NULL)
  # The last QC step before mapping, e.g. host removal after trimming
  count_columns.v <- grep("_count$", names(stats.df), value = TRUE)
  qc_columns.v <- setdiff(count_columns.v[!grepl("^(Reads|Bases)_mapped_", count_columns.v)], "Raw_count")
  sets.v <- mapping_set_labels(names(stats.df))
  group.s <- primary_group(project.l)
  groups.l <- list(All = stats.df$Sample_ID)
  if (!is.null(group.s)){
    levels.v <- as.character(unique(stats::na.omit(analysed.df[[group.s]])))
    if (is.factor(analysed.df[[group.s]])) levels.v <- intersect(levels(analysed.df[[group.s]]), levels.v)
    for (level.s in levels.v) groups.l[[level.s]] <- analysed.df$Sample_ID[analysed.df[[group.s]] %in% level.s]
  }
  rows.l <- list(Samples = function(x.df) nrow(x.df),
                 `Raw reads (median)` = function(x.df) format_count(stats::median(x.df$Raw_count, na.rm = TRUE)))
  if (length(qc_columns.v) > 0){
    qc.s <- utils::tail(qc_columns.v, 1)
    rows.l$`Kept after QC (median %)` <- function(x.df) {
      round(stats::median(100 * x.df[[qc.s]] / x.df$Raw_count, na.rm = TRUE), 1)
    }
  }
  median_percent.f <- function(column.s){
    force(column.s)
    function(x.df) round(stats::median(x.df[[column.s]], na.rm = TRUE), 1)
  }
  for (set.s in sets.v){
    row.s <- paste0("Mapped to ", gsub("_", " ", set.s), " (median %)")
    rows.l[[row.s]] <- median_percent.f(paste0("Reads_mapped_", set.s, "_percent"))
  }
  table.df <- data.frame(Measure = names(rows.l), stringsAsFactors = FALSE)
  for (name.s in names(groups.l)){
    x.df <- stats.df[stats.df$Sample_ID %in% groups.l[[name.s]], , drop = FALSE]
    table.df[[name.s]] <- vapply(rows.l, function(f) as.character(f(x.df)), character(1))
  }
  notes.v <- if (length(sets.v) > 0) "Mapping percentages are of raw reads, as in the pipeline read statistics."
  list(section = "reads", title = "Samples and reads",
       tables = stats::setNames(list(table.df), if (is.null(group.s)) "Reads" else paste("Reads by", variable_label(group.s))),
       notes = notes.v, files = project_path(project.l$config, project.l$config$outputs$tables, "Metadata_and_read_stats.xlsx"))
}

report_mags <- function(project.l){
  bins.df <- project.l$tables$bin_summary
  if (is.null(bins.df)) return(NULL)
  tiers.v <- intersect(c("High", "Medium", "Low"), unique(bins.df$MIMAG_tier))
  hq.v <- bins.df$High_quality %in% TRUE
  tiers.df <- data.frame(MIMAG_tier = tiers.v,
                         Bins = vapply(tiers.v, function(t) sum(bins.df$MIMAG_tier %in% t), integer(1)),
                         High_quality = vapply(tiers.v, function(t) sum(bins.df$MIMAG_tier %in% t & hq.v), integer(1)),
                         stringsAsFactors = FALSE)
  tables.l <- list(`Bins by MIMAG tier` = tiers.df)
  hq.df <- bins.df[bins.df$High_quality %in% TRUE, , drop = FALSE]
  if (nrow(hq.df) > 0 && "Phylum" %in% names(hq.df)){
    phyla.v <- clean_taxon_name(hq.df$Phylum)
    phyla.v <- sort(table(ifelse(is.na(phyla.v) | phyla.v == "", "Unclassified", phyla.v)), decreasing = TRUE)
    tables.l$`High-quality bins by phylum` <- data.frame(Phylum = names(phyla.v), Bins = as.integer(phyla.v),
                                                          stringsAsFactors = FALSE)
  }
  config.l <- project.l$config
  notes.v <- paste0("High quality: completeness - ", config.l$mags$quality_weight, " x contamination >= ",
                    config.l$mags$quality_threshold, " (", config.l$mags$quality_source, ").")
  list(section = "mags", title = "Bins and MAGs", tables = tables.l, notes = notes.v,
       files = project_path(config.l, config.l$outputs$tables, "MAG_summary.xlsx"))
}

report_genomes <- function(project.l){
  genomes.df <- project.l$tables$genome_summary
  if (is.null(genomes.df)) return(NULL)
  columns.v <- intersect(c("Sample_label", "Species", "Completeness", "Contamination", "Genome_size", "Contigs",
                           "Contig_N50", "ST", "Mean_coverage", "AMR_elements"), names(genomes.df))
  shown.df <- genomes.df[, columns.v, drop = FALSE]
  if ("Species" %in% columns.v) shown.df$Species <- clean_taxon_name(shown.df$Species)
  tables.l <- list(Genomes = shown.df)
  if ("Species" %in% names(genomes.df)){
    species.v <- sort(table(clean_taxon_name(genomes.df$Species)), decreasing = TRUE)
    tables.l <- c(list(`Genomes by species` = data.frame(Species = names(species.v), Genomes = as.integer(species.v),
                                                          stringsAsFactors = FALSE)), tables.l)
  }
  list(section = "genomes", title = "Genomes", tables = tables.l,
       files = project_path(project.l$config, project.l$config$outputs$tables, "Genome_summary.xlsx"))
}

format_count <- function(x){
  if (is.na(x)) return("")
  if (x >= 1e6) return(paste0(format(round(x / 1e6, 1), nsmall = 1), "M"))
  if (x >= 1e3) return(paste0(format(round(x / 1e3, 1), nsmall = 1), "k"))
  format(x)
}

# Values for display: p values with two significant digits, other numbers with three. Test_label
# (the p value as drawn on figures) repeats P_value, so it is left out.
format_report_table <- function(table.df){
  table.df <- as.data.frame(table.df, stringsAsFactors = FALSE)
  table.df <- table.df[setdiff(names(table.df), "Test_label")]
  for (column.s in names(table.df)){
    x <- table.df[[column.s]]
    if (is.logical(x)){
      x <- ifelse(is.na(x), "", ifelse(x, "yes", "no"))
    } else if (is.numeric(x)){
      p_value.b <- grepl("(^|_)(P|p|q)(_value|_adjusted|_adj)?$|P_value|P_adjusted|q_value", column.s)
      x <- vapply(x, function(v){
        if (is.na(v)) return("")
        if (p_value.b) return(if (v < 0.001) "< 0.001" else formatC(signif(v, 2), format = "fg", digits = 2))
        if (v == round(v) && abs(v) < 1e15) return(formatC(v, format = "d", big.mark = ","))
        formatC(signif(v, 3), format = "fg", digits = 3, big.mark = ",")
      }, character(1))
    } else {
      x <- ifelse(is.na(x), "", as.character(x))
    }
    table.df[[column.s]] <- x
  }
  names(table.df) <- gsub("_", " ", names(table.df))
  table.df
}

# Relative link from the report directory, so the links work when the project folder is moved or shared
relative_path <- function(path, from_dir){
  resolve.f <- function(p){
    p <- path.expand(p)
    if (file.exists(p)) return(normalizePath(p))
    parent.s <- dirname(p)
    if (dir.exists(parent.s)) file.path(normalizePath(parent.s), basename(p)) else p
  }
  target.v <- strsplit(resolve.f(path), "/", fixed = TRUE)[[1]]
  from.v <- strsplit(normalizePath(from_dir), "/", fixed = TRUE)[[1]]
  shared.n <- 0
  while (shared.n < min(length(target.v), length(from.v)) && target.v[shared.n + 1] == from.v[shared.n + 1]){
    shared.n <- shared.n + 1
  }
  if (shared.n <= 1) return(paste0("file://", resolve.f(path)))
  paste(c(rep("..", length(from.v) - shared.n), target.v[-seq_len(shared.n)]), collapse = "/")
}

# Width over height of a png, from its IHDR header
png_aspect_ratio <- function(png.r){
  size.f <- function(bytes.r) sum(as.integer(bytes.r) * 256^(3:0))
  size.f(png.r[17:20]) / size.f(png.r[21:24])
}

html_escape <- function(x){
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

# The report body as pandoc markdown with raw HTML blocks: a level 2 heading per section (for the table of
# contents) followed by its cards, tables, figures, files and notes
report_body_html <- function(sections.l, project.l, report_dir, max_rows){
  config.l <- project.l$config
  lines.v <- c(sprintf("<p class=\"gm-meta\">Generated %s with gmaio %s \u00b7 %s project</p>",
                       format(Sys.time(), "%Y-%m-%d %H:%M"), utils::packageVersion("gmaio"), config.l$mode), "")
  for (section.l in sections.l){
    lines.v <- c(lines.v, paste("##", section.l$title), "")
    recorded.s <- if (!is.null(section.l$created)){
      paste0("From ", if_missing(section.l$script, "an analysis script"), ", ", format(section.l$created, "%Y-%m-%d %H:%M"))
    } else "From the processed project (main.R)"
    lines.v <- c(lines.v, sprintf("<p class=\"gm-meta\">%s</p>", html_escape(recorded.s)), "")
    if (length(section.l$values) > 0){
      cards.v <- vapply(names(section.l$values), function(name.s){
        value <- section.l$values[[name.s]]
        value.s <- if (is.numeric(value)) format_report_table(data.frame(v = value))$v else as.character(value)
        sprintf("<div class=\"gm-card\"><div class=\"gm-card-label\">%s</div><div class=\"gm-card-value\">%s</div></div>",
                html_escape(name.s), html_escape(value.s))
      }, character(1))
      lines.v <- c(lines.v, "<div class=\"gm-cards\">", cards.v, "</div>", "")
    }
    for (name.s in names(section.l$tables)){
      table.df <- section.l$tables[[name.s]]
      shown.df <- format_report_table(utils::head(table.df, max_rows))
      table.s <- knitr::kable(shown.df, format = "html", align = "c", row.names = FALSE,
                              table.attr = "class=\"table table-condensed gm-table\"")
      lines.v <- c(lines.v, sprintf("<p class=\"gm-table-title\">%s</p>", html_escape(name.s)),
                   "<div class=\"gm-table-wrap\">", gsub("\n", "", table.s), "</div>")
      if (nrow(table.df) > max_rows){
        lines.v <- c(lines.v, sprintf("<p class=\"gm-meta\">First %d of %d rows; the full table is in the files below.</p>",
                                      max_rows, nrow(table.df)))
      }
      lines.v <- c(lines.v, "")
    }
    figures.l <- Filter(function(f) file.exists(f$thumbnail), section.l$figures)
    if (length(figures.l) > 0){
      figures.v <- vapply(figures.l, function(f){
        png.r <- readBin(f$thumbnail, "raw", file.info(f$thumbnail)$size)
        image.s <- jsonlite::base64_enc(png.r)
        link.s <- html_escape(relative_path(f$path, report_dir))
        # Figures more than twice as wide as tall would be too small in half the page width
        class.s <- if (png_aspect_ratio(png.r) > 2) "gm-figure gm-figure-wide" else "gm-figure"
        sprintf(paste0("<figure class=\"%s\"><a href=\"%s\"><img src=\"data:image/png;base64,%s\" alt=\"%s\"></a>",
                       "<figcaption>%s \u00b7 <a href=\"%s\">full figure</a></figcaption></figure>"),
                class.s, link.s, gsub("\n", "", image.s), html_escape(f$caption), html_escape(f$caption), link.s)
      }, character(1))
      lines.v <- c(lines.v, "<div class=\"gm-figures\">", figures.v, "</div>", "")
    }
    files.v <- as.character(section.l$files)
    files.v <- unique(files.v[file.exists(files.v)])
    if (length(files.v) > 0){
      links.v <- sprintf("<a href=\"%s\">%s</a>", html_escape(vapply(files.v, relative_path, character(1), report_dir)),
                         html_escape(basename(files.v)))
      lines.v <- c(lines.v, sprintf("<p class=\"gm-files\">Full results: %s</p>", paste(links.v, collapse = " \u00b7 ")), "")
    }
    for (note.s in section.l$notes) lines.v <- c(lines.v, sprintf("<p class=\"gm-note\">%s</p>", html_escape(note.s)), "")
  }
  raw_html_blocks(c(lines.v, report_outputs_html(sections.l, config.l, report_dir)))
}

# Pandoc reads markdown inside plain HTML blocks (a table cell holding "-", such as an unknown MLST type, becomes a
# list and breaks the page), so everything but the section headings goes in raw HTML blocks that pandoc passes on as is
raw_html_blocks <- function(lines.v){
  lines.v <- lines.v[nzchar(lines.v)]
  heading.v <- startsWith(lines.v, "## ")
  block.v <- cumsum(heading.v | c(TRUE, heading.v[-length(heading.v)]))
  unlist(lapply(split(lines.v, block.v), function(x.v){
    if (startsWith(x.v[1], "## ")) c(x.v, "") else c("```{=html}", x.v, "```", "")
  }), use.names = FALSE)
}

report_outputs_html <- function(sections.l, config.l, report_dir){
  lines.v <- c("## Outputs", "")
  indexes.v <- c(Tables = project_path(config.l, config.l$outputs$tables, "README.md"),
                 Figures = project_path(config.l, config.l$outputs$figures, "README.md"))
  indexes.v <- indexes.v[file.exists(indexes.v)]
  if (length(indexes.v) > 0){
    links.v <- sprintf("<a href=\"%s\">%s index</a>", html_escape(vapply(indexes.v, relative_path, character(1), report_dir)),
                       names(indexes.v))
    lines.v <- c(lines.v, sprintf("<p class=\"gm-files\">%s</p>", paste(links.v, collapse = " \u00b7 ")), "")
  }
  # A script counts as not run when its section has no record (built-in sections exist without one)
  sections.df <- report_sections(config.l$mode)
  sections.df <- sections.df[!is.na(sections.df$Script), , drop = FALSE]
  recorded.v <- names(Filter(function(x) !is.null(x$created), sections.l))
  missing.v <- sections.df$Script[!sections.df$Section %in% recorded.v]
  if (length(missing.v) > 0){
    lines.v <- c(lines.v, sprintf(paste("<p class=\"gm-note\">No results recorded from %s: not run yet, skipped",
                                        "for lack of inputs, or run before <code>outputs: report: true</code> was set.</p>"),
                                  paste(missing.v, collapse = ", ")), "")
  }
  lines.v
}

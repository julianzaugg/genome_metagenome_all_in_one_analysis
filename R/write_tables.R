#' Write data frames to a styled Excel workbook
#'
#' One sheet per list element: bold header, frozen header row (and first column when
#' `freeze_first_column`), auto column widths and filters. Empty data frames are
#' written as a header-only sheet. An empty list writes nothing.
#'
#' @param tables.l Named list of data frames (names become sheet names, max 31 characters).
#' @param path Output `.xlsx` path; the directory is created.
#' @param freeze_first_column Freeze the first column as well as the header.
#' @param number_format Optional Excel number format for numeric columns, e.g. `"0.0000"`.
#' @return The path, or `NULL` when there was nothing to write, invisibly.
#' @export
write_xlsx_tables <- function(tables.l, path, freeze_first_column = TRUE, number_format = NULL){
  if (is.data.frame(tables.l)) tables.l <- list(Sheet1 = tables.l)
  if (length(tables.l) == 0){
    cli::cli_inform(c("i" = "No tables to write to {.path {path}}"))
    return(invisible(NULL))
  }
  if (is.null(names(tables.l)) || any(names(tables.l) == "")) cli::cli_abort("{.arg tables.l} must be a named list")
  sheet_names.v <- substr(gsub("[][*?:/\\\\]", "_", names(tables.l)), 1, 31)
  if (anyDuplicated(sheet_names.v)) cli::cli_abort("Sheet names are not unique after truncation: {.val {sheet_names.v}}")

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  workbook <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "Bottom",
                                        halign = "center", valign = "center", wrapText = TRUE)
  number_style <- if (!is.null(number_format)) openxlsx::createStyle(numFmt = number_format) else NULL

  for (i in seq_along(tables.l)){
    sheet.s <- sheet_names.v[i]
    table.df <- as.data.frame(tables.l[[i]], check.names = FALSE, stringsAsFactors = FALSE)
    factor_columns.v <- vapply(table.df, is.factor, logical(1))
    table.df[factor_columns.v] <- lapply(table.df[factor_columns.v], as.character)

    openxlsx::addWorksheet(workbook, sheet.s)
    openxlsx::writeData(workbook, sheet.s, table.df, headerStyle = header_style,
                        withFilter = ncol(table.df) > 0 && nrow(table.df) > 0)
    openxlsx::freezePane(workbook, sheet.s, firstRow = TRUE, firstCol = freeze_first_column && ncol(table.df) > 1)
    if (ncol(table.df) > 0){
      widths.v <- vapply(seq_along(table.df), function(j){
        values.v <- c(names(table.df)[j], utils::head(format(table.df[[j]]), 500))
        min(max(nchar(values.v, keepNA = FALSE), 8) + 2, 60)
      }, numeric(1))
      openxlsx::setColWidths(workbook, sheet.s, cols = seq_along(table.df), widths = widths.v)
    }
    numeric_columns.v <- which(vapply(table.df, is.double, logical(1)))
    if (!is.null(number_style) && length(numeric_columns.v) > 0 && nrow(table.df) > 0){
      openxlsx::addStyle(workbook, sheet.s, number_style, rows = seq_len(nrow(table.df)) + 1,
                         cols = numeric_columns.v, gridExpand = TRUE)
    }
  }
  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
  invisible(path)
}

#' Write one profile per sheet
#'
#' @param profiles.l Named list of `gm_profile` objects.
#' @param path Output `.xlsx` path.
#' @param metadata.df Optional metadata; when given, sample columns use `Sample_label`.
#' @param annotation_columns Feature columns written before the samples.
#' @return The path, invisibly.
#' @export
write_profiles_xlsx <- function(profiles.l, path, metadata.df = NULL,
                                annotation_columns = c("Feature_ID", "Label")){
  sample_labels.v <- if (!is.null(metadata.df)) stats::setNames(metadata.df$Sample_label, metadata.df$Sample_ID) else NULL
  tables.l <- lapply(profiles.l, profile_to_wide_df, annotation_columns = annotation_columns,
                     sample_labels = sample_labels.v)
  write_xlsx_tables(tables.l, path, number_format = "0.0000")
}

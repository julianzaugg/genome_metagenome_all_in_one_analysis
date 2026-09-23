#' Convert a data frame to a matrix using one column as row names
#'
#' @param input.df Data frame.
#' @param rownames_column Column (name or index) holding the row names.
#' @return Numeric matrix.
#' @export
df2m <- function(input.df, rownames_column = 1){
  input.df <- as.data.frame(input.df, check.names = FALSE)
  if (is.numeric(rownames_column)) rownames_column <- names(input.df)[rownames_column]
  row_names.v <- as.character(input.df[[rownames_column]])
  output.m <- as.matrix(input.df[, setdiff(names(input.df), rownames_column), drop = FALSE])
  rownames(output.m) <- row_names.v
  output.m
}

#' Convert a matrix to a data frame, moving row names into the first column
#'
#' @param input.m Matrix.
#' @param column_name Name of the new first column.
#' @return Data frame.
#' @export
m2df <- function(input.m, column_name = "Feature"){
  output.df <- data.frame(rownames(input.m), as.data.frame(input.m, check.names = FALSE),
                          check.names = FALSE, stringsAsFactors = FALSE)
  names(output.df)[1] <- column_name
  rownames(output.df) <- NULL
  output.df
}

read_tsv_fast <- function(path, na.strings = c("", "NA", "N/A"), ...){
  data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE, na.strings = na.strings, ...)
}

require_pkg <- function(pkg, reason = NULL){
  rlang::check_installed(pkg, reason = reason)
}

require_columns <- function(input.df, columns.v, what = "table"){
  missing.v <- setdiff(columns.v, names(input.df))
  if (length(missing.v) > 0){
    # Column names are case sensitive; point out near misses such as "label" for "Label"
    near.v <- names(input.df)[tolower(names(input.df)) %in% tolower(missing.v)]
    hint.v <- if (length(near.v) > 0) c("i" = "Column names are case sensitive; did you mean {.val {near.v}}?") else
      c("i" = "Available: {.val {names(input.df)}}")
    cli::cli_abort(c("{what} is missing column{?s} {.val {missing.v}}", hint.v))
  }
  invisible(TRUE)
}

as_yes_no_logical <- function(x){
  if (is.logical(x)) return(x)
  x.v <- tolower(trimws(as.character(x)))
  out.v <- rep(NA, length(x.v))
  out.v[x.v %in% c("yes", "y", "true", "t", "1", "exclude", "excluded")] <- TRUE
  out.v[x.v %in% c("no", "n", "false", "f", "0", "")] <- FALSE
  out.v
}

safe_file_name <- function(x){
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

significance_stars <- function(p.v){
  as.character(cut(p.v, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns")))
}

#' Run an analysis step, turning an error into a warning
#'
#' Used in template loops so one failing dataset is reported without stopping the others.
#'
#' @param expr Expression to evaluate.
#' @param label Label for the warning.
#' @return The value of `expr`, or `NULL` on error.
#' @export
try_step <- function(expr, label){
  tryCatch(expr, error = function(e){
    cli::cli_warn("{label} skipped: {conditionMessage(e)}")
    NULL
  })
}

# Join metadata onto a per-sample table by Sample_ID. Metadata columns that share a name with
# a column of x.df (e.g. a metadata "Label" column next to the feature Label) are left out,
# so the data being plotted keeps its column names instead of becoming .x/.y pairs.
join_metadata <- function(x.df, metadata.df, type = c("left", "inner")){
  type <- match.arg(type)
  metadata.df <- metadata.df[, c("Sample_ID", setdiff(names(metadata.df), names(x.df))), drop = FALSE]
  if (type == "left") dplyr::left_join(x.df, metadata.df, by = "Sample_ID") else
    dplyr::inner_join(x.df, metadata.df, by = "Sample_ID")
}

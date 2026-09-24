# Which gmaio build is running, so load_processed() can tell when the R session or
# processed.rds predates the installed package.
gm_build <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname){
  gm_build$loaded <- installed_build(pkgname, libname)
  gm_build$library <- libname
}

# Version plus the git commit (install_github) or build/install time (local installs); NA when unknown
installed_build <- function(pkgname = "gmaio", lib.loc = NULL){
  description.l <- tryCatch(suppressWarnings(utils::packageDescription(pkgname, lib.loc = lib.loc)),
                            error = function(e) NA)
  if (!is.list(description.l)) return(NA_character_)
  # Built is "R version; platform; install time; OS"
  built.s <- if (!is.null(description.l$Built)) strsplit(description.l$Built, "; ", fixed = TRUE)[[1]][3] else NULL
  stamp.s <- if (!is.null(description.l$RemoteSha)) substr(description.l$RemoteSha, 1, 7) else
    sub(";.*", "", description.l$Packaged %||% built.s)
  trimws(paste(description.l$Version, stamp.s %||% ""))
}

# devtools::load_all() has no stable installed build to compare against
is_dev_session <- function(){
  exists(".__DEVTOOLS__", envir = asNamespace("gmaio"), inherits = FALSE)
}

# Warn when gmaio was reinstalled after this session loaded it, since the session keeps running
# the old code until R restarts, and when processed.rds was made by a different build.
check_build <- function(project.l){
  loaded.s <- gm_build$loaded
  if (is_dev_session() || is.null(loaded.s) || is.na(loaded.s)) return(invisible(NULL))
  on_disk.s <- installed_build("gmaio", gm_build$library)
  if (!is.na(on_disk.s) && on_disk.s != loaded.s){
    cli::cli_warn(c("gmaio {.val {on_disk.s}} is installed, but this R session is still running {.val {loaded.s}}",
                    "i" = "Restart R (RStudio: Session > Restart R) so scripts use the installed version"))
  }
  processed.s <- project.l$gmaio_build
  if (!identical(processed.s, loaded.s)){
    made.s <- if (is.null(processed.s)) "an older gmaio" else "gmaio {.val {processed.s}}"
    cli::cli_warn(c(paste0("{.file processed.rds} was made with ", made.s, "; this session runs {.val {loaded.s}}"),
                    "i" = "Rerun {.file code/main.R} so the analyses see the current processing"))
  }
  invisible(NULL)
}

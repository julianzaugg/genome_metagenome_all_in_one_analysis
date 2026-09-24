palette_families <- function(){
  kelly.v <- c("#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#008856", "#E68FAC",
               "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600",
               "#654522", "#E25822", "#2B3D26", "#00A6A6")
  list(
    okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
    tableau10 = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1",
                  "#FF9DA7", "#9C755F", "#8CD17D"),
    distinct20 = kelly.v,
    distinct30 = c(kelly.v, "#1F78B4", "#FB9A99", "#6A3D9A", "#B2DF8A", "#FF7F00", "#CAB2D6", "#B15928",
                   "#33A02C", "#FDBF6F", "#E31A1C")
  )
}

#' Fixed colours for special categories
#'
#' @return Named character vector for `"Other"`, `"Unassigned"`, missing values and reference genomes.
#' @export
special_colours <- function(){
  c(Other = "#D9D9D9", Unassigned = "#9E9E9E", Missing = "#EFEFEF", Reference = "#4D4D4D")
}

#' Generate a qualitative palette of n colours
#'
#' Uses Okabe-Ito (up to 7), Tableau 10, Kelly's colours of maximum contrast (up to 20),
#' an extended 30, and beyond 30 those 30 followed by an HCL wheel with alternating lightness.
#' None contain greys, which are reserved for `Other` and `Unassigned`.
#'
#' @param n Number of colours.
#' @return Character vector of hex colours.
#' @export
qualitative_palette <- function(n){
  families.l <- palette_families()
  if (n <= 0) return(character())
  if (n <= 7) return(families.l$okabe_ito[seq_len(n)])
  if (n <= 10) return(families.l$tableau10[seq_len(n)])
  if (n <= 20) return(families.l$distinct20[seq_len(n)])
  if (n <= 30) return(families.l$distinct30[seq_len(n)])
  extra.n <- n - 30
  hues.v <- seq(15, 375, length.out = extra.n + 1)[seq_len(extra.n)]
  c(families.l$distinct30,
    grDevices::hcl(h = hues.v, c = rep(c(70, 50), length.out = extra.n), l = rep(c(55, 75, 40), length.out = extra.n)))
}

#' Assign colours to the values of a variable
#'
#' Values keep colours from `fixed.v` (user config) first, then `existing.v` (saved
#' palette); `Other`, `Unassigned` and `Reference` get their [special_colours()], and any
#' remaining values get unused colours from [qualitative_palette()] in level order, so the
#' result is deterministic.
#'
#' @param values.v Vector (factor levels are respected, otherwise sorted unique values).
#' @param existing.v Named vector of previously assigned colours.
#' @param fixed.v Named vector of user-fixed colours.
#' @param avoid.v Colours to leave out when generating new ones (e.g. those of another
#'   variable), as long as enough other colours remain.
#' @return Named character vector of colours.
#' @export
assign_colours <- function(values.v, existing.v = NULL, fixed.v = NULL, avoid.v = NULL){
  levels.v <- if (is.factor(values.v)) levels(droplevels(values.v)) else sort(unique(as.character(stats::na.omit(values.v))))
  special.v <- special_colours()
  colours.v <- stats::setNames(rep(NA_character_, length(levels.v)), levels.v)
  for (source.v in list(fixed.v, existing.v, special.v[c("Other", "Unassigned", "Reference")])){
    if (length(source.v) == 0) next
    hit.v <- levels.v[is.na(colours.v[levels.v]) & levels.v %in% names(source.v)]
    colours.v[hit.v] <- unlist(source.v)[hit.v]
  }
  remaining.v <- names(colours.v)[is.na(colours.v)]
  if (length(remaining.v) > 0){
    candidates.v <- setdiff(toupper(qualitative_palette(length(levels.v))), toupper(colours.v))
    if (length(avoid.v) > 0){
      spread.v <- setdiff(toupper(qualitative_palette(length(levels.v) + length(avoid.v))), toupper(c(colours.v, avoid.v)))
      if (length(spread.v) >= length(remaining.v)) candidates.v <- spread.v
    }
    if (length(candidates.v) < length(remaining.v)){
      candidates.v <- c(candidates.v, qualitative_palette(length(remaining.v) + 30))
      candidates.v <- setdiff(unique(toupper(candidates.v)), toupper(colours.v))
    }
    colours.v[remaining.v] <- candidates.v[seq_along(remaining.v)]
  }
  colours.v
}

#' Assign taxon colours
#'
#' Phyla get qualitative colours in order of abundance. With `scheme = "distinct"`,
#' taxa below phylum get their own qualitative colours in order of overall abundance, so
#' the most abundant taxa get the most distinct colours. With `scheme = "hierarchical"`,
#' they get shades of their phylum colour (clearer grouping, but similar shades when a
#' phylum has many taxa). Unassigned taxa are grey.
#'
#' @param taxa.df Data frame with `Label`, `Phylum` and `Abundance` (used for ordering).
#' @param existing.v Named vector of saved colours (by label).
#' @param fixed.v Named vector of user-fixed colours (by label).
#' @param scheme `"distinct"` or `"hierarchical"`.
#' @return Named character vector of colours by label.
#' @export
assign_taxa_colours <- function(taxa.df, existing.v = NULL, fixed.v = NULL, scheme = c("distinct", "hierarchical")){
  scheme <- match.arg(scheme)
  require_columns(taxa.df, c("Label", "Phylum", "Abundance"), "Taxa table")
  taxa.df <- taxa.df[!is.na(taxa.df$Label), , drop = FALSE]
  taxa.df$Phylum[is.na(taxa.df$Phylum)] <- "Unassigned"
  taxa.df <- stats::aggregate(taxa.df["Abundance"], by = taxa.df[c("Label", "Phylum")], FUN = max)
  phylum_abundance.v <- sort(tapply(taxa.df$Abundance, taxa.df$Phylum, sum), decreasing = TRUE)
  phyla.v <- setdiff(names(phylum_abundance.v), "Unassigned")

  phylum_colours.v <- assign_colours(factor(phyla.v, levels = phyla.v), existing.v, fixed.v)
  colours.v <- c(phylum_colours.v, special_colours()[c("Other", "Unassigned")])

  lower.df <- taxa.df[!taxa.df$Label %in% names(colours.v), , drop = FALSE]
  unassigned.v <- is_unassigned_label(lower.df$Label)
  grey.v <- rep(special_colours()[["Unassigned"]], sum(unassigned.v))
  colours.v <- c(colours.v, stats::setNames(grey.v, lower.df$Label[unassigned.v]))
  lower.df <- lower.df[!unassigned.v, , drop = FALSE]
  if (scheme == "distinct"){
    lower.df <- lower.df[order(-lower.df$Abundance, lower.df$Label), , drop = FALSE]
    lower_existing.v <- unlist(existing.v)[intersect(names(unlist(existing.v)), lower.df$Label)]
    colours.v <- c(colours.v, assign_colours(factor(lower.df$Label, levels = lower.df$Label), lower_existing.v, fixed.v))
    lower.df <- lower.df[0, , drop = FALSE]
  }
  lower.df <- lower.df[order(lower.df$Phylum, -lower.df$Abundance, lower.df$Label), , drop = FALSE]
  for (phylum.s in unique(lower.df$Phylum)){
    labels.v <- lower.df$Label[lower.df$Phylum == phylum.s]
    base.s <- if (phylum.s %in% names(colours.v)) colours.v[[phylum.s]] else special_colours()[["Unassigned"]]
    colours.v <- c(colours.v, stats::setNames(shade_colours(base.s, length(labels.v)), labels.v))
  }
  for (source.v in list(existing.v, fixed.v)){
    shared.v <- intersect(names(source.v), names(colours.v))
    colours.v[shared.v] <- unlist(source.v)[shared.v]
  }
  colours.v
}

shade_colours <- function(base.s, n){
  if (n == 0) return(character())
  base.m <- farver::decode_colour(base.s, to = "hcl")
  index.v <- seq_len(n) - 1
  lightness_offset.v <- ((index.v * 0.618034) %% 1) * 2 - 1
  hue_offset.v <- ((index.v * 0.414214) %% 1) * 2 - 1
  lightness_offset.v[1] <- 0
  hue_offset.v[1] <- 0
  hcl.m <- cbind(h = (base.m[, "h"] + 18 * hue_offset.v) %% 360,
                 c = pmax(base.m[, "c"] - 10 * abs(lightness_offset.v), 15),
                 l = pmin(pmax(base.m[, "l"] + 28 * lightness_offset.v, 25), 88))
  farver::encode_colour(hcl.m, from = "hcl")
}

#' Assign point shapes to the values of a variable
#'
#' @param values.v Vector (factor levels respected).
#' @return Named integer vector of filled shapes (21 to 25, then open shapes).
#' @export
assign_shapes <- function(values.v){
  levels.v <- if (is.factor(values.v)) levels(droplevels(values.v)) else sort(unique(as.character(stats::na.omit(values.v))))
  pool.v <- c(21, 22, 24, 23, 25, 0:14)
  stats::setNames(pool.v[(seq_along(levels.v) - 1) %% length(pool.v) + 1], levels.v)
}

#' Build all project palettes
#'
#' @param metadata.df Metadata from [read_metadata()].
#' @param variables.v Metadata columns to colour.
#' @param taxa.df Optional taxa table for [assign_taxa_colours()].
#' @param existing.l Previously saved palettes (see [read_palettes()]).
#' @param fixed.l User-fixed colours from the config `colours:` section.
#' @param taxa_scheme Passed to [assign_taxa_colours()].
#' @return List with one named colour vector per variable plus `taxa`.
#' @export
build_palettes <- function(metadata.df, variables.v, taxa.df = NULL, existing.l = list(), fixed.l = list(),
                           taxa_scheme = "distinct"){
  variables.v <- unique(c("Sample_ID", variables.v))
  require_columns(metadata.df, variables.v, "Metadata")
  # Each grouping variable avoids the colours of the ones before it, so e.g. Treatment and
  # Time annotations on one figure do not share colours
  palettes.l <- list()
  used.v <- character()
  for (variable.s in variables.v){
    values.v <- metadata.df[[variable.s]]
    if (variable.s == "Sample_ID") values.v <- factor(values.v, levels = metadata.df$Sample_ID)
    palettes.l[[variable.s]] <- assign_colours(values.v, unlist(existing.l[[variable.s]]), unlist(fixed.l[[variable.s]]),
                                               avoid.v = if (variable.s != "Sample_ID") used.v)
    if (variable.s != "Sample_ID") used.v <- unique(c(used.v, palettes.l[[variable.s]]))
  }
  palettes.l$Sample_label <- stats::setNames(palettes.l$Sample_ID[metadata.df$Sample_ID], metadata.df$Sample_label)
  if (!is.null(taxa.df)){
    palettes.l$taxa <- assign_taxa_colours(taxa.df, unlist(existing.l$taxa), unlist(fixed.l$taxa), scheme = taxa_scheme)
  }
  palettes.l
}

#' Read and write the project palette file
#'
#' @param path Path to `palettes.yml`.
#' @param palettes.l Palettes from [build_palettes()].
#' @return `read_palettes()` returns a list of named character vectors (empty when the file is absent).
#' @export
read_palettes <- function(path){
  if (!file.exists(path)) return(list())
  lapply(yaml::read_yaml(path) %||% list(), unlist)
}

#' @rdname read_palettes
#' @export
write_palettes <- function(palettes.l, path){
  store.l <- palettes.l[setdiff(names(palettes.l), "Sample_label")]
  header.s <- paste("# Colour assignments written by gmaio. Edit any colour; it is kept on the next run.",
                    "# Colours set in config.yml (colours:) take precedence over this file.", sep = "\n")
  writeLines(c(header.s, yaml::as.yaml(lapply(store.l, as.list))), path)
  invisible(path)
}

#' Get a palette for a variable
#'
#' @param palettes.l Project palettes, or a processed project list containing `palettes`.
#' @param variable Variable name, or `"taxa"`.
#' @param values.v Optional values to return colours for. Taxa without a saved colour
#'   get an unused colour; other variables must already be covered.
#' @return Named character vector.
#' @export
get_palette <- function(palettes.l, variable, values.v = NULL){
  if (!is.null(palettes.l$palettes)) palettes.l <- palettes.l$palettes
  colours.v <- palettes.l[[variable]]
  if (is.null(colours.v) && is.null(values.v)) cli::cli_abort("No palette for {.val {variable}}")
  if (is.null(values.v)) return(colours.v)
  values.v <- unique(as.character(stats::na.omit(values.v)))
  missing.v <- setdiff(values.v, names(colours.v))
  if (length(missing.v) > 0){
    if (variable != "taxa"){
      cli::cli_abort(c("No {.val {variable}} colour for {.val {missing.v}}",
                       "i" = "Rerun {.file main.R} after changing metadata"))
    }
    colours.v <- c(colours.v, assign_colours(c(names(colours.v), missing.v), existing.v = colours.v)[missing.v])
  }
  colours.v[values.v]
}

da_result_columns <- function(){
  c("Method", "Model", "Feature_ID", "Label", "Variable", "Contrast", "Reference", "Effect", "Effect_type",
    "SE", "P_value", "Q_value", "Stability", "Enriched_in")
}

prepare_da_inputs <- function(profile, metadata.df, variables.v){
  check_profile(profile)
  require_columns(metadata.df, variables.v, "Metadata")
  metadata.df <- metadata.df[metadata.df$Sample_ID %in% profile_samples(profile), , drop = FALSE]
  metadata.df <- metadata.df[stats::complete.cases(metadata.df[, variables.v, drop = FALSE]), , drop = FALSE]
  metadata.df[variables.v] <- lapply(metadata.df[variables.v], function(x){
    if (is.numeric(x)) x else droplevels(as.factor(x))
  })
  values.m <- profile$values[, metadata.df$Sample_ID, drop = FALSE]
  values.m <- values.m[rowSums(values.m != 0) > 0, , drop = FALSE]
  safe_ids.v <- stats::setNames(sprintf("F%05d", seq_len(nrow(values.m))), rownames(values.m))
  safe.m <- values.m
  rownames(safe.m) <- safe_ids.v
  model.df <- data.frame(metadata.df[, variables.v, drop = FALSE], row.names = metadata.df$Sample_ID)
  list(values = safe.m, metadata = model.df, feature_ids = stats::setNames(names(safe_ids.v), safe_ids.v),
       labels = profile_labels(profile))
}

finish_da_table <- function(result.df, inputs.l){
  result.df$Feature_ID <- unname(inputs.l$feature_ids[result.df$Feature_ID])
  result.df$Label <- unname(inputs.l$labels[result.df$Feature_ID])
  for (column.s in setdiff(da_result_columns(), names(result.df))) result.df[[column.s]] <- NA
  result.df[, da_result_columns()]
}

da_formula <- function(variable, covariates, random_effects){
  random.v <- if (length(random_effects) > 0) paste0("(1|", random_effects, ")") else character()
  paste("~", paste(c(variable, covariates, random.v), collapse = " + "))
}

#' Differential abundance with MaAsLin3
#'
#' Features are renamed to safe IDs before fitting and mapped back afterwards, so
#' taxonomy strings are never mangled. Both the abundance and prevalence models are returned.
#'
#' @param profile A `gm_profile` (relative abundance, counts or normalised RPKM).
#' @param metadata.df Metadata.
#' @param variable Variable of interest.
#' @param covariates Additional fixed effects.
#' @param random_effects Optional random effect columns (e.g. subject for repeated measures).
#' @param normalization,transform Passed to [maaslin3::maaslin3()]; the default normalisation is
#'   `"TSS"` for relative abundance and counts and `"NONE"` for already normalised values.
#' @param min_prevalence Minimum fraction of samples with the feature.
#' @param output_dir Directory for MaAsLin3's own output.
#' @param cores Number of cores.
#' @return Tidy data frame (see [run_differential_abundance()]).
#' @export
run_maaslin3 <- function(profile, metadata.df, variable, covariates = character(), random_effects = character(),
                         normalization = NULL, transform = "LOG", min_prevalence = 0.1, output_dir = tempfile("maaslin3_"),
                         cores = 1){
  require_pkg("maaslin3", "for MaAsLin3 differential abundance")
  inputs.l <- prepare_da_inputs(profile, metadata.df, c(variable, covariates, random_effects))
  normalization <- normalization %||% if (profile$value_type %in% c("normalised_rpkm", "rpkm", "copy_number")) "NONE" else "TSS"
  formula.s <- da_formula(variable, covariates, random_effects)
  fit.l <- suppressWarnings(maaslin3::maaslin3(
    input_data = as.data.frame(t(inputs.l$values)), input_metadata = inputs.l$metadata, output = output_dir,
    formula = formula.s, normalization = normalization, transform = transform, min_prevalence = min_prevalence,
    plot_summary_plot = FALSE, plot_associations = FALSE, save_models = FALSE, verbosity = "ERROR", cores = cores
  ))
  results.l <- lapply(c(abundance = "fit_data_abundance", prevalence = "fit_data_prevalence"), function(slot.s){
    x.df <- fit.l[[slot.s]]$results
    x.df <- x.df[x.df$metadata == variable, , drop = FALSE]
    data.frame(Method = "MaAsLin3", Model = sub("fit_data_", "", slot.s), Feature_ID = x.df$feature, Variable = variable,
               Contrast = x.df$value, Effect = x.df$coef, Effect_type = "log2 coefficient", SE = x.df$stderr,
               P_value = x.df$pval_individual, Q_value = x.df$qval_individual, stringsAsFactors = FALSE)
  })
  result.df <- do.call(rbind, results.l)
  result.df$Reference <- if (is.factor(inputs.l$metadata[[variable]])) levels(inputs.l$metadata[[variable]])[1] else NA
  result.df$Enriched_in <- enriched_group(result.df$Effect, result.df$Contrast, result.df$Reference)
  finish_da_table(result.df, inputs.l)
}

#' Differential abundance with LinDA
#'
#' @inheritParams run_maaslin3
#' @param data_type `"count"` for read counts, otherwise `"proportion"` (default from the value type).
#' @param min_prevalence Minimum fraction of samples with the feature.
#' @param min_nonzero Features nonzero in fewer samples are left out before fitting; LinDA
#'   flags features below 3 as having virtually no statistical power.
#' @return Tidy data frame (see [run_differential_abundance()]).
#' @export
run_linda <- function(profile, metadata.df, variable, covariates = character(), random_effects = character(),
                      data_type = NULL, min_prevalence = 0.1, min_nonzero = 3){
  require_pkg("MicrobiomeStat", "for LinDA differential abundance")
  inputs.l <- prepare_da_inputs(profile, metadata.df, c(variable, covariates, random_effects))
  inputs.l$values <- inputs.l$values[rowSums(inputs.l$values != 0) >= min_nonzero, , drop = FALSE]
  if (nrow(inputs.l$values) == 0) cli::cli_abort("No features are nonzero in at least {min_nonzero} samples")
  data_type <- data_type %||% if (profile$value_type == "read_count") "count" else "proportion"
  formula.s <- da_formula(variable, covariates, random_effects)
  # With random effects lme4 reports every feature with no between-subject variance; that is expected, not a problem
  fit.l <- withCallingHandlers(
    MicrobiomeStat::linda(inputs.l$values, inputs.l$metadata, formula = formula.s, feature.dat.type = data_type,
                          prev.filter = min_prevalence, verbose = FALSE),
    message = function(m) if (grepl("boundary \\(singular\\) fit", conditionMessage(m))) invokeRestart("muffleMessage"))
  terms.v <- names(fit.l$output)
  terms.v <- terms.v[startsWith(terms.v, variable)]
  result.df <- do.call(rbind, lapply(terms.v, function(term.s){
    x.df <- fit.l$output[[term.s]]
    data.frame(Method = "LinDA", Model = "abundance", Feature_ID = rownames(x.df), Variable = variable,
               Contrast = if (is.factor(inputs.l$metadata[[variable]])) substring(term.s, nchar(variable) + 1) else variable,
               Effect = x.df$log2FoldChange, Effect_type = "log2 fold change", SE = x.df$lfcSE,
               P_value = x.df$pvalue, Q_value = x.df$padj, stringsAsFactors = FALSE)
  }))
  result.df$Reference <- if (is.factor(inputs.l$metadata[[variable]])) levels(inputs.l$metadata[[variable]])[1] else NA
  result.df$Enriched_in <- enriched_group(result.df$Effect, result.df$Contrast, result.df$Reference)
  finish_da_table(result.df, inputs.l)
}

enriched_group <- function(effect.v, contrast.v, reference.v){
  ifelse(is.na(reference.v) | is.na(effect.v), NA_character_, ifelse(effect.v > 0, contrast.v, reference.v))
}

#' Discriminant features with sparse PLS-DA (mixOmics)
#'
#' Values get a pseudocount of half the smallest non-zero value and a CLR transform.
#' The number of features per component is tuned by repeated cross-validation, and
#' selection stability is taken from the per-component `perf()` output.
#'
#' @inheritParams run_maaslin3
#' @param n_components Maximum number of components.
#' @param folds,n_repeats Cross-validation folds (capped at the smallest group size) and repeats.
#' @param test_keep_x Candidate numbers of features per component.
#' @param seed Random seed.
#' @return List with `results` (tidy table of selected features: loading as effect,
#'   stability, and the group with the highest mean), `model` and `performance`.
#' @export
run_splsda <- function(profile, metadata.df, variable, n_components = 3, folds = 5, n_repeats = 50,
                       test_keep_x = c(5, 10, 15, 20, 30, 50), seed = 1234){
  require_pkg("mixOmics", "for sPLS-DA")
  inputs.l <- prepare_da_inputs(profile, metadata.df, variable)
  outcome.v <- inputs.l$metadata[[variable]]
  if (!is.factor(outcome.v) || nlevels(outcome.v) < 2){
    cli::cli_abort("sPLS-DA needs a categorical {.field {variable}} with 2+ levels")
  }
  folds <- max(2, min(folds, min(table(outcome.v))))
  values.m <- inputs.l$values[apply(inputs.l$values, 1, stats::var) > 0, , drop = FALSE]
  pseudocount.n <- min(values.m[values.m > 0]) / 2
  x.m <- t(values.m + pseudocount.n)
  test_keep_x <- test_keep_x[test_keep_x < ncol(x.m)]
  if (length(test_keep_x) == 0) test_keep_x <- max(1, ncol(x.m) - 1)
  n_components <- max(2, min(n_components, nlevels(outcome.v) + 1, nrow(x.m) - 1))

  set.seed(seed)
  tuning <- tryCatch(mixOmics::tune.splsda(x.m, outcome.v, ncomp = n_components, logratio = "CLR", test.keepX = test_keep_x,
                                           validation = "Mfold", folds = folds, nrepeat = n_repeats, dist = "max.dist",
                                           measure = "BER", progressBar = FALSE),
                     error = function(e){
                       cli::cli_inform(c("!" = "sPLS-DA tuning failed ({conditionMessage(e)}); using keepX = {min(test_keep_x)}"))
                       NULL
                     })
  if (is.null(tuning)){
    n_components <- 2
    keep_x.v <- rep(min(test_keep_x), n_components)
  } else {
    n_components <- max(2, tuning$choice.ncomp$ncomp %||% 2, na.rm = TRUE)
    keep_x.v <- tuning$choice.keepX[seq_len(n_components)]
  }
  model <- mixOmics::splsda(x.m, outcome.v, ncomp = n_components, keepX = keep_x.v, logratio = "CLR")
  set.seed(seed)
  performance <- tryCatch(mixOmics::perf(model, validation = "Mfold", folds = folds, nrepeat = n_repeats, dist = "max.dist",
                                         progressBar = FALSE),
                          error = function(e){
                            cli::cli_inform(c("!" = "sPLS-DA performance estimation failed: {conditionMessage(e)}"))
                            NULL
                          })

  clr.m <- model$X
  group_means.m <- apply(clr.m, 2, function(x) tapply(x, outcome.v, mean))
  result.df <- do.call(rbind, lapply(seq_len(n_components), function(component.n){
    selected.v <- mixOmics::selectVar(model, comp = component.n)$name
    loadings.v <- model$loadings$X[selected.v, component.n]
    stability.v <- performance$features$stable[[paste0("comp", component.n)]]
    stability.v <- if (is.null(stability.v)) rep(NA_real_, length(selected.v)) else
      stability.v[match(selected.v, names(stability.v))]
    data.frame(Method = "sPLS-DA", Model = paste0("component_", component.n), Feature_ID = selected.v,
               Variable = variable, Contrast = NA_character_, Effect = loadings.v, Effect_type = "loading",
               Stability = as.numeric(stability.v), stringsAsFactors = FALSE)
  }))
  result.df$Enriched_in <- rownames(group_means.m)[apply(group_means.m[, result.df$Feature_ID, drop = FALSE], 2, which.max)]
  result.df$Reference <- levels(outcome.v)[1]
  list(results = finish_da_table(result.df, inputs.l), model = model, performance = performance, tuning = tuning)
}

#' Run differential abundance with several methods
#'
#' Every method returns the same columns: `Method`, `Model`, `Feature_ID`, `Label`,
#' `Variable`, `Contrast` (level compared with `Reference`), `Effect`, `Effect_type`, `SE`,
#' `P_value`, `Q_value`, `Stability` (sPLS-DA) and `Enriched_in` (the group in which the
#' feature is higher). All results are kept, not only significant ones.
#'
#' @inheritParams run_maaslin3
#' @param methods Any of `"maaslin3"`, `"linda"`, `"splsda"`.
#' @param splsda_repeats Cross-validation repeats for sPLS-DA.
#' @param seed Random seed.
#' @return List with `results` (combined table) and `splsda` (model, when run).
#' @export
run_differential_abundance <- function(profile, metadata.df, variable, covariates = character(),
                                       random_effects = character(), methods = c("maaslin3", "linda", "splsda"),
                                       min_prevalence = 0.1, splsda_repeats = 50, output_dir = tempfile("maaslin3_"),
                                       seed = 1234){
  methods <- match.arg(methods, several.ok = TRUE)
  warn_within_sample(profile, "Differential abundance")
  results.l <- list()
  splsda.l <- NULL
  if ("maaslin3" %in% methods){
    results.l$maaslin3 <- try_step(run_maaslin3(profile, metadata.df, variable, covariates, random_effects,
                                                min_prevalence = min_prevalence, output_dir = output_dir), "MaAsLin3")
  }
  if ("linda" %in% methods){
    results.l$linda <- try_step(run_linda(profile, metadata.df, variable, covariates, random_effects,
                                          min_prevalence = min_prevalence), "LinDA")
  }
  if ("splsda" %in% methods){
    filtered.p <- filter_profile(profile, min_prevalence = min_prevalence)
    splsda.l <- try_step(run_splsda(filtered.p, metadata.df, variable, n_repeats = splsda_repeats, seed = seed), "sPLS-DA")
    results.l$splsda <- splsda.l$results
  }
  list(results = do.call(rbind, Filter(Negate(is.null), results.l)), splsda = splsda.l)
}

#' Combine differential abundance results into a consensus
#'
#' A feature counts as significant for MaAsLin3 (abundance model) and LinDA when
#' `Q_value <= alpha`, and for sPLS-DA when selected with `Stability >= min_stability`.
#' Methods agree when they call the feature enriched in the same group, which works for
#' any number of groups (the dump scripts compared incompatible group labels).
#'
#' @param results.df `results` from [run_differential_abundance()].
#' @param alpha Q value threshold.
#' @param min_stability sPLS-DA selection stability threshold.
#' @return Data frame, one row per feature and enriched group, with per-method effects,
#'   q values and `N_methods` (number of methods calling it).
#' @export
da_consensus <- function(results.df, alpha = 0.05, min_stability = 0.7){
  if (is.null(results.df) || nrow(results.df) == 0) return(data.frame())
  significant.df <- results.df[
    (results.df$Method %in% c("MaAsLin3", "LinDA") & results.df$Model == "abundance" &
       !is.na(results.df$Q_value) & results.df$Q_value <= alpha) |
      (results.df$Method == "sPLS-DA" & !is.na(results.df$Stability) & results.df$Stability >= min_stability), , drop = FALSE]
  if (nrow(significant.df) == 0) return(data.frame())
  significant.df <- significant.df[order(significant.df$Method, -abs(significant.df$Effect)), , drop = FALSE]
  significant.df <- significant.df[!duplicated(significant.df[, c("Method", "Feature_ID", "Enriched_in")]), , drop = FALSE]

  keys.df <- unique(significant.df[, c("Feature_ID", "Label", "Variable", "Enriched_in")])
  for (method.s in c("MaAsLin3", "LinDA", "sPLS-DA")){
    method.df <- significant.df[significant.df$Method == method.s, , drop = FALSE]
    index.v <- match(paste(keys.df$Feature_ID, keys.df$Enriched_in), paste(method.df$Feature_ID, method.df$Enriched_in))
    prefix.s <- gsub("[^A-Za-z0-9]", "", method.s)
    keys.df[[paste0(prefix.s, "_effect")]] <- method.df$Effect[index.v]
    keys.df[[paste0(prefix.s, if (method.s == "sPLS-DA") "_stability" else "_q")]] <-
      if (method.s == "sPLS-DA") method.df$Stability[index.v] else method.df$Q_value[index.v]
  }
  keys.df$N_methods <- rowSums(!is.na(keys.df[, grep("_effect$", names(keys.df)), drop = FALSE]))
  conflict.v <- keys.df$Feature_ID[duplicated(keys.df$Feature_ID)]
  keys.df$Conflicting_direction <- keys.df$Feature_ID %in% conflict.v
  keys.df <- keys.df[order(-keys.df$N_methods, keys.df$Enriched_in, keys.df$Label), , drop = FALSE]
  rownames(keys.df) <- NULL
  keys.df
}

#' Plot consensus differential abundance effects
#'
#' One bar per method and feature, coloured by the group the feature is higher in.
#'
#' @param consensus.df Result of [da_consensus()].
#' @param colours.v Named colours for the enriched groups.
#' @param min_methods Minimum number of agreeing methods to show a feature.
#' @param max_features Maximum number of features shown (most methods first).
#' @param methods Methods shown, as panels in this order.
#' @param bar_width Bar width.
#' @param label_width Wrap feature labels longer than this many characters; `NULL` for no wrapping.
#' @param x_label X axis title.
#' @param legend_title Legend title.
#' @param base_size Base font size, in points.
#' @param legend_position Legend position.
#' @return A ggplot, or `NULL` when nothing passes.
#' @export
plot_da_consensus <- function(consensus.df, colours.v = NULL, min_methods = 2, max_features = 40,
                              methods = c("MaAsLin3", "LinDA", "sPLS-DA"), bar_width = 0.7, label_width = NULL,
                              x_label = "Effect (log2 coefficient, log2 fold change or loading)", legend_title = "Higher in",
                              base_size = 9, legend_position = "right"){
  if (nrow(consensus.df) == 0) return(NULL)
  shown.df <- utils::head(consensus.df[consensus.df$N_methods >= min_methods, , drop = FALSE], max_features)
  if (nrow(shown.df) == 0) return(NULL)
  effect_columns.v <- da_effect_columns()[methods]
  long.df <- do.call(rbind, lapply(names(effect_columns.v), function(method.s){
    data.frame(Label = shown.df$Label, Enriched_in = shown.df$Enriched_in, Method = method.s,
               Effect = shown.df[[effect_columns.v[[method.s]]]] %||% NA_real_)
  }))
  long.df <- long.df[!is.na(long.df$Effect), , drop = FALSE]
  if (nrow(long.df) == 0) return(NULL)
  if (!is.null(label_width)) long.df$Label <- wrap_labels(long.df$Label, label_width)
  shown_labels.v <- if (is.null(label_width)) shown.df$Label else wrap_labels(shown.df$Label, label_width)
  long.df$Label <- factor(long.df$Label, levels = rev(unique(shown_labels.v)))
  long.df$Method <- factor(long.df$Method, levels = methods)
  if (is.null(colours.v)) colours.v <- assign_colours(long.df$Enriched_in)
  consensus.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data$Effect, y = .data$Label, fill = .data$Enriched_in)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_col(width = bar_width) +
    ggplot2::facet_grid(~Method, scales = "free_x") +
    ggplot2::scale_fill_manual(values = colours.v, name = legend_title) +
    ggplot2::labs(x = x_label, y = NULL) +
    theme_gmaio(base_size = base_size, legend_position = legend_position)
  with_size(consensus.gg, 8 + 5 * length(unique(long.df$Method)), 3 + 0.4 * length(unique(long.df$Label)))
}

da_effect_columns <- function(){
  c(MaAsLin3 = "MaAsLin3_effect", LinDA = "LinDA_effect", `sPLS-DA` = "sPLSDA_effect")
}

wrap_labels <- function(labels.v, width){
  vapply(labels.v, function(x) paste(strwrap(x, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
}

#' Heatmap of differentially abundant features
#'
#' Abundance of every feature called by at least `min_methods` methods, with row tracks for
#' the group the feature is higher in, the group each method called it for (so disagreement
#' between methods shows as a colour mismatch) and features whose methods conflict in direction.
#' Built on [plot_heatmap()]; any of its options can be passed to change the figure.
#'
#' @param profile The `gm_profile` the differential abundance was run on.
#' @param consensus.df Result of [da_consensus()].
#' @param metadata.df Metadata.
#' @param group Metadata column that was tested; used for the column split and the track colours.
#' @param palettes.l Project palettes (`project.l$palettes`).
#' @param min_methods Minimum number of methods calling a feature.
#' @param max_features Maximum number of features shown (most methods first).
#' @param methods Methods drawn as tracks.
#' @param rank_annotation Feature column drawn as the first track (e.g. `"Phylum"`), when present.
#' @param not_called_colour Track colour where a method did not call the feature.
#' @param conflict_colour Track colour for features called higher in different groups.
#' @param ... Options passed to [plot_heatmap()], overriding the defaults set here
#'   (`column_split = group`, `row_split = "Higher_in"`).
#' @return A `Heatmap`, or `NULL` when fewer than two features pass.
#' @export
plot_da_heatmap <- function(profile, consensus.df, metadata.df, group, palettes.l = list(), min_methods = 1,
                            max_features = 60, methods = c("MaAsLin3", "LinDA", "sPLS-DA"), rank_annotation = "Phylum",
                            not_called_colour = "#F2F2F2", conflict_colour = "#D7191C", ...){
  if (nrow(consensus.df) == 0) return(NULL)
  shown.df <- consensus.df[consensus.df$N_methods >= min_methods, , drop = FALSE]
  shown.df <- shown.df[order(-shown.df$N_methods), , drop = FALSE]
  features.v <- utils::head(unique(shown.df$Feature_ID), max_features)
  features.v <- intersect(features.v, rownames(profile$values))
  if (length(features.v) < 2) return(NULL)

  subset.p <- subset_features(profile, features.v)
  features.df <- subset.p$features
  best.df <- shown.df[match(features.df$Feature_ID, shown.df$Feature_ID), , drop = FALSE]
  features.df$Higher_in <- best.df$Enriched_in
  effect_columns.v <- da_effect_columns()[methods]
  effect_columns.v <- effect_columns.v[effect_columns.v %in% names(consensus.df)]
  effect_columns.v <- effect_columns.v[vapply(effect_columns.v, function(x) any(!is.na(consensus.df[[x]])), logical(1))]
  for (method.s in names(effect_columns.v)){
    called.df <- consensus.df[!is.na(consensus.df[[effect_columns.v[[method.s]]]]), , drop = FALSE]
    called.v <- called.df$Enriched_in[match(features.df$Feature_ID, called.df$Feature_ID)]
    features.df[[method.s]] <- ifelse(is.na(called.v), "Not called", called.v)
  }
  conflict.v <- features.df$Feature_ID %in% consensus.df$Feature_ID[consensus.df$Conflicting_direction]
  features.df$Conflict <- ifelse(conflict.v, "Yes", "No")
  subset.p$features <- features.df

  group_levels.v <- unique(stats::na.omit(c(as.character(metadata.df[[group]]), features.df$Higher_in)))
  group_colours.v <- if (is.null(palettes.l[[group]])) assign_colours(group_levels.v) else
    get_palette(palettes.l, group, group_levels.v)
  track_colours.l <- c(list(Higher_in = group_colours.v, Conflict = c(Yes = conflict_colour, No = not_called_colour)),
                       stats::setNames(rep(list(c(group_colours.v, `Not called` = not_called_colour)),
                                           length(effect_columns.v)), names(effect_columns.v)))
  rank.v <- intersect(rank_annotation, names(features.df))
  tracks.v <- c(rank.v, "Higher_in", names(effect_columns.v), "Conflict")
  defaults.l <- list(profile = subset.p, metadata.df = metadata.df, palettes.l = palettes.l, column_split = group,
                     row_split = "Higher_in", row_annotation = tracks.v, row_annotation_colours = track_colours.l,
                     row_annotation_legends = c(rank.v, "Conflict"),
                     legend_title = sub(" (", "\n(", value_type_label(profile$value_type), fixed = TRUE))
  do.call(plot_heatmap, utils::modifyList(defaults.l, list(...)))
}

#' Boxplots of feature abundance by group
#'
#' One panel per feature, e.g. the features that several differential abundance methods agree on.
#'
#' @param profile A `gm_profile`.
#' @param metadata.df Metadata.
#' @param feature_ids.v Features to show, in panel order.
#' @param group Metadata column on the x axis.
#' @param colours.v Named colours for `group`.
#' @param log_scale Use a log10 y axis (zeros are shown at half the smallest non-zero value).
#' @param ncol Facet columns.
#' @param shape_by Optional metadata column shown as point shape.
#' @param shapes Named point shapes for `shape_by`; filled shapes (21-25) take the group colour.
#' @param box_width,box_alpha Box width and fill transparency.
#' @param point_size,jitter_width Point size and horizontal jitter.
#' @param label_width Wrap panel titles longer than this many characters.
#' @param free_y Give every panel its own y axis.
#' @param x_label,y_label Axis titles; `y_label` defaults from the profile value type.
#' @param x_text_angle Angle of the group labels.
#' @param base_size Base font size, in points.
#' @param legend_position Legend position when `shape_by` is used.
#' @return A ggplot.
#' @export
plot_feature_boxplots <- function(profile, metadata.df, feature_ids.v, group, colours.v = NULL, log_scale = TRUE, ncol = 4,
                                  shape_by = NULL, shapes = NULL, box_width = 0.6, box_alpha = 0.5, point_size = 1.5,
                                  jitter_width = 0.12, label_width = 22, free_y = TRUE, x_label = NULL, y_label = NULL,
                                  x_text_angle = 30, base_size = 8, legend_position = "right"){
  subset.p <- subset_features(profile, feature_ids.v)
  long.df <- profile_to_long(subset.p, metadata.df[metadata.df$Sample_ID %in% profile_samples(subset.p), , drop = FALSE])
  long.df$Label <- factor(long.df$Label, levels = subset.p$features$Label[match(feature_ids.v, subset.p$features$Feature_ID)])
  if (log_scale){
    floor.n <- min(long.df$Value[long.df$Value > 0]) / 2
    long.df$Value[long.df$Value == 0] <- floor.n
  }
  if (is.null(colours.v)) colours.v <- assign_colours(long.df[[group]])
  point_aes <- if (is.null(shape_by)) ggplot2::aes(fill = .data[[group]]) else
    ggplot2::aes(fill = .data[[group]], shape = .data[[shape_by]])
  boxplots.gg <- ggplot2::ggplot(long.df, ggplot2::aes(x = .data[[group]], y = .data$Value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[group]]), outlier.shape = NA, width = box_width, alpha = box_alpha,
                          linewidth = 0.3) +
    do.call(ggplot2::geom_point, c(list(mapping = point_aes, size = point_size, colour = "grey15", stroke = 0.3,
                                        position = ggplot2::position_jitter(width = jitter_width, height = 0, seed = 1)),
                                   if (is.null(shape_by)) list(shape = 21))) +
    ggplot2::facet_wrap(~Label, scales = if (free_y) "free_y" else "fixed", ncol = ncol,
                        labeller = ggplot2::label_wrap_gen(label_width)) +
    ggplot2::scale_fill_manual(values = colours.v, guide = "none") +
    ggplot2::labs(x = x_label, y = y_label %||% value_type_label(profile$value_type)) +
    theme_gmaio(base_size = base_size, legend_position = legend_position) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = if (x_text_angle == 0) 0.5 else 1))
  if (!is.null(shape_by)){
    shape_levels.v <- levels(as.factor(long.df[[shape_by]]))
    shapes.v <- shapes %||% stats::setNames(rep_len(c(21, 24, 22, 23, 25), length(shape_levels.v)), shape_levels.v)
    boxplots.gg <- boxplots.gg + ggplot2::scale_shape_manual(values = shapes.v, name = variable_label(shape_by)) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey60")))
  }
  if (log_scale) boxplots.gg <- boxplots.gg + ggplot2::scale_y_log10(labels = plain_number)
  n_panels.n <- length(feature_ids.v)
  n_columns.n <- min(ncol, n_panels.n)
  with_size(boxplots.gg, 4.5 * n_columns.n + 1 + if (is.null(shape_by)) 0 else 2.5, 4.5 * ceiling(n_panels.n / ncol) + 1)
}

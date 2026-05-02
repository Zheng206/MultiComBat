.get_target_residuals <- function(data, bat, covar,
                                  model = lm, formula = NULL,
                                  ref.batch = NULL) {
  m <- length(data)
  data <- lapply(seq_len(m), function(i) data.frame(data[[i]]))

  batch_result  <- lapply(seq_len(m), function(i)
    batch_matrix(bat[[i]], ref.batch = ref.batch))

  fitted_model  <- lapply(seq_len(m), function(i)
    model_fitting(data[[i]], batch_result[[i]]$batch_matrix,
                  covar[[i]], model, formula))

  stand_result  <- lapply(seq_len(m), function(i)
    standardize_data(fitted_model[[i]], batch_result[[i]]))

  resid_list <- lapply(seq_len(m), function(i) stand_result[[i]]$data_stand)

  list(resid_list = resid_list,
       batch_result = batch_result)
}

.feature_metric_cov <- function(resid_list, g, idx = NULL) {
  mat <- do.call(cbind, lapply(resid_list, function(D) {
    if (is.null(idx)) D[, g] else D[idx, g]
  }))
  cov(mat)
}

.mat_to_long <- function(S, label, m) {
  metric_labels <- paste0("M", seq_len(m))
  expand.grid(row = metric_labels, col = metric_labels) |>
    dplyr::mutate(value = as.vector(S), source = label)
}

.compute_avg_cov <- function(resid_list, batch_result, use_correlation = TRUE) {
  m  <- length(resid_list)
  G  <- ncol(resid_list[[1]])

  batch_vec    <- batch_result[[1]]$batch_index
  batch_levels <- names(batch_vec)

  S_pool <- Reduce(`+`, lapply(seq_len(G), function(g) {
    S <- .feature_metric_cov(resid_list, g, idx = NULL)
    if (use_correlation) cov2cor(S) else S
  })) / G

  S_by_batch <- lapply(batch_levels, function(b) {
    idx <- batch_vec[[b]]
    Reduce(`+`, lapply(seq_len(G), function(g) {
      S <- .feature_metric_cov(resid_list, g, idx = idx)
      if (use_correlation) cov2cor(S) else S
    })) / G
  })
  names(S_by_batch) <- batch_levels

  list(pool = S_pool, by_batch = S_by_batch)
}

.compute_feature_cov <- function(resid_list, batch_result, g,
                                 use_correlation = TRUE) {
  batch_vec    <- batch_result[[1]]$batch_index
  batch_levels <- names(batch_vec)

  S_pool <- .feature_metric_cov(resid_list, g, idx = NULL)
  if (use_correlation) S_pool <- cov2cor(S_pool)

  S_by_batch <- lapply(batch_levels, function(b) {
    S <- .feature_metric_cov(resid_list, g, idx = batch_vec[[b]])
    if (use_correlation) cov2cor(S) else S
  })
  names(S_by_batch) <- batch_levels

  list(pool = S_pool, by_batch = S_by_batch)
}

.cov_heatmap <- function(df, title, use_correlation) {
  lim <- max(abs(df$value), na.rm = TRUE)

  ggplot2::ggplot(df, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::facet_wrap(~ source, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low      = "#2471A3", mid = "white", high = "#C0392B",
      midpoint = 0, limits = c(-lim, lim),
      name     = if (use_correlation) "Correlation" else "Covariance"
    ) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "#F0F0F0"),
      strip.text       = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = title, x = "Modality", y = "Modality")
}



#' Plot the pooled cross-modality covariance or correlation structure
#'
#' This diagnostic visualizes the average cross-modality covariance structure
#' after removing fixed effects and, optionally, batch effects relative to a
#' reference batch. It pools information across all features to display a
#' single summary heatmap, which can help assess whether the residual
#' cross-modality structure is strong and potentially worth preserving.
#'
#' If the pooled structure shows clear and stable off-diagonal patterns,
#' `robust_cov = TRUE` may be worth considering, as aggressive whitening could
#' remove covariance structure shared across modalities.
#'
#' @param data      list of m matrices (observations x features)
#' @param bat       list of m batch factor vectors
#' @param covar     list of m covariate data.frames
#' @param model     model function passed to model_fitting (default lm)
#' @param formula   optional formula
#' @param ref.batch optional reference batch label
#' @param use_correlation if TRUE (default), display correlations not covariances
#'
#' @return ggplot object: single pooled heatmap averaged across features
diag_plot_pooled_cov <- function(data, bat, covar,
                                 model         = lm,
                                 formula       = NULL,
                                 ref.batch     = NULL,
                                 use_correlation = TRUE) {
  prep <- .get_target_residuals(data, bat, covar, model, formula, ref.batch)
  obj  <- .compute_avg_cov(prep$resid_list, prep$batch_result, use_correlation)

  m  <- length(data)
  df <- .mat_to_long(obj$pool, "Pooled (avg across features)", m)

  type_label <- if (use_correlation) "correlation" else "covariance"
  .cov_heatmap(df,
               title = paste("Average cross-modality", type_label,
                             "(pooled) \u2014 all features"),
               use_correlation = use_correlation)
}


#' Plot pooled and batch-specific cross-modality covariance structure
#'
#' This diagnostic visualizes the cross-modality covariance or correlation
#' structure after removing fixed effects and, optionally, batch effects
#' relative to a reference batch. It displays the pooled structure together
#' with one heatmap per batch, making it easier to assess heterogeneity in
#' residual covariance patterns across batches.
#'
#' Large differences between the pooled panel and the batch-specific panels
#' may indicate substantial batch-related covariance heterogeneity, in which
#' case `robust_cov = TRUE` may be worth considering.
#'
#' If `feature_idx = NULL`, the function averages covariance structure across
#' all features. If a feature index is provided, the heatmap is constructed
#' using only that feature.
#'
#' @param data A list of length `m`, where each element is a matrix of
#'   observations by features for one modality.
#' @param bat A list of length `m` containing batch labels for each modality.
#' @param covar A list of length `m` containing covariate data frames.
#' @param model A modeling function passed to `model_fitting()`. Default is `lm`.
#' @param formula An optional model formula.
#' @param ref.batch An optional reference batch label.
#' @param use_correlation Logical; if `TRUE` (default), plot correlations.
#'   Otherwise, plot covariances.
#' @param feature_idx Integer; if `NULL` (default), the covariance structure is
#'   averaged across all features. If a feature index is provided, only that
#'   feature is shown.
#'
#' @return A `ggplot2` faceted heatmap showing the pooled covariance or
#'   correlation structure together with one panel per batch.
diag_plot_batch_cov <- function(data, bat, covar,
                                model           = lm,
                                formula         = NULL,
                                ref.batch       = NULL,
                                use_correlation = TRUE,
                                feature_idx     = NULL) {
  prep <- .get_target_residuals(data, bat, covar, model, formula, ref.batch)
  m    <- length(data)

  if (is.null(feature_idx)) {
    obj        <- .compute_avg_cov(prep$resid_list, prep$batch_result,
                                   use_correlation)
    title_stem <- "Average cross-modality"
    feat_label <- "all features"
  } else {
    obj        <- .compute_feature_cov(prep$resid_list, prep$batch_result,
                                       g = feature_idx, use_correlation)
    title_stem <- "Cross-modality"
    feat_label <- paste("feature", feature_idx)
  }

  type_label <- if (use_correlation) "correlation" else "covariance"

  df <- dplyr::bind_rows(
    .mat_to_long(obj$pool, "Pooled", m),
    dplyr::bind_rows(lapply(names(obj$by_batch), function(b)
      .mat_to_long(obj$by_batch[[b]], paste0("Batch ", b), m)))
  )

  df$source <- factor(df$source,
                      levels = c("Pooled", paste0("Batch ", names(obj$by_batch))))

  .cov_heatmap(df,
               title = paste(title_stem, type_label, "\u2014", feat_label),
               use_correlation = use_correlation)
}





#' Diagnostics: summarize fitted models for harmonization checks (S3)
#'
#' @param diag_model An object returned by \code{diag_model_gen()} whose class
#'   has been tagged to \code{"lm_model"}, \code{"gam_model"}, or \code{"lmerMod_model"}.
#' @param ... Passed to methods.
#' @export
diag_model_summary <- function(diag_model, ...) {
  UseMethod("diag_model_summary")
}

#' @describeIn diag_model_summary Summary for linear models
#' @export
diag_model_summary.lm_model <- function(diag_model, ...){
  fits <- diag_model$fits
  coefs <- sapply(fits, coef)
  data <- diag_model$data
  covar_name <- names(diag_model$mod)
  covar_name <- covar_name[which(!grepl("batch", covar_name))]
  covar <- diag_model$mod[, covar_name] |> as.data.frame()
  colnames(covar) <- covar_name
  covar <- model.matrix(as.formula(paste0("~ ", paste0(covar_name, collapse = "+"), "-1")), data = covar)
  resid_add <- data - covar %*% coefs[colnames(covar), ]
  resid_mul <- sapply(fits, resid)
  diag_summary <- list(coef = coefs, resid_add = resid_add, resid_mul = resid_mul)
  class(diag_summary) <- "diag_summary"
  return(diag_summary)
}


#' @describeIn diag_model_summary Summary for GAMs (mgcv)
#' @export
diag_model_summary.gam_model <- function(diag_model, ...){
  fits <- diag_model$fits
  coefs <- sapply(fits, coef)
  data <- diag_model$data
  fix_effects <- sapply(fits, function(model){
    covar_matrix <- model.matrix(model)
    covar_matrix <- covar_matrix[, which(!grepl("batch", colnames(covar_matrix)))]
    coef_value <- coef(model)[which(!grepl("batch", names(coef(model))))]
    covar_matrix %*% coef_value
  }) |> data.frame()
  resid_add <- data - fix_effects
  resid_mul <- sapply(fits, resid)
  diag_summary <- list(coef = coefs, resid_add = resid_add, resid_mul = resid_mul)
  class(diag_summary) <- "diag_summary"
  return(diag_summary)
}


#' @describeIn diag_model_summary Summary for mixed models (lme4)
#' @param random Character; name of the grouping factor for random effects (e.g., subject id).
#' @export
diag_model_summary.lmerMod_model <- function(diag_model, random, ...){
  fits <- diag_model$fits
  coefs <- sapply(fits, function(x) {
    coef_result <- coef(x)[[1]]
    coef_result <- coef_result[, which(!grepl("(Intercept)", names(coef_result)))] %>% distinct()
    return(coef_result)
  }) %>% as.matrix()
  rand_effects <- sapply(fits, function(x) {
    coef_result <- coef(x)[[1]]
    coef_result <- coef_result[, which(grepl("(Intercept)", names(coef_result)))]
    return(coef_result)
  }) %>% as.matrix()
  sub_id <- rownames(coef(fits[[1]])[[1]])
  rand_df <- cbind(sub_id, rand_effects) |> data.frame() %>% mutate(sub_id = as.factor(sub_id), across(-1, ~as.numeric(.)))
  rand_df <- data.frame(sub_id = as.factor(diag_model$mod[[random]])) %>% left_join(rand_df, by = "sub_id")
  data <- diag_model$data
  covar_name <- names(diag_model$mod)
  covar_name <- covar_name[which(!grepl("batch", covar_name))]
  covar_name <- covar_name[which(!grepl(random, covar_name))]
  covar <- diag_model$mod[, covar_name] |> as.data.frame()
  colnames(covar) <- covar_name
  covar <- model.matrix(as.formula(paste0("~ ", paste0(covar_name, collapse = "+"), "-1")), data = covar)
  coef_mat <- as.matrix(coefs[colnames(covar), ])
  coef_mat <- apply(coef_mat, 2, as.numeric)
  resid_add <- data - (covar %*% coef_mat + as.matrix(rand_df[, -1]))
  resid_mul <- sapply(fits, resid)
  diag_summary <- list(coef = coefs, resid_add = resid_add, resid_mul = resid_mul)
  class(diag_summary) <- "diag_summary"
  return(diag_summary)
}


#' Generate a diagnostics-ready model bundle
#'
#' Fits per-feature models with or without batch dummies and tags the result
#' with a class used by \code{diag_model_summary()}.
#'
#' @param bat Factor (or coercible) of batch labels (length n).
#' @param data Numeric n x G matrix/data frame of responses.
#' @param covar Data frame of covariates (or NULL for a null model).
#' @param model Modeling function accepting \code{formula} and \code{data}.
#' @param formula RHS formula for the covariates (exclude batch; added internally if \code{bat_adjust = TRUE}).
#' @param ref.batch Optional reference batch level.
#' @param bat_adjust Logical; if TRUE (default) include batch in fitting via \code{model_fitting()}.
#' @param ... Passed to the modeling function.
#'
#' @return A list containing at least \code{$fits}, \code{$data}, \code{$mod}, tagged with class
#'   \code{"lm_model"}, \code{"gam_model"}, or \code{"lmerMod_model"} accordingly.
#' @export
diag_model_gen <- function(bat, data, covar, model, formula = NULL, ref.batch = NULL, bat_adjust = TRUE, ...){
  if (bat_adjust){
    batch_result <- batch_matrix(bat, ref.batch = ref.batch)
    fitted_model <- model_fitting(data, batch_result$batch_matrix, covar, model, formula, ...)
    class_type <- class(fitted_model$fits[[1]])
  }else{
    fitted_model <- apply(data, 2, function(y) {
      if(is.null(covar)){
        dat <- data.frame(y = y)}else{
          dat <- data.frame(y = y, covar)
        }
      do.call(model, list(formula = formula, data = dat, ...))
    })
    class_type <- class(fitted_model[[1]])
  }
  class(fitted_model) <- paste0(class_type[1], "_model")
  return(fitted_model)
}

# ---- Univariate tests (means/variances) ------------------------------------

#' One-way ANOVA across batches (per feature)
#'
#' Tests mean differences across batches feature-by-feature using one-way
#' ANOVA on additive residuals. P-values are Bonferroni-adjusted across
#' features.
#'
#' @param resid_add Numeric matrix (n × G). Additive residuals; rows = samples,
#'   columns = features.
#' @param batch Factor of length n with batch labels.
#'
#' @return A list with:
#' \describe{
#'   \item{test_table}{data.frame with columns \code{feature}, \code{p_value}.}
#'   \item{perc.sig}{Character string with percentage of features with adjusted
#'                   p-value < 0.05.}
#' }
#' @examples
#' # anova_test(resid_add, batch)
#' @family univariate-tests
#' @export
anova_test <- function(resid_add, batch){
  G <- ncol(resid_add)
  p_value <- sapply(1:G, function(g) summary(aov(resid_add[,g] ~ batch))[[1]]$"Pr(>F)"[1]) |> p.adjust(method = "bonferroni") |> round(3)
  test_result <- data.frame(cbind(feature = colnames(resid_add), p_value))
  percent_below_0.05 <- mean(p_value < 0.05) * 100
  num_sig <- paste0(sprintf("%.2f", percent_below_0.05), "%")
  return(list(test_table = test_result, perc.sig = num_sig))
}


#' Kruskal–Wallis test across batches (per feature)
#'
#' Non-parametric alternative to one-way ANOVA on additive residuals.
#' P-values are Bonferroni-adjusted across features.
#'
#' @param resid_add Numeric matrix (n × G). Additive residuals; rows = samples,
#'   columns = features.
#' @param batch Factor of length n with batch labels.
#'
#' @return See \code{\link{anova_test}}.
#' @examples
#' # kruskal_test(resid_add, batch)
#' @family univariate-tests
#' @export
kruskal_test <- function(resid_add, batch){
  G <- ncol(resid_add)
  p_value <- sapply(1:G, function(g) kruskal.test(resid_add[,g] ~ batch) %>% tidy() %>% pull(.data[["p.value"]])) |> p.adjust(method = "bonferroni") |> round(3)
  test_result <- data.frame(cbind(feature = colnames(resid_add), p_value))
  percent_below_0.05 <- mean(p_value < 0.05) * 100
  num_sig <- paste0(sprintf("%.2f", percent_below_0.05), "%")
  return(list(test_table = test_result, perc.sig = num_sig))
}

#' Levene's test for equal variances (per feature)
#'
#' Tests variance homogeneity across batches on multiplicative residuals.
#' Requires \pkg{car}. P-values are Bonferroni-adjusted across features.
#'
#' @param resid_mul Numeric matrix (n × G). Multiplicative residuals; rows = samples,
#'   columns = features.
#' @param batch Factor of length n with batch labels.
#'
#' @return See \code{\link{anova_test}}.
#' @examples
#' # lv_test(resid_mul, batch)
#' @family univariate-tests
#' @export


lv_test <- function(resid_mul, batch){
  G <- ncol(resid_mul)
  p_value <- sapply(1:G, function(g) car::leveneTest(resid_mul[,g] ~ batch) %>% tidy() %>% pull(.data[["p.value"]])) |> p.adjust(method = "bonferroni") |> round(3)
  test_result <- data.frame(cbind(feature = colnames(resid_mul), p_value))
  percent_below_0.05 <- mean(p_value < 0.05) * 100
  num_sig <- paste0(sprintf("%.2f", percent_below_0.05), "%")
  return(list(test_table = test_result, perc.sig = num_sig))
}

#' Bartlett's test for equal variances (per feature)
#'
#' Tests equality of variances across batches on multiplicative residuals.
#' Sensitive to non-normality. P-values are Bonferroni-adjusted across features.
#'
#' @param resid_mul Numeric matrix (n × G). Multiplicative residuals; rows = samples,
#'   columns = features.
#' @param batch Factor of length n with batch labels.
#'
#' @return See \code{\link{anova_test}}.
#' @examples
#' # bl_test(resid_mul, batch)
#' @family univariate-tests
#' @export
bl_test <- function(resid_mul, batch){
  G <- ncol(resid_mul)
  p_value <- sapply(1:G, function(g) bartlett.test(resid_mul[,g] ~ batch) %>% tidy() %>% pull(.data[["p.value"]])) |> p.adjust(method = "bonferroni") |> round(3)
  test_result <- data.frame(cbind(feature = colnames(resid_mul), p_value))
  percent_below_0.05 <- mean(p_value < 0.05) * 100
  num_sig <- paste0(sprintf("%.2f", percent_below_0.05), "%")
  return(list(test_table = test_result, perc.sig = num_sig))
}

#' Fligner–Killeen test for equal variances (per feature)
#'
#' Rank-based, robust alternative to Bartlett's test on multiplicative residuals.
#' P-values are Bonferroni-adjusted across features.
#'
#' @param resid_mul Numeric matrix (n × G). Multiplicative residuals; rows = samples,
#'   columns = features.
#' @param batch Factor of length n with batch labels.
#'
#' @return See \code{\link{anova_test}}.
#' @examples
#' # fk_test(resid_mul, batch)
#' @family univariate-tests
#' @export
fk_test <- function(resid_mul, batch){
  G <- ncol(resid_mul)
  p_value <- sapply(1:G, function(g) fligner.test(resid_mul[,g] ~ batch) %>% tidy() %>% pull(.data[["p.value"]])) |> p.adjust(method = "bonferroni") |> round(3)
  test_result <- data.frame(cbind(feature = colnames(resid_mul), p_value))
  percent_below_0.05 <- mean(p_value < 0.05) * 100
  num_sig <- paste0(sprintf("%.2f", percent_below_0.05), "%")
  return(list(test_table = test_result, perc.sig = num_sig))
}


# --- Multivariate tests -----------------------------------------------------

#' Reshape multivariate residuals for MANOVA / Box's M
#'
#' Given diagnostic summaries from multiple measurements, stacks the
#' measurement-specific residuals column-wise by feature to produce inputs
#' for multivariate tests.
#'
#' @param diag_summary List of length m; each element is a \code{"diag_summary"}
#'   with components \code{resid_add} (n × G) and \code{resid_mul} (n × G).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{resid_add}{List of length G; each element is an n × m matrix of
#'                    additive residuals (measurements in columns).}
#'   \item{resid_mul}{List of length G; each element is an n × m matrix of
#'                    multiplicative residuals (measurements in columns).}
#' }
#' @examples
#' # mm <- multi_test_reshape(list(diag1, diag2))
#' @family multivariate-tests
#' @export
multi_test_reshape <- function(diag_summary){
  m <- length(diag_summary)
  G <- ncol(diag_summary[[1]]$resid_add)
  resid_add <- lapply(1:G, function(g){
    do.call(cbind, lapply(1:m, function(i) diag_summary[[i]]$resid_add[, g]))
  })
  resid_mul <- lapply(1:G, function(g){
    do.call(cbind, lapply(1:m, function(i) diag_summary[[i]]$resid_mul[, g]))
  })
  names(resid_add) <- names(resid_mul) <- colnames(diag_summary[[1]]$resid_add)
  return(list(resid_add = resid_add, resid_mul = resid_mul))
}

#' MANOVA across batches (per feature block)
#'
#' Tests mean differences across batches using MANOVA on each feature's
#' stacked measurement matrix from \code{multi_test_reshape()}. P-values are
#' Bonferroni-adjusted across features.
#'
#' @param resid_add List of length G; each element an n × m matrix of additive
#'   residuals (measurements in columns).
#' @param batch Factor of length n with batch labels.
#' @param test Character scalar; which MANOVA statistic to report
#'   (e.g., \code{"Pillai"}, \code{"Wilks"}, \code{"Hotelling-Lawley"}, \code{"Roy"}).
#'
#' @return A list with:
#' \describe{
#'   \item{test_table}{data.frame with columns \code{feature}, \code{p.value}.}
#'   \item{perc.sig}{Character string with percentage of features with adjusted
#'                   p-value ≤ 0.05.}
#' }
#' @examples
#' # mv <- manova_test(mm$resid_add, batch)
#' @family multivariate-tests
#' @export
manova_test <- function(resid_add, batch, test = "Pillai"){
  G <- length(resid_add)
  g_names <- names(resid_add)
  test_result <- lapply(1:G, function(g) {
    manova_result <- manova(resid_add[[g]] ~ batch)
    summary_output <- summary(manova_result, test = test)
    p_value <- summary_output$stats["batch", "Pr(>F)"]
  })
  test_result_mat <- do.call(rbind, test_result) |> p.adjust(method = "bonferroni")
  test_result_df <- data.frame(cbind(g_names, test_result_mat))
  colnames(test_result_df) <- c("feature", "p.value")
  test_result_df <- test_result_df %>% mutate(p.value = round(as.numeric(p.value), 3)) %>% arrange(p.value)
  num_sig <- paste0(sprintf("%.2f", 100*mean(test_result_mat <= 0.05)), "%")
  return(list(test_table = test_result_df, perc.sig = num_sig))
}

#' Box's M test
#'
#' Tests equality of covariance matrices across groups using the
#' chi-square approximation with Box's small-sample correction.
#'
#' @param Y Numeric matrix \eqn{n \times p} (rows = samples, cols = variables).
#' @param group Factor (or coercible) of length n with group labels.
#' @param ridge Small diagonal ridge added to each group covariance for
#'   numerical stability (default \code{1e-8}).
#'
#' @return A list with elements \code{statistic},
#'   \code{parameter} (df), \code{p.value}, and \code{method}.
#' @examples
#' set.seed(123)
#' n <- 60; p <- 3
#' Z1 <- matrix(rnorm(n * p), n, p)
#' Z2 <- matrix(rnorm(n * p), n, p)
#' Sigma2 <- matrix(c(1, 0.6, 0.3,
#'                    0.6, 1, 0.4,
#'                    0.3, 0.4, 1), p, p)
#' Y1 <- Z1
#' Y2 <- Z2 %*% chol(Sigma2)
#' Y  <- rbind(Y1, Y2)
#' g  <- factor(rep(c("A","B"), each = n))
#' out_diff <- boxM_simple(Y, g)
#' out_diff$p.value
#'
#'
#' @export
boxM_simple <- function(Y, group, ridge = 1e-8) {
  Y <- as.matrix(Y)
  if (!is.factor(group)) group <- factor(group)

  ok <- stats::complete.cases(Y, group)
  Y <- Y[ok, , drop = FALSE]
  group <- droplevels(group[ok])

  p <- ncol(Y); k <- nlevels(group)
  ni <- as.integer(table(group)); N <- sum(ni)

  cov_i <- lapply(split(as.data.frame(Y), group), function(df) {
    stats::cov(df) + diag(ridge, p)
  })
  Sp <- Reduce(`+`, Map(function(S, n) (n - 1) * S, cov_i, ni)) / (N - k)

  logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
  Mstat <- (N - k) * logdet(Sp) - sum((ni - 1) * vapply(cov_i, logdet, numeric(1)))
  cfac  <- ((2*p^2 + 3*p - 1) / (6*(p + 1)*(k - 1))) * (sum(1/(ni - 1)) - 1/(N - k))
  X2    <- Mstat * (1 - cfac)
  df    <- (k - 1) * p * (p + 1) / 2
  pval  <- stats::pchisq(X2, df = df, lower.tail = FALSE)

  structure(list(statistic = unname(X2), parameter = unname(df),
                 p.value = unname(pval), method = "Box's M (complete-case)"))
}


#' Box's M test for equality of covariance matrices (per feature block)
#'
#' Tests whether covariance matrices are equal across batches for each feature's
#' stacked measurement matrix (from \code{multi_test_reshape()}). Requires \pkg{heplots}.
#' P-values are Bonferroni-adjusted across features.
#'
#' @param resid_mul List of length G; each element an \eqn{n \times m} matrix of multiplicative
#'   residuals (measurements in columns).
#' @param batch Factor of length n with batch labels.
#'
#' @return See \code{\link{manova_test}}.
#' @examples
#' # bm <- boxM_test(mm$resid_mul, batch)
#' @family multivariate-tests
#' @export
boxM_test <- function(resid_mul, batch){
  G <- length(resid_mul)
  g_names <- names(resid_mul)
  test_result <- lapply(1:G, function(g) boxM_simple(resid_mul[[g]], batch)$p.value)
  test_result_mat <- do.call(rbind, test_result) |> p.adjust(method = "bonferroni")
  test_result_df <- data.frame(cbind(g_names, test_result_mat))
  colnames(test_result_df) <- c("feature", "p.value")
  test_result_df <- test_result_df %>% mutate(p.value = round(as.numeric(p.value), 3)) %>% arrange(p.value)
  num_sig <- paste0(sprintf("%.2f", 100*mean(test_result_mat <= 0.05)), "%")
  return(list(test_table = test_result_df, perc.sig = num_sig))
}


# --- Covariance recovery metrics -------------------------------------------

#' Frobenius norm error between covariance matrices
#'
#' @param Sigma True covariance matrix.
#' @param Sigma_hat Estimated covariance matrix.
#'
#' @return Numeric scalar: \eqn{\| \Sigma - \hat\Sigma \|_F}.
#' @examples
#' # frobenius_norm(S, S_hat)
#' @family covariance-metrics
#' @export
frobenius_norm <- function(Sigma, Sigma_hat) {
  sqrt(sum((Sigma - Sigma_hat)^2))
}

#' Mean squared entrywise error (covariance matrices)
#' @inheritParams frobenius_norm
#' @return Numeric scalar: mean of squared entrywise errors.
#' @examples
#' # mse_cov(S, S_hat)
#' @family covariance-metrics
#' @export
mse_cov <- function(Sigma, Sigma_hat) {
  mean((Sigma - Sigma_hat)^2)
}

#' Spectral norm (operator 2-norm) of error
#' @inheritParams frobenius_norm
#' @details For symmetric matrices, equals the largest absolute eigenvalue
#'   of \code{Sigma - Sigma_hat}.
#' @return Numeric scalar.
#' @examples
#' # spectral_norm(S, S_hat)
#' @family covariance-metrics
#' @export
spectral_norm <- function(Sigma, Sigma_hat) {
  max(eigen(Sigma - Sigma_hat, symmetric = TRUE, only.values = TRUE)$values)
}


#' Mean absolute eigenvalue error
#' @inheritParams frobenius_norm
#' @return Numeric scalar: mean absolute difference between sorted eigenvalues.
#' @examples
#' # eigen_error(S, S_hat)
#' @family covariance-metrics
#' @export
eigen_error <- function(Sigma, Sigma_hat) {
  true_eigen <- eigen(Sigma, symmetric = TRUE)$values
  est_eigen <- eigen(Sigma_hat, symmetric = TRUE)$values
  mean(abs(true_eigen - est_eigen))
}

#' KL divergence between zero-mean Gaussians
#'
#' Computes \deqn{ \mathrm{KL}\big(N(0,\Sigma)\,\|\,N(0,\hat\Sigma)\big)
#'   = \tfrac12\{\mathrm{tr}(\hat\Sigma^{-1}\Sigma) - p + \log|\hat\Sigma|
#'   - \log|\Sigma|\}. }
#' Adds a small ridge to \eqn{\hat\Sigma} for stability.
#'
#' @inheritParams frobenius_norm
#' @return Numeric scalar.
#' @examples
#' # kl_divergence(S, S_hat)
#' @family covariance-metrics
#' @export
kl_divergence <- function(Sigma, Sigma_hat) {
  p <- ncol(Sigma)
  epsilon <- 1e-6
  inv_Sigma_hat <- solve(Sigma_hat + diag(epsilon, nrow(Sigma_hat)))
  0.5 * (sum(diag(inv_Sigma_hat %*% Sigma)) - p + log(det(Sigma_hat) - log(det(Sigma))))
}

#' Evaluate covariance recovery with multiple metrics
#' @inheritParams frobenius_norm
#' @return A named list with elements \code{Frobenius}, \code{MSE},
#'   \code{Spectral}, \code{EigenError}, \code{KL}.
#' @examples
#' # evaluate_cov_recovery(S, S_hat)
#' @family covariance-metrics
#' @export
evaluate_cov_recovery <- function(Sigma, Sigma_hat) {
  list(
    Frobenius = frobenius_norm(Sigma, Sigma_hat),
    MSE = mse_cov(Sigma, Sigma_hat),
    Spectral = spectral_norm(Sigma, Sigma_hat),
    EigenError = eigen_error(Sigma, Sigma_hat),
    KL = kl_divergence(Sigma, Sigma_hat)
  )
}

#' PCA prep for batch diagnostics (single or multiple measurements)
#'
#' Fits per-feature models (via [diag_model_gen()] + [diag_model_summary()]),
#' extracts **additive residuals**, and runs PCA to obtain score matrices for
#' visualization/diagnostics:
#'
#' - **Multivariate input** (`data` is a list of length *m*): run PCA within
#'   each measurement to get `F_list` (per-measurement scores), then concatenate
#'   scores and run a second PCA to get the shared score space `G`.
#' - **Univariate input** (`data` is a matrix/data.frame): run a single PCA on
#'   the residual matrix to get `F_t`.
#'
#' @param bat Either a factor of length *n* (univariate case), or a list of *m*
#'   factors (one per measurement) of length *n* (multivariate).
#' @param data Either an numeric matrix/data frame, or a list of *m*
#'   such matrices/data frames (matching `bat`).
#' @param covar Data frame of covariates. `NULL` for a null model,
#'   or a list of data frames if `data` is a list.
#' @param model Modeling function that accepts `formula` and `data`
#'   (e.g., [stats::lm()], [mgcv::gam()], [lme4::lmer()]).
#' @param formula RHS formula for covariates (do **not** include batch; it is
#'   added internally when batch adjustment is enabled).
#' @param ref.batch Optional reference batch level (forwarded to internal
#'   modeling, if applicable).
#' @param bat_adjust Logical; if `TRUE` (default) include batch dummies in the
#'   per-feature model via `model_fitting()` (same behavior as in harmonization).
#' @param ... Additional arguments forwarded to `model`.
#'
#' @return
#' - **If `data` is a list (multivariate)**: a list with
#'   \describe{
#'     \item{F_list}{List of length *m*; each element is an \eqn{n \times r_t} PCA score
#'       matrix from the within-measurement PCA on additive residuals.}
#'     \item{G}{\eqn{n \times r_concat} matrix of scores from a second PCA on the
#'       column-bound `F_list` (shared subject score space).}
#'     \item{bat}{The input list of batch factors (returned for convenience).}
#'   }
#' - **If `data` is a matrix/data.frame (univariate)**: a list with
#'   \describe{
#'     \item{F_t}{\eqn{n \times r} PCA score matrix from the residual PCA.}
#'     \item{bat}{The input batch factor (returned for convenience).}
#'   }
#'
#' @details
#' Residuals are the **additive** residuals produced by
#' [diag_model_summary()], i.e., the data after removing non-batch fixed effects
#' (and batch dummies if `bat_adjust = TRUE`).
#'
#' @examples
#' set.seed(1)
#' n <- 40; G <- 10
#' X <- matrix(rnorm(n*G), n, G); colnames(X) <- paste0("g",1:G)
#' b <- factor(rep(LETTERS[1:2], each = n/2))
#' # Univariate (single measurement)
#' pp1 <- pca_prep(b, X, covar = NULL, model = stats::lm, formula = y ~ 1)
#' str(pp1$F_t)
#'
#' # Multivariate (two measurements)
#' Y1 <- X + matrix(rnorm(n*G, sd=.3), n, G)
#' Y2 <- X + matrix(rnorm(n*G, sd=.5), n, G)
#' bats <- list(b, b)
#' pp2 <- pca_prep(bats, list(Y1, Y2), covar = list(NULL,NULL),
#'                 model = stats::lm, formula = y ~ 1)
#' lapply(pp2$F_list, dim); dim(pp2$G)
#'
#' @seealso [diag_model_gen()], [diag_model_summary()], [pca_plot()]
#' @export

pca_prep <- function(bat, data, covar, model, formula = NULL, ref.batch = NULL, bat_adjust = TRUE, ...){
  if(is.list(data)&& !is.data.frame(data) && !is.matrix(data)){
    m <- length(data)
    diag_model_list <- lapply(1:m, function(i) diag_model_gen(bat[[i]], data[[i]], covar[[i]], model, formula))
    resid_add_list <- lapply(1:m, function(i) diag_model_summary(diag_model_list[[i]])$resid_add)
    F_list <- lapply(1:m, function(i){
      F_t <- prcomp(resid_add_list[[i]], center = TRUE, scale. = TRUE)$x
      return(F_t)
    })
    F_concat <- do.call(cbind, F_list)
    pc_concat <- prcomp(F_concat, center = TRUE, scale. = FALSE)
    G <- pc_concat$x
    return(list(F_list = F_list, G = G, bat = bat))
  }else if(is.data.frame(data) || is.matrix(data)){
    diag_model_result <- diag_model_gen(bat, data, covar, model, formula)
    resid_add_df <- diag_model_summary(diag_model_result)$resid_add
    F_t <- prcomp(resid_add_df, center = TRUE, scale. = TRUE)$x
    return(list(F_t = F_t, bat = bat))
  }else{
    stop("`data` must be a data.frame, matrix, or a list of those.")
  }

}


#' PCA scatterplots colored by batch/site
#'
#' Convenience plotting for the objects returned by [pca_prep()].
#' Produces 2D PC scatterplots colored by `bat`, with optional 68% ellipses.
#'
#' - **When `pca_prep_result` came from multivariate input**:
#'   * `type = "within"` (default) plots per-measurement PC scores from `F_list`
#'     faceted by measurement.
#'   * `type = "shared"` plots the shared subject score space `G` (second-stage PCA).
#' - **When `pca_prep_result` came from univariate input**:
#'   plots the single `F_t` score matrix.
#'
#' @param pca_prep_result The list returned by [pca_prep()].
#' @param type For multivariate results, either `"within"` (per-measurement
#'   PCs) or `"shared"` (shared score space `G`). Ignored for univariate results.
#' @param pc_1,pc_2 Integers giving which principal components to plot on
#'   the x/y axes (default `1` and `2`).
#' @param ellipse Logical; draw 68% normal ellipses per batch (default `TRUE`).
#'
#' @return A `ggplot` object you can print or add layers to.
#'
#' @examples
#' \dontrun{
#' set.seed(2)
#' n <- 40; G <- 8
#' X <- matrix(rnorm(n*G), n, G); colnames(X) <- paste0("g",1:G)
#' b <- factor(rep(LETTERS[1:2], each=n/2))
#' ppu <- pca_prep(b, X, covar = NULL, model = stats::lm, formula = y ~ 1)
#' p <- pca_plot(ppu, pc_1 = 1, pc_2 = 2, ellipse = TRUE)
#' # print(p)
#'
#' Y1 <- X + matrix(rnorm(n*G, sd=.3), n, G)
#' Y2 <- X + matrix(rnorm(n*G, sd=.5), n, G)
#' bats <- list(b, b)
#' ppm <- pca_prep(bats, list(Y1, Y2), covar = list(NULL,NULL),
#'                 model = stats::lm, formula = y ~ 1)
#' p_within <- pca_plot(ppm, type = "within")
#' p_shared <- pca_plot(ppm, type = "shared")
#' }
#'
#' @seealso [pca_prep()], [ggplot2::ggplot()], [ggplot2::stat_ellipse()]
#' @importFrom ggplot2 scale_color_brewer facet_wrap element_text element_rect
#' @importFrom ggplot2 element_blank
#' @export
pca_plot <- function(pca_prep_result, type = "within", pc_1 = 1, pc_2 = 2, ellipse = TRUE){
  pc_pair <- paste0("PC", c(pc_1, pc_2))
  element_names <- names(pca_prep_result)
  if("F_list" %in% element_names){
    if(type == "within"){
      m <- length(pca_prep_result$F_list); p <- ncol(pca_prep_result$F_list[[1]])
      F_df_con <- lapply(1:m, function(i) pca_prep_result$F_list[[i]] %>% data.frame() %>%
                           mutate(measurement = paste0("M_", i), bat = pca_prep_result$bat[[i]])) %>%
        bind_rows() %>% dplyr::select(dplyr::all_of(pc_pair), measurement, bat)
      plt <- ggplot(
        F_df_con,
        aes(x = .data[[pc_pair[1]]], y = .data[[pc_pair[2]]], color = bat)
      ) +
        geom_point(alpha = 0.7, size = 1.6, stroke = 0) +
        facet_wrap(~ measurement, ncol = 3) +
        coord_equal() +
        labs(
          x = pc_pair[1],
          y = pc_pair[2],
          color = "Batch / Site"
        ) +
        scale_color_brewer(palette = "Set2") +   # or: scale_color_viridis_d()
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          plot.title.position = "plot",
          panel.spacing = unit(6, "pt"),
          axis.title = element_text(face = "bold")
        )
      if(ellipse){
        plt + stat_ellipse(level = 0.68, linewidth = 0.6, show.legend = FALSE)
      }else{
        plt
      }

    }else{
      G_df <- pca_prep_result$G %>% data.frame() %>% mutate(bat = pca_prep_result$bat[[1]])
      plt <- ggplot(
        G_df,
        aes(x = .data[[pc_pair[1]]], y = .data[[pc_pair[2]]], color = bat)
      ) +
        geom_point(alpha = 0.7, size = 1.6, stroke = 0) +
        coord_equal() +
        labs(
          x = pc_pair[1],
          y = pc_pair[2],
          color = "Batch / Site"
        ) +
        scale_color_brewer(palette = "Set2") +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          plot.title.position = "plot",
          panel.spacing = unit(6, "pt"),
          axis.title = element_text(face = "bold")
        )

      if(ellipse){
        plt + stat_ellipse(level = 0.68, linewidth = 0.6, show.legend = FALSE)
      }else{
        plt
      }
    }
  }else{
    F_df <- pca_prep_result$F_t %>% data.frame() %>% mutate(bat = pca_prep_result$bat)
    plt <- ggplot(
      F_df,
      aes(x = .data[[pc_pair[1]]], y = .data[[pc_pair[2]]], color = bat)
    ) +
      geom_point(alpha = 0.7, size = 1.6, stroke = 0) +
      stat_ellipse(level = 0.68, linewidth = 0.6, show.legend = FALSE) +
      coord_equal() +
      labs(
        x = pc_pair[1],
        y = pc_pair[2],
        color = "Batch / Site"
      ) +
      scale_color_brewer(palette = "Set2") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title.position = "plot",
        panel.spacing = unit(6, "pt"),
        axis.title = element_text(face = "bold")
      )

    if(ellipse){
      plt + stat_ellipse(level = 0.68, linewidth = 0.6, show.legend = FALSE)
    }else{
      plt
    }
  }

}

utils::globalVariables(c(".data", "p.value", "bat"))


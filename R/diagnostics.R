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

utils::globalVariables(c(".data", "p.value"))


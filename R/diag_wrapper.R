#' Univariate batch-effect diagnostics (percent significant)
#'
#' Fits a per-feature model (via \code{\link{diag_model_gen}}) and summarizes
#' univariate batch tests on the resulting residuals:
#' \itemize{
#'   \item Mean-shift: one-way ANOVA and Kruskal–Wallis on additive residuals.
#'   \item Variance-shift: Levene, Bartlett, Fligner–Killeen on multiplicative residuals.
#' }
#' For each test, the function reports the percentage of features with adjusted
#' p-value < 0.05 (the underlying test helpers use Bonferroni by default).
#'
#' @param bat Either a factor of length \eqn{n} giving batch labels (univariate
#'   case), or a list of \eqn{m} factors (one per measurement) for the multivariate
#'   input variant.
#' @param data Either a numeric \eqn{n \times G} matrix/data frame of responses
#'   (rows = samples, columns = features), or a list of \eqn{m} such matrices
#'   (one per measurement) matching \code{bat} when multivariate inputs are supplied.
#' @param covar Data frame of covariates (\eqn{n \times q}), \code{NULL} for a
#'   null model, or a list of data frames if \code{data} is a list.
#' @param model A modeling function that accepts \code{formula} and \code{data}
#'   (e.g., \code{stats::lm}, \code{mgcv::gam}, \code{lme4::lmer}).
#' @param formula RHS formula for covariates (do not include batch; it is handled
#'   internally when batch adjustment is enabled).
#' @param ref.batch Optional reference batch level (forwarded to internal modeling
#'   if applicable).
#' @param test_param List of extra arguments passed to
#'   \code{\link{diag_model_summary}} (e.g., \code{list(random = "id")} for mixed models).
#' @param ... Additional arguments forwarded to the modeling function.
#'
#' @return A \code{data.frame} with five columns containing the percentage (as a
#'   formatted character, e.g. \code{"12.50%"}) of features called significant for:
#'   ANOVA, Kruskal–Wallis, Levene, Bartlett, and Fligner–Killeen tests.
#'   If \code{data} is a list, each column contains one percentage per measurement.
#'
#' @seealso \code{\link{diag_model_gen}}, \code{\link{diag_model_summary}},
#'   \code{\link{anova_test}}, \code{\link{kruskal_test}},
#'   \code{\link{lv_test}}, \code{\link{bl_test}}, \code{\link{fk_test}}
#'
#' @examples
#' set.seed(1)
#' n <- 60; G <- 10
#' bat <- factor(rep(letters[1:3], each = n/3))
#' X   <- matrix(rnorm(n * G), n, G); colnames(X) <- paste0("g", 1:G)
#'
#' # Simple null model per feature (lm with intercept only)
#' uni_test(bat, data = X, covar = NULL, model = lm)

#'
#' @export

uni_test <- function(bat, data, covar, model, formula = NULL, ref.batch = NULL, test_param = list(random = NULL), ...){
  if(is.matrix(data) || is.data.frame(data)){
    m <- 1
    diag_model <- diag_model_gen(bat, data, covar, model = model, formula = formula, ...)
    diag_summary <- do.call(diag_model_summary, c(list(diag_model = diag_model), test_param))
    anova_result <- anova_test(diag_summary$resid_add, bat)$perc.sig
    kruskal_result <- kruskal_test(diag_summary$resid_add, bat)$perc.sig
    lv_result <- lv_test(diag_summary$resid_mul, bat)$perc.sig
    bl_result <- bl_test(diag_summary$resid_mul, bat)$perc.sig
    fk_result <- fk_test(diag_summary$resid_mul, bat)$perc.sig
  }else{
    m <- length(data)
    diag_model <- lapply(1:m, function(i) {
      diag_model_gen(bat[[i]], data[[i]], covar[[i]], model = model, formula = formula, ...)
    })
    diag_summary <- lapply(1:m, function(i) do.call(diag_model_summary, c(list(diag_model = diag_model[[i]]), test_param)))
    anova_result <- sapply(1:m, function(i) anova_test(diag_summary[[i]]$resid_add, bat[[i]])$perc.sig)
    kruskal_result <- sapply(1:m, function(i) kruskal_test(diag_summary[[i]]$resid_add, bat[[i]])$perc.sig)
    lv_result <- sapply(1:m, function(i) lv_test(diag_summary[[i]]$resid_mul, bat[[i]])$perc.sig)
    bl_result <- sapply(1:m, function(i) bl_test(diag_summary[[i]]$resid_mul, bat[[i]])$perc.sig)
    fk_result <- sapply(1:m, function(i) fk_test(diag_summary[[i]]$resid_mul, bat[[i]])$perc.sig)
  }
  result_table <- data.frame(cbind(anova_result, kruskal_result, lv_result, bl_result, fk_result))
  return(result_table)
}


#' Multivariate batch-effect diagnostics (percent significant)
#'
#' For multiple measurements per feature, fits per-feature models for each
#' measurement, reshapes residuals with \code{\link{multi_test_reshape}}, and
#' runs multivariate batch tests per feature block:
#' \itemize{
#'   \item \strong{MANOVA} on additive residual matrices.
#'   \item \strong{Box's M} on multiplicative residual matrices.
#' }
#' Reports the percentage of features with adjusted p-value < 0.05 (Bonferroni
#' by default).
#'
#' @param bat List of \eqn{m} factors (one per measurement) of length \eqn{n}
#'   giving batch labels.
#' @param data List of \eqn{m} numeric \eqn{n \times G} matrices/data frames
#'   (rows = samples, columns = features), one per measurement.
#' @param covar List of \eqn{m} data frames of covariates (or \code{NULL} entries).
#' @param model A modeling function that accepts \code{formula} and \code{data}
#'   (e.g., \code{stats::lm}, \code{mgcv::gam}, \code{lme4::lmer}).
#' @param formula RHS formula for covariates (do not include batch; added internally).
#' @param ref.batch Optional reference batch level (forwarded to internal modeling
#'   if applicable).
#' @param ... Additional arguments forwarded to the modeling function.
#'
#' @return A \code{data.frame} with two columns (characters formatted as percentages):
#'   \code{manova_result} and \code{boxM_result}, giving the percentage of features
#'   with adjusted p-value < 0.05 for MANOVA and Box's M, respectively.
#'
#' @seealso \code{\link{diag_model_gen}}, \code{\link{diag_model_summary}},
#'   \code{\link{multi_test_reshape}}, \code{\link{manova_test}},
#'   \code{\link{boxM_test}}
#'
#' @examples
#' set.seed(2)
#' n <- 60; G <- 8; m <- 2
#' bat <- replicate(m, factor(rep(letters[1:3], each = n/3)), simplify = FALSE)
#' data <- replicate(m, {
#'   X <- matrix(rnorm(n * G), n, G); colnames(X) <- paste0("g", 1:G); X
#' }, simplify = FALSE)
#'
#' mul_test(bat, data, covar = replicate(m, NULL, simplify = FALSE), model = lm)
#'
#' @export
mul_test <- function(bat, data, covar, model, formula = NULL, ref.batch = NULL, ...){
  m <- length(data)
  diag_model <- lapply(1:m, function(i) {
    diag_model_gen(bat[[i]], data[[i]], covar[[i]], model = model, formula = formula, ...)
  })
  diag_summary <- lapply(1:m, function(i) diag_model_summary(diag_model[[i]]))
  multi_result <- multi_test_reshape(diag_summary)
  resid_add <- multi_result$resid_add
  resid_mul <- multi_result$resid_mul
  boxM_result <- boxM_test(resid_mul, bat[[1]])$perc.sig
  manova_result <- manova_test(resid_add, bat[[1]])$perc.sig
  result_table <- data.frame(cbind(manova_result, boxM_result))
  return(result_table)
}

#' Heatmap of univariate diagnostics (percent significant)
#'
#' Builds a tile heatmap summarizing the percentage of significant features
#' per test across measurements.
#'
#' @param uni_result An \eqn{M \times 5} data frame or matrix where rows are
#'   measurements and columns are the univariate tests.
#' @param ms Optional character vector of length \eqn{M} with measurement labels.
#'   If \code{NULL} (default) uses \code{"M_1"}, \code{"M_2"}, …, \code{"M_M"}.
#' @param angle Numeric; x-axis label rotation in degrees (default \code{45}).
#' @param hjust Numeric; horizontal justification for x-axis labels (default \code{1}).
#' @param vjust Numeric or \code{NULL}; vertical justification for x-axis labels
#'   (passed to \code{ggplot2::element_text()}).
#'
#' @return A \code{ggplot} object (tile heatmap) that you can print or add layers to.
#'
#' @details
#' Univariate tests are displayed in the fixed order:
#' ANOVA, Kruskal–Wallis, Levene, Bartlett, Fligner–Killeen.
#'
#' @examples
#' set.seed(1)
#' M <- 4
#' mock <- data.frame(
#'   anova_result   = paste0(round(runif(M,  5, 40), 1), "%"),
#'   kruskal_result = paste0(round(runif(M, 10, 50), 1), "%"),
#'   lv_result      = paste0(round(runif(M,  3, 35), 1), "%"),
#'   bl_result      = paste0(round(runif(M,  2, 30), 1), "%"),
#'   fk_result      = paste0(round(runif(M,  1, 25), 1), "%")
#' )
#' p <- uni_plot(mock, ms = paste0("Modality_", seq_len(M)))
#' # print(p)
#'
#' @export
uni_plot <- function(uni_result, ms = NULL, angle = 45, hjust = 1, vjust = NULL){
  uni_result <- apply(uni_result, 2, function(x) as.numeric(sub("%", "", x))) %>% data.frame()
  if(is.null(ms)){ms <- paste0("M_", 1:nrow(uni_result))}
  test_long <- uni_result %>% mutate(measurement = ms) %>% pivot_longer(c(1:5), names_to = "test", values_to = "sig")
  test_long <- test_long %>%
    mutate(
      test = gsub("_result", "", test),
      # Force desired order
      test = factor(test, levels = c("anova", "kruskal", "lv", "bl", "fk"))
    )

  ggplot(test_long, aes(x = test, y = measurement, fill = sig)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(sig, 1), "%")), size = 2.8, color = "black") +
    scale_fill_gradientn(
      colours = c("#f0f9e8", "#bae4bc", "#7bccc4", "#0868ac"),
      trans = "sqrt"
    ) +
    labs(
      x = "Test Type",
      y = "Measurement",
      fill = "%"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

utils::globalVariables(c("test", "sig"))

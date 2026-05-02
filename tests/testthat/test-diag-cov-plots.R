# tests/testthat/test-diag-cov-plots.R

testthat::test_that("diag_plot_pooled_cov returns ggplot for correlation", {
  testthat::skip_if_not_installed("ggplot2")

  data <- list(
    matrix(rnorm(20), nrow = 5, ncol = 4),
    matrix(rnorm(20), nrow = 5, ncol = 4)
  )
  bat <- list(
    factor(c("a", "a", "b", "b", "a")),
    factor(c("a", "a", "b", "b", "a"))
  )
  covar <- list(
    data.frame(age = 1:5),
    data.frame(age = 1:5)
  )

  fake_get_target_residuals <- function(data, bat, covar,
                                        model = lm, formula = NULL,
                                        ref.batch = NULL) {
    list(
      resid_list = list(
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5)),
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5))
      ),
      batch_result = list(
        list(batch_index = list(a = c(1, 2, 5), b = c(3, 4)))
      )
    )
  }

  ns <- asNamespace("MultiComBat")
  old_fun <- get(".get_target_residuals", envir = ns)
  unlockBinding(".get_target_residuals", ns)
  assign(".get_target_residuals", fake_get_target_residuals, envir = ns)
  lockBinding(".get_target_residuals", ns)
  withr::defer({
    unlockBinding(".get_target_residuals", ns)
    assign(".get_target_residuals", old_fun, envir = ns)
    lockBinding(".get_target_residuals", ns)
  })

  p <- diag_plot_pooled_cov(
    data = data,
    bat = bat,
    covar = covar,
    formula = y ~ age,
    use_correlation = TRUE
  )

  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_match(p$labels$title, "Average cross-modality correlation")
  testthat::expect_equal(p$labels$x, "Modality")
  testthat::expect_equal(p$labels$y, "Modality")
})

testthat::test_that("diag_plot_pooled_cov returns ggplot for covariance", {
  testthat::skip_if_not_installed("ggplot2")

  data <- list(
    matrix(rnorm(20), nrow = 5, ncol = 4),
    matrix(rnorm(20), nrow = 5, ncol = 4)
  )
  bat <- list(
    factor(c("a", "a", "b", "b", "a")),
    factor(c("a", "a", "b", "b", "a"))
  )
  covar <- list(
    data.frame(age = 1:5),
    data.frame(age = 1:5)
  )

  fake_get_target_residuals <- function(data, bat, covar,
                                        model = lm, formula = NULL,
                                        ref.batch = NULL) {
    list(
      resid_list = list(
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5)),
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5))
      ),
      batch_result = list(
        list(batch_index = list(a = c(1, 2, 5), b = c(3, 4)))
      )
    )
  }

  ns <- asNamespace("MultiComBat")
  old_fun <- get(".get_target_residuals", envir = ns)
  unlockBinding(".get_target_residuals", ns)
  assign(".get_target_residuals", fake_get_target_residuals, envir = ns)
  lockBinding(".get_target_residuals", ns)
  withr::defer({
    unlockBinding(".get_target_residuals", ns)
    assign(".get_target_residuals", old_fun, envir = ns)
    lockBinding(".get_target_residuals", ns)
  })

  p <- diag_plot_pooled_cov(
    data = data,
    bat = bat,
    covar = covar,
    formula = y ~ age,
    use_correlation = FALSE
  )

  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_match(p$labels$title, "Average cross-modality covariance")
})

testthat::test_that("diag_plot_batch_cov returns ggplot for averaged features", {
  testthat::skip_if_not_installed("ggplot2")

  data <- list(
    matrix(rnorm(20), nrow = 5, ncol = 4),
    matrix(rnorm(20), nrow = 5, ncol = 4)
  )
  bat <- list(
    factor(c("a", "a", "b", "b", "a")),
    factor(c("a", "a", "b", "b", "a"))
  )
  covar <- list(
    data.frame(age = 1:5),
    data.frame(age = 1:5)
  )

  fake_get_target_residuals <- function(data, bat, covar,
                                        model = lm, formula = NULL,
                                        ref.batch = NULL) {
    list(
      resid_list = list(
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5)),
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5))
      ),
      batch_result = list(
        list(batch_index = list(a = c(1, 2, 5), b = c(3, 4)))
      )
    )
  }

  ns <- asNamespace("MultiComBat")
  old_fun <- get(".get_target_residuals", envir = ns)
  unlockBinding(".get_target_residuals", ns)
  assign(".get_target_residuals", fake_get_target_residuals, envir = ns)
  lockBinding(".get_target_residuals", ns)
  withr::defer({
    unlockBinding(".get_target_residuals", ns)
    assign(".get_target_residuals", old_fun, envir = ns)
    lockBinding(".get_target_residuals", ns)
  })

  p <- diag_plot_batch_cov(
    data = data,
    bat = bat,
    covar = covar,
    formula = y ~ age,
    use_correlation = TRUE,
    feature_idx = NULL
  )

  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_match(p$labels$title, "Average cross-modality correlation")
})

testthat::test_that("diag_plot_batch_cov returns ggplot for a single feature", {
  testthat::skip_if_not_installed("ggplot2")

  data <- list(
    matrix(rnorm(20), nrow = 5, ncol = 4),
    matrix(rnorm(20), nrow = 5, ncol = 4)
  )
  bat <- list(
    factor(c("a", "a", "b", "b", "a")),
    factor(c("a", "a", "b", "b", "a"))
  )
  covar <- list(
    data.frame(age = 1:5),
    data.frame(age = 1:5)
  )

  fake_get_target_residuals <- function(data, bat, covar,
                                        model = lm, formula = NULL,
                                        ref.batch = NULL) {
    list(
      resid_list = list(
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5)),
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5))
      ),
      batch_result = list(
        list(batch_index = list(a = c(1, 2, 5), b = c(3, 4)))
      )
    )
  }

  ns <- asNamespace("MultiComBat")
  old_fun <- get(".get_target_residuals", envir = ns)
  unlockBinding(".get_target_residuals", ns)
  assign(".get_target_residuals", fake_get_target_residuals, envir = ns)
  lockBinding(".get_target_residuals", ns)
  withr::defer({
    unlockBinding(".get_target_residuals", ns)
    assign(".get_target_residuals", old_fun, envir = ns)
    lockBinding(".get_target_residuals", ns)
  })

  p <- diag_plot_batch_cov(
    data = data,
    bat = bat,
    covar = covar,
    formula = y ~ age,
    use_correlation = FALSE,
    feature_idx = 2
  )

  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_match(p$labels$title, "Cross-modality covariance")
})

testthat::test_that("diag_plot_batch_cov includes pooled and batch panels", {
  testthat::skip_if_not_installed("ggplot2")

  data <- list(
    matrix(rnorm(20), nrow = 5, ncol = 4),
    matrix(rnorm(20), nrow = 5, ncol = 4)
  )
  bat <- list(
    factor(c("a", "a", "b", "b", "a")),
    factor(c("a", "a", "b", "b", "a"))
  )
  covar <- list(
    data.frame(age = 1:5),
    data.frame(age = 1:5)
  )

  fake_get_target_residuals <- function(data, bat, covar,
                                        model = lm, formula = NULL,
                                        ref.batch = NULL) {
    list(
      resid_list = list(
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5)),
        data.frame(V1 = rnorm(5), V2 = rnorm(5), V3 = rnorm(5), V4 = rnorm(5))
      ),
      batch_result = list(
        list(batch_index = list(a = c(1, 2, 5), b = c(3, 4)))
      )
    )
  }

  ns <- asNamespace("MultiComBat")
  old_fun <- get(".get_target_residuals", envir = ns)
  unlockBinding(".get_target_residuals", ns)
  assign(".get_target_residuals", fake_get_target_residuals, envir = ns)
  lockBinding(".get_target_residuals", ns)
  withr::defer({
    unlockBinding(".get_target_residuals", ns)
    assign(".get_target_residuals", old_fun, envir = ns)
    lockBinding(".get_target_residuals", ns)
  })

  p <- diag_plot_batch_cov(data = data, bat = bat, covar = covar, formula = y ~ age)

  plot_data <- ggplot2::ggplot_build(p)$plot$data

  testthat::expect_true("source" %in% names(plot_data))
  testthat::expect_true(any(grepl("Pooled", plot_data$source)))
  testthat::expect_true(any(grepl("Batch a", plot_data$source)))
  testthat::expect_true(any(grepl("Batch b", plot_data$source)))
})

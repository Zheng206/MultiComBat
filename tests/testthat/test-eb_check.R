test_that("eb_check.univariate returns tidy prior/empirical df and eb_plot works", {
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 90, p = 18, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    add_outlier = TRUE, outlier_size = 6,
    prior_type = "iw", seed = 1004
  )
  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]
  batches <- levels(bat)

  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm(
    bat    = bat,
    data   = Y,
    covar  = cov,
    model  = lm,
    formula = form,
    ref.batch = NULL,
    eb     = TRUE,
    stan   = FALSE,
    robust.LS = FALSE
  )


  chk <- eb_check(out$eb_result)
  expect_s3_class(chk, "univariate")
  expect_type(chk, "list")
  expect_true(is.data.frame(chk$eb_df))

  # columns and values
  expect_true(all(c("batch", "parm", "type", "dist") %in% names(chk$eb_df)))
  expect_true(all(unique(chk$eb_df$parm) %in% c("gamma", "delta")))
  expect_true(all(unique(chk$eb_df$type) %in% c("prior", "emp.dist")))
  expect_true(all(unique(chk$eb_df$batch) %in% batches))

  # plotting should return a ggplot object
  p <- eb_plot(chk)
  expect_s3_class(p, "ggplot")
})

test_that("eb_check.multivariate returns eb_gamma/eb_delta and eb_plot works (gamma)", {
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 90, p = 18, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    add_outlier = TRUE, outlier_size = 6,
    prior_type = "iw", seed = 1004
  )
  Y   <- sim$data
  bat <- sim$batch
  cov <- sim$covariates
  batches <- levels(bat[[1]])

  rhs  <- paste(colnames(cov[[1]]), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm.multivariate(
    bat    = bat,
    data   = Y,
    covar  = cov,
    model  = lm,
    formula = form,
    ref.batch = NULL,
    eb     = TRUE,
    stan   = FALSE,
    robust.LS = FALSE
  )


  chk <- eb_check(out$eb_result)
  expect_s3_class(chk, "multivariate")
  expect_type(chk, "list")
  expect_true(is.data.frame(chk$eb_gamma))
  expect_true(is.data.frame(chk$eb_delta))

  # columns and values
  expect_true(all(c("batch", "measurement", "type", "dist") %in% names(chk$eb_gamma)))
  expect_true(all(c("method", "value", "type", "batch") %in% names(chk$eb_delta)))
  expect_true(all(unique(chk$eb_gamma$type) %in% c("prior", "emp.dist")))
  expect_true(all(unique(chk$eb_delta$type) %in% c("prior", "emp.dist")))
  expect_true(all(unique(chk$eb_gamma$batch) %in% batches))
  expect_true(all(unique(chk$eb_delta$batch) %in% batches))

  # plotting should return a ggplot object
  p <- eb_plot(chk)
  expect_s3_class(p, "ggplot")
  p1 <- eb_plot(chk, param = "delta", bat = "a")
  expect_s3_class(p1, "ggplot")

})


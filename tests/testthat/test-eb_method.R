test_that("combat_type sets class correctly and returns the tag", {
  u <- combat_type("univariate")
  m <- combat_type("multivariate")

  expect_type(u, "character")
  expect_s3_class(u, "univariate")
  expect_identical(unclass(u), "univariate")

  expect_type(m, "character")
  expect_s3_class(m, "multivariate")
  expect_identical(unclass(m), "multivariate")
})

test_that("EB (univariate) runs and returns well-formed results", {
  set.seed(123)

  # Simulate one measurement -> univariate
  sim <- simulate_data.m(
    m = 3, n = 60, p = 20, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 123
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]

  # Build formula from available covariates
  rhs <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  bres <- batch_matrix(bat)
  fit  <- model_fitting(Y, batch = bres$batch_matrix, covar = cov, model = lm, formula = form)
  std  <- standardize_data(fit, bres, robust.LS = FALSE)

  type <- combat_type("univariate")
  eb   <- eb_algorithm(type, std, bres, eb = TRUE, robust.LS = FALSE)

  # Structure
  expect_true(is.list(eb))
  expect_true(all(c("gamma_star", "delta_star", "gamma_hat", "delta_hat") %in% names(eb)))

  # Dimensions
  B <- nlevels(bres$batch_vector)
  G <- ncol(Y)
  expect_equal(dim(eb$gamma_star), c(B, G))
  expect_equal(dim(eb$delta_star), c(B, G))
  expect_equal(rownames(eb$gamma_star), levels(bres$batch_vector))

  # Finite, positive where expected
  expect_false(any(!is.finite(eb$gamma_star)))
  expect_false(any(!is.finite(eb$delta_star)))
  expect_true(all(eb$delta_star > 0))

  var_hat  <- apply(eb$gamma_hat,  2, var)
  var_star <- apply(eb$gamma_star, 2, var)
  expect_lte(mean(var_star, na.rm = TRUE), mean(var_hat, na.rm = TRUE) + 1e-8)

  eb0 <- eb_algorithm(type, std, bres, eb = FALSE, robust.LS = FALSE)
  expect_equal(eb0$gamma_star, eb0$gamma_hat)
  expect_equal(eb0$delta_star, eb0$delta_hat)
})

test_that("EB (multivariate) runs and returns PD covariance by feature/batch", {
  set.seed(456)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 80, p = 15, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 456
  )

  bres_list <- vector("list", m_meas)
  fit_list  <- vector("list", m_meas)
  std_list  <- vector("list", m_meas)

  for (i in seq_len(m_meas)) {
    Y   <- as.matrix(sim$data[[i]])
    bat <- sim$batch[[i]]
    cov <- sim$covariates[[i]]

    rhs  <- paste(colnames(cov), collapse = " + ")
    form <- stats::as.formula(paste("y ~", rhs))

    bres_list[[i]] <- batch_matrix(bat)
    fit_list[[i]]  <- model_fitting(Y, batch = bres_list[[i]]$batch_matrix, covar = cov, model = stats::lm, formula = form)
    std_list[[i]]  <- standardize_data(fit_list[[i]], bres_list[[i]], robust.LS = FALSE)
  }

  type <- combat_type("multivariate")
  eb   <- eb_algorithm(type, std_list, bres_list, eb = TRUE, robust.LS = FALSE)

  # Structure sanity
  expect_true(is.list(eb$gamma_star))
  expect_true(is.list(eb$delta_star))
  batches <- levels(bres_list[[1]]$batch_vector)
  expect_equal(sort(names(eb$gamma_star)), sort(batches))
  expect_equal(sort(names(eb$delta_star)), sort(batches))

  # Dim checks and PD checks
  G <- ncol(sim$data[[1]])
  for (b in batches) {
    expect_equal(dim(eb$gamma_star[[b]]), c(m_meas, G))
    expect_false(any(!is.finite(eb$gamma_star[[b]])))

    expect_equal(length(eb$delta_star[[b]]), G)
    for (g in seq_len(G)) {
      S <- eb$delta_star[[b]][[g]]
      expect_equal(dim(S), c(m_meas, m_meas))
      expect_true(isTRUE(all(abs(S - t(S)) < 1e-8)))
      ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
      expect_true(all(is.finite(ev)))
      expect_true(min(ev) > 0)
    }
  }

  changed <- sum(abs(eb$gamma_star[[batches[1]]] - eb$gamma_hat[[batches[1]]])) > 0
  expect_true(changed)
})

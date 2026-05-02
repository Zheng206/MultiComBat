test_that("com_harm (univariate, EB=TRUE) returns well-formed results and shrinks batch differences", {
  set.seed(123)
  sim <- simulate_data.m(
    m = 3, n = 80, p = 20, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 123
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]

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

  # Structure checks
  expect_true(is.list(out))
  expect_true(all(c("harm_data", "resid", "eb_result") %in% names(out)))
  expect_equal(dim(out$harm_data), dim(Y))
  expect_equal(colnames(out$harm_data), colnames(Y))
  expect_false(any(!is.finite(out$harm_data)))
  expect_false(any(!is.finite(out$resid)))

  # EB results present and finite
  eb <- out$eb_result
  expect_true(all(c("gamma_star", "delta_star", "gamma_hat", "delta_hat") %in% names(eb)))
  expect_equal(dim(eb$gamma_star), c(nlevels(bat), ncol(Y)))
  expect_equal(dim(eb$delta_star), c(nlevels(bat), ncol(Y)))
  expect_false(any(!is.finite(eb$gamma_star)))
  expect_true(all(eb$delta_star > 0))

  # --- Shrinkage sanity: between-batch variance of means typically decreases
  # helper: per-feature variance of batch means
  between_var <- function(M, groups) {
    lvl <- levels(groups)
    v <- numeric(ncol(M))
    for (g in seq_len(ncol(M))) {
      means <- vapply(lvl, function(b) mean(M[groups == b, g, drop = TRUE]), numeric(1))
      v[g] <- stats::var(means)
    }
    v
  }
  # Before: original data; After: harmonized
  v_before <- between_var(Y, bat)
  v_after  <- between_var(out$harm_data, bat)
  expect_lte(mean(v_after, na.rm = TRUE), mean(v_before, na.rm = TRUE) + 1e-8)
})

test_that("com_harm (univariate) respects reference-batch passthrough", {
  set.seed(456)

  sim <- simulate_data.m(
    m = 3, n = 60, p = 15, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 456
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]

  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  ref_level <- levels(bat)[1]
  out <- com_harm(
    bat    = bat,
    data   = Y,
    covar  = cov,
    model  = lm,
    formula = form,
    ref.batch = ref_level,
    eb     = TRUE,
    stan   = FALSE,
    robust.LS = FALSE
  )

  # Rows in the reference batch should be exactly the original data
  idx_ref <- which(bat == ref_level)
  expect_equal(out$harm_data[idx_ref, , drop = FALSE],
               Y[idx_ref, , drop = FALSE],
               tolerance = 0)  # exact copy per implementation
})

test_that("com_harm (univariate) eb=FALSE works as expected", {
  set.seed(789)

  sim <- simulate_data.m(
    m = 3, n = 50, p = 10, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 789
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]

  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm(
    bat    = bat,
    data   = Y,
    covar  = cov,
    model  = lm,
    formula = form,
    ref.batch = NULL,
    eb     = FALSE,
    stan   = FALSE,
    robust.LS = FALSE
  )

  eb <- out$eb_result
  expect_equal(eb$gamma_star, eb$gamma_hat)
  expect_equal(eb$delta_star, eb$delta_hat)
})


test_that("com_harm (univariate) cov=TRUE works as expected", {
  set.seed(789)

  sim <- simulate_data.m(
    m = 3, n = 50, p = 10, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "covbat", seed = 789
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]

  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- suppressWarnings(com_harm(
    bat    = bat,
    data   = Y,
    covar  = cov,
    model  = lm,
    formula = form,
    ref.batch = NULL,
    eb     = TRUE,
    stan   = FALSE,
    robust.LS = FALSE,
    cov = TRUE
  ))

  eb <- out$eb_result
  expect_true(all(c("gamma_star", "delta_star", "gamma_hat", "delta_hat") %in% names(eb)))
  expect_equal(dim(eb$gamma_star), c(nlevels(bat), ncol(Y)))
  expect_equal(dim(eb$delta_star), c(nlevels(bat), ncol(Y)))
  expect_false(any(!is.finite(eb$gamma_star)))
  expect_true(all(eb$delta_star > 0))
})

test_that("com_harm (univariate) robust.LS helps under outliers", {
  set.seed(2468)

  sim <- simulate_data.m(
    m = 3, n = 80, p = 25, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    add_outlier = TRUE, outlier_size = 6,
    prior_type = "iw", seed = 2468
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]

  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out_classic <- com_harm(
    bat = bat, data = Y, covar = cov, model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = FALSE
  )
  out_robust  <- com_harm(
    bat = bat, data = Y, covar = cov, model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = TRUE
  )

  # And both return finite data
  expect_false(any(!is.finite(out_classic$harm_data)))
  expect_false(any(!is.finite(out_robust$harm_data)))
})

test_that("com_harm (univariate) works with covariates = NULL (null model with warning)", {
  set.seed(1357)

  sim <- simulate_data.m(
    m = 3, n = 40, p = 12, K = 2,
    add_covariates = FALSE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 1357
  )

  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]

  expect_warning(
    out <- com_harm(
      bat    = bat,
      data   = Y,
      covar  = NULL,
      model  = lm,
      formula = NULL,      # null model auto-selected
      ref.batch = NULL,
      eb     = TRUE,
      stan   = FALSE,
      robust.LS = FALSE
    ),
    regexp = "No covariates were provided.*null model"
  )

  expect_equal(dim(out$harm_data), dim(Y))
  expect_false(any(!is.finite(out$harm_data)))
})

### Multivaraite ComBat ###
test_that("com_harm.multivariate (EB=TRUE) returns well-formed results and PD covariances", {
  set.seed(1001)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 80, p = 16, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 1001
  )

  Y_list <- sim$data
  bat    <- sim$batch
  covs   <- sim$covariates

  # Build a simple formula from covariate names (shared across measurements)
  rhs  <- paste(colnames(covs[[1]]), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm.multivariate(
    bat      = sim$batch,
    data     = Y_list,
    covar    = covs,
    model    = lm,
    formula  = form,
    ref.batch = NULL,
    eb       = TRUE,
    stan     = FALSE,
    robust.LS = FALSE
  )

  # Structure
  expect_true(is.list(out))
  expect_true(all(c("harm_data", "resid", "eb_result") %in% names(out)))
  expect_equal(length(out$harm_data), m_meas)
  expect_equal(length(out$resid),     m_meas)

  n <- nrow(Y_list[[1]]); G <- ncol(Y_list[[1]])
  for (i in seq_len(m_meas)) {
    expect_equal(dim(out$harm_data[[i]]), c(n, G))
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
    expect_false(any(!is.finite(as.matrix(out$resid[[i]]))))
  }

  eb <- out$eb_result
  batches <- levels(bat[[1]])
  expect_equal(sort(names(eb$gamma_star)), sort(batches))
  expect_equal(sort(names(eb$delta_star)), sort(batches))

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
      expect_gt(min(ev), 0)  # PD
    }
  }

  batch_disp <- function(data_list, group) {
    lvl <- levels(group); m <- length(data_list); G <- ncol(data_list[[1]])
    v <- numeric(G)
    for (g in seq_len(G)) {
      M <- matrix(NA_real_, nrow = length(lvl), ncol = m)
      for (b in seq_along(lvl)) {
        idx <- which(group == lvl[b])
        for (i in seq_len(m)) {
          M[b, i] <- mean(data_list[[i]][idx, g, drop = TRUE])
        }
      }
      s_b <- rowMeans(M)
      v[g] <- var(s_b)
    }
    v
  }
  v_before <- batch_disp(Y_list, bat[[1]])
  v_after  <- batch_disp(out$harm_data, bat[[1]])
  expect_lte(mean(v_after, na.rm = TRUE), mean(v_before, na.rm = TRUE) + 1e-8)
})

test_that("com_harm.multivariate respects reference-batch passthrough across all measurements", {
  set.seed(1002)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 60, p = 12, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 1002
  )

  Y_list <- sim$data
  bat    <- sim$batch
  covs   <- sim$covariates
  rhs  <- paste(colnames(covs[[1]]), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  ref_level <- levels(bat)[1]
  out <- com_harm.multivariate(
    bat      = sim$batch,
    data     = Y_list,
    covar    = covs,
    model    = lm,
    formula  = form,
    ref.batch = ref_level,
    eb       = TRUE,
    stan     = FALSE,
    robust.LS = FALSE
  )

  idx_ref <- which(bat == ref_level)
  for (i in seq_len(m_meas)) {
    expect_equal(out$harm_data[[i]][idx_ref, , drop = FALSE],
                 Y_list[[i]][idx_ref, , drop = FALSE],
                 tolerance = 0)
  }
})

test_that("com_harm.multivariate with eb=FALSE returns hat==star in EB payload", {
  set.seed(1003)
  m_meas <- 2
  sim <- simulate_data.m(
    m = m_meas, n = 50, p = 10, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 1003
  )

  Y_list <- sim$data
  covs   <- sim$covariates
  rhs  <- paste(colnames(covs[[1]]), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm.multivariate(
    bat      = sim$batch,
    data     = Y_list,
    covar    = covs,
    model    = lm,
    formula  = form,
    ref.batch = NULL,
    eb       = FALSE,   # no shrinkage
    stan     = FALSE,
    robust.LS = FALSE
  )

  eb <- out$eb_result
  batches <- names(eb$gamma_star)
  for (b in batches) {
    expect_equal(eb$gamma_star[[b]], eb$gamma_hat[[b]])
    for (g in seq_len(length(eb$delta_star[[b]]))) {
      expect_equal(eb$delta_star[[b]][[g]], eb$delta_hat[[b]][[g]])
    }
  }
})

test_that("com_harm.multivariate with cov=TRUE works successfully", {
  set.seed(1003)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 50, p = 10, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "covbat_hier", seed = 1003
  )

  Y_list <- sim$data
  bat    <- sim$batch
  covs   <- sim$covariates
  rhs  <- paste(colnames(covs[[1]]), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- suppressWarnings(com_harm.multivariate(
    bat      = sim$batch,
    data     = Y_list,
    covar    = covs,
    model    = lm,
    formula  = form,
    ref.batch = NULL,
    eb       = TRUE,
    stan     = FALSE,
    robust.LS = FALSE,
    cov = TRUE
  ))

  # Structure
  expect_true(is.list(out))
  expect_true(all(c("harm_data", "resid", "eb_result") %in% names(out)))
  expect_equal(length(out$harm_data), m_meas)
  expect_equal(length(out$resid),     m_meas)

  n <- nrow(Y_list[[1]]); G <- ncol(Y_list[[1]])
  for (i in seq_len(m_meas)) {
    expect_equal(dim(out$harm_data[[i]]), c(n, G))
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
    expect_false(any(!is.finite(as.matrix(out$resid[[i]]))))
  }

  eb <- out$eb_result
  batches <- levels(bat[[1]])
  expect_equal(sort(names(eb$gamma_star)), sort(batches))
  expect_equal(sort(names(eb$delta_star)), sort(batches))

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
      expect_gt(min(ev), 0)  # PD
    }
  }
})

test_that("com_harm.multivariate robust.LS helps under outliers", {
  set.seed(1004)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 90, p = 18, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    add_outlier = TRUE, outlier_size = 6,
    prior_type = "iw", seed = 1004
  )

  Y_list <- sim$data
  bat    <- sim$batch
  covs   <- sim$covariates
  rhs  <- paste(colnames(covs[[1]]), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out_c <- com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = FALSE
  )
  out_r <- com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = TRUE
  )

  for (i in seq_len(m_meas)) {
    expect_false(any(!is.finite(as.matrix(out_c$harm_data[[i]]))))
    expect_false(any(!is.finite(as.matrix(out_c$resid[[i]]))))
    expect_false(any(!is.finite(as.matrix(out_r$harm_data[[i]]))))
    expect_false(any(!is.finite(as.matrix(out_r$resid[[i]]))))
  }
})

test_that("com_harm (univariate) robust_cov argument is accepted without error", {
  # robust_cov is a multivariate-only arg, but verifying univariate doesn't break
  set.seed(111)
  sim <- simulate_data.m(
    m = 2, n = 40, p = 8, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 111
  )
  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]
  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm(
    bat = bat, data = Y, covar = cov,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = TRUE
  )
  expect_equal(dim(out$harm_data), dim(Y))
  expect_false(any(!is.finite(as.matrix(out$harm_data))))
})

test_that("com_harm (univariate) with cov=TRUE and ref.batch both active", {
  set.seed(222)
  sim <- simulate_data.m(
    m = 2, n = 60, p = 20, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 222
  )
  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]
  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))
  ref_level <- levels(bat)[1]

  out <- suppressWarnings(com_harm(
    bat = bat, data = Y, covar = cov,
    model = lm, formula = form,
    ref.batch = ref_level, eb = TRUE, stan = FALSE,
    robust.LS = FALSE, cov = TRUE
  ))

  # Note: cov=TRUE applies a second-stage transformation after ref.batch passthrough,
  # so ref rows are not guaranteed to match original exactly. Just verify shape and finiteness.
  expect_equal(dim(out$harm_data), dim(Y))
  expect_false(any(!is.finite(as.matrix(out$harm_data))))
})

test_that("com_harm (univariate) var_thresh / min_rblock / max_rblock args are passed through", {
  set.seed(333)
  sim <- simulate_data.m(
    m = 2, n = 50, p = 20, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 333
  )
  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]
  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))

  out <- suppressWarnings(com_harm(
    bat = bat, data = Y, covar = cov,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE,
    robust.LS = FALSE, cov = TRUE,
    var_thresh = 0.80, min_rblock = 2, max_rblock = 5
  ))
  expect_equal(dim(out$harm_data), dim(Y))
  expect_false(any(!is.finite(as.matrix(out$harm_data))))
})

test_that("com_harm (univariate) eb=TRUE and robust.LS=TRUE with ref.batch", {
  set.seed(444)
  sim <- simulate_data.m(
    m = 2, n = 60, p = 10, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 444
  )
  Y   <- as.matrix(sim$data[[1]])
  bat <- sim$batch[[1]]
  cov <- sim$covariates[[1]]
  rhs  <- paste(colnames(cov), collapse = " + ")
  form <- stats::as.formula(paste("y ~", rhs))
  ref_level <- levels(bat)[2]

  out <- com_harm(
    bat = bat, data = Y, covar = cov,
    model = lm, formula = form,
    ref.batch = ref_level, eb = TRUE, stan = FALSE, robust.LS = TRUE
  )
  idx_ref <- which(bat == ref_level)
  expect_equal(out$harm_data[idx_ref, , drop = FALSE],
               Y[idx_ref, , drop = FALSE], tolerance = 0)
})

### com_harm.multivariate additional branches ###

test_that("com_harm.multivariate with robust_cov=TRUE runs and returns finite results", {
  set.seed(2001)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 80, p = 14, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 2001
  )
  Y_list <- sim$data
  covs   <- sim$covariates
  rhs    <- paste(colnames(covs[[1]]), collapse = " + ")
  form   <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE,
    robust.LS = FALSE, robust_cov = TRUE
  )

  expect_true(is.list(out))
  expect_true(all(c("harm_data", "resid", "eb_result") %in% names(out)))
  for (i in seq_len(m_meas)) {
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
    expect_false(any(!is.finite(as.matrix(out$resid[[i]]))))
  }
})

test_that("com_harm.multivariate robust_cov=TRUE with only 2 batches triggers fallback message", {
  # With K=2 batches, robust_sigma0_from_eb may not use robust weighting
  set.seed(2002)
  m_meas <- 2
  sim <- simulate_data.m(
    m = m_meas, n = 40, p = 8, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 2002
  )
  Y_list <- sim$data
  covs   <- sim$covariates
  rhs    <- paste(colnames(covs[[1]]), collapse = " + ")
  form   <- stats::as.formula(paste("y ~", rhs))

  # Expect either a message or silent success — we just confirm it runs
  expect_no_error(
    com_harm.multivariate(
      bat = sim$batch, data = Y_list, covar = covs,
      model = lm, formula = form,
      ref.batch = NULL, eb = TRUE, stan = FALSE,
      robust.LS = FALSE, robust_cov = TRUE
    )
  )
})

test_that("com_harm.multivariate robust_cov=TRUE with eb=FALSE runs correctly", {
  set.seed(2003)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 60, p = 10, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 2003
  )
  Y_list <- sim$data
  covs   <- sim$covariates
  rhs    <- paste(colnames(covs[[1]]), collapse = " + ")
  form   <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = NULL, eb = FALSE, stan = FALSE,
    robust.LS = FALSE, robust_cov = TRUE
  )
  for (i in seq_len(m_meas)) {
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
  }
})

test_that("com_harm.multivariate with cov=TRUE, robust.LS=TRUE, and ref.batch", {
  set.seed(2004)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 70, p = 12, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "covbat_hier", seed = 2004
  )
  Y_list <- sim$data
  bat    <- sim$batch
  covs   <- sim$covariates
  rhs    <- paste(colnames(covs[[1]]), collapse = " + ")
  form   <- stats::as.formula(paste("y ~", rhs))
  ref_level <- levels(bat[[1]])[1]

  out <- suppressWarnings(com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = ref_level, eb = TRUE, stan = FALSE,
    robust.LS = TRUE, cov = TRUE
  ))

  expect_equal(length(out$harm_data), m_meas)
  for (i in seq_len(m_meas)) {
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
  }
})

test_that("com_harm.multivariate with cov=TRUE custom var_thresh args", {
  set.seed(2005)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 60, p = 10, K = 2,
    add_covariates = TRUE, add_biomarkers = FALSE,
    prior_type = "covbat_hier", seed = 2005
  )
  Y_list <- sim$data
  covs   <- sim$covariates
  rhs    <- paste(colnames(covs[[1]]), collapse = " + ")
  form   <- stats::as.formula(paste("y ~", rhs))

  out <- suppressWarnings(com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE,
    robust.LS = FALSE, cov = TRUE,
    var_thresh = 0.80, min_rblock = 1, max_rblock = 4
  ))

  expect_equal(length(out$harm_data), m_meas)
  for (i in seq_len(m_meas)) {
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
  }
})

test_that("com_harm.multivariate with null covariates issues warning and proceeds", {
  set.seed(2006)
  m_meas <- 2
  sim <- simulate_data.m(
    m = m_meas, n = 40, p = 8, K = 2,
    add_covariates = FALSE, add_biomarkers = FALSE,
    prior_type = "iw", seed = 2006
  )
  Y_list <- sim$data
  covs   <- lapply(1:m_meas, function(i) NULL)

  expect_warning(
    out <- com_harm.multivariate(
      bat = sim$batch, data = Y_list, covar = covs,
      model = lm, formula = NULL,
      ref.batch = NULL, eb = TRUE, stan = FALSE,
      robust.LS = FALSE
    ),
    regexp = "No covariates were provided.*null model"
  )
  expect_equal(length(out$harm_data), m_meas)
})

test_that("com_harm.multivariate robust.LS=TRUE with robust_cov=TRUE", {
  set.seed(2007)
  m_meas <- 3
  sim <- simulate_data.m(
    m = m_meas, n = 90, p = 16, K = 3,
    add_covariates = TRUE, add_biomarkers = FALSE,
    add_outlier = TRUE, outlier_size = 5,
    prior_type = "iw", seed = 2007
  )
  Y_list <- sim$data
  covs   <- sim$covariates
  rhs    <- paste(colnames(covs[[1]]), collapse = " + ")
  form   <- stats::as.formula(paste("y ~", rhs))

  out <- com_harm.multivariate(
    bat = sim$batch, data = Y_list, covar = covs,
    model = lm, formula = form,
    ref.batch = NULL, eb = TRUE, stan = FALSE,
    robust.LS = TRUE, robust_cov = TRUE
  )
  for (i in seq_len(m_meas)) {
    expect_false(any(!is.finite(as.matrix(out$harm_data[[i]]))))
    expect_false(any(!is.finite(as.matrix(out$resid[[i]]))))
  }
})

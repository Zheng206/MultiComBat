test_that("biweight_midvar returns a finite numeric for typical data", {
  set.seed(1)
  x <- rnorm(200, sd = 2)
  v <- biweight_midvar(x)
  expect_type(v, "double")
  expect_length(v, 1L)
  expect_true(is.finite(v))
})

test_that("biweight_midvar is (approximately) translation invariant with default center", {
  set.seed(2)
  x  <- rnorm(300, sd = 1.5)
  y  <- x + 5
  v1 <- biweight_midvar(x)
  v2 <- biweight_midvar(y)
  expect_equal(v1, v2, tolerance = 1e-8)
})

test_that("biweight_midvar respects the 'center' argument", {
  set.seed(3)
  x <- rnorm(250, mean = 0.2, sd = 1)
  v_def   <- biweight_midvar(x)  # uses median(x)
  v_median <- biweight_midvar(x, center = median(x))
  v_zero   <- biweight_midvar(x, center = 0)
  expect_equal(v_def, v_median, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(v_def, v_zero)))
})

test_that("biweight_midvar is robust to a single large outlier compared to var()", {
  set.seed(4)
  x_base <- rnorm(80, sd = 1)
  x_out  <- c(x_base, 20)  # big outlier
  v_base <- var(x_base)
  v_classic <- var(x_out)
  v_robust  <- biweight_midvar(x_out)

  # Robust estimate should be much closer to baseline than classical variance
  expect_lt(abs(v_robust - v_base), abs(v_classic - v_base))
  # And typically smaller than the inflated classical variance with outlier
  expect_lt(v_robust, v_classic)
})

test_that("biweight_midvar responds to norm.unbiased flag", {
  set.seed(5)
  x <- rnorm(150)
  v_true  <- biweight_midvar(x, norm.unbiased = TRUE)
  v_false <- biweight_midvar(x, norm.unbiased = FALSE)
  # Usually different because c changes
  expect_false(isTRUE(all.equal(v_true, v_false)))
})

test_that("biweight_midvar returns NaN/Inf when MAD=0 (edge case)", {
  x <- rep(1, 50)
  v <- biweight_midvar(x)
  expect_true(is.nan(v) || is.infinite(v) || is.na(v))
})


############# Biweight_mid_val_mul #######################
test_that("biweight_midvar_mul returns symmetric PSD matrix with correct dims", {
  set.seed(10)
  m <- 6; n <- 80
  X <- matrix(rnorm(m * n), nrow = m, ncol = n)

  S <- biweight_midvar_mul(X)
  expect_equal(dim(S), c(m, m))
  # Symmetry
  expect_equal(S, t(S), tolerance = 1e-10)
  # PSD (allow tiny negative eigs due to numeric)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-8)
  # Finite
  expect_false(any(!is.finite(S)))
})

test_that("biweight_midvar_mul is close to classic covariance without outliers", {
  set.seed(11)
  m <- 5; n <- 200
  X <- matrix(rnorm(m * n), nrow = m, ncol = n)
  S_classic <- tcrossprod(X) / n
  S_robust  <- biweight_midvar_mul(X)
  # Should be reasonably close for Gaussian data
  diff <- sqrt(sum((S_classic - S_robust)^2))
  ref  <- sqrt(sum(S_classic^2))
  expect_lt(diff / ref, 0.15)  # within 15% Frobenius; adjust if needed
})

test_that("biweight_midvar_mul downweights a single extreme outlier column", {
  set.seed(12)
  m <- 5; n <- 150
  X_base <- matrix(rnorm(m * n, sd = 1), nrow = m, ncol = n)
  S_base <- tcrossprod(X_base) / n

  # Inject a strong multivariate outlier in the last column
  X_out <- X_base
  X_out[, n] <- X_out[, n] + 10

  S_classic_out <- tcrossprod(X_out) / n
  S_robust_out  <- biweight_midvar_mul(X_out)

  # Robust covariance should be closer to the baseline (no outlier) than classical
  d_classic <- sqrt(sum((S_classic_out - S_base)^2))
  d_robust  <- sqrt(sum((S_robust_out  - S_base)^2))
  expect_lt(d_robust, d_classic)

  # Also expect no NAs/Infs
  expect_false(any(!is.finite(S_robust_out)))
})

test_that("biweight_midvar_mul responds to norm.unbiased flag", {
  set.seed(13)
  m <- 4; n <- 120
  X <- matrix(rnorm(m * n), nrow = m)
  S1 <- biweight_midvar_mul(X, norm.unbiased = TRUE)
  S2 <- biweight_midvar_mul(X, norm.unbiased = FALSE)
  # Usually different because scaling constant c changes
  expect_false(isTRUE(all.equal(S1, S2)))
})

test_that("biweight_midvar_mul handles rows with tiny MADs (no division-by-zero)", {
  set.seed(14)
  m <- 5; n <- 100
  X <- matrix(rnorm(m * n, sd = 1), nrow = m)
  # Make one row almost constant (very small MAD, but not zero)
  X[1, ] <- rnorm(n, mean = 0, sd = 1e-6)

  S <- biweight_midvar_mul(X)
  expect_equal(dim(S), c(m, m))
  expect_false(any(!is.finite(S)))
})


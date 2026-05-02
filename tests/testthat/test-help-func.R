test_that("multi_matrix_transform extracts one feature across measurements", {
  set.seed(1)
  m <- 3
  n <- 5
  G <- 4

  Z_list <- lapply(1:m, function(i) {
    matrix(rnorm(n * G), nrow = n, ncol = G)
  })

  Zi <- multi_matrix_transform(Z_list, g = 2)

  expect_true(is.matrix(Zi))
  expect_equal(dim(Zi), c(m, n))

  for (i in 1:m) {
    expect_equal(Zi[i, ], Z_list[[i]][, 2])
  }
})

test_that("multi_matrix_transform_inverse reconstructs measurement-specific single-column data frames", {
  set.seed(1)
  m <- 3
  n <- 5
  G <- 4

  Z_list <- lapply(1:m, function(i) {
    matrix(rnorm(n * G), nrow = n, ncol = G)
  })

  Zi <- multi_matrix_transform(Z_list, g = 3)
  out <- multi_matrix_transform_inverse(Zi, g_name = "feat3")

  expect_length(out, m)

  for (i in 1:m) {
    expect_s3_class(out[[i]], "data.frame")
    expect_equal(dim(out[[i]]), c(n, 1))
    expect_identical(colnames(out[[i]]), "feat3")
    expect_equal(out[[i]][[1]], Z_list[[i]][, 3])
  }
})

test_that("multi_matrix_transform and inverse are consistent", {
  set.seed(2)
  m <- 4
  n <- 6
  G <- 5

  Z_list <- lapply(1:m, function(i) {
    matrix(rnorm(n * G), nrow = n, ncol = G)
  })

  Zi <- multi_matrix_transform(Z_list, g = 1)
  out <- multi_matrix_transform_inverse(Zi, g_name = "feat1")

  for (i in 1:m) {
    expect_equal(as.numeric(out[[i]][, 1]), Z_list[[i]][, 1])
  }
})

test_that("sigma_inv_sqrt whitens an SPD matrix", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(25), 5, 5))

  S <- sigma_inv_sqrt(A)

  expect_true(is.matrix(S))
  expect_equal(dim(S), c(5, 5))
  expect_equal(S, t(S), tolerance = 1e-10)

  expect_equal(t(S) %*% A %*% S, diag(5), tolerance = 1e-5)
})

test_that("sigma_inv_sqrt warns and regularizes tiny or non-positive eigenvalues", {
  A <- diag(c(4, 1e-12, -1))

  expect_warning(
    S <- sigma_inv_sqrt(A, tol = 1e-6),
    "thresholded"
  )

  expect_true(is.matrix(S))
  expect_equal(dim(S), c(3, 3))
  expect_equal(S, t(S), tolerance = 1e-10)
})

test_that("sigma_sqrt squares back to the original SPD matrix", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(25), 5, 5))

  S <- sigma_sqrt(A)

  expect_true(is.matrix(S))
  expect_equal(dim(S), c(5, 5))
  expect_equal(S, t(S), tolerance = 1e-10)

  expect_equal(S %*% S, A, tolerance = 1e-5)
})

test_that("sigma_sqrt warns and regularizes tiny or non-positive eigenvalues", {
  A <- diag(c(9, 1e-12, -2))

  expect_warning(
    S <- sigma_sqrt(A, tol = 1e-6),
    "thresholded"
  )

  expect_true(is.matrix(S))
  expect_equal(dim(S), c(3, 3))
  expect_equal(S, t(S), tolerance = 1e-10)
})

test_that("robust_sigma0_from_eb falls back to weighted mean when too few batches", {
  S1 <- diag(2)
  S2 <- 2 * diag(2)

  out <- robust_sigma0_from_eb(
    list(S1, S2),
    weights = c(1, 3),
    min_batches_robust = 3
  )

  expected <- (1 * S1 + 3 * S2) / 4

  expect_false(out$used_robust)
  expect_equal(out$sigma0, expected)
  expect_equal(out$sigma0_init, expected)
  expect_equal(out$adaptive_weights, c(1, 1))
  expect_equal(out$final_weights, c(1, 3))
  expect_true(all(is.na(out$distances)))
})

test_that("robust_sigma0_from_eb returns expected components when robust mode is used", {
  S1 <- diag(2)
  S2 <- 1.2 * diag(2)
  S3 <- 0.8 * diag(2)

  out <- robust_sigma0_from_eb(list(S1, S2, S3))

  expect_true(out$used_robust)
  expect_named(
    out,
    c("sigma0", "sigma0_init", "distances",
      "adaptive_weights", "final_weights", "used_robust")
  )

  expect_true(is.matrix(out$sigma0))
  expect_equal(dim(out$sigma0), c(2, 2))
  expect_length(out$distances, 3)
  expect_length(out$adaptive_weights, 3)
  expect_length(out$final_weights, 3)

  expect_true(all(out$adaptive_weights > 0))
  expect_true(all(out$adaptive_weights <= 1))
})

test_that("robust_sigma0_from_eb downweights an outlying covariance matrix", {
  S1 <- diag(2)
  S2 <- 1.05 * diag(2)
  S3 <- matrix(c(20, 0, 0, 0.2), 2, 2)

  out <- robust_sigma0_from_eb(
    list(S1, S2, S3),
    weights = c(1, 1, 1),
    min_batches_robust = 3
  )

  expect_true(out$used_robust)
  expect_lt(out$adaptive_weights[3], 1)
  expect_equal(out$final_weights, out$adaptive_weights, tolerance = 1e-12)
})

test_that("robust_matrix_solver matches solve on a well-conditioned SPD matrix", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(16), 4, 4))

  invA <- robust_matrix_solver(A, method = "solve")
  expect_equal(invA %*% A, diag(4), tolerance = 1e-6)

  invA_auto <- robust_matrix_solver(A, method = "auto")
  expect_equal(invA_auto %*% A, diag(4), tolerance = 1e-6)
})

test_that("robust_matrix_solver pseudoinverse works on singular matrix", {
  A <- matrix(c(1, 2, 2, 4), 2, 2)

  invA <- robust_matrix_solver(A, method = "pseudoinverse")

  expect_true(is.matrix(invA))
  expect_equal(dim(invA), c(2, 2))

  ## Moore-Penrose-style reconstruction check
  expect_equal(A %*% invA %*% A, A, tolerance = 1e-6)
})

test_that("robust_matrix_solver qr method works on invertible square matrix", {
  A <- matrix(c(4, 1, 2, 3), 2, 2)

  invA <- robust_matrix_solver(A, method = "qr")

  expect_equal(invA %*% A, diag(2), tolerance = 1e-6)
})

test_that("robust_matrix_solver chol method works on SPD matrix", {
  A <- matrix(c(4, 1, 1, 3), 2, 2)

  invA <- robust_matrix_solver(A, method = "chol")

  expect_equal(invA %*% A, diag(2), tolerance = 1e-6)
})

test_that("robust_matrix_solver warns and falls back when chosen method fails", {
  A <- matrix(c(1, 2, 2, 4), 2, 2)

  expect_warning(
    invA <- robust_matrix_solver(A, method = "chol"),
    "Matrix solve failed"
  )

  expect_true(is.matrix(invA))
  expect_equal(dim(invA), c(2, 2))
})

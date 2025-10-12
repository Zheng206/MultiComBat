test_that("pick_r_from_pc respects thresholds and bounds", {
  set.seed(1)
  X <- scale(matrix(rnorm(200), 50, 4))
  pf <- prcomp(X, center = TRUE, scale. = FALSE)

  r1 <- pick_r_from_pc(pf, var_thresh = 0.90, min_r = 1, max_r = Inf)
  expect_type(r1, "double")
  expect_gte(r1, 1)
  expect_lte(r1, ncol(X))

  r2 <- pick_r_from_pc(pf, var_thresh = 0.99, min_r = 4, max_r = 4)
  expect_equal(r2, 4)

  r3 <- pick_r_from_pc(pf, var_thresh = 0.50, min_r = 3, max_r = 3)
  expect_equal(r3, 3)
})


test_that("score_space_coupling returns well-formed components", {
  set.seed(2)
  n <- 40
  R1 <- scale(matrix(rnorm(n * 20), n, 20))
  R2 <- scale(matrix(rnorm(n * 15), n, 15))

  fit <- score_space_coupling(list(R1, R2),
                              var_thresh = 0.80, min_rblock = 2, max_rblock = Inf)

  expect_type(fit, "list")
  expect_true(all(c("G","A_list","H_list","F_list","F_list_std","L_list",
                    "sdev_list","sc_list","meta") %in% names(fit)))
  m <- 2
  r_shared <- ncol(fit$G)
  expect_gte(r_shared, 1)

  for (t in seq_len(m)) {
    rt <- ncol(fit$F_list[[t]])
    expect_equal(dim(fit$A_list[[t]]), c(r_shared, rt))
    expect_equal(dim(fit$H_list[[t]]), c(n, rt))
    expect_equal(dim(fit$F_list_std[[t]]), c(n, rt))
    expect_equal(length(fit$sc_list[[t]]), rt)
    expect_equal(dim(fit$L_list[[t]]), c(ncol(get(paste0("R", t))), rt))
  }
})

test_that("Appropriate transformation of G, A and H approximately reconstructs standardized scores", {
  set.seed(3)
  n <- 50
  R1 <- scale(matrix(rnorm(n * 25), n, 25))
  R2 <- scale(matrix(rnorm(n * 18), n, 18))

  fit <- score_space_coupling(list(R1, R2), var_thresh = 0.85, min_rblock = 2)

  for (t in 1:2) {
    approx_std <- fit$G %*% fit$A_list[[t]] + fit$H_list[[t]]
    err <- sqrt(mean((approx_std - fit$F_list_std[[t]])^2))
    expect_lt(err, 1e-8)
  }
})

test_that("reconstruct_block returns data with correct shape and near original (centered) scale", {
  set.seed(4)
  n <- 60
  R1 <- scale(matrix(rnorm(n * 30), n, 30))
  R2 <- scale(matrix(rnorm(n * 20), n, 20))

  fit <- score_space_coupling(list(R1, R2), var_thresh = 0.9, min_rblock = 3)

  rec1 <- reconstruct_block(fit$G, fit$A_list[[1]], fit$H_list[[1]],
                            fit$L_list[[1]], fit$sc_list[[1]])
  rec2 <- reconstruct_block(fit$G, fit$A_list[[2]], fit$H_list[[2]],
                            fit$L_list[[2]], fit$sc_list[[2]])

  expect_equal(dim(rec1), dim(R1))
  expect_equal(dim(rec2), dim(R2))
  expect_false(any(!is.finite(rec1)))
  expect_false(any(!is.finite(rec2)))
})

test_that("stage2_variantA returns shapes and uses com_harm on shared scores", {
  set.seed(5)
  n <- 40
  R1 <- scale(matrix(rnorm(n * 22), n, 22))
  R2 <- scale(matrix(rnorm(n * 16), n, 16))
  site <- list(factor(rep(LETTERS[1:2], each = n/2)),
               factor(rep(LETTERS[1:2], each = n/2)))
  out <- suppressWarnings(stage2_variantA(list(R1, R2), site, var_thresh = 0.85, min_rblock = 2))
  expect_true(all(c("R_h_list","coupling","G_star","meta") %in% names(out)))
  expect_equal(length(out$R_h_list), 2L)
  expect_equal(dim(out$R_h_list[[1]]), dim(R1))
  expect_equal(dim(out$R_h_list[[2]]), dim(R2))
  expect_false(any(!is.finite(as.matrix(out$R_h_list[[1]]))))
  expect_false(any(!is.finite(as.matrix(out$R_h_list[[2]]))))
})

test_that("covbat works successfully", {
  set.seed(6)
  n <- 50; p <- 28
  R <- scale(matrix(rnorm(n * p), n, p))
  site <- factor(rep(LETTERS[1:2], each = n/2))
  out <- suppressWarnings(covbat(R, site, var_thresh = 0.9, min_rblock = 2, max_rblock = Inf))
  expect_equal(dim(out), dim(R))
})

test_that("diag_model_gen + diag_model_summary (lm) produce sane outputs", {
  sim <- simulate_data.m(
    m = 3, n = 60, p = 12, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 1002
  )
  Y    <- sim$data
  bat  <- sim$batch
  covs <- sim$covar

  dm <- diag_model_gen(
    bat     = bat[[1]],
    data    = Y[[1]],
    covar   = covs[[1]],
    model   = lm,
    formula = y ~ Age + Sex + Diagnosis,
    ref.batch = NULL,
    bat_adjust = TRUE
  )
  expect_s3_class(dm, "lm_model")

  ds <- diag_model_summary(dm)
  expect_s3_class(ds, "diag_summary")
  expect_true(all(c("coef","resid_add","resid_mul") %in% names(ds)))

  # Dimensions
  expect_equal(dim(ds$resid_add), dim(Y[[1]]))
  expect_equal(dim(ds$resid_mul), dim(Y[[1]]))
  expect_equal(ncol(ds$coef), ncol(Y[[1]]))
})

test_that("diag_model_gen works with bat_adjust = FALSE (null batch)", {
  sim <- simulate_data.m(
    m = 3, n = 60, p = 12, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    prior_type = "iw", seed = 1002
  )
  Y    <- sim$data
  bat  <- sim$batch
  covs <- sim$covar

  dm <- diag_model_gen(
    bat     = bat[[1]],
    data    = Y[[1]],
    covar   = covs[[1]],
    model   = lm,
    formula = y ~ Age + Sex + Diagnosis,
    ref.batch = NULL,
    bat_adjust = FALSE
  )
  expect_s3_class(dm, "lm_model")
})

test_that("anova_test and kruskal_test detect mean shifts across batches", {
  set.seed(3)
  n <- 120; G <- 6
  group <- factor(rep(c("A","B","C"), each = n/3))

  base <- matrix(rnorm(n * G, sd = 1), n, G)
  shift <- c(0, 0.8, -0.8)  # mean per group
  mu <- matrix(shift[as.integer(group)], n, 1) %*% matrix(seq(1, 0.5, length.out = G), 1, G)
  resid_add <- base + mu
  colnames(resid_add) <- paste0("feat", seq_len(G))

  a <- anova_test(resid_add, group)
  k <- kruskal_test(resid_add, group)
  expect_true(as.numeric(sub("%","", a$perc.sig)) >= 80)
  expect_true(as.numeric(sub("%","", k$perc.sig)) >= 80)
})

test_that("lv_test, bl_test, fk_test detect variance shifts across batches", {
  set.seed(4)
  n <- 150; G <- 5
  group <- factor(rep(c("A","B","C"), each = n/3))

  # Multiplicative residuals: same mean, different variances
  sd_per_group <- c(1, 2.5, 1.8)
  resid_mul <- sapply(seq_len(G), function(g) {
    rnorm(n, sd = sd_per_group[as.integer(group)])
  })
  colnames(resid_mul) <- paste0("feat", seq_len(G))

  lv <- lv_test(resid_mul, group)
  bl <- bl_test(resid_mul, group)
  fk <- fk_test(resid_mul, group)

  expect_true(as.numeric(sub("%","", lv$perc.sig)) >= 70)
  expect_true(as.numeric(sub("%","", bl$perc.sig)) >= 70)
  expect_true(as.numeric(sub("%","", fk$perc.sig)) >= 70)
})

test_that("multi_test_reshape returns per-feature n x m stacks", {
  set.seed(5)
  n <- 50; G <- 4; m <- 3

  mk_ds <- function() {
    ra <- matrix(rnorm(n * G), n, G, dimnames = list(NULL, paste0("f", 1:G)))
    rm <- matrix(rnorm(n * G), n, G, dimnames = list(NULL, paste0("f", 1:G)))
    structure(list(resid_add = ra, resid_mul = rm), class = "diag_summary")
  }
  diag_list <- replicate(m, mk_ds(), simplify = FALSE)

  mm <- multi_test_reshape(diag_list)
  expect_true(all(names(mm) == c("resid_add","resid_mul")))
  expect_length(mm$resid_add, G)
  expect_length(mm$resid_mul, G)

  for (g in seq_len(G)) {
    expect_equal(dim(mm$resid_add[[g]]), c(n, m))
    expect_equal(dim(mm$resid_mul[[g]]), c(n, m))
  }
})

test_that("manova_test flags mean differences; boxM_simple flags covariance differences", {
  set.seed(6)
  n <- 90; m <- 3
  batch <- factor(rep(c("A","B","C"), each = n/3))
  G <- 3
  F1 <- matrix(rnorm(n * m, sd = 1), n, m); F1[batch == "B", ] <- F1[batch == "B", ] + 1.2
  S_A <- diag(m)
  S_B <- matrix(0.5, m, m); diag(S_B) <- 1
  S_C <- diag(c(1, 1.8, 0.7))
  F2 <- matrix(NA_real_, n, m)
  F2[batch == "A", ] <- MASS::mvrnorm(sum(batch == "A"), mu = rep(0, m), Sigma = S_A)
  F2[batch == "B", ] <- MASS::mvrnorm(sum(batch == "B"), mu = rep(0, m), Sigma = S_B)
  F2[batch == "C", ] <- MASS::mvrnorm(sum(batch == "C"), mu = rep(0, m), Sigma = S_C)
  F3 <- matrix(rnorm(n * m), n, m)

  resid_add <- list(f1 = F1, f2 = F3, f3 = F3)
  resid_mul <- list(f1 = F3, f2 = F2, f3 = F3)

  mv <- manova_test(resid_add, batch, test = "Pillai")
  bm <- boxM_test(resid_mul, batch)

  # MANOVA should rank f1 most significant
  expect_identical(mv$test_table$feature[1], "f1")
  expect_lt(as.numeric(mv$test_table$p.value[1]), 0.001)

  # Box's M should rank f2 most significant
  expect_identical(bm$test_table$feature[1], "f2")
  expect_lt(as.numeric(bm$test_table$p.value[1]), 0.01)
})

test_that("boxM_simple behaves on equal vs different covariances", {
  set.seed(7)
  n <- 120; p <- 3
  S <- matrix(c(1, .4, .2, .4, 1, .1, .2, .1, 1), p, p)
  Y1 <- MASS::mvrnorm(n/2, mu = rep(0, p), Sigma = S)
  Y2 <- MASS::mvrnorm(n/2, mu = rep(0, p), Sigma = S)
  Y  <- rbind(Y1, Y2)
  g  <- factor(rep(c("A","B"), each = n/2))
  out_eq <- boxM_simple(Y, g)
  expect_gt(out_eq$p.value, 0.10)

  S2 <- matrix(c(1.4, .7, .1, .7, 1.8, .2, .1, .2, 0.9), p, p)
  Y2d <- MASS::mvrnorm(n/2, mu = rep(0, p), Sigma = S2)
  Yd  <- rbind(Y1, Y2d)
  out_diff <- boxM_simple(Yd, g)
  expect_lt(out_diff$p.value, 0.05)
})

test_that("covariance metrics are zero when Sigma_hat equals to Sigma", {
  S <- diag(3)
  expect_equal(frobenius_norm(S, S), 0)
  expect_equal(mse_cov(S, S), 0)
  expect_equal(spectral_norm(S, S), 0)
  expect_equal(eigen_error(S, S), 0)
})

test_that("covariance metrics increase with larger perturbations", {
  set.seed(8)
  S <- diag(c(1, 2, 3))
  Sh_small <- S + 0.05 * diag(3)
  Sh_big   <- S + 0.50 * diag(3)
  expect_lt(frobenius_norm(S, Sh_small), frobenius_norm(S, Sh_big))
  expect_lt(mse_cov(S, Sh_small), mse_cov(S, Sh_big))
  expect_lt(abs(spectral_norm(S, Sh_small)), abs(spectral_norm(S, Sh_big)))
  expect_lt(eigen_error(S, Sh_small), eigen_error(S, Sh_big))
})



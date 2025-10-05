test_that("uni_test and mul_test detect mean and variance shifts across batches", {
  set.seed(1101)
  sim <- simulate_data.m(
    m = 3, n = 90, p = 18, K = 3,
    add_covariates = TRUE, add_biomarkers = TRUE,
    add_outlier = TRUE, outlier_size = 6,
    prior_type = "iw", seed = 1004
  )
  out <- uni_test(bat = sim$batch, data = sim$data, covar = sim$covariates, model = lm, formula = y ~ Age + Sex + Diagnosis)
  mul_out <- mul_test(bat = sim$batch, data = sim$data, covar = sim$covariates, model = lm, formula = y ~ Age + Sex + Diagnosis)
  # Structure: 1x5 with the % strings
  expect_true(is.data.frame(out))
  expect_true(is.data.frame(mul_out))
  expect_equal(ncol(out), 5)
  expect_equal(ncol(mul_out), 2)
  expect_setequal(
    colnames(out),
    c("anova_result","kruskal_result","lv_result","bl_result","fk_result")
  )
  expect_setequal(
    colnames(mul_out),
    c("manova_result","boxM_result")
  )
})

test_that("uni_plot returns a ggplot object and accepts custom labels", {
  M <- 4
  mock <- data.frame(
    anova_result   = paste0(round(runif(M,  5, 40), 1), "%"),
    kruskal_result = paste0(round(runif(M, 10, 50), 1), "%"),
    lv_result      = paste0(round(runif(M,  3, 35), 1), "%"),
    bl_result      = paste0(round(runif(M,  2, 30), 1), "%"),
    fk_result      = paste0(round(runif(M,  1, 25), 1), "%")
  )

  p <- uni_plot(mock, ms = paste0("Mod_", seq_len(M)), angle = 30, hjust = 0.9, vjust = 1)
  expect_s3_class(p, "ggplot")
})




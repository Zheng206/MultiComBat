skip_if_no_stan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    testthat::skip("cmdstanr not installed")
  }
  ver <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  if (inherits(ver, "try-error") || is.na(ver)) {
    testthat::skip("CmdStan not installed (run cmdstanr::install_cmdstan())")
  }
}

skip_if_no_stan()
test_that("get_stan_model errors cleanly when model file is missing", {
  expect_error(
    get_stan_model("__this_model_does_not_exist__"),
    "Stan model not found"
  )
})

test_that("get_stan_model can be called with a temp compile_dir", {
  stan_path <- system.file("stan", "univariate_model.stan", package = "MultiComBat")
  if (!nzchar(stan_path)) skip("No installed stan file to point at (univariate_model.stan).")
  tmp <- tempfile("mc-cache-")
  dir.create(tmp)
  mdl <- get_stan_model("univariate_model", compile_dir = tmp, force_recompile = TRUE)
  expect_s3_class(mdl, "CmdStanModel")
})

test_that("data_reshape builds the correct long table", {
  set.seed(1)
  n <- 4; G <- 3; M <- 2

  # Minimal standardized data_stand_result list
  dsr <- lapply(seq_len(M), function(i) {
    list(data_stand = matrix(rnorm(n * G), n, G,
                             dimnames = list(NULL, paste0("g", 1:G))))
  })

  # Two measurement batch vectors
  b1 <- factor(rep(c("A","B"), each = n/2))
  b2 <- factor(rep(c("A","B"), each = n/2))
  y <- data_reshape(dsr, list(b1, b2))

  expect_equal(nrow(y), n * G * M)
  expect_true(all(c("value","measurement","feature","observation",
                    "batch","batch_id","feature_id") %in% names(y)))
  expect_equal(max(y$batch_id), nlevels(b1))
  expect_equal(max(y$feature_id), G)
})

test_that("stan_data_prep.univariate produces vectors/indices", {
  set.seed(2)
  n <- 5; G <- 3
  dsr <- list(data_stand = matrix(rnorm(n * G), n, G))
  bat <- factor(rep(c("A","B"), length.out = n))

  sd <- stan_data_prep(combat_type("univariate"), dsr, bat)

  expect_named(sd, c("N","G","I","y","g","i","nu_psi","nu_tau","nu_delta","A_psi","A_tau","A_delta"))
  expect_equal(sd$N, n * G)
  expect_equal(sd$G, G)
  expect_equal(sd$I, nlevels(bat))
  expect_length(sd$y, n * G)
  expect_length(sd$g, n * G)
  expect_length(sd$i, n * G)
})

test_that("stan_data_prep.multivariate produces y[N,G,M] and indices", {
  set.seed(3)
  n <- 6; G <- 2; M <- 2
  dsr <- lapply(seq_len(M), function(i) list(data_stand = matrix(rnorm(n*G), n, G)))
  b1  <- factor(rep(c("A","B"), each = n/2))
  b2  <- factor(rep(c("A","B"), each = n/2))
  sd  <- stan_data_prep(combat_type("multivariate"), dsr, list(b1, b2))

  expect_named(sd, c("N","G","M","I","i_id","y"))
  expect_equal(sd$N, n)
  expect_equal(sd$G, G)
  expect_equal(sd$M, M)
  expect_equal(sd$I, nlevels(b1))
  expect_equal(length(sd$i_id), n)
  expect_equal(dim(sd$y), c(n, G, M))
})

test_that("stan_algorithm.univariate parses summaries into matrices (mocked)", {
  fake_fit <- list(
    summary = function(variable, ...) {
      if (identical(variable, "mu_ig")) {
        # feature, batch
        data.frame(
          variable = c("mu_ig[1,1]", "mu_ig[2,1]", "mu_ig[1,2]", "mu_ig[2,2]"),
          mean     = c(0.1, 0.2, 0.3, 0.4)
        )
      } else if (identical(variable, "sigma_delta")) {
        data.frame(
          variable = c("sigma_delta[1,1]", "sigma_delta[2,1]",
                       "sigma_delta[1,2]", "sigma_delta[2,2]"),
          mean     = c(1.0, 1.1, 1.2, 1.3)  # will be squared by code
        )
      } else {
        stop("Unexpected variable in fake summary: ", variable)
      }
    }
  )

  fake_model <- structure(list(
    compile = function(dir, force_recompile = FALSE) invisible(NULL),
    sample  = function(...) fake_fit
  ), class = "CmdStanModel")

  # Mock our package's get_stan_model() to return the fake model:
  testthat::local_mocked_bindings(
    get_stan_model = function(name, compile_dir = tempfile(), force_recompile = FALSE) fake_model,
    .package = "MultiComBat"
  )

  # Minimal stan_data (values unused by mocked fit) and names
  sd <- list(N = 10, G = 2, I = 2, y = numeric(), g = integer(), i = integer())
  res <- stan_algorithm(combat_type("univariate"), sd,
                        chains = 1, parallel_chains = 1,
                        batch_names = c("A","B"))

  expect_true(is.matrix(res$gamma_star))
  expect_true(is.matrix(res$delta_star))
  expect_equal(dim(res$gamma_star), c(2, 2))
  expect_equal(dim(res$delta_star), c(2, 2))
  expect_identical(rownames(res$gamma_star), c("A","B"))
  expect_identical(rownames(res$delta_star), c("A","B"))

  # Check numeric content (delta_star should be mean^2 from the fake)
  expect_equal(sort(as.numeric(res$gamma_star)), c(0.1,0.2,0.3,0.4), tolerance = 1e-8)
  expect_equal(sort(as.numeric(res$delta_star)),
               c(1.0,1.1,1.2,1.3)^2, tolerance = 1e-8)
})

test_that("stan_algorithm.multivariate parses summaries into lists of matrices (mocked)", {
  I <- 2; G <- 2; M <- 2; N <- 5

  mu_df <- expand.grid(batch = seq_len(I),
                       feature = seq_len(G),
                       measurement = seq_len(M))
  sig_df <- expand.grid(batch = seq_len(I),
                        feature = seq_len(G),
                        m1 = seq_len(M), m2 = seq_len(M))

  fake_fit <- list(
    summary = function(variable, ...) {
      if (identical(variable, "mu_ig")) {
        data.frame(
          variable = sprintf("mu_ig[%d,%d,%d]",
                             mu_df$batch, mu_df$feature, mu_df$measurement),
          mean     = seq_len(nrow(mu_df)) / 10
        )
      } else if (identical(variable, "Sigma_ig")) {
        data.frame(
          variable = sprintf("Sigma_ig[%d,%d,%d,%d]",
                             sig_df$batch, sig_df$feature, sig_df$m1, sig_df$m2),
          mean     = 0.5
        )
      } else {
        stop("Unexpected variable in fake summary: ", variable)
      }
    }
  )

  fake_model <- structure(list(
    compile = function(dir, force_recompile = FALSE) invisible(NULL),
    sample  = function(...) fake_fit
  ), class = "CmdStanModel")

  testthat::local_mocked_bindings(
    get_stan_model = function(name, compile_dir = tempfile(), force_recompile = FALSE) fake_model,
    .package = "MultiComBat"
  )

  sd <- list(
    N = N, G = G, M = M, I = I,
    i_id = c(1L, 1L, 2L, 2L, 2L),
    y = array(0, dim = c(N, G, M))
  )

  out <- stan_algorithm(combat_type("multivariate"), sd,
                        chains = 1, parallel_chains = 1,
                        batch_names = c("A","B"))

  expect_true(is.list(out$gamma_star))
  expect_identical(names(out$gamma_star), c("A","B"))
  expect_equal(dim(out$gamma_star[[1]]), c(M, G))

  expect_true(is.list(out$delta_star))
  expect_identical(names(out$delta_star), c("A","B"))
  expect_true(is.list(out$delta_star[[1]]))
  expect_length(out$delta_star[[1]], G)
  expect_equal(dim(out$delta_star[[1]][[1]]), c(M, M))
})

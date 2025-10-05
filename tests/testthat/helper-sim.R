simulate_data.m <- function(
    m = 8, n = 100, p = 80, K = 2,
    sd_p = 2, sd_m = 3, sd_b = 2, sigma = 1,
    bt_type = "AR", w_type = "AR",
    bt_params = list(rho = 0.9, scale = 5),
    w_params = list(rho = 0.9, scale = 5),
    batch_cov = 0.5,
    seed = 123,
    global_batch_strength = 30,
    add_covariates = TRUE,
    add_biomarkers = TRUE,
    biomarker_prop = 0.6,
    add_outlier = FALSE,
    outlier_size = 6,
    add_noise = FALSE,
    prior_type = "iw",
    covbat_r = 3,
    covbat_rotate = c("none","householder","random"),
    covbat_severity = c("mild","moderate","severe"),
    covbat_jitter_sd = 0.02,
    covbat_nu0 = NULL,
    ...) {

  set.seed(seed)
  covbat_rotate   <- match.arg(covbat_rotate)
  covbat_severity <- match.arg(covbat_severity)
  if (is.null(covbat_nu0)) covbat_nu0 <- m + 6

  # 1. Generate Base Data Structure ------------------------------------------
  D_p <- diag(sd_p, p)
  D_m <- diag(sd_m, m)
  Sigma_p <- generate_cov.AR(m = p, rho = 0.7)
  Sigma_p <- generate_cov.ind(m = p)
  Sigma_m <- generate_cov.random_pd(m = m)
  Sigma_m <- D_m %*% Sigma_m %*% D_m
  mu <- rep(0, m * p)
  Sigma_total <- kronecker(Sigma_m, Sigma_p)
  simulated_data <- MASS::mvrnorm(n, mu, Sigma_total)
  individual_noise <- matrix(rnorm(n * m * p, mean = 0, sd = sigma), nrow = n)
  simulated_data <- simulated_data + individual_noise
  data_array <- array(simulated_data, dim = c(n, m, p))
  data_list <- lapply(1:m, function(i) {
    mat <- data_array[, i, ]
    colnames(mat) <- paste0("Feature_", 1:p)
    return(mat)
  })

  # 2. Add Covariates --------------------------------------------
  covariates <- NULL
  beta_age <- NULL
  beta_sex <- NULL
  beta_diag <- NULL

  if (add_covariates) {
    Age <- rnorm(n, mean = 50, sd = 10)
    Sex <- rbinom(n, size = 1, prob = 0.5)
    covariates <- data.frame(Age = Age, Sex = Sex)

    # Add covariate effects while preserving correlation structure
    beta_age <- matrix(rnorm(m * p, mean = 2, sd = 0.3), m, p)
    beta_sex <- matrix(rnorm(m * p, mean = 1.8, sd = 0.3), m, p)

    for (i in 1:m) {
      for (j in 1:p) {
        data_list[[i]][, j] <- data_list[[i]][, j] +
          beta_age[i, j] * Age +
          beta_sex[i, j] * Sex
      }
    }
  }

  # 3. Add Biomarkers --------------------------------------------
  Diagnosis <- NULL
  biomarkers <- NULL
  if (add_biomarkers) {
    Diagnosis <- rbinom(n, size = 1, prob = 0.3)
    n_biomarkers <- round(biomarker_prop * p)
    biomarkers <- sample(1:p, n_biomarkers)

    # Add biological effect
    bio_effect <- lapply(1:m, function(i) rnorm(n_biomarkers, mean = ifelse(i %in% c(1, 2), 0.05, 0.9), sd = 0.1))
    beta_diag <- sapply(1:m, function(i) {
      beta_diag_seq <- rep(0, p)
      for (j in seq_along(biomarkers)) {
        col <- biomarkers[j]
        beta_diag_seq[col] <- bio_effect[[i]][j]
      }
      return(beta_diag_seq)
    }) |> t()

    for (i in 1:m) {
      for (j in seq_along(biomarkers)) {
        col <- biomarkers[j]
        data_list[[i]][, col] <- data_list[[i]][, col] + bio_effect[[i]][j] * Diagnosis
      }
    }

    if (!is.null(covariates)) {
      covariates$Diagnosis <- Diagnosis
    }else{
      covariates <- lapply(1:m, function(i) data.frame(Diagnosis = Diagnosis))
    }
  }

  true_data_list <- data_list

  # 4. Add Batch Effects ---------------------------------------------------
  batch <- sample(letters[1:K], n, replace = TRUE)
  batch <- factor(batch, levels = letters[1:K])

  true_params <- list(
    gamma_i = lapply(1:K, function(i) rnorm(m, mean = 0, sd = global_batch_strength)),
    T_i = lapply(1:K, function(i) {
      R <- do.call(generate_cov, c(list(type = bt_type, m = m), bt_params))
      D <- diag(rnorm(m, mean = 1, sd = sd_b))
      D %*% R %*% D + batch_cov * i
    }),
    Psi_i = lapply(1:K, function(i) {
      R <- do.call(generate_cov, c(list(type = w_type, m = m), w_params))
      D <- diag(rnorm(m, mean = 1, sd = sd_b))
      D %*% R %*% D + batch_cov * i
    }),
    nu_i = rep(m + 2, K)
  )

  gamma_ig <- lapply(1:K, function(k) {
    t(MASS::mvrnorm(p, true_params$gamma_i[[k]], true_params$T_i[[k]]))
  })

  if (prior_type == "covbat") {
    covbat_obj <- make_covbat_sigmas(
      Sigma0 = Sigma_m, K = K, r = covbat_r,
      rotate = covbat_rotate, severity = covbat_severity,
      jitter_sd = covbat_jitter_sd, seed = seed + 999
    )
    Sigma_ig <- lapply(1:K, function(k) {
      base_Sk <- covbat_obj$Sigmas_k[[k]]
      replicate(p, covbat_obj$jitter_fun(base_Sk, sd = covbat_jitter_sd), simplify = FALSE)
    })

  } else if (prior_type == "covbat_hier") {
    covbat_obj <- make_covbat_sigmas(
      Sigma0 = Sigma_m, K = K, r = covbat_r,
      rotate = covbat_rotate, severity = covbat_severity,
      jitter_sd = covbat_jitter_sd, seed = seed + 999
    )
    Sigmas_k <- covbat_obj$Sigmas_k
    Sigma_ig <- lapply(1:K, function(k) {
      Sk   <- Sigmas_k[[k]]
      PsiK <- (covbat_nu0 - m - 1) * Sk
      replicate(p, MCMCpack::riwish(covbat_nu0, PsiK), simplify = FALSE)
    })

  } else if(prior_type == "iw"){
    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, MCMCpack::riwish(true_params$nu_i[i], true_params$Psi_i[[i]]),
                simplify = FALSE)
    })
  }else if(prior_type == "fa"){
    cov_mix <- c(iw = 0, lkj_half_t = 0, fa = 1, ar1_scales = 0, cs_scales = 0)
    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, {
        dist_i <- sample(names(cov_mix), size = 1, prob = cov_mix)
        draw_cov(
          m,
          dist = dist_i,
          iw_df   = true_params$nu_i[i],
          iw_scale= true_params$Psi_i[[i]],
          lkj_eta = 2, half_t_df = 4, half_t_scale = 2,
          fa_rank = min(2, m), ar1_rho = 0.6, cs_rho = 0.3
        )
      }, simplify = FALSE)
    })
  }else if (prior_type == "AR"){
    cov_mix <- c(iw = 0, lkj_half_t = 0, fa = 0, ar1_scales = 1, cs_scales = 0)
    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, {
        dist_i <- sample(names(cov_mix), size = 1, prob = cov_mix)
        draw_cov(
          m,
          dist = dist_i,
          iw_df   = true_params$nu_i[i],
          iw_scale= true_params$Psi_i[[i]],
          lkj_eta = 2, half_t_df = 4, half_t_scale = 2,
          fa_rank = min(2, m), ar1_rho = 0.8, cs_rho = 0.3
        )
      }, simplify = FALSE)
    })
  }else if (prior_type == "cs") {
    cov_mix <- c(iw = 0, lkj_half_t = 0, fa = 0, ar1_scales = 0, cs_scales = 1)
    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, {
        dist_i <- sample(names(cov_mix), size = 1, prob = cov_mix)
        draw_cov(
          m,
          dist = dist_i,
          iw_df   = true_params$nu_i[i],
          iw_scale= true_params$Psi_i[[i]],
          lkj_eta = 2, half_t_df = 4, half_t_scale = 2,
          fa_rank = min(2, m), ar1_rho = 0.6, cs_rho = 0.8
        )
      }, simplify = FALSE)
    })
  }else if (prior_type == "LKJ"){
    cov_mix <- c(iw = 0.5, lkj_half_t = 0.5, fa = 0, ar1_scales = 0, cs_scales = 0)
    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, {
        dist_i <- sample(names(cov_mix), size = 1, prob = cov_mix)
        draw_cov(
          m,
          dist = dist_i,
          iw_df   = true_params$nu_i[i],
          iw_scale= true_params$Psi_i[[i]],
          lkj_eta = 2, half_t_df = 4, half_t_scale = 2,
          fa_rank = min(2, m), ar1_rho = 0.6, cs_rho = 0.8
        )
      }, simplify = FALSE)
    })
  }else if (prior_type == "random"){
    cov_mix <- c(iw = 0, lkj_half_t = 0, fa = 0, ar1_scales = 0, cs_scales = 0, random = 1)
    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, {
        dist_i <- sample(names(cov_mix), size = 1, prob = cov_mix)
        draw_cov(
          m,
          dist = dist_i,
          iw_df   = true_params$nu_i[i],
          iw_scale= true_params$Psi_i[[i]],
          lkj_eta = 2, half_t_df = 4, half_t_scale = 2,
          fa_rank = min(2, m), ar1_rho = 0.6, cs_rho = 0.8
        )
      }, simplify = FALSE)
    })
  }else if (prior_type == "mix"){
    cov_mix <- c(iw = 0.2, lkj_half_t = 0.3, fa = 0.2, ar1_scales = 0.2, cs_scales = 0.1)

    Sigma_ig <- lapply(1:K, function(i) {
      replicate(p, {
        dist_i <- sample(names(cov_mix), size = 1, prob = cov_mix)
        draw_cov(
          m,
          dist = dist_i,
          iw_df   = true_params$nu_i[i],
          iw_scale= true_params$Psi_i[[i]],
          lkj_eta = 2, half_t_df = 4, half_t_scale = 1,
          fa_rank = min(2, m), ar1_rho = 0.6, cs_rho = 0.3
        )
      }, simplify = FALSE)
    })
  }

  # Apply batch effects
  for (g in 1:p) {
    for (i in 1:K) {
      obs_in_batch <- which(batch == letters[i])
      gamma_ig_current <- gamma_ig[[i]][, g]
      Sigma_ig_current <- Sigma_ig[[i]][[g]]

      for (meas in 1:m) {
        data_list[[meas]][obs_in_batch, g] <-
          data_list[[meas]][obs_in_batch, g] + gamma_ig_current[meas]
      }

      noise <- MASS::mvrnorm(length(obs_in_batch), rep(0, m), Sigma_ig_current)
      for (meas in 1:m) {
        data_list[[meas]][obs_in_batch, g] <-
          data_list[[meas]][obs_in_batch, g] + noise[, meas]
      }

      for (meas in 1:m){
        if(add_outlier){
          outlier_idx <- lapply(1:2, function(j) sample(obs_in_batch, size = outlier_size))
          if (meas %in% c(1, 2)) {  # Make meas 1 have stronger effects
            data_list[[meas]][outlier_idx[[meas]], g] <- data_list[[meas]][outlier_idx[[meas]], g] + rnorm(outlier_size, mean = 0, sd = 2)
          }
        }
        if(add_noise && meas %in% c(1,2)){
          data_list[[meas]][obs_in_batch[[meas]], g] <- data_list[[meas]][obs_in_batch[[meas]], g] + rnorm(1, 0, 3)
        }
      }
    }
  }


  # 5. Return Results ------------------------------------------------------
  structure(
    list(
      data = lapply(1:m, function(i) data.frame(data_list[[i]])),
      true_data_list = lapply(1:m, function(i) data.frame(true_data_list[[i]])),
      batch = lapply(1:m, function(i) batch),
      covariates = lapply(1:m, function(i) covariates),
      beta = list(Age = beta_age, Sex = beta_sex, Diagnosis = beta_diag),
      biomarkers = biomarkers,
      diagnosis = lapply(1:m, function(i) Diagnosis),
      true_params = true_params,
      settings = list(
        m = m, n = n, p = p, K = K,
        global_strength = global_batch_strength,
        bt_type = bt_type, w_type = w_type,
        seed = seed
      )
    ),
    class = "simulated_multivariate_data"
  )
}


generate_cov <- function(type, m, ..., scale = 1) {
  # Input validation
  if (!is.character(type)) {
    stop("'type' must be a character string specifying covariance structure")
  }
  if (m <= 0 || !is.numeric(m) || length(m) > 1) {
    stop("'m' must be a positive integer")
  }

  # Dispatch to appropriate method
  UseMethod("generate_cov", object = structure(list(), class = type))
}

generate_cov.default <- function(type, m, ..., scale = 1) {
  stop("Unknown covariance type: ", type,
       "\nSupported types are: AR, compound_symmetry, banded, factor, random_pd, block_diag, toeplitz, wishart, ind")
}

generate_cov.ind <- function(type, m, ..., scale = 1) {
  diag(scale, m)
}

generate_cov.AR <- function(type, m, rho, ..., scale = 1) {
  if (rho <= 0 || rho >= 1) stop("rho must be between 0 and 1")
  cov_mat <- rho^abs(outer(1:m, 1:m, "-"))
  diag(cov_mat) <- 1
  cov_mat * scale
}

generate_cov.compound_symmetry <- function(type, m, rho, ..., scale = 1) {
  if (abs(rho) >= 1) stop("rho must be between -1 and 1")

  cov_mat <- matrix(rho, m, m)
  diag(cov_mat) <- 1
  cov_mat * scale
}

generate_cov.banded <- function(type, m, k = 1, decay = 0.5, ..., scale = 1) {
  if (k < 0 || k >= m) stop("k must be between 0 and m-1")

  cov_mat <- diag(m)
  for (i in 1:m) {
    for (j in max(1, i-k):min(m, i+k)) {
      if (i != j) cov_mat[i, j] <- decay^abs(i-j)
    }
  }
  cov_mat * scale
}

generate_cov.factored <- function(type, m, k = 2, noise_level = 0.1, ..., scale = 1) {
  if (k <= 0 || k > m) stop("k must be between 1 and m")

  F_mat <- matrix(rnorm(m * k), m, k)
  D <- diag(runif(m, min = 0.1, max = noise_level))
  cov_mat <- tcrossprod(F_mat) + D
  cov_mat * scale
}

generate_cov.random_pd <- function(type, m, scale = 1, ...) {
  # Input validation
  if (m <= 0 || !is.numeric(m) || length(m) > 1) {
    stop("'m' must be a positive integer")
  }
  if (scale <= 0) stop("'scale' must be positive")

  A <- matrix(rnorm(m^2, mean = 0, sd = 1/sqrt(m)), m, m)
  cov_mat <- tcrossprod(A)  # A %*% t(A) is positive definite
  cov_mat * scale
}

generate_cov.block_diag <- function(type, m, block_size = 2, rho, scale = 1, ...) {
  if (m %% block_size != 0) {
    stop("'block_size' must exactly divide 'm'")
  }
  if (abs(rho) >= 1) stop("'rho' must be between -1 and 1")
  if (scale <= 0) stop("'scale' must be positive")

  n_blocks <- m / block_size
  cov_mat <- matrix(0, m, m)

  for (b in 1:n_blocks) {
    idx <- ((b-1)*block_size + 1):(b*block_size)
    block <- matrix(rho, block_size, block_size)
    diag(block) <- 1
    cov_mat[idx, idx] <- block
  }

  cov_mat * scale
}

generate_cov.toeplitz <- function(type, m, rho, scale = 1, ...) {
  # Input validation
  if (rho <= 0 || rho >= 1) stop("'rho' must be between 0 and 1")
  if (scale <= 0) stop("'scale' must be positive")

  # Create Toeplitz matrix with more stable computation
  exponents <- abs(outer(1:m, 1:m, "-"))
  cov_mat <- rho^exponents
  cov_mat * scale
}

generate_cov.wishart <- function(type, m, df = 15, scale = 1, ...) {
  # Input validation
  if (df <= m - 1) stop("'df' must be greater than m-1")
  if (scale <= 0) stop("'scale' must be positive")

  # More numerically stable Wishart generation
  A <- matrix(rnorm(m * ceiling(df)), ceiling(df), m)
  cov_mat <- crossprod(A) / df  # A'A/df ~ Wishart(I, df)
  cov_mat * scale
}

draw_lkj_corr <- function(m, eta = 2) {
  A <- matrix(rnorm(m*m), m, m); S <- crossprod(A) + diag(1e-3, m)
  sd <- sqrt(diag(S)); S / (sd %*% t(sd))
}

draw_cov <- function(m,
                     dist = c("iw", "lkj_half_t", "fa", "ar1_scales", "cs_scales", "random"),
                     iw_df = m + 2, iw_scale = diag(m),
                     lkj_eta = 2, half_t_df = 4, half_t_scale = 1,
                     fa_rank = min(2, m), fa_loading_sd = 0.7,
                     fa_uni_min = 0.2, fa_uni_max = 1.2,
                     ar1_rho = 0.6, cs_rho = 0.3) {
  dist <- match.arg(dist)

  if (dist == "iw") {
    return(MCMCpack::riwish(iw_df, iw_scale))
  }

  if (dist == "lkj_half_t") {
    C <- draw_lkj_corr(m, eta = lkj_eta)
    s <- abs(stats::rt(m, df = half_t_df)) * half_t_scale  # half-t scales
    S <- diag(s, m)
    return(S %*% C %*% S)
  }

  if (dist == "fa") {
    r <- fa_rank
    Lambda <- matrix(rnorm(m * r, 0, fa_loading_sd), m, r)
    D <- diag(runif(m, fa_uni_min, fa_uni_max), m)
    return(Lambda %*% t(Lambda) + D)
  }

  if (dist == "ar1_scales") {
    C <- generate_cov.AR(m = m, rho = ar1_rho)
    s <- abs(stats::rt(m, df = half_t_df)) * half_t_scale
    S <- diag(s, m)
    return(S %*% C %*% S)
  }

  if (dist == "cs_scales") {
    C <- generate_cov.compound_symmetry(m = m, rho = cs_rho)
    s <- abs(stats::rt(m, df = half_t_df)) * half_t_scale
    S <- diag(s, m)
    return(S %*% C %*% S)
  }

  if (dist == "random") {
    C <- generate_cov.random_pd(m = m)
    s <- abs(stats::rt(m, df = half_t_df)) * half_t_scale
    S <- diag(s, m)
    return(S %*% C %*% S)
  }
}

random_orthogonal <- function(m) {
  Z <- matrix(rnorm(m*m), m)
  qr.Q(qr(Z))
}

householder_rot <- function(m) {
  v <- rnorm(m); v <- v / sqrt(sum(v^2))
  diag(m) - 2 * (v %*% t(v))
}

covbat_severity_gen <- function(level = c("mild","moderate","severe"), r, m) {
  level <- match.arg(level)
  mult <- rep(1, m)
  bump <- switch(level,
                 mild     = c(1.2, 0.9),
                 moderate = c(1.8, 0.75),
                 severe   = c(3.0, 0.6))
  mult[1:r] <- seq(bump[1], (bump[1]+1)/2, length.out = r)
  list(mult = mult,
       delta = switch(level, mild=0.10, moderate=0.25, severe=0.50))
}

make_covbat_sigmas <- function(Sigma0, K, r = 3,
                               rotate = c("none","householder","random"),
                               severity = c("mild","moderate","severe"),
                               jitter_sd = 0.02, seed = 1) {
  set.seed(seed)
  rotate   <- match.arg(rotate)
  severity <- match.arg(severity)

  eig0 <- eigen(Sigma0, symmetric = TRUE)
  Q0   <- eig0$vectors
  L0   <- eig0$values
  sev  <- covbat_severity_gen(severity, r = r, m = length(L0))

  Sigmas_k <- vector("list", K)
  for (k in 1:K) {
    Qk <- switch(rotate,
                 none        = Q0,
                 householder = Q0 %*% householder_rot(ncol(Q0)),
                 random      = Q0 %*% random_orthogonal(ncol(Q0)))
    Lk <- L0 * sev$mult
    Uk <- Qk[, 1:r, drop = FALSE]
    bump <- sev$delta * (Uk %*% t(Uk))
    Sigmak <- Qk %*% diag(Lk) %*% t(Qk) + bump

    # ensure PD
    ev <- eigen(Sigmak, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-8) Sigmak <- Sigmak + (abs(min(ev)) + 1e-3) * diag(ncol(Sigmak))
    Sigmas_k[[k]] <- Sigmak
  }

  jitter_fun <- function(S, sd = jitter_sd) S + diag(abs(rnorm(nrow(S), 0, sd)))
  list(Sigmas_k = Sigmas_k, jitter_fun = jitter_fun)
}

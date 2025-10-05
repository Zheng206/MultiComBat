.need_cmdstanr <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "This feature requires the 'cmdstanr' package.\n",
      "Install via:\n",
      "install.packages('cmdstanr', ",
      "repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n",
      call. = FALSE
    )
  }
  if (is.na(cmdstanr::cmdstan_version())) {
    stop("CmdStan is not installed. Run: cmdstanr::install_cmdstan()", call. = FALSE)
  }
  invisible(TRUE)
}


#' Locate and (lazily) compile a Stan model
#'
#' Finds \code{inst/stan/<name>.stan} inside the package and compiles it to a
#' user-writable cache directory (reused across sessions), returning a
#' \code{cmdstanr::CmdStanModel}. This avoids compiling into the read-only
#' package directory and avoids heavy work at install/load time.
#'
#' @param name Basename (without \code{.stan}) of the model
#'   (e.g., \code{"univariate_model"} or \code{"multivariate_model"}).
#' @param compile_dir Directory for the compiled executable. Defaults to
#'   \code{tools::R_user_dir("MultiComBat", "cache")}.
#' @param force_recompile Logical; recompile even if an existing binary is found.
#' @return A \code{cmdstanr::CmdStanModel}.
#' @seealso \code{cmdstanr::cmdstan_model}
#' @examples
#' \dontrun{
#' mdl <- get_stan_model("univariate_model")
#' fit <- mdl$sample(data = list(...), chains = 2)
#' }
#' @export
get_stan_model <- function(name = "my_model",
                           compile_dir = tools::R_user_dir("MultiComBat", "cache"),
                           force_recompile = FALSE) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required. See https://mc-stan.org/cmdstanr/")
  }
  stan_path <- system.file("stan", paste0(name, ".stan"), package = "MultiComBat")
  if (!file.exists(stan_path)) stop("Stan model not found: ", name)
  dir.create(compile_dir, showWarnings = FALSE, recursive = TRUE)

  mod <- cmdstanr::cmdstan_model(stan_path, compile = FALSE)
  mod$compile(dir = compile_dir, force_recompile = isTRUE(force_recompile))
  mod
}

#' Stan algorithm (S3 generic)
#' @param type S3 tag (e.g., \code{combat_type("univariate")} or \code{"multivariate"}).
#' @param ... Passed to methods.
#' @export
stan_algorithm <- function(type, ...) {
  UseMethod("stan_algorithm")
}

#' Build Stan data (S3 generic)
#' @param type S3 tag (e.g., \code{combat_type("univariate")} or \code{"multivariate"}).
#' @param ... Passed to methods.
#' @export
stan_data_prep <- function(type, ...) {
  UseMethod("stan_data_prep")
}

#' Long-format stacking across measurements for Stan
#'
#' Convert standardized per-measurement matrices into a long data frame with
#' observation, feature, measurement, and batch labels.
#'
#' @param data_stand_result For the multivariate method, a \emph{list} of length
#'   \eqn{m} (measurements/modalities). Each element is the output of
#'   \code{standardize_data()} for that measurement.
#' @param batch_list List of length \eqn{m} of length-\eqn{n} batch labels.
#' @return A data frame with columns: \code{value}, \code{measurement}, \code{feature},
#'   \code{observation}, \code{batch}, \code{batch_id}, \code{feature_id}.
#' @export

data_reshape <- function(data_stand_result, batch_list){
  m <- length(data_stand_result)
  data_list <- lapply(1:m, function(i) data_stand_result[[i]]$data_stand)
  y_long <- map_dfr(1:length(data_list), function(meas) {
    data.frame(
      value = as.vector(as.matrix(data_list[[meas]])),  # Flatten all values
      measurement = meas,                              # Measurement ID
      feature = rep(1:ncol(data_list[[meas]]), each = nrow(data_list[[meas]])),
      observation = rep(1:nrow(data_list[[meas]]), times = ncol(data_list[[meas]]))
    ) %>%
      mutate(batch = batch_list[[meas]][observation])
  }) %>%
    mutate(
      batch_id = as.integer(factor(batch)),
      feature_id = as.integer(factor(feature))
    )
  return(y_long)
}


#' Prepare Stan data (multivariate)
#'
#' Builds a 3D array \eqn{y[N, G, M]} and batch indices for multivariate Stan models.
#'
#' @param type Ignored; S3 signature alignment.
#' @param data_stand_result Output from \code{standardize_data()}.
#' @param batch_list List of batch vectors (one per measurement).
#' @param ... Unused; for S3 compatibility.
#' @return A list with elements \code{N, G, M, I, i_id, y}.
#' @method stan_data_prep multivariate
#' @export

stan_data_prep.multivariate <- function(type, data_stand_result, batch_list, ...){
  M <- length(data_stand_result)
  G <- ncol(data_stand_result[[1]]$data_stand)
  N <- nrow(data_stand_result[[1]]$data_stand)
  y_long <- data_reshape(data_stand_result, batch_list)
  y_array <- y_long %>%
    dplyr::select(observation, feature_id, measurement, value) %>%
    pivot_wider(names_from = measurement, values_from = value) %>%
    arrange(feature_id, observation) %>%
    dplyr::select(-observation, -feature_id) %>%
    as.matrix() %>%
    array(dim = c(N, G, M))

  i_id <- y_long$batch_id[1:max(y_long$observation)]
  N <- max(y_long$observation)
  M <- max(y_long$measurement)
  I <- max(y_long$batch_id)
  G <- max(y_long$feature_id)

  stan_data <- list(
    N = N,
    G = G,
    M = M,
    I = I,
    i_id = i_id,
    y = y_array
  )
  return(stan_data)
}

#' Prepare Stan data (univariate)
#'
#' Builds long-format vectors/indices for a univariate Stan model over all
#' observations and features.
#'
#' @param type Ignored; S3 signature alignment.
#' @param data_stand_result Output from \code{standardize_data()}.
#' @param batch_list Vector of batch labels (length \eqn{n}).
#' @param ... Unused; for S3 compatibility.
#' @return A list with elements \code{N, G, I, y, g, i} and weakly-informative hyperparameters.
#' @method stan_data_prep univariate
#' @export
stan_data_prep.univariate <- function(type, data_stand_result, batch_list, ...){
  data_list <- data_stand_result$data_stand
  y_long <- data.frame(
    value = as.vector(as.matrix(data_list)),  # Flatten all values
    feature = rep(1:ncol(data_list), each = nrow(data_list)),
    observation = rep(1:nrow(data_list), times = ncol(data_list))
  ) %>% mutate(batch = batch_list[observation]) %>%
    mutate(
      batch_id = as.integer(factor(batch)),
      feature_id = as.integer(factor(feature))
    )
  stan_data <- list(
    N = nrow(y_long),
    G = max(y_long$feature_id),
    I = max(y_long$batch_id),
    y = y_long$value,
    g = y_long$feature_id,
    i = y_long$batch_id,
    nu_psi = 3,
    nu_tau = 4,
    nu_delta = 3,
    A_psi = 2.5,
    A_tau = 2.5,
    A_delta = 2.5
  )
  return(stan_data)
}

# Internal helpers (not exported)
gamma_star_reshape.multivariate <- function(gamma_star, batch_ind){
  result <- gamma_star %>% filter(batch == batch_ind) %>% pivot_wider(names_from = feature, values_from = mean) %>% dplyr::select(-(1:2)) %>% as.matrix()
  return(result)
}

# Internal helpers (not exported)
gamma_star_reshape.univariate <- function(gamma_star){
  result <- gamma_star %>% pivot_wider(names_from = feature, values_from = mean) %>% dplyr::select(-1) %>% as.matrix()
  return(result)
}

# Internal helpers (not exported)
delta_star_reshape <- function(delta_star, batch_ind){
  result <- delta_star %>% filter(batch == batch_ind) %>% split(list(.$batch, .$feature)) %>%
    map(~ {
      m <- max(.x$m1)
      matrix(.x$mean, nrow = m, ncol = m, byrow = FALSE)
    })
  return(result)
}

#' @describeIn stan_algorithm Multivariate Stan routine
#' @param stan_data List created by \code{stan_data_prep.multivariate()}.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param batch_names Character vector of batch level names (used to name outputs).
#' @param ... Unused; for S3 compatibility.
#' @return For multivariate: a list with \code{gamma_star} (named list of \eqn{m \times G} matrices),
#'   \code{delta_star} (named list of length-\eqn{G} lists of \eqn{m \times m} matrices), and mcmc fitted model.
#' @method stan_algorithm multivariate
#' @export
stan_algorithm.multivariate <- function(type, stan_data, chains = 3, parallel_chains = 3, batch_names, ...){
  message("Starting multivariate Stan algorithm...")
  start_time <- Sys.time()
  message("[1/4] Loading Stan model...")
  stan_model <- get_stan_model("multivariate_model")
  message("[2/4] Running MCMC sampling (this may take a while)...")
  fit <- stan_model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    adapt_delta = 0.98
  )
  message("[3/4] Processing gamma_star results...")
  gamma_star <- fit$summary(variable = "mu_ig", "mean") %>%
    tidyr::extract(
      variable,
      into = c("batch", "feature", "measurement"),
      regex = "mu_ig\\[([0-9]+),([0-9]+),([0-9]+)\\]",
      convert = TRUE  # Converts to integers
    )
  message("[4/4] Processing delta_star results...")
  delta_star <- fit$summary(variable = "Sigma_ig", "mean") %>%
    tidyr::extract(
      variable,
      into = c("batch", "feature", "m1", "m2"),
      regex = "Sigma_ig\\[([0-9]+),([0-9]+),([0-9]+),([0-9]+)\\]",
      convert = TRUE
    )
  message("Reshaping results...")
  gamma_star <- lapply(1:max(gamma_star$batch), function(i) gamma_star_reshape.multivariate(gamma_star, batch_ind = i))
  names(gamma_star) <- batch_names

  delta_star <- lapply(1:max(delta_star$batch), function(i) delta_star_reshape(delta_star, batch_ind = i))
  names(delta_star) <- batch_names

  end_time <- Sys.time()
  message(sprintf("\nCompleted in %s", format(end_time - start_time)))
  return(list(gamma_star = gamma_star, delta_star = delta_star, fit = fit))
}

#' @describeIn stan_algorithm Univariate Stan routine
#' @param stan_data List created by \code{stan_data_prep.univariate()}.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param batch_names Character vector of batch level names (row names for outputs).
#' @param ... Unused; for S3 compatibility.
#' @return For univariate: a list with matrices \code{gamma_star} and \code{delta_star}
#'   (rows = batches, columns = features).
#' @method stan_algorithm univariate
#' @export
stan_algorithm.univariate <- function(type, stan_data, chains = 3, parallel_chains = 3, batch_names, ...){
  message("Starting univariate Stan algorithm...")
  start_time <- Sys.time()
  message("[1/4] Loading Stan model...")
  stan_model <- get_stan_model("univariate_model")
  message("[2/4] Running MCMC sampling (this may take a while)...")
  fit <- stan_model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    adapt_delta = 0.98,
    show_messages = TRUE
  )
  message("[3/4] Processing gamma_star results...")
  gamma_star <- fit$summary(variable = "mu_ig", "mean") %>%
    tidyr::extract(
      variable,
      into = c("feature", "batch"),
      regex = "mu_ig\\[([0-9]+),([0-9]+)\\]",
      convert = TRUE  # Converts to integers
    )
  message("[4/4] Processing delta_star results...")
  delta_star <- fit$summary(variable = "sigma_delta", "mean") %>%
    dplyr::mutate(mean = mean^2) %>%
    tidyr::extract(
      variable,
      into = c("feature", "batch"),
      regex = "sigma_delta\\[([0-9]+),([0-9]+)\\]",
      convert = TRUE  # Converts to integers
    )
  message("Reshaping results...")
  gamma_star <- gamma_star_reshape.univariate(gamma_star)
  rownames(gamma_star) <- batch_names


  delta_star <- gamma_star_reshape.univariate(delta_star)
  rownames(delta_star) <- batch_names

  end_time <- Sys.time()
  message(sprintf("\nCompleted in %s", format(end_time - start_time)))
  return(list(gamma_star = gamma_star, delta_star = delta_star))
}

utils::globalVariables(c(
  "variable","feature","batch","measurement","m1","m2","mean","observation", ".", "feature_id", "value"
))

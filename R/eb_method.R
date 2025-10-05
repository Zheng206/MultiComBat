#' Create a ComBat type tag for S3 dispatch
#'
#' Convenience constructor for tagging an object with a class used to choose
#' univariate vs. multivariate ComBat methods.
#'
#' @param type Character scalar, typically `"univariate"` or `"multivariate"`.
#'
#' @return An object whose class is set to \code{type}; primarily used for S3 dispatch.
#'
#' @examples
#' class(combat_type("multivariate"))
#'
#' @export

combat_type <- function(type){
  result <- type
  class(result) <- type
  return(result)
}

#' Empirical Bayes shrinkage for ComBat parameters
#'
#' S3 generic for running the empirical Bayes (EB) step that shrinks per-batch
#' location (\eqn{\gamma}) and scale/covariance (\eqn{\delta}) parameters toward
#' pooled priors, either feature-wise (univariate) or measurement-wise
#' (multivariate).
#'
#' @param type An object created by \code{\link{combat_type}} or a character
#'   scalar equal to its class (\code{"univariate"} or \code{"multivariate"}).
#' @param ... Passed to methods.
#'
#' @return See method-specific return values below.
#'
#' @examples
#' \dontrun{
#' type <- combat_type("univariate")
#' # eb_algorithm(type, data_stand_result, batch_result)
#' }
#'
#' @export
eb_algorithm <- function(type, ...) {
  UseMethod("eb_algorithm")
}


#' Compute batch-wise means of standardized data (S3 generic)
#'
#' @param type An object created by \code{\link{combat_type}} or a character
#'   scalar equal to its class.
#' @param ... Passed to methods.
#'
#' @return Method-specific structure summarizing \eqn{\gamma} (batch means).
#'
#' @export

gamm_hat_gen <- function(type, ...){
  UseMethod("gamm_hat_gen")
}

#' Compute batch-wise dispersion (variance/covariance) estimates (S3 generic)
#'
#' @param type An object created by \code{\link{combat_type}} or a character scalar.
#' @param ... Passed to methods.
#'
#' @return Method-specific structure summarizing \eqn{\delta}
#'   (variances or covariance matrices).
#'
#' @export

delta_hat_gen <- function(type, ...){
  UseMethod("delta_hat_gen")
}

#' Method-of-moments hyperparameters for EB (S3 generic)
#'
#' Estimate prior hyperparameters for the EB step from batch-wise summaries.
#'
#' @param type An object created by \code{\link{combat_type}} or a character scalar.
#' @param ... Passed to methods.
#'
#' @return Method-specific list of hyperparameters used by \code{\link{eb_one_iteration}}
#'   and \code{\link{eb_algorithm}}.
#'
#' @export

mom_calculation <- function(type, ...){
  UseMethod("mom_calculation")
}

#' One EB fixed-point iteration (S3 generic)
#'
#' Perform a single iteration updating batch mean(s) \eqn{\gamma} and
#' variance/covariance \eqn{\delta}, given current values and EB hyperparameters.
#'
#' @param type An object created by \code{\link{combat_type}} or a character scalar.
#' @param ... Passed to methods.
#'
#' @return A list with updated parameters and a scalar convergence measure
#'   \code{change}.
#'
#' @export

eb_one_iteration <- function(type, ...){
  UseMethod("eb_one_iteration")
}


#' @rdname eb_one_iteration
#' @method eb_one_iteration univariate
#'
#' @param bdat Subset of standardized data for a single batch (\eqn{n_b \times G}).
#' @param g_orig Vector of original per-feature means for this batch (length \eqn{G}).
#' @param g_old Current estimate of \eqn{\gamma} (length \eqn{G}).
#' @param d_old Current estimate of \eqn{\delta} (length \eqn{G}).
#' @param mom Hyperparameters list from \code{\link{mom_calculation.univariate}}.
#'
#' @return List with \code{g_new} (length \eqn{G}), \code{d_new} (length \eqn{G}),
#'   and \code{change} (scalar).
#'
#' @export

eb_one_iteration.univariate <- function(type, bdat, g_orig, g_old, d_old, mom, ...){
  n_b <- nrow(bdat)
  g_new <- (n_b*mom$g_var*g_orig + d_old*mom$g_bar)/(n_b*mom$g_var + d_old)
  sum2   <- colSums(sweep(bdat, 2, g_new)^2)
  d_new <- (sum2/2 + mom$d_b)/(n_b/2 + mom$d_a - 1)
  change <- max(abs(g_new - g_old)/g_old, abs(d_new - d_old)/d_old)
  return(list("g_new" = g_new, "d_new" = d_new, "change" = change))
}

#' @rdname eb_one_iteration
#' @method eb_one_iteration multivariate
#'
#' @param bdat List of length \eqn{m}; each element is the standardized data for
#'   one measurement in the same batch (\eqn{n_b \times G}).
#' @param g_orig Matrix \eqn{m \times G} of original per-feature mean vectors.
#' @param g_old Matrix \eqn{m \times G} of current mean estimates.
#' @param d_old Length-\eqn{G} list of \eqn{m \times m} covariance estimates.
#' @param mom Hyperparameters list from \code{\link{mom_calculation.multivariate}}.
#'
#' @details
#' Uses helper functions \code{robust_matrix_solver()} and
#' \code{multi_matrix_transform()} (must be available in the package namespace)
#' to stabilize inversions and assemble feature-specific matrices.
#'
#' @return List with \code{g_new} (\eqn{m \times G}), \code{d_new} (list of
#'   \eqn{m \times m}), and \code{change} (scalar).
#'
#' @export

eb_one_iteration.multivariate <- function(type, bdat, g_orig, g_old, d_old, mom, ...){
  n_b <- nrow(bdat[[1]])
  G <- length(d_old)
  g_new_list <- list()
  d_new_list <- list()

  for (g in 1:G) {
    T_i_inv <- robust_matrix_solver(mom$g_var)
    d_old_inv <- robust_matrix_solver(d_old[[g]])

    # Compute new gamma estimate
    Omega_ig <- robust_matrix_solver(n_b * d_old_inv + T_i_inv)
    mu_ig <- Omega_ig %*% (n_b * d_old_inv %*% g_orig[, g] + T_i_inv %*% mom$g_bar)
    g_new_list[[g]] <- mu_ig

    # Compute new covariance estimate
    bdat_matrix <- multi_matrix_transform(bdat, g)
    centered_data <- sweep(bdat_matrix, 1, mu_ig, FUN = "-")
    sum2 <- tcrossprod(centered_data)
    d_new_list[[g]] <- (sum2 + mom$d_b) / (n_b + mom$d_a - length(mu_ig) - 1)
  }

  # Convert lists to matrices
  g_new <- do.call(cbind, g_new_list)
  colnames(g_new) <- colnames(g_old)
  d_new <- d_new_list
  names(d_new) <- names(d_old)

  # Compute convergence criterion
  change_g <- max(sapply(1:G, function(g) max(abs(g_new[, g] - g_old[, g]) / abs(g_old[, g]))))
  change_d <- max(sapply(1:G, function(g) max(abs(d_new[[g]] - d_old[[g]]) / abs(d_old[[g]]))))
  change <- max(change_g, change_d)
  return(list("g_new" = g_new, "d_new" = d_new, "change" = change))
}


#' @rdname gamm_hat_gen
#' @method gamm_hat_gen univariate
#'
#' @param data_stand_result Output from \code{standardize_data()}.
#' @param batch_result Output from \code{batch_matrix()}.
#'
#' @return For the univariate method, a matrix \code{gamma_hat} of size
#'   (batches \eqn{\times} features) with per-batch feature means.
#'
#' @export

gamm_hat_gen.univariate <- function(type, data_stand_result, batch_result, ...){
  data_stand <- data_stand_result$data_stand
  batch_vector <- batch_result$batch_vector
  gamma_hat <- Reduce(rbind, by(data_stand, batch_vector, function(x) apply(x, 2, mean)))
  rownames(gamma_hat) <- levels(batch_vector)
  return(gamma_hat)
}

#' @rdname delta_hat_gen
#' @method delta_hat_gen univariate
#'
#' @param data_stand_result Output from \code{standardize_data()}.
#' @param batch_result Output from \code{batch_matrix()}.
#' @param robust.LS Logical; if \code{TRUE}, use \code{\link{biweight_midvar}}
#'   instead of \code{\link[stats]{var}}.
#'
#' @return For the univariate method, a matrix \code{delta_hat} of size
#'   (batches \eqn{\times} features) with per-batch feature variances.
#'
#' @export

delta_hat_gen.univariate <- function(type, data_stand_result, batch_result, robust.LS = FALSE, ...){
  data_stand <- data_stand_result$data_stand
  batch_vector <- batch_result$batch_vector
  if(robust.LS){
    delta_hat <- Reduce(rbind, by(data_stand, batch_vector, function(x) apply(x, 2, biweight_midvar)))
  }else{
    delta_hat <- Reduce(rbind, by(data_stand, batch_vector, function(x) apply(x, 2, var)))
  }
  rownames(delta_hat) <- levels(batch_vector)
  return(delta_hat)
}

#' @rdname mom_calculation
#' @method mom_calculation univariate
#'
#' @param gamma_hat Matrix of per-batch means (batches \eqn{\times} features).
#' @param delta_hat Matrix of per-batch variances (batches \eqn{\times} features).
#'
#' @details
#' Univariate method matches moments to a normal prior for \eqn{\gamma} with
#' mean \code{g_bar} and variance \code{g_var}, and an inverse-gamma prior for
#' \eqn{\delta} parameterized by \code{d_a}, \code{d_b}.
#'
#' @return List with \code{g_bar}, \code{g_var}, \code{d_bar}, \code{d_var},
#'   \code{d_a}, \code{d_b}.
#'
#' @export

mom_calculation.univariate <- function(type, gamma_hat, delta_hat, ...){
  g_bar <- mean(gamma_hat)
  g_var <- var(gamma_hat)
  d_bar <- mean(delta_hat)
  d_var <- var(delta_hat)
  d_a <- (2 * d_var + d_bar^2)/d_var
  d_b <- (d_bar * d_var + d_bar^3)/d_var
  return(list("g_bar" = g_bar, "g_var" = g_var, "d_bar" = d_bar, "d_var" = d_var, "d_a" = d_a, "d_b" = d_b))
}

#' @rdname gamm_hat_gen
#' @method gamm_hat_gen multivariate
#'
#' @param data_stand_result List of \code{standardize_data()} outputs (one per measurement).
#' @param batch_result List of \code{batch_matrix()} outputs (one per measurement).
#'
#' @return For the multivariate method, a named list over batches where each
#'   element is an \eqn{m \times G} matrix: rows = measurements/modalities,
#'   columns = features.
#'
#' @export

gamm_hat_gen.multivariate <- function(type, data_stand_result, batch_result, ...){
  m = length(batch_result)
  gamma_hat_list <- lapply(1:m, function(i) {
    gamma_hat_df <- Reduce(rbind, by(data_stand_result[[i]]$data_stand, batch_result[[i]]$batch_vector, function(x) apply(x, 2, mean)))
    rownames(gamma_hat_df) <- levels(batch_result[[i]]$batch_vector)
    return(gamma_hat_df)
  })
  gamma_hat_batch_list <- lapply(levels(batch_result[[1]]$batch_vector), function(b) {
    gamma_hat_batch <- lapply(1:m, function(i) gamma_hat_list[[i]][b, ]) %>% bind_rows() %>% as.matrix()
  })
  names(gamma_hat_batch_list) <- levels(batch_result[[1]]$batch_vector)
  return(gamma_hat_batch_list)
}

#' @rdname delta_hat_gen
#' @method delta_hat_gen multivariate
#'
#' @param data_stand_result List of \code{standardize_data()} outputs (one per measurement).
#' @param batch_result List of \code{batch_matrix()} outputs (one per measurement).
#' @param gamma_hat_list Output of \code{gamm_hat_gen.multivariate()} (used for centering).
#' @param robust.LS Logical; if \code{TRUE}, use \code{\link{biweight_midvar_mul}}
#'   to form robust \eqn{m \times m} covariances; otherwise use classical covariance.
#'
#' @return For the multivariate method, a named list over batches; each element is
#'   a length-\eqn{G} list of \eqn{m \times m} covariance matrices (one per feature).
#'
#' @export

delta_hat_gen.multivariate <- function(type, data_stand_result, batch_result, gamma_hat_list, robust.LS = FALSE, ...){
  m <- length(batch_result)
  delta_hat_batch_list <- lapply(levels(batch_result[[1]]$batch_vector), function(b) {
    Z_i <- lapply(1:m, function(i) data_stand_result[[i]]$data_stand[batch_result[[i]]$batch_index[[b]], ])
    G <- ncol(Z_i[[1]])
    delta_hat_batch <- lapply(1:G, function(g) {
      Zi_matrix <- multi_matrix_transform(Z_i, g)
      gamma_hat_ig <- gamma_hat_list[[b]][, g]
      centered_Zi <- sweep(Zi_matrix, 1, gamma_hat_ig, FUN = "-")
      ni <- ncol(centered_Zi)
      if(robust.LS){
        Sigma_hat_ig <- biweight_midvar_mul(centered_Zi)
      }else{
        Sigma_hat_ig <- (centered_Zi %*% t(centered_Zi)) / (ni - 1)
      }
      return(Sigma_hat_ig) }
    )
    names(delta_hat_batch) <- colnames(Z_i[[1]])
    return(delta_hat_batch)
  })
  names(delta_hat_batch_list) <- levels(batch_result[[1]]$batch_vector)
  return(delta_hat_batch_list)
}

#' @rdname mom_calculation
#' @method mom_calculation multivariate
#'
#' @param gamma_hat_batch An \eqn{m \times G} matrix of batch-wise means
#'   across measurements (rows) and features (columns) for a single batch.
#' @param delta_hat_batch Length-\eqn{G} list of \eqn{m \times m} covariance
#'   matrices for that batch.
#'
#' @details
#' Multivariate method matches moments to a matrix-normal / inverse-Wishart–like
#' prior, producing mean vector \code{g_bar} (\eqn{m}), covariance \code{g_var}
#' (\eqn{m \times m}) for \eqn{\gamma}, and IW parameters \code{d_a} (df)
#' and \code{d_b} (\eqn{m \times m} scale) for \eqn{\Delta}. Diagnostics
#' \code{d_var} provide dispersion summaries across features.
#'
#' @return List with \code{g_bar} (\eqn{m}), \code{g_var} (\eqn{m \times m}),
#'   \code{d_bar} (\eqn{m \times m}), \code{d_var} (\eqn{m \times m}),
#'   \code{d_a} (scalar df), \code{d_b} (\eqn{m \times m} scale).
#'
#' @export

mom_calculation.multivariate <- function(type, gamma_hat_batch, delta_hat_batch, ...) {
  if (length(delta_hat_batch) == 0) stop("No delta matrices provided.")
  p <- nrow(delta_hat_batch[[1]])
  G <- length(delta_hat_batch)
  if (G < p + 3) warning("Sample size too small for reliable estimates.")

  g_bar <- rowMeans(gamma_hat_batch)
  centered_gamma <- sweep(gamma_hat_batch, 1, g_bar, FUN = "-")
  T_i <- tcrossprod(centered_gamma) / (G - 1)

  V_bar_i <- Reduce(`+`, delta_hat_batch) / G

  # Compute element-wise variances
  S_diff <- lapply(delta_hat_batch, function(S) S - V_bar_i)
  S_var <- Reduce(`+`, lapply(S_diff, function(S) S^2)) / (G - 1)
  avg_S_var <- mean(diag(S_var))
  avg_V_bar_sq <- mean(diag(V_bar_i^2))
  nu_i <- p + 3 + (2 * avg_V_bar_sq) / avg_S_var
  nu_i <- max(nu_i, p + 1 + 1e-6)

  # Estimate Psi
  Psi_i <- V_bar_i * (nu_i - p - 1)

  # Optional diagnostics
  S_bar_i <- Reduce(`+`, lapply(S_diff, function(S) tcrossprod(S))) / (G - 1)

  return(list(
    g_bar = g_bar,
    g_var = T_i,
    d_bar = V_bar_i,
    d_var = S_bar_i,
    d_a = nu_i,
    d_b = Psi_i
  ))
}


#' @rdname eb_algorithm
#' @method eb_algorithm univariate
#'
#' @param data_stand_result List returned by \code{standardize_data()} with
#'   components \code{$data_stand}, \code{$stand_mean}, \code{$sd_mat}.
#' @param batch_result List returned by \code{batch_matrix()} with components
#'   including \code{$batch_vector}, \code{$batch_index}, \code{$n_batches}.
#' @param eb Logical; if \code{TRUE} (default) perform EB shrinkage. If
#'   \code{FALSE}, return method-of-moments estimates without shrinkage.
#' @param robust.LS Logical; when \code{TRUE}, use robust variance estimates in
#'   \code{\link{delta_hat_gen}} via \code{\link{biweight_midvar}}.
#'
#' @return For the univariate method, a list with:
#' \describe{
#'   \item{\code{gamma_star}}{Batch-by-feature matrix of shrunken means.}
#'   \item{\code{delta_star}}{Batch-by-feature matrix of shrunken variances.}
#'   \item{\code{gamma_hat}}{Batch-by-feature sample means before EB.}
#'   \item{\code{delta_hat}}{Batch-by-feature sample/robust variances before EB.}
#'   \item{\code{mom}}{List of hyperparameters (\code{g_bar}, \code{g_var},
#'     \code{d_bar}, \code{d_var}, \code{d_a}, \code{d_b}) per batch.}
#' }
#'
#' @export

eb_algorithm.univariate <- function(type, data_stand_result, batch_result, eb = TRUE, robust.LS = FALSE, ...){
  data_stand <- data_stand_result$data_stand
  gamma_hat <- gamm_hat_gen(type, data_stand_result, batch_result)
  delta_hat <- delta_hat_gen(type, data_stand_result, batch_result, robust.LS = robust.LS)
  if(eb){
    gamma_star <- NULL
    delta_star <- NULL
    batch_level <- levels(batch_result$batch_vector)
    eb_result <- lapply(batch_level, function(b){
      n_b <- batch_result$n_batches[b]
      mom <- mom_calculation.univariate(type, gamma_hat[b, ], delta_hat[b, ])
      # adjust within batch
      bdat <- data_stand[batch_result$batch_index[[b]],]
      g_orig <- gamma_hat[b,]
      g_old  <- gamma_hat[b,]
      d_old  <- delta_hat[b,]
      change_old <- 1
      change <- 1
      count  <- 0
      while(change > 10e-5){
        eb_one_result <- eb_one_iteration(type, bdat, g_orig, g_old, d_old, mom)
        change <- eb_one_result$change
        if (count > 30) {
          if (change > change_old) {
            warning("Empirical Bayes step failed to converge after 30 iterations,
    	            using estimate before change between iterations increases.")
            break
          }
        }
        g_old <- eb_one_result$g_new
        d_old <- eb_one_result$d_new
        change_old <- eb_one_result$change
        count <- count+1
      }
      gamma_star <- data.frame(t(eb_one_result$g_new))
      delta_star <- data.frame(t(eb_one_result$d_new))
      return(list("gamma_star" = gamma_star, "delta_star" = delta_star, "mom" = mom))
    })
    gamma_star <- lapply(1:length(batch_level), function(i) eb_result[[i]]$gamma_star) %>% bind_rows() %>% as.matrix()
    delta_star <- lapply(1:length(batch_level), function(i) eb_result[[i]]$delta_star) %>% bind_rows() %>% as.matrix()
    mom <- lapply(1:length(batch_level), function(i) eb_result[[i]]$mom)
    rownames(gamma_star) <- rownames(delta_star) <- batch_level
  }else{
    gamma_star <- gamma_hat
    delta_star <- delta_hat
    mom <- NULL
  }
  return(list("gamma_star" = gamma_star, "delta_star" = delta_star, "gamma_hat" = gamma_hat, "delta_hat" = delta_hat, "mom" = mom))
}

#' @rdname eb_algorithm
#' @method eb_algorithm multivariate
#'
#' @param data_stand_result For the multivariate method, a \emph{list} of length
#'   \eqn{m} (measurements/modalities). Each element is the output of
#'   \code{standardize_data()} for that measurement.
#' @param batch_result For the multivariate method, a \emph{list} of length
#'   \eqn{m}, each element the output of \code{batch_matrix()} for that
#'   measurement. All must share the same batch levels.
#'
#' @return For the multivariate method, a list with:
#' \describe{
#'   \item{\code{gamma_star}}{Named list over batches; each element is an
#'     \eqn{m \times G} matrix of shrunken means across measurements (rows)
#'     for each feature \eqn{g=1,\dots,G} (columns).}
#'   \item{\code{delta_star}}{Named list over batches; each element is a
#'     length-\eqn{G} list of \eqn{m \times m} covariance matrices (one per feature).}
#'   \item{\code{gamma_hat}, \code{delta_hat}, \code{mom}}{Pre-EB estimates and
#'     EB hyperparameters per batch (see Details of \code{\link{mom_calculation}}).}
#' }
#'
#' @details
#' The multivariate EB step updates, for each feature \eqn{g}, the mean vector
#' \eqn{\gamma_{\cdot g}} and covariance \eqn{\Delta_g} across measurements via
#' normal–inverse-Wishart-like shrinkage with hyperparameters estimated by
#' method of moments.
#'
#' @export

eb_algorithm.multivariate <- function(type, data_stand_result, batch_result, eb = TRUE, robust.LS = FALSE, ...){
  gamma_hat_list <- gamm_hat_gen(type, data_stand_result, batch_result)
  delta_hat_list <- delta_hat_gen(type, data_stand_result, batch_result, gamma_hat_list, robust.LS = robust.LS)
  m <- length(batch_result)
  if(eb){
    batch_level <- levels(batch_result[[1]]$batch_vector)
    eb_result <- lapply(batch_level, function(b){
      n_b <- batch_result[[1]]$n_batches[b]
      mom <- mom_calculation(type, gamma_hat_list[[b]], delta_hat_list[[b]])

      # adjust within batch
      bdat <- lapply(1:m, function(i) data_stand_result[[i]]$data_stand[batch_result[[i]]$batch_index[[b]],])
      g_orig <- gamma_hat_list[[b]]
      g_old  <- gamma_hat_list[[b]]
      d_old  <- delta_hat_list[[b]]
      change_old <- 1
      change <- 1
      count  <- 0
      while(change > 10e-5){
        eb_one_result <- eb_one_iteration(type, bdat, g_orig, g_old, d_old, mom)
        change <- eb_one_result$change
        if (count > 100) {
          if (change > change_old) {
            warning("Empirical Bayes step failed to converge after 100 iterations,
    	            using estimate before change between iterations increases.")
            break
          }
        }
        g_old <- eb_one_result$g_new
        d_old <- eb_one_result$d_new
        change_old <- eb_one_result$change
        count <- count+1
      }
      gamma_star <- eb_one_result$g_new
      delta_star <- eb_one_result$d_new
      return(list("gamma_star" = gamma_star, "delta_star" = delta_star, "mom" = mom))
    })
    gamma_star <- lapply(1:length(batch_level), function(i) eb_result[[i]]$gamma_star)
    delta_star <- lapply(1:length(batch_level), function(i) eb_result[[i]]$delta_star)
    mom <- lapply(1:length(batch_level), function(i) eb_result[[i]]$mom)
    names(gamma_star) <- names(delta_star) <- batch_level
  }else{
    gamma_star <- gamma_hat_list
    delta_star <- delta_hat_list
    mom <- NULL
  }
  return(list("gamma_star" = gamma_star, "delta_star" = delta_star, "gamma_hat" = gamma_hat_list, "delta_hat" = delta_hat_list, "mom" = mom))
}

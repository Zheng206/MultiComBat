#' Empirical-Bayes diagnostic distributions & plots
#'
#' Tools to compare the fitted **empirical distributions** of batch parameters
#' with their **prior** distributions used by the EB procedure.
#'
#' \describe{
#'   \item{\code{eb_check()}}{Builds tidy data frames of prior samples and empirical
#'   summaries for plotting.}
#'   \item{\code{eb_plot()}}{Convenience plots (densities/histograms) comparing
#'   prior vs empirical by batch, for either mean effects (\eqn{\gamma}) or
#'   dispersion/covariance (\eqn{\delta}).}
#' }
#'
#' @section Univariate EB:
#' For each batch \eqn{i}, \code{eb_check.univariate()} draws
#' \eqn{\gamma_{ig} \sim \mathcal{N}(g\_bar,\; g\_var)} and
#' \eqn{\delta_{ig} \sim \mathrm{InvGamma}(d\_a,\; d\_b)} (using
#' \code{MCMCpack::rinvgamma}'s \emph{shape/scale} parameterization),
#' and collates them with the empirical \code{gamma_hat} / \code{delta_hat}.
#'
#' @section Multivariate EB:
#' For each batch \eqn{i}, \code{eb_check.multivariate()} draws
#' \eqn{\boldsymbol{\gamma}_{ig} \sim \mathcal{N}_m(g\_bar,\; g\_var)}
#' via \code{MASS::mvrnorm} and
#' \eqn{\boldsymbol{\Sigma}_{ig} \sim \mathcal{IW}(\nu, \Psi)} via
#' \code{MCMCpack::riwish}. It then compares empirical \code{delta_hat} against
#' the prior by several distances (Frobenius, MSE, spectral, eigen-error, KL).
#'
#' @name eb_diagnostics
#' @keywords diagnostics visualization
NULL


#' @rdname eb_diagnostics
#' @param eb_result An EB result list returned by \code{eb_algorithm()}.
#' @param ... Passed to methods.
#' @export
eb_check <- function(eb_result, ...) {
  UseMethod("eb_check")
}

#' @rdname eb_diagnostics
#' @param eb_df Output produced by \code{eb_check()} for the corresponding mode.
#' @export
eb_plot <- function(eb_df, ...) {
  UseMethod("eb_plot")
}



#' @rdname eb_diagnostics
#' @return For \code{eb_check.univariate}, a list with element \code{eb_df}
#'   (tidy data frame with columns \code{batch}, \code{parm} in \{gamma,delta\},
#'   \code{type} in \{prior,emp.dist\}, and \code{dist} = sampled/empirical values).
#' @examples
#' \dontrun{
#' # suppose eb_uni is the output of eb_algorithm(type = combat_type("univariate"), ...)
#' chk <- eb_check(eb_uni)
#' p   <- eb_plot(chk)   # densities for gamma & delta
#' }
#' @method eb_check univariate
#' @export
eb_check.univariate <- function(eb_result, ...){
  batch_level <- rownames(eb_result$gamma_star)
  g_prior <- sapply(1:length(batch_level), function(i) rnorm(length(eb_result$gamma_hat[i,]), eb_result$mom[[i]]$g_bar, sqrt(eb_result$mom[[i]]$g_var)))
  colnames(g_prior) <- batch_level
  g_prior <- data.frame(g_prior) %>% pivot_longer(cols = everything(), names_to = "batch", values_to = "dist") %>% mutate(type = "prior", parm = "gamma")
  d_prior <- sapply(1:length(batch_level), function(i) rinvgamma(length(eb_result$gamma_hat[i,]), eb_result$mom[[i]]$d_a, scale = 1/eb_result$mom[[i]]$d_b))
  colnames(d_prior) <- batch_level
  d_prior <- data.frame(d_prior) %>% pivot_longer(cols = everything(), names_to = "batch", values_to = "dist") %>% mutate(type = "prior", parm = "delta")
  g_hat <- sapply(1:length(batch_level), function(i) eb_result$gamma_hat[i,])
  colnames(g_hat) <- batch_level
  g_hat <- data.frame(g_hat) %>% pivot_longer(cols = everything(), names_to = "batch", values_to = "dist") %>% mutate(type = "emp.dist", parm = "gamma")
  d_hat <- sapply(1:length(batch_level), function(i) eb_result$delta_hat[i,])
  colnames(d_hat) <- batch_level
  d_hat <- data.frame(d_hat) %>% pivot_longer(cols = everything(), names_to = "batch", values_to = "dist") %>% mutate(type = "emp.dist", parm = "delta")
  eb_df <- list(eb_df = rbind(g_prior, g_hat, d_prior, d_hat))
  class(eb_df) <- class(eb_result)
  return(eb_df)
}




#' @rdname eb_diagnostics
#' @return For \code{eb_plot.univariate}, a \code{ggplot} object with density
#'   overlays of prior vs empirical per batch, faceted by parameter
#'   (\code{gamma}, \code{delta}).
#' @method eb_plot univariate
#' @export
eb_plot.univariate <- function(eb_df, ...){
  ggplot(eb_df$eb_df, aes(x = dist, color = batch, linetype = type)) +
    geom_density(linewidth = 1.1, alpha = 0.8) +
    facet_wrap(~parm, scales = "free") +
    labs(
      title = "Empirical Bayes Prior vs Posterior by Batch",
      x = "Value",
      y = "Density",
      color = "Batch",
      linetype = "Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
}




#' @rdname eb_diagnostics
#' @return For \code{eb_check.multivariate}, a list with
#'   \code{eb_gamma} (tidy df of prior/empirical mean draws by measurement) and
#'   \code{eb_delta} (tidy df of distance metrics comparing covariance priors vs
#'   empirical estimates).
#' @method eb_check multivariate
#' @export
eb_check.multivariate <- function(eb_result, ...){
  m <- nrow(eb_result$gamma_hat[[1]])
  p <- ncol(eb_result$gamma_hat[[1]])
  batch_level <- names(eb_result$gamma_hat)
  g_prior <- lapply(1:length(batch_level), function(i){
    g_prior <- MASS::mvrnorm(n = p, eb_result$mom[[i]]$g_bar, eb_result$mom[[i]]$g_var)
    colnames(g_prior) <- paste0("M_", 1:m)
    g_prior <- g_prior %>% data.frame() %>% pivot_longer(cols = everything(), names_to = "measurement", values_to = "dist") %>% mutate(batch = batch_level[i], type = "prior")
  }) %>% bind_rows()

  g_hat <- lapply(1:length(batch_level), function(i){
    g_hat <- eb_result$gamma_hat[[i]] %>% t()
    colnames(g_hat) <- paste0("M_", 1:m)
    g_hat <- g_hat %>% data.frame() %>% pivot_longer(cols = everything(), names_to = "measurement", values_to = "dist") %>% mutate(batch = batch_level[i], type = "emp.dist")
  }) %>% bind_rows()

  eb_gamma <- rbind(g_prior, g_hat)

  eb_delta <- lapply(1:length(batch_level), function(i){
    Psi <- eb_result$mom[[i]]$d_b
    nu <- eb_result$mom[[i]]$d_a
    n_sim <- 1000
    d_prior <- array(NA, dim = c(m, m, n_sim))
    for (s in 1:n_sim) {
      d_prior[, , s] <- MCMCpack::riwish(nu, Psi)
    }

    prior_mean <- Psi * nu / (nu - m - 1)
    results <- lapply(eb_result$delta_hat[[i]], function(S_hat) {
      evaluate_cov_recovery(S_hat, prior_mean) %>% data.frame()
    }) %>% bind_rows() %>% pivot_longer(cols = everything(), names_to = "method", values_to = "value") %>% mutate(type = "emp.dist", batch = batch_level[i])

    prior_results <- lapply(1:n_sim, function(s) {
      evaluate_cov_recovery(d_prior[, , s], prior_mean) %>% data.frame()
    }) %>% bind_rows() %>% pivot_longer(cols = everything(), names_to = "method", values_to = "value") %>% mutate(type = "prior", batch = batch_level[i])

    eb_delta <- rbind(prior_results, results)
  }) %>% bind_rows() %>% mutate(batch = factor(batch))
  eb_df <- list(eb_gamma = eb_gamma, eb_delta = eb_delta)
  class(eb_df) <- class(eb_result)
  return(eb_df)
}



#' @rdname eb_diagnostics
#' @param param Character, either \code{"gamma"} (mean effects) or
#'   \code{"delta"} (covariance/variance diagnostics).
#' @param bat Character scalar; when \code{param = "delta"}, choose which batch
#'   to display (required if multiple batches are present).
#' @return For \code{eb_plot.multivariate}, a \code{ggplot} object:
#'   densities for \code{param = "gamma"} (faceted by measurement), or
#'   histograms of distance metrics for \code{param = "delta"} (faceted by metric).
#' @method eb_plot multivariate
#' @export
eb_plot.multivariate <- function(eb_df, param = "gamma", bat = NULL, ...){
  if(param == "gamma"){
    ggplot(eb_df$eb_gamma, aes(x = dist, color = batch, linetype = type)) +
      geom_density(linewidth = 1.1, alpha = 0.8) +
      facet_wrap(~measurement, scales = "free") +
      labs(
        title = "Empirical Bayes Prior vs Posterior by Batch (Gamma)",
        x = "Value",
        y = "Density",
        color = "Batch",
        linetype = "Type"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank()
      )
  }else{
    eb_delta <- eb_df$eb_delta %>% filter(batch == bat)
    ggplot(eb_delta, aes(x = value, fill = type)) +
      geom_histogram(
        data = filter(eb_delta, type == "prior"),
        aes(y = after_stat(density)),  # Use density for fair comparison
        bins = 80, alpha = 0.5, color = "gray30"
      ) +
      geom_histogram(
        data = filter(eb_delta, type == "emp.dist"),
        aes(y = after_stat(density)),
        bins = 80, alpha = 0.7, color = "navy"
      ) +
      scale_fill_manual(
        values = c("prior" = "gray80", "emp.dist" = "steelblue"),
        labels = c("Prior Samples", "Empirical")
      ) +
      facet_wrap(~ method, scales = "free", ncol = 2) +
      #facet_grid(batch ~ method, scales = "free"       # Allow axis scales to vary
      #) +
      labs(
        title = "Distribution of Distance Metrics: Empirical vs Prior",
        x = "Distance Value",
        y = "Density",
        fill = "Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank()
      )
  }
}

utils::globalVariables(c("dist", "type", "density"))

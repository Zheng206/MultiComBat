#' ComBat harmonization (univariate)
#'
#' End-to-end \emph{univariate} ComBat-style harmonization for an
#' \eqn{n \times G} matrix (rows = samples, columns = features), with optional
#' empirical Bayes (EB) shrinkage and robust scale estimation. A \emph{fully
#' Bayesian (Stan)} variant is also available via \code{stan = TRUE}.
#'
#' @section Full Bayesian (Stan) variant:
#' With \code{stan = TRUE}, the batch parameters \eqn{\gamma} (means) and
#' \eqn{\delta} (variances) are inferred from their joint posterior using MCMC.
#' Internally, \code{stan_data_prep()} builds the data list and
#' \code{stan_algorithm()} runs Stan and returns posterior summaries used for
#' harmonization (e.g., posterior means of \code{gamma_star} and
#' \code{delta_star}).
#'
#' @param bat Factor (or coercible) of length \code{n} giving batch labels.
#' @param data Numeric \eqn{n \times G} matrix/data frame of features (columns must be named).
#' @param covar \code{NULL} or an \eqn{n \times q} data frame of covariates. If \code{NULL},
#'   a null model \code{y ~ 1} is used (warning issued).
#' @param model Modeling function that accepts \code{formula} and \code{data}
#'   (e.g., \code{stats::lm}, \code{stats::glm}, \code{mgcv::gam}). Default \code{lm}.
#' @param formula RHS formula for covariates (excluding batch), e.g., \code{y ~ age + sex}.
#'   Required when \code{covar} is not \code{NULL}.
#' @param ref.batch Optional reference batch level. Predictions for standardization
#'   are taken under this batch; final output replaces those rows with the original values.
#' @param eb Logical; if \code{TRUE} (default) perform EB shrinkage of batch means/variances.
#' @param stan Logical; if \code{TRUE}, use the full Bayesian (Stan) routine
#'   (\code{stan_data_prep()}, \code{stan_algorithm()}) instead of EB.
#' @param robust.LS Logical; if \code{TRUE}, use robust scale estimators
#'   (\code{\link{biweight_midvar}}) in standardization/dispersion.
#' @param cov Logical; if \code{TRUE}, apply an additional CovBat-style covariance
#'   harmonization in PCA score space after ComBat mean/variance adjustment (default \code{FALSE}).
#' @param var_thresh Numeric in (0,1]; cumulative explained-variance threshold used
#'   to choose how many PC scores are harmonized when \code{cov = TRUE}
#'   (passed to \code{\link{pick_r_from_pc}}; default \code{0.95}).
#' @param min_rblock Integer; minimum number of PC scores to harmonize when \code{cov = TRUE}
#'   (default \code{1}).
#' @param max_rblock Integer; maximum number of PC scores to harmonize when \code{cov = TRUE}
#'   (default \code{Inf}).
#' @param ... Additional arguments forwarded to \code{model}.
#'
#' @return A list with:
#' \describe{
#'   \item{harm_data}{\eqn{n \times G} harmonized matrix on the original scale.}
#'   \item{resid}{Standardized, batch-adjusted residuals (\eqn{n \times G}).}
#'   \item{eb_result}{(when \code{stan = FALSE}) EB outputs: \code{gamma_star}, \code{delta_star},
#'                    and pre-EB summaries.}
#'   \item{stan_result}{(when \code{stan = TRUE}) Output of \code{stan_algorithm()}, typically
#'                      including posterior summaries for \code{gamma_star}/\code{delta_star}
#'                      and MCMC diagnostics.}
#' }
#'
#'
#' @examples
#' set.seed(1)
#' n <- 30; G <- 5
#' Y <- matrix(rnorm(n*G, sd = 2), n, G, dimnames = list(NULL, paste0("feat", 1:G)))
#' bat <- factor(sample(c("A","B","C"), n, TRUE))
#' X   <- data.frame(age = rnorm(n, 60, 8))
#' out <- com_harm(bat = bat, data = Y, covar = X, model = stats::lm,
#'                 formula = y ~ age, ref.batch = "A", eb = TRUE)
#' str(out$harm_data)
#'
#' @export
com_harm <- function(bat, data, covar, model = lm, formula = NULL, ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = FALSE, cov = FALSE, var_thresh = 0.95, min_rblock = 1, max_rblock = Inf, ...){
  feature <- colnames(data)
  batch_result <- batch_matrix(bat, ref.batch = ref.batch)
  fitted_model <- model_fitting(data, batch_result$batch_matrix, covar, model, formula, ...)
  data_stand_result <- standardize_data(fitted_model, batch_result, robust.LS = robust.LS)
  type <- combat_type(type = "univariate")
  if(!stan){
    eb_result <- eb_algorithm(type, data_stand_result, batch_result, eb = eb, robust.LS = robust.LS)
    class(eb_result) <- type
    data_nb <- data_stand_result$data_stand
    batch_levels <- rownames(eb_result$gamma_star)
    for (b in batch_levels){
      data_nb[batch_result$batch_index[[b]],] <- sweep(data_nb[batch_result$batch_index[[b]],, drop = FALSE], 2,
                                                       eb_result$gamma_star[b,], "-")
      data_nb[batch_result$batch_index[[b]],] <- sweep(data_nb[batch_result$batch_index[[b]],, drop = FALSE], 2,
                                                       sqrt(eb_result$delta_star[b,]), "/")
    }
  }else{
    stan_data <- stan_data_prep(type, data_stand_result, bat)
    batch_levels <- levels(bat)
    stan_result <- stan_algorithm(type, stan_data, batch_names = batch_levels)
    data_nb <- data_stand_result$data_stand
    G <- length(feature)
    for (b in batch_levels){
      data_nb[batch_result$batch_index[[b]],] <- sweep(data_nb[batch_result$batch_index[[b]],, drop = FALSE], 2,
                                                       stan_result$gamma_star[b,], "-")
      data_nb[batch_result$batch_index[[b]],] <- sweep(data_nb[batch_result$batch_index[[b]],, drop = FALSE], 2,
                                                       sqrt(stan_result$delta_star[b,]), "/")
    }
  }
  data_combat <- data_nb*data_stand_result$sd_mat + data_stand_result$stand_mean
  if (!is.null(batch_result$ref)) {
    data_combat[batch_result$ref,] <- data[batch_result$ref,]
  }

  if(cov){
    R <- data_combat - data_stand_result$stand_mean
    R_star <- covbat(R, site = bat, center = TRUE, scale_scores = TRUE, var_thresh = var_thresh, min_rblock = min_rblock,
                     max_rblock = max_rblock,
                     ref.batch = ref.batch)
    data_combat <- R_star + data_stand_result$stand_mean
  }

  if(!stan){
    return(list(harm_data = data_combat, eb_result = eb_result, resid = data_nb))
  }else{
    return(list(harm_data = data_combat, stan_result = stan_result, resid = data_nb))
  }
}

#' ComBat harmonization across measurements (multivariate)
#'
#' End-to-end \emph{multivariate} ComBat-style harmonization for a list of
#' \eqn{m} measurements/modalities. For each feature, the batch effect is a mean
#' vector across measurements and an \eqn{m \times m} covariance; these are
#' shrunken via EB (or inferred via a \emph{fully Bayesian (Stan)} model),
#' residuals are whitened with \code{\link{sigma_inv_sqrt}}, and then mapped
#' back to the original scale.
#'
#'
#' @section Full Bayesian (Stan) variant:
#' With \code{stan = TRUE}, the per-feature mean vectors and covariance matrices
#' across measurements are sampled from their joint posterior (e.g., multivariate normal /
#' LKJ + half-t) using MCMC. \code{stan_data_prep()} prepares the
#' data; \code{stan_algorithm()} returns posterior summaries for the
#' harmonization step.
#'
#' @param bat List of length \eqn{m}; each element is a factor (or coercible) of
#'   length \code{n} with batch labels for that measurement. All must share the
#'   same batch levels.
#' @param data List of length \eqn{m}; each element is an \eqn{n \times G}
#'   matrix/data frame. All share identical feature columns (names and order) and
#'   the same sample order (rows).
#' @param covar List of length \eqn{m} of covariate data frames (each \eqn{n \times q_i}),
#'   or \code{NULL} to fit null models.
#' @param model Modeling function accepting \code{formula} and \code{data}
#'   (e.g., \code{stats::lm}, \code{mgcv::gam}). Default \code{lm}.
#' @param formula RHS formula for covariates (excluding batch), e.g., \code{y ~ age + sex}.
#' @param ref.batch Optional reference batch level shared across measurements.
#' @param eb Logical; if \code{TRUE} (default) run multivariate EB shrinkage.
#' @param stan Logical; if \code{TRUE}, use the full Bayesian (Stan) multivariate routine.
#' @param robust.LS Logical; if \code{TRUE}, use robust estimators
#'   (\code{\link{biweight_midvar}}, \code{\link{biweight_midvar_mul}}).
#' @param cov Logical; if \code{TRUE}, apply an additional CovBat-style PCA
#'   score-space coupling across measurements after the main multivariate step
#'   (default \code{FALSE}).
#' @param var_thresh Numeric in (0,1]; cumulative explained-variance threshold used
#'   to choose how many PC scores per measurement are used in the stage-2 coupling
#'   (default \code{0.95}).
#' @param min_rblock Integer; minimum number of PC scores per measurement used in
#'   stage-2 coupling (default \code{1}).
#' @param max_rblock Integer; maximum number of PC scores per measurement used in
#'   stage-2 coupling (default \code{Inf}).
#' @param ... Additional arguments forwarded to \code{model}.
#'
#' @return A list with:
#' \describe{
#'   \item{harm_data}{List of length \eqn{m}; each element is an \eqn{n \times G}
#'                    harmonized matrix on the original scale.}
#'   \item{resid}{List of length \eqn{m} of standardized, batch-adjusted residuals.}
#'   \item{eb_result}{(when \code{stan = FALSE}) Multivariate EB outputs:
#'                    \code{gamma_star} (per-batch \eqn{m \times G} means) and
#'                    \code{delta_star} (per-batch length-\eqn{G} list of \eqn{m \times m} covariances).}
#'   \item{stan_result}{(when \code{stan = TRUE}) Output of \code{stan_algorithm()}, typically
#'                      including posterior summaries for \code{gamma_star}/\code{delta_star}
#'                      and MCMC diagnostics.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 30; G <- 5
#' Y1 <- matrix(rnorm(n*G), n, G, dimnames = list(NULL, paste0("feat", 1:G)))
#' Y2 <- Y1 + matrix(rnorm(n*G, sd = 0.5), n, G, dimnames = dimnames(Y1))
#' bat <- factor(sample(c("A","B","C"), n, TRUE))
#' bats <- lapply(1:2, function(i) bat)
#' covs <- list(data.frame(age = rnorm(n)), data.frame(age = rnorm(n)))
#' out_m <- com_harm.multivariate(
#'   bat = bats, data = list(Y1, Y2), covar = covs,
#'   model = lm, formula = y ~ age, ref.batch = "A", eb = FALSE
#' )
#' length(out_m$harm_data)
#' }
#'
#' @export
com_harm.multivariate <- function(bat, data, covar, model = lm, formula = NULL, ref.batch = NULL, eb = TRUE, stan = FALSE, robust.LS = FALSE, cov = FALSE, var_thresh = 0.95, min_rblock = 1, max_rblock = Inf, ...){
  m <- length(data)
  data <- lapply(1:m, function(i) data.frame(data[[i]]))
  feature <- colnames(data[[1]])
  batch_result <- lapply(1:m, function(i) batch_matrix(bat[[i]], ref.batch = ref.batch))
  fitted_model <- lapply(1:m, function(i) model_fitting(data[[i]], batch_result[[i]]$batch_matrix, covar[[i]], model, formula, ...))
  data_stand_result <- lapply(1:m, function(i) standardize_data(fitted_model[[i]], batch_result[[i]], robust.LS = robust.LS))
  type <- combat_type(type = "multivariate")
  if(!stan){
    eb_result <- eb_algorithm(type, data_stand_result, batch_result, eb = eb, robust.LS = robust.LS)
    class(eb_result) <- type
    data_nb <- lapply(1:m, function(i) data_stand_result[[i]]$data_stand)
    batch_levels <- names(eb_result$gamma_star)
    G <- length(feature)
    for (b in batch_levels){
      for(i in 1:m){
        data_nb[[i]][batch_result[[i]]$batch_index[[b]],] <- sweep(data_nb[[i]][batch_result[[i]]$batch_index[[b]],, drop = FALSE], 2,
                                                                   eb_result$gamma_star[[b]][i,], "-")
      }

      gene_list <- lapply(1:G, function(g){
        data_nb_list <- lapply(1:m, function(i) data_nb[[i]][batch_result[[i]]$batch_index[[b]],])
        g_matrix <- multi_matrix_transform(data_nb_list, g)
        chol_inv_delta <- sigma_inv_sqrt(eb_result$delta_star[[b]][[g]])
        data_batch <- chol_inv_delta %*% g_matrix
        batch_corrected_resi <- multi_matrix_transform_inverse(data_batch, feature[g])
        return(batch_corrected_resi)
      })

      resid_df_list <- lapply(1:m, function(i) {
        sub_list <- lapply(1:G, function(g) gene_list[[g]][[i]])
        sub_df <- do.call(cbind, sub_list)
        return(sub_df)
      })

      for(i in 1:m){
        data_nb[[i]][batch_result[[i]]$batch_index[[b]],] <- resid_df_list[[i]]
      }
    }
  }else{
    stan_data <- stan_data_prep(type, data_stand_result, bat)
    batch_levels <- levels(bat[[1]])
    stan_result <- stan_algorithm(type, stan_data, batch_names = batch_levels)
    data_nb <- lapply(1:m, function(i) data_stand_result[[i]]$data_stand)
    G <- length(feature)
    for (b in batch_levels){
      for(i in 1:m){
        data_nb[[i]][batch_result[[i]]$batch_index[[b]],] <- sweep(data_nb[[i]][batch_result[[i]]$batch_index[[b]],, drop = FALSE], 2,
                                                                   stan_result$gamma_star[[b]][i,], "-")
      }

      gene_list <- lapply(1:G, function(g){
        data_nb_list <- lapply(1:m, function(i) data_nb[[i]][batch_result[[i]]$batch_index[[b]],])
        g_matrix <- multi_matrix_transform(data_nb_list, g)
        chol_inv_delta <- sigma_inv_sqrt(stan_result$delta_star[[b]][[g]])
        data_batch <- chol_inv_delta %*% g_matrix
        batch_corrected_resi <- multi_matrix_transform_inverse(data_batch, feature[g])
        return(batch_corrected_resi)
      })

      resid_df_list <- lapply(1:m, function(i) {
        sub_list <- lapply(1:G, function(g) gene_list[[g]][[i]])
        sub_df <- do.call(cbind, sub_list)
        return(sub_df)
      })

      for(i in 1:m){
        data_nb[[i]][batch_result[[i]]$batch_index[[b]],] <- resid_df_list[[i]]
      }
    }
  }
  data_combat <- lapply(1:m, function(i) data_nb[[i]]*data_stand_result[[i]]$sd_mat + data_stand_result[[i]]$stand_mean)

  if (!is.null(batch_result[[1]]$ref)) {
    lapply(1:m, function(i) data_combat[[i]][batch_result[[i]]$ref,] <- data[[i]][batch_result[[i]]$ref,])
  }

  if(cov){
    R_list <- lapply(1:m, function(i) data_combat[[i]] - data_stand_result[[i]]$stand_mean)
    R_h_list <- stage2_variantA(R_list, site = bat, var_thresh = var_thresh, min_rblock = min_rblock, max_rblock = max_rblock)$R_h_list
    data_combat <- lapply(1:m, function(i) R_h_list[[i]] + data_stand_result[[i]]$stand_mean)
  }

  if(!stan){
    return(list(harm_data = data_combat, eb_result = eb_result, resid = data_nb))
  }else{
    return(list(harm_data = data_combat, stan_result = stan_result, resid = data_nb))
  }
}

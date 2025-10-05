#' Build batch design matrix and index helpers
#'
#' Constructs a one-hot batch design matrix, per-batch indices, counts, and an
#' optional logical indicator for a reference batch.
#'
#' @param bat A **factor** vector of batch labels. Unused levels are dropped.
#' @param ref.batch Optional length-1 character or factor level naming the
#'   reference batch. Must be one of \code{levels(bat)}; otherwise an error is thrown.
#'
#' @return A list with components:
#' \describe{
#'   \item{batch_vector}{The (dropped-levels) batch factor.}
#'   \item{batch_matrix}{An \eqn{n \times B} model matrix of batch dummies
#'     with no intercept (\code{model.matrix(~ -1 + bat)}).}
#'   \item{batch_index}{A named list of integer indices for each batch level.}
#'   \item{n_batches}{A named integer vector: counts per batch level.}
#'   \item{ref}{A logical vector of length \eqn{n} marking the reference batch,
#'     or \code{NULL} if \code{ref.batch} is not supplied.}
#' }
#'
#' @examples
#' set.seed(1)
#' bat <- factor(sample(c("siteA", "siteB", "siteC"), 10, replace = TRUE))
#'
#' # With a reference batch
#' out <- batch_matrix(bat, ref.batch = "siteA")
#' names(out)
#' head(out$batch_matrix)
#' out$n_batches
#' table(out$ref)
#'
#' # Without a reference batch
#' out2 <- batch_matrix(bat)
#' str(out2$batch_index)
#'
#' @export
batch_matrix <- function(bat, ref.batch = NULL){
  bat <- as.factor(bat)
  bat <- droplevels(bat)
  batch <- model.matrix(~ -1 + bat)
  batch_index <- lapply(levels(bat), function(x) which(bat == x))
  names(batch_index) <- levels(bat)
  n_batches <- sapply(batch_index, length)

  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(bat))) {
      stop("Reference batch must be in the batch levels")
    }
    ref <- bat == ref.batch
  }else{ref <- NULL}
  return(list("batch_vector" = bat, "batch_matrix" = batch, "batch_index" = batch_index, "n_batches" = n_batches, "ref" = ref))
}

#' Fit a per-defined model with batch effects
#'
#' Fits the user-specified model to each column of \code{data}, augmenting the
#' model with batch indicators (no intercept) and optional covariates.
#'
#' @param data A numeric matrix or data frame of responses with observations in rows
#'   and features (e.g., ROIs, genes) in columns. Each column is modeled separately.
#' @param batch A factor (recommended) or vector coercible to factor giving the batch
#'   label for each observation (length must match \code{nrow(data)}).
#' @param covar Optional data frame (or list/data.frame-coercible object)
#'   of covariates with \code{nrow(covar) == nrow(data)}. If \code{NULL}, a null
#'   model (intercept only) is fitted and a warning is issued.
#' @param model A modeling function (e.g., \code{\link[stats]{lm}},
#'   \code{\link[lme4]{lmer}}, \code{mgcv::gam}). Must accept a
#'   \code{formula} and \code{data} argument.
#' @param formula A right-hand-side formula for the covariate part of the model
#'   (excluding batch). For example, \code{y ~ age + sex}. If \code{covar} is
#'   provided and \code{formula} is \code{NULL}, the function stops with an error.
#'   When \code{covar} is \code{NULL}, a null model \code{y ~ 1} is used.
#' @param ... Additional arguments passed on to \code{model}.
#'
#'
#' @return A list with components:
#' \describe{
#'   \item{data}{The input \code{data}.}
#'   \item{mod}{The constructed data frame of predictors (covariates and batch).}
#'   \item{formula}{The base formula used before adding batch terms.}
#'   \item{fits}{A named list (one per column of \code{data}) of fitted model objects.}
#' }
#'
#' @section Warnings/Errors:
#' \itemize{
#'   \item If \code{covar} is \code{NULL}, a null model is used and a warning is issued.
#'   \item If \code{covar} is provided but \code{formula} is \code{NULL}, the function errors.
#'   \item \code{batch} is expected to be a valid factor aligning with \code{nrow(data)}.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 5
#' Y <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("feat", 1:p)))
#' bat <- factor(sample(c("siteA", "siteB", "siteC"), n, TRUE))
#'
#' # Null model (no covariates)
#' fit0 <- model_fitting(Y, batch = bat, covar = NULL, model = lm)
#' names(fit0$fits)
#' summary(fit0$fits[[1]])
#'
#' # With covariates and a formula
#' covars <- data.frame(age = rnorm(n, 60, 10), sex = factor(sample(c("F","M"), n, TRUE)))
#' out <- batch_matrix(bat, ref.batch = "siteA")
#' fit1 <- model_fitting(
#'   data   = Y,
#'   batch  = out$batch_matrix,
#'   covar  = covars,
#'   model  = lm,
#'   formula = y ~ age + sex
#' )
#' coef(fit1$fits[[1]])
#'
#' @export

model_fitting <- function(data, batch, covar, model, formula = NULL, ...){

  if (is.null(covar)) {
    mod <- data.frame(I(batch))
    formula <- y ~ 1
    warning("No covariates were provided. A null model has been fitted.")
  } else {
    mod <- data.frame(covar, I(batch))
    if(is.null(formula)){
      stop("Please provide a formula to fit the model!")
    }
  }

  fits <- apply(data, 2, function(y) {
    dat <- data.frame(y = y, mod)
    bat_formula <- update(formula, ~ . + batch + -1)
    do.call(model, list(formula = bat_formula, data = dat, ...))
  })
  names(fits) <- colnames(data)
  return(list("data" = data, "mod" = mod, "formula" = formula, "fits" = fits))
}


#' Standardize features by removing fitted batch/covariate means and pooling variance
#'
#' Uses fitted models (from \code{model_fitting()}) and batch metadata (from
#' \code{batch_matrix()}) to compute per-feature standardized values:
#' \deqn{Z = (Y - \hat{\mu}_{\mathrm{stand}}) / \hat{\sigma}}.
#' The standardization mean \eqn{\hat{\mu}_{\mathrm{stand}}} is obtained by predicting
#' from the fitted models under either (i) a mixture batch profile given by batch
#' proportions (default), or (ii) a designated reference batch if provided.
#' Variance is pooled across observations (optionally using a robust scale).
#'
#' @param model A list returned by \code{\link{model_fitting}} containing
#'   \code{$data}, \code{$mod}, \code{$formula}, and \code{$fits}.
#' @param batch_result A list returned by \code{\link{batch_matrix}} containing
#'   \code{$batch_vector}, \code{$n_batches}, and optionally \code{$ref}.
#' @param robust.LS Logical; if \code{TRUE}, use a robust scale estimator
#'   (\code{biweight_midvar}) for the pooled variance; otherwise use \code{\link[stats]{var}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{data_stand}{An \eqn{n \times p} matrix of standardized values.}
#'   \item{stand_mean}{An \eqn{n \times p} matrix of standardization means used in the numerator.}
#'   \item{sd_mat}{An \eqn{n \times p} matrix of per-feature pooled standard deviations.}
#' }
#'
#'
#' @examples
#' set.seed(1)
#' n <- 80; p <- 6
#' Y <- matrix(rnorm(n * p, sd = 2), n, p,
#'             dimnames = list(NULL, paste0("feat", 1:p)))
#' bat <- factor(sample(c("A","B","C"), n, TRUE, prob = c(.4,.4,.2)))
#' X   <- data.frame(age = rnorm(n, 60, 8),
#'                   sex = factor(sample(c("F","M"), n, TRUE)))
#'
#' # Build batch helpers and fit per-feature models
#' bres <- batch_matrix(bat)
#' mf   <- model_fitting(
#'   data    = Y,
#'   batch   = bres$batch_matrix,
#'   covar   = X,
#'   model   = lm,
#'   formula = y ~ age + sex
#' )
#'
#' # Standardize using mixture-of-batches mean and pooled variance
#' out1 <- standardize_data(mf, bres, robust.LS = FALSE)
#' dim(out1$data_stand); colMeans(out1$stand_mean)
#'
#' @export

standardize_data <- function(model, batch_result, robust.LS = FALSE){
  pmod <- model$mod
  n = length(batch_result$batch_vector)
  pmod$batch[] <- matrix(batch_result$n_batches/n, n, nlevels(batch_result$batch_vector), byrow = TRUE)
  if (!is.null(batch_result$ref)) {
    pmod$batch[] <- 0
    pmod$batch[,batch_result$batch_vector[which(batch_result$ref == 1)] |> unique()] <- 1
  }
  stand_mean <- sapply(model$fits, predict, newdata = pmod, type = "response")
  resid_mean <- sapply(model$fits, predict, newdata = model$mod, type = "response")

  if(robust.LS){
    scale = biweight_midvar
  }else{
    scale = var
  }

  if (!is.null(batch_result$ref)) {
    nref <- sum(batch_result$ref)
    var_pooled <- apply((model$data - resid_mean)[batch_result$ref, , drop = FALSE], 2, scale) *
      (nref - 1)/nref
  } else {
    var_pooled <- apply(model$data - resid_mean, 2, scale) * (n - 1)/n
  }
  sd_mat <- sapply(sqrt(var_pooled), rep, n)
  data_stand <- (model$data-stand_mean)/sd_mat
  return(list("data_stand" = data_stand, "stand_mean" = stand_mean, "sd_mat" = sd_mat))
}

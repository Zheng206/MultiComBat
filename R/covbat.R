#' Pick PC count by cumulative explained variance
#'
#' Chooses the number of principal components whose cumulative explained
#' variance first reaches (or exceeds) `var_thresh`, then clamps the choice
#' to `[min_r, max_r]` and to `length(pca_fit$sdev)`.
#'
#' @param pca_fit A PCA object (e.g., from [stats::prcomp()] or any object
#'   with a numeric vector `sdev` of singular values / SDs).
#' @param var_thresh Numeric in (0,1]; target cumulative variance threshold
#'   (default `0.95`).
#' @param min_r Integer minimum allowed PCs (default `5`).
#' @param max_r Integer maximum allowed PCs (default `Inf`).
#'
#' @return Integer number of PCs to retain.
#'
#' @examples
#' set.seed(1)
#' X <- scale(matrix(rnorm(200), 50, 4))
#' pf <- prcomp(X, center = TRUE, scale. = FALSE)
#' pick_r_from_pc(pf, var_thresh = 0.9, min_r = 1)
#'
#' @export
pick_r_from_pc <- function(pca_fit, var_thresh = 0.95, min_r = 5, max_r = Inf) {
  cumve <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  r <- which(cumve >= var_thresh)[1]
  r <- max(min_r, r)
  r <- min(r, length(pca_fit$sdev))
  r <- min(r, max_r)
  return(r)
}


#' Score-space coupling across multiple blocks
#'
#' For a list of blocks (same subjects/rows), performs:
#' 1) per-block PCA (keeping a data-driven number of PCs per block),
#' 2) concatenates standardized block scores and runs a second PCA to extract
#'    *shared* subject scores `G`,
#' 3) projects each block's scores onto the shared space to get block
#'    loadings in score-space `A_t` and residual scores `H_t`.
#'
#' This is a simple, PCA-based coupling of blocks that exposes common
#' variation across measurements while allowing block-specific residuals.
#'
#' @param R_list List of length \eqn{m}; each element an \eqn{n \times p_t}
#'   numeric matrix/data frame (rows = subjects, columns = variables for block *t*).
#'   All elements **must** have the same number of rows.
#' @param center Logical; center columns before PCA in each block (default `TRUE`).
#' @param scale_scores Logical; rescale per-block retained scores to unit SD
#'   before concatenation (default `TRUE`).
#' @param var_thresh Numeric in (0,1]; per-PCA cumulative variance threshold used
#'   by [pick_r_from_pc()] to choose retained dimensions (default `0.95`).
#' @param min_rblock Integer; minimum PCs per block (default `5`).
#' @param max_rblock Integer; maximum PCs per block (default `Inf`).
#'
#' @return A list with components:
#' \describe{
#'   \item{G}{\eqn{n \times r_\mathrm{shared}} matrix of shared subject scores from
#'            the second-stage PCA on concatenated scores.}
#'   \item{A_list}{List of block-specific projection matrices \eqn{A_t} of size
#'            \eqn{r_\mathrm{shared} \times r_t} such that
#'            \eqn{F^{(t)}_\mathrm{std} \approx G A_t + H_t}.}
#'   \item{H_list}{List of residual score matrices \eqn{H_t} (\eqn{n \times r_t}).}
#'   \item{F_list}{List of per-block retained scores (\eqn{n \times r_t}).}
#'   \item{F_list_std}{Same as `F_list` but column-standardized if `scale_scores=TRUE`.}
#'   \item{L_list}{List of per-block PCA loadings (\eqn{p_t \times r_t}).}
#'   \item{sdev_list}{List of per-block retained singular values (length `r_t`).}
#'   \item{sc_list}{List of per-block score standardization factors (length `r_t`).}
#'   \item{meta}{List with `n`, `m`, `r_shared`, and flags.}
#' }
#'
#' @seealso [stats::prcomp()], [pick_r_from_pc()], [reconstruct_block()]
#'
#' @examples
#' set.seed(2)
#' n <- 60
#' R1 <- scale(matrix(rnorm(n*30), n, 30))
#' R2 <- scale(matrix(rnorm(n*20), n, 20))
#' fit <- score_space_coupling(list(R1, R2), var_thresh = 0.9, min_rblock = 3)
#' str(fit$meta)
#'
#' @export
score_space_coupling <- function(R_list,
                                 center = TRUE,
                                 scale_scores = TRUE, var_thresh = 0.95, min_rblock = 5, max_rblock = Inf) {
  stopifnot(is.list(R_list), length(R_list) >= 1)
  n_rows <- vapply(R_list, nrow, 1L)
  if (length(unique(n_rows)) != 1L) stop("All blocks must have same number of rows (subjects).")
  n <- n_rows[1]; m <- length(R_list)

  r_block_vec <- integer(m)
  prefits <- vector("list", m)
  for (t in seq_len(m)) {
    pf <- prcomp(R_list[[t]], center = center, scale. = FALSE)
    prefits[[t]] <- pf
    r_block_vec[t] <- pick_r_from_pc(pf, var_thresh, min_r = min_rblock, max_r = max_rblock)
  }

  F_list <- vector("list", m)
  L_list <- vector("list", m)
  sdev_list <- vector("list", m)

  for (t in seq_len(m)) {
    keep <- seq_len(r_block_vec[t])
    F_list[[t]]   <- prefits[[t]]$x[, keep, drop = FALSE]
    L_list[[t]]   <- prefits[[t]]$rotation[, keep, drop = FALSE]
    sdev_list[[t]]<- prefits[[t]]$sdev[keep]
  }

  F_list_std <- F_list
  sc_list <- vector("list", m)
  if (scale_scores) {
    for (t in seq_len(m)) {
      sc <- apply(F_list[[t]], 2, sd)
      sc[sc == 0] <- 1
      sc_list[[t]] <- sc
      F_list_std[[t]] <- sweep(F_list[[t]], 2, sc, "/")
    }
  }else {
    for (t in seq_len(m)) sc_list[[t]] <- rep(1, ncol(F_list[[t]]))
  }

  # Concatenate scores and get shared subject scores G (top r_shared PCs of F_concat)
  F_concat <- do.call(cbind, F_list_std)
  pc_concat <- prcomp(F_concat, center = TRUE, scale. = FALSE)
  r_shared <- pick_r_from_pc(pc_concat, var_thresh, min_r = 1)
  G <- pc_concat$x[, seq_len(r_shared), drop = FALSE]

  # Projection matrices per block: A_t and residual scores H_t
  GTG_inv <- solve(crossprod(G))
  A_list <- vector("list", m)
  H_list <- vector("list", m)
  for (t in seq_len(m)) {
    A_t <- GTG_inv %*% crossprod(G, F_list_std[[t]])
    H_t <- F_list_std[[t]] - G %*% A_t
    A_list[[t]] <- A_t
    H_list[[t]] <- H_t
  }

  list(
    G = G,
    A_list = A_list,
    H_list = H_list,
    F_list = F_list,
    F_list_std = F_list_std,
    L_list = L_list,
    sdev_list = sdev_list,
    sc_list = sc_list,
    meta = list(n = n, m = m, r_shared = r_shared,
                center = center, scale_scores = scale_scores)
  )
}


#' Reconstruct a block back to variable space
#'
#' Given shared scores `G`, a block's score-space projection `A_t`, residual
#' scores `H_t`, and the block loadings `L_t`, reconstructs the block in the
#' original variable space. If score standardization was applied, pass `sc_t`
#' to undo it.
#'
#' @param G \eqn{n \times r_\mathrm{shared}} matrix of shared scores.
#' @param A_t \eqn{r_\mathrm{shared} \times r_t} projection matrix for block *t*.
#' @param H_t \eqn{n \times r_t} residual scores for block *t*.
#' @param L_t \eqn{p_t \times r_t} PCA loadings (variables by retained PCs) for block *t*.
#' @param sc_t Optional length-\eqn{r_t} vector of score scaling factors used
#'   during coupling (if scores were standardized).
#'
#' @return \eqn{n \times p_t} numeric matrix: reconstructed block *t*.
#'
#' @seealso [score_space_coupling()]
#'
#' @examples
#' set.seed(2)
#' n <- 60
#' R1 <- scale(matrix(rnorm(n*30), n, 30))
#' R2 <- scale(matrix(rnorm(n*20), n, 20))
#' fit <- score_space_coupling(list(R1, R2), var_thresh = 0.9, min_rblock = 3)
#' result <- reconstruct_block(fit$G, fit$A_list[[1]], fit$H_list[[1]], fit$L_list[[1]],
#' fit$sc_list[[1]])
#' @export
reconstruct_block <- function(G, A_t, H_t, L_t, sc_t = NULL) {
  F_t_star_std <- G %*% A_t + H_t
  if (!is.null(sc_t)) {
    F_t_star <- F_t_star_std %*% diag(sc_t)
  } else {
    F_t_star <- F_t_star_std
  }
  F_t_star %*% t(L_t)
}

#' Stage-2 variant A: harmonize shared scores, then reconstruct blocks
#'
#' A two-stage pipeline for multiple blocks:
#' \enumerate{
#'   \item Couple blocks in score space via [score_space_coupling()] to obtain
#'         shared subject scores \eqn{G}.
#'   \item Harmonize \eqn{G} across sites with [com_harm()] (univariate),
#'         yielding \eqn{G^\*}.
#'   \item Reconstruct each block using \eqn{G^\*}, the block projection
#'         \eqn{A_t}, residual scores \eqn{H_t}, and loadings \eqn{L_t}.
#' }
#'
#' @param R_list List of \eqn{m} blocks (\eqn{n \times p_t} each), same rows/subjects.
#' @param site List of \eqn{m} site/batch factors (one per block). Only the
#'   first element is used for harmonizing `G` in this variant.
#' @param center,scale_scores,var_thresh,min_rblock,max_rblock
#'   Passed to [score_space_coupling()] (see its documentation).
#'
#' @return A list with:
#' \describe{
#'   \item{R_h_list}{List of harmonized blocks (same shapes as inputs).}
#'   \item{coupling}{The object returned by [score_space_coupling()].}
#'   \item{G_star}{Harmonized shared scores \eqn{G^\*}.}
#'   \item{meta}{List with `r_shared`, `center`, `scale_scores`.}
#' }
#'
#' @seealso [score_space_coupling()], [reconstruct_block()], [com_harm()]
#'
#' @examples
#' \dontrun{
#' set.seed(3)
#' n <- 50
#' R1 <- scale(matrix(rnorm(n*25), n, 25))
#' R2 <- scale(matrix(rnorm(n*15), n, 15))
#' site <- list(factor(rep(LETTERS[1:2], each = n/2)),
#'              factor(rep(LETTERS[1:2], each = n/2)))
#' out <- stage2_variantA(list(R1, R2), site, var_thresh = 0.9, min_rblock = 2)
#' lapply(out$R_h_list, dim)
#' }
#'
#' @export
stage2_variantA <- function(R_list, site,
                            center = TRUE,
                            scale_scores = TRUE, var_thresh = 0.8, min_rblock = 1, max_rblock = Inf) {
  stopifnot(is.list(R_list), length(R_list) >= 1)
  n <- nrow(R_list[[1]]); m <- length(R_list)
  if (length(unique(vapply(R_list, nrow, 1L))) != 1L) {
    stop("All blocks must have same number of rows (subjects).")
  }
  site <- lapply(1:length(site), function(i) as.factor(site[[i]]))
  fit <- score_space_coupling(R_list, center = center, scale_scores = scale_scores, var_thresh = var_thresh, min_rblock = min_rblock, max_rblock = max_rblock)
  G_star <- com_harm(bat = site[[1]], data = fit$G, covar = NULL, cov = FALSE)$harm_data
  R_h_list <- lapply(seq_len(m), function(t){
    result <- reconstruct_block(G = G_star, A_t = fit$A_list[[t]], H_t = fit$H_list[[t]], L_t = fit$L_list[[t]], sc_t = fit$sc_list[[t]])
    colnames(result) <- colnames(R_list[[t]])
    rownames(result) <- rownames(R_list[[t]])
    return(data.frame(result))
  })
  list(
    R_h_list = R_h_list,
    coupling = fit,
    G_star = G_star,
    meta = list(r_shared = fit$meta$r_shared,
                center = center, scale_scores = scale_scores)
  )
}

#' CovBat-style harmonization of a single block via PCA scores
#'
#' Implements a simple CovBat pipeline for one block:
#' \enumerate{
#'   \item PCA on `R`, keep `r` PCs so that cumulative variance \eqn{\ge} `var_thresh`.
#'   \item Harmonize the first `r` PC scores across `site` with [com_harm()] (univariate EB off).
#'   \item Reconstruct the data using harmonized first `r` scores and original
#'         remaining scores, then invert the PCA transform.
#' }
#'
#' @param R Numeric matrix/data frame \eqn{n \times p} (rows = subjects).
#' @param site Factor (or coercible) of length \eqn{n} with site/batch labels.
#' @param center Logical; pass to [stats::prcomp()] (default `TRUE`).
#' @param scale_scores Logical; if `TRUE`, scale variables before PCA score
#'   harmonization (passed as `scale.` to [stats::prcomp()]) (default `TRUE`).
#' @param var_thresh,min_rblock,max_rblock See [pick_r_from_pc()]; determine
#'   how many PC scores are harmonized.
#' @param ref.batch Optional reference batch for [com_harm()] (univariate).
#'
#' @return \eqn{n \times p} numeric matrix: harmonized data on the original scale.
#'
#' @seealso [pick_r_from_pc()], [com_harm()], [stats::prcomp()]
#'
#' @examples
#' set.seed(4)
#' n <- 40; p <- 30
#' R <- scale(matrix(rnorm(n*p), n, p))
#' site <- factor(rep(LETTERS[1:2], each = n/2))
#' R_h <- covbat(R, site, var_thresh = 0.9, min_rblock = 2)
#' dim(R_h)
#'
#' @export
covbat <- function(R, site, center = TRUE,
                   scale_scores = TRUE, var_thresh = 0.95, min_rblock = 1, max_rblock = Inf,
                   ref.batch = NULL){
  R <- as.matrix(R)
  n <- nrow(R)
  p <- ncol(R)
  site <- as.factor(site)
  pf <- prcomp(R, center = center, scale. = scale_scores)
  r <- pick_r_from_pc(pf, var_thresh, min_r = min_rblock, max_r = max_rblock)
  scores <- pf$x[,1:r]
  scores_com <- com_harm(bat = site, data = scores, covar = NULL, eb = FALSE, ref.batch = ref.batch, cov = FALSE)
  full_scores <- pf$x
  full_scores[,1:r] <- scores_com$harm_data
  if (scale_scores) {
    data_covbat <- full_scores %*% t(pf$rotation) *
      matrix(pf$scale, n, p, byrow = TRUE) +
      matrix(pf$center, n, p, byrow = TRUE)
  } else {
    data_covbat <- full_scores %*% t(pf$rotation) +
      matrix(pf$center, n, p, byrow = TRUE)
  }
  return(data_covbat)
}


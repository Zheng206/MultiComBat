#' Reshape a single feature across measurements (and back)
#'
#' Helpers to (a) extract one feature/column \eqn{g} from a list of matrices
#' \eqn{Z_1,\dots,Z_m} and stack it into an \eqn{m \times n} matrix
#' (rows = measurements, cols = samples), and (b) invert that operation.
#'
#' @param Z_list A list of length \eqn{m}. Each element is an \eqn{n \times G}
#'   numeric matrix/data frame (rows = samples, columns = features), all with the
#'   same \eqn{n} and \eqn{G}, and columns ordered consistently.
#' @param g Integer column index of the feature to extract (1…\eqn{G}).
#'
#' @return \code{multi_matrix_transform()} returns an \eqn{m \times n} numeric
#'   matrix whose \eqn{i}-th row is the \eqn{g}-th column of \eqn{Z_i}.
#'
#' @examples
#' set.seed(1)
#' m <- 3; n <- 5; G <- 4
#' Z_list <- lapply(1:m, function(i) matrix(rnorm(n*G), n, G,
#'                    dimnames = list(NULL, paste0("feat", 1:G))))
#' Zi <- multi_matrix_transform(Z_list, g = 2)   # m x n
#' dim(Zi)
#'
#' @seealso \code{\link{multi_matrix_transform_inverse}}
#' @export

multi_matrix_transform <- function(Z_list, g){
  m <- length(Z_list); n <- nrow(Z_list[[1]])
  Zi_matrix <- matrix(NA, nrow = m, ncol = n)
  for (i in 1:m) {
    Zi_matrix[i, ] <- Z_list[[i]][, g]
  }
  return(Zi_matrix)
}

#' @rdname multi_matrix_transform
#'
#' @param Zi_matrix An \eqn{m \times n} numeric matrix produced by
#'   \code{multi_matrix_transform()}.
#' @param g_name Character scalar used as the single column name for each
#'   reconstructed \eqn{n \times 1} data frame.
#'
#' @return \code{multi_matrix_transform_inverse()} returns a list of length
#'   \eqn{m}; each element is an \eqn{n \times 1} data frame containing the
#'   corresponding row of \code{Zi_matrix}, with column name \code{g_name}.
#'
#' @examples
#' g_list <- multi_matrix_transform_inverse(Zi, g_name = "feat2")
#' length(g_list); nrow(g_list[[1]])
#'
#' @export
multi_matrix_transform_inverse <- function(Zi_matrix, g_name){
  m <- nrow(Zi_matrix); n <- ncol(Zi_matrix)
  Z_list <- vector("list", m)
  for (i in 1:m) {
    Z_list[[i]] <- data.frame(Zi_matrix[i, ])
    colnames(Z_list[[i]]) <- g_name
  }
  return(Z_list)
}


#' Numerically stabilized inverse square root of a symmetric matrix
#'
#' Computes a stabilized matrix inverse square root \eqn{\Sigma^{-1/2}} via
#' eigen-decomposition. For numerical stability, eigenvalues smaller than
#' \code{tol} are replaced by \code{tol} before inversion.
#'
#' @param sigma A symmetric \eqn{p \times p} matrix, ideally positive
#'   semi-definite.
#' @param tol Numeric tolerance. Eigenvalues smaller than \code{tol} are
#'   replaced by \code{tol}. Default is \code{1e-6}.
#'
#' @return A \eqn{p \times p} symmetric matrix giving the inverse square root
#'   based on the thresholded eigenvalues of \code{sigma}.
#'
#' @details
#' If any eigenvalues of \code{sigma} are smaller than \code{tol}, they are
#' replaced by \code{tol}, and a warning is issued.
#'
#' @examples
#' set.seed(1)
#' A <- crossprod(matrix(rnorm(25), 5, 5))  # SPD
#' S <- sigma_inv_sqrt(A)
#' round(max(abs(t(S) %*% A %*% S - diag(5))), 6)
#'
#' @seealso \code{\link{robust_matrix_solver}}
#' @export

sigma_inv_sqrt <- function(sigma, tol = 1e-6){
  eig <- eigen(sigma, symmetric = TRUE)
  orig_values <- eig$values
  eig_values <- orig_values

  eig_values[eig_values < tol] <- tol

  if (any(orig_values < tol)) {
    warning("Some eigenvalues were very small or non-positive and were thresholded.")
  }

  Q <- eig$vectors
  Lambda_inv_sqrt <- diag(1 / sqrt(eig_values), nrow = length(eig_values))
  inv_sqrt <- Q %*% Lambda_inv_sqrt %*% t(Q)
  return(inv_sqrt)
}


#' Numerically stabilized square root of a symmetric matrix
#'
#' Computes a stabilized matrix square root \eqn{\Sigma^{1/2}} via
#' eigen-decomposition. For numerical stability, eigenvalues smaller than
#' \code{tol} are replaced by \code{tol} before taking square roots.
#'
#' @param sigma A symmetric \eqn{p \times p} matrix, ideally positive
#'   semi-definite.
#' @param tol Numeric tolerance. Eigenvalues smaller than \code{tol} are
#'   replaced by \code{tol}. Default is \code{1e-6}.
#'
#' @return A \eqn{p \times p} symmetric matrix giving the square root
#'   based on the thresholded eigenvalues of \code{sigma}.
#'
#' @details
#' If any eigenvalues of \code{sigma} are smaller than \code{tol}, they are
#' replaced by \code{tol}, and a warning is issued. As a result, the returned
#' matrix is the square root of a regularized version of \code{sigma}, rather
#' than necessarily the exact square root of the original matrix.
#'
#' @examples
#' set.seed(1)
#' A <- crossprod(matrix(rnorm(25), 5, 5))  # SPD
#' S <- sigma_sqrt(A)
#' round(max(abs(S %*% S - A)), 6)
#'
#' @export

sigma_sqrt <- function(sigma, tol = 1e-6){
  eig <- eigen(sigma, symmetric = TRUE)
  orig_values <- eig$values
  eig_values <- orig_values

  eig_values[eig_values < tol] <- tol

  if (any(orig_values < tol)) {
    warning("Some eigenvalues were very small or non-positive and were thresholded.")
  }

  Q <- eig$vectors
  Lambda_sqrt <- diag(sqrt(eig_values), nrow = length(eig_values))
  sqrt_mat <- Q %*% Lambda_sqrt %*% t(Q)
  return(sqrt_mat)
}


#' Robust target covariance from batch-specific covariance estimates
#'
#' Computes a pooled target covariance matrix from a list of batch-specific
#' covariance estimates, with optional robust down-weighting of batches whose
#' covariance matrices are far from an initial weighted mean in Frobenius norm.
#'
#' @param Sigma_list A list of symmetric covariance matrices, typically one
#'   per batch.
#' @param weights Optional nonnegative numeric vector of batch weights of the
#'   same length as \code{Sigma_list}. If \code{NULL}, equal weights are used.
#' @param eps Small positive constant added for numerical stability when
#'   computing robust scale estimates and adaptive weights. Default is
#'   \code{1e-8}.
#' @param min_batches_robust Minimum number of batch levels required to apply
#'   robust reweighting. If the number of matrices in \code{Sigma_list} is less
#'   than this threshold, the function returns the default weighted mean without
#'   robust adjustment. Default is \code{3}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{sigma0}}{The final pooled target covariance matrix.}
#'   \item{\code{sigma0_init}}{The initial weighted mean covariance matrix
#'   before robust reweighting.}
#'   \item{\code{distances}}{Frobenius distances from each batch covariance
#'   matrix to \code{sigma0_init}.}
#'   \item{\code{adaptive_weights}}{Robust adaptive weights derived from the
#'   distances. These are all 1 when robust reweighting is not used.}
#'   \item{\code{final_weights}}{The final weights used to pool the covariance
#'   matrices, equal to the product of \code{weights} and
#'   \code{adaptive_weights}.}
#'   \item{\code{used_robust}}{Logical indicator of whether robust reweighting
#'   was applied.}
#' }
#'
#' @details
#' The function first computes an initial weighted mean covariance matrix.
#' When enough batch levels are available, it then measures the Frobenius
#' distance from each batch-specific covariance matrix to this initial center.
#' These distances are summarized using the median and median absolute
#' deviation (MAD), and a Huber-like down-weighting rule is used to reduce the
#' influence of outlying batch covariance matrices. The final target covariance
#' is the weighted average under these robustified weights.
#'
#' If fewer than \code{min_batches_robust} batch levels are provided, the
#' function falls back to the default weighted mean and returns
#' \code{used_robust = FALSE}.
#'
#' @examples
#' set.seed(1)
#' S1 <- crossprod(matrix(rnorm(9), 3, 3))
#' S2 <- crossprod(matrix(rnorm(9), 3, 3))
#' S3 <- crossprod(matrix(rnorm(9), 3, 3))
#'
#' out <- robust_sigma0_from_eb(list(S1, S2, S3))
#' out$sigma0
#' out$adaptive_weights
#'
#' @seealso \code{\link{frobenius_norm}}
#' @export
robust_sigma0_from_eb <- function(Sigma_list, weights = NULL, eps = 1e-8,
                                  min_batches_robust = 3) {
  B <- length(Sigma_list)
  if (is.null(weights)) weights <- rep(1, B)

  # initial weighted mean
  wsum <- sum(weights)
  Sigma_init <- Reduce(`+`, Map(function(S, w) w * S, Sigma_list, weights)) / wsum

  # fallback to default weighted mean when too few batch levels
  if (B < min_batches_robust) {
    return(list(
      sigma0 = Sigma_init,
      sigma0_init = Sigma_init,
      distances = rep(NA_real_, B),
      adaptive_weights = rep(1, B),
      final_weights = weights,
      used_robust = FALSE
    ))
  }

  # distances to initial center
  d <- sapply(Sigma_list, function(S) frobenius_norm(S, Sigma_init))

  # robust scale of distances
  med_d <- median(d, na.rm = TRUE)
  mad_d <- mad(d, center = med_d, constant = 1, na.rm = TRUE) + eps

  # Huber-like adaptive weights
  cutoff <- med_d + 1.5 * mad_d
  a <- pmin(1, cutoff / (d + eps))

  # combined weights
  w_rob <- weights * a
  w_rob_sum <- sum(w_rob)

  Sigma_rob <- Reduce(`+`, Map(function(S, w) w * S, Sigma_list, w_rob)) / w_rob_sum

  list(
    sigma0 = Sigma_rob,
    sigma0_init = Sigma_init,
    distances = d,
    adaptive_weights = a,
    final_weights = w_rob,
    used_robust = TRUE
  )
}

#' Robust matrix inversion/solve with fallbacks
#'
#' Attempts to invert a matrix (or solve \eqn{A X = I}) using several methods,
#' with an automatic choice based on conditioning and a safe pseudoinverse
#' fallback. Useful when covariance/precision matrices may be ill-conditioned.
#'
#' @param mat Numeric matrix (square for a true inverse).
#' @param method One of \code{"auto"}, \code{"solve"}, \code{"pseudoinverse"},
#'   \code{"qr"}, or \code{"chol"}. With \code{"auto"} (default), uses
#'   \code{MASS::ginv} when \code{kappa(mat)} is large; otherwise \code{solve()}.
#'
#' @details
#' If the chosen method errors, the function warns and returns
#' \code{MASS::ginv(mat + diag(1e-8, nrow(mat)))} as a regularized pseudoinverse.
#' The \code{"qr"} method computes an inverse by solving \eqn{A X = I} via QR;
#' the \code{"chol"} method uses \code{chol2inv(chol(A))} and requires SPD input.
#'
#' @return A matrix approximating \eqn{A^{-1}}. For rectangular \code{mat} and
#'   methods that allow it, this is a Moore–Penrose pseudoinverse.
#'
#' @examples
#' set.seed(1)
#' A <- crossprod(matrix(rnorm(16), 4, 4))
#' invA <- robust_matrix_solver(A, method = "auto")
#' max(abs(invA %*% A - diag(4)))
#'
#' # Ill-conditioned example falls back to pseudoinverse
#' B <- crossprod(cbind(1, 1 + 1e-10, rnorm(4)))
#' invB <- robust_matrix_solver(B)
#'
#' @export

robust_matrix_solver <- function(mat, method = c("auto", "solve", "pseudoinverse", "qr", "chol")) {
  method <- match.arg(method)

  if (method == "auto") {
    method <- if (kappa(mat) > 1/.Machine$double.eps) "pseudoinverse" else "solve"
  }

  tryCatch({
    switch(method,
           "solve" = solve(mat),
           "pseudoinverse" = MASS::ginv(mat),
           "qr" = qr.solve(mat),
           "chol" = chol2inv(chol(mat)))
  }, error = function(e) {
    warning("Matrix solve failed with ", method, " method: ", e$message)
    # Fallback to regularized pseudoinverse
    MASS::ginv(mat + diag(1e-8, nrow(mat)))
  })
}


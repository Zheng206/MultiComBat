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


#' Inverse square root of a (near) SPD matrix via eigen-decomposition
#'
#' Computes \eqn{\Sigma^{-1/2}} from an eigen factorization, thresholding small
#' or non-positive eigenvalues to \code{tol} for numerical stability.
#'
#' @param sigma Symmetric (ideally positive semi-definite) \eqn{p \times p} matrix.
#' @param tol Numeric tolerance; eigenvalues smaller than \code{tol} are replaced
#'   by \code{tol}. Default \code{1e-6}.
#'
#' @return A \eqn{p \times p} matrix \eqn{\Sigma^{-1/2}}.
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
  eig <- eigen(sigma)
  eig_values <- eig$values

  # Check for non-positive eigenvalues (close to zero) and threshold them
  eig_values[eig_values < tol] <- tol

  # If any eigenvalue is zero or negative, print a warning
  if (any(eig_values <= 0)) {
    warning("Some eigenvalues are zero or negative. Small eigenvalues have been thresholded.")
  }

  Q <- eig$vectors
  Lambda_inv_sqrt <- diag(1 / sqrt(eig_values))
  inv_sqrt <- Q %*% Lambda_inv_sqrt %*% t(Q)
  return(inv_sqrt)
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


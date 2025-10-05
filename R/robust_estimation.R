#' Tukey biweight mid-variance (univariate robust variance)
#'
#' Computes a robust variance estimate using Tukey's biweight with
#' median/MAD centering and scaling. Observations with \eqn{|u_i| \ge 1} receive
#' zero weight, where \eqn{u_i = (x_i - \mathrm{center})/(c \cdot \mathrm{MAD})}.
#'
#' @param data Numeric vector.
#' @param center Optional numeric center. Defaults to \code{median(data)}.
#' @param norm.unbiased Logical; if \code{TRUE} (default), uses
#'   \code{c = 9 / qnorm(0.75)} under
#'   normality. If \code{FALSE}, uses \code{c = 9}.
#'
#' @details
#' Let \eqn{d_i = x_i - \mathrm{center}}, \eqn{\mathrm{MAD} = \mathrm{median}(|d_i|)},
#' and \eqn{u_i = d_i/(c\,\mathrm{MAD})}. The estimator is
#' \deqn{\hat{\sigma}^2_{\mathrm{BI}} =
#' n \cdot \frac{\sum_i \mathbb{1}(|u_i|<1)\, d_i^2 (1-u_i^2)^4}
#' {\left\{\sum_i \mathbb{1}(|u_i|<1)\, (1-u_i^2)(1-5u_i^2)\right\}^2}.}
#'
#' @return A single numeric: robust variance estimate.
#'
#' @note
#' If \code{MAD == 0} or the denominator is zero (e.g., all points trimmed),
#' the result may be \code{NaN}/\code{Inf}. Consider falling back to
#' \code{var(data)} or adding a small jitter in such edge cases.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 0, 1), 10)  # one outlier
#' var(x)
#' biweight_midvar(x)           # robust to the outlier
#'
#' @seealso \code{\link[stats]{var}}
#' @export

biweight_midvar <- function(data, center=NULL, norm.unbiased = TRUE) {
  if (is.null(center)) {
    center <- median(data)
  }

  mad <- median(abs(data - center))
  d <- data - center
  c <- ifelse(norm.unbiased, 9/qnorm(0.75), 9)
  u <- d/(c*mad)

  n <- length(data)
  indic <- abs(u) < 1

  num <- sum(indic * d^2 * (1 - u^2)^4)
  dem <- sum(indic * (1 - u^2) * (1 - 5*u^2))^2

  n * num/dem
}

#' Robust covariance via Tukey biweight column aggregation
#'
#' Builds a robust covariance matrix from a (features \eqn{\times} samples)
#' matrix by downweighting entire \emph{columns} (samples) that are extreme on
#' many rows using Tukey's biweight. For each row, standardized residuals
#' \eqn{u_{ij}} are formed with a rowwise MAD; per-cell weights are
#' \eqn{w_{ij} = (1-u_{ij}^2)^2 \mathbb{1}(|u_{ij}|<1)}. Column weights
#' \eqn{w_j \propto \sum_i w_{ij}} are normalized to sum to 1, and the
#' covariance is \eqn{\sum_j w_j\, x_{\cdot j} x_{\cdot j}^\top}.
#'
#' @param centered_Zi Numeric matrix of size \eqn{m \times n} (rows = features,
#'   columns = samples). Ideally row-centered; the function still computes row
#'   medians and MADs internally.
#' @param norm.unbiased Logical; if \code{TRUE} (default), uses
#'   \code{c = 9 / qnorm(0.75)} in the standardization under normality. If \code{FALSE}, uses \code{c = 9}.
#'
#' @details
#' For each row \eqn{i}, compute \eqn{\mathrm{MAD}_i} from that row and
#' \eqn{u_{ij} = x_{ij}/(c\,\mathrm{MAD}_i)} (with the row median as center).
#' Aggregate per-cell weights across rows to obtain column weights \eqn{w_j},
#' normalize them to sum to 1, then form the weighted sum of outer products.
#' This targets column-wise outliers (samples) rather than row-wise ones.
#'
#' @return An \eqn{m \times m} numeric covariance matrix \code{Sigma_hat_ig}.
#'
#' @note
#' If any row has \code{MAD == 0} or if all \eqn{w_{ij}} are zero,
#' \code{sum(w_vec)} will be zero and the normalized column weights become
#' undefined. You may wish to guard against this and fall back to an unweighted
#' covariance or add a small epsilon to MADs.
#'
#' @examples
#' set.seed(1)
#' m <- 5; n <- 60
#' X <- matrix(rnorm(m * n), m, n)
#' X[, n] <- X[, n] + 8  # make one column a multivariate outlier
#'
#' S_classic <- tcrossprod(X) / n
#' S_robust  <- biweight_midvar_mul(X)
#' dim(S_robust); all.equal(dim(S_robust), c(m, m))
#'
#' @seealso \code{robustbase::covMcd} for a high-breakdown alternative.
#' @export

biweight_midvar_mul <- function(centered_Zi, norm.unbiased = TRUE) {

  mad_vec <- apply(centered_Zi, 1, function(x) median(abs(x - median(x))))
  c_val <- if (norm.unbiased) 9 / qnorm(0.75) else 9
  u_matrix <- sweep(centered_Zi, 1, c_val * mad_vec, "/")

  w_matrix <- (1 - u_matrix^2)^2 * (abs(u_matrix) < 1)
  w_vec <- colSums(w_matrix)
  w_vec <- w_vec / sum(w_vec)

  Sigma_hat_ig <- matrix(0, nrow = nrow(centered_Zi), ncol = nrow(centered_Zi))
  ni <- ncol(centered_Zi)
  for (j in 1:ni) {
    x <- centered_Zi[, j, drop = FALSE]
    Sigma_hat_ig <- Sigma_hat_ig + w_vec[j] * (x %*% t(x))
  }

  return(Sigma_hat_ig)
}

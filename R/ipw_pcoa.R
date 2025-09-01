#' IPW-adjusted PCoA (ordination with inverse-probability weights)
#'
#' Computes an inverse-probability-weighted principal coordinates analysis.
#' Starting from an (unweighted) Gower double-centered Gram matrix \eqn{G},
#' we apply *weighted centering* using propensity-based weights \eqn{w}:
#' \deqn{G_w = (I - \mathbf{v}\mathbf{1}^\top)\, G\, (I - \mathbf{1}\mathbf{v}^\top),
#' \quad \mathbf{v} = \frac{w}{\sum_i w_i}.}
#' An eigendecomposition of \eqn{G_w} yields IPW-adjusted principal coordinates.
#'
#' If \code{z = NULL} and \code{weights = NULL}, the method reduces to classical PCoA.
#' If \code{weights} are supplied, they are used directly (rescaled to sum to 1).
#'
#' @param D A square distance object (\code{dist}) or numeric matrix of size \eqn{n\times n}.
#' @param x Optional trait (length \eqn{n}) used only to *estimate* IPW from \code{z}
#'   when \code{weights} are not supplied. Can be numeric (continuous) or factor
#'   (binary or multiclass). If \code{x=NULL} and \code{weights=NULL}, all weights are 1.
#' @param z Optional confounder matrix or \code{data.frame} with \eqn{n} rows,
#'   used to fit the propensity model for \code{x} and construct IPW.
#' @param weights Optional vector of positive inverse-probability weights of length \eqn{n}.
#'   If supplied, \code{x}, \code{z} and \code{family} are ignored for weighting.
#'   Internally rescaled to sum 1 (centering operator depends on \eqn{\mathbf{v}} only).
#' @param family Character; one of \code{"auto"}, \code{"binomial"}, or \code{"gaussian"}.
#'   With \code{"auto"}, logistic is used for binary \code{x}, Gaussian working model for
#'   continuous, and multinomial (via \pkg{nnet}) for factors with \eqn{>2} levels.
#' @param trim Numeric in \([0,1)\); symmetric trim proportion for estimated weights.
#' @param stabilise Logical; return stabilised IPW (marginal over conditional).
#' @param k Optional integer; number of principal coordinates to return (default: all
#'   positive-eigenvalue axes).
#' @param weight_warn_cutoff Numeric; warn if \code{max(weights)} exceeds this threshold.
#'
#' @details
#' \strong{Why this works.} In IPW, the reweighted pseudo-population balances \code{x}
#' w.r.t. \code{z}. Applying weights in the centering operators (not in the distance itself)
#' produces an embedding whose group structure visualises the \emph{adjusted} separation
#' (i.e., what you would expect in the balanced pseudo-population). This is the ordination-level
#' companion to IPW-PERMANOVA.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{scores}: \eqn{n \times k} matrix of IPW principal coordinates.
#'   \item \code{eigenvalues}: vector of retained (nonnegative) eigenvalues.
#'   \item \code{prop_explained}: fraction of total positive eigenvalue sum explained by each axis.
#'   \item \code{cumulative}: cumulative explained proportions.
#'   \item \code{weights}: the final (rescaled-to-sum-1) weights used for centering.
#'   \item \code{gram}: the unweighted, Gower double-centered Gram matrix \eqn{G}.
#'   \item \code{gram_weighted}: the IPW-centered Gram matrix \eqn{G_w}.
#' }
#'
#' @examples
#' ## Confounded binary example: visual contrast (ordinary PCoA vs IPW-PCoA)
#' set.seed(1)
#' n <- 90
#' z <- data.frame(batch = sample(1:3, n, TRUE))
#' x <- rbinom(n, 1, c(0.2, 0.5, 0.8)[z$batch])  # confounded with batch
#' Y <- matrix(rnorm(n * 40, mean = 0.7 * x), n)
#' D <- dist(Y)
#' ipw_out <- ipw_pcoa(D, x, z)
#' head(ipw_out$scores[,1:2]); ipw_out$prop_explained[1:3]
#'
#' @importFrom stats dist glm predict lm fitted residuals dnorm var quantile model.matrix
#' @export
ipw_pcoa <- function(
    D, x = NULL, z = NULL, weights = NULL,
    family = c("auto","binomial","gaussian"),
    trim = 0.01, stabilise = TRUE, k = NULL,
    weight_warn_cutoff = 10
){
  ## ---- helpers --------------------------------------------------------------
  as_matrix <- function(D) {
    if (inherits(D, "dist")) return(as.matrix(D))
    D <- as.matrix(D)
    if (nrow(D) != ncol(D)) stop("Distance matrix must be square.")
    if (!isTRUE(all.equal(D, t(D)))) warning("D is not exactly symmetric; symmetrising.")
    (D + t(D)) / 2
  }
  est_weights <- function(x, z, fam, trim, stabilise) {
    z <- as.data.frame(z)
    if (identical(fam, "auto")) {
      if (is.factor(x) && nlevels(x) > 2) fam <- "multinomial" else
        if (is.factor(x) || all(x %in% c(0,1))) fam <- "binomial" else fam <- "gaussian"
    }
    if (fam == "binomial") {
      x_bin <- if (is.factor(x)) as.numeric(x == levels(x)[2]) else as.numeric(x)
      fit <- stats::glm(x_bin ~ ., data = data.frame(z), family = stats::binomial())
      pi <- stats::predict(fit, type = "response")
      p  <- mean(x_bin == 1)
      w  <- ifelse(x_bin == 1, p / pmax(pi, .Machine$double.eps),
                   (1 - p) / pmax(1 - pi, .Machine$double.eps))
      # if not stabilised, remove the marginal factor
      if (!stabilise) w <- w / c(p, 1 - p)[x_bin + 1L]
    } else if (fam == "multinomial") {
      if (!requireNamespace("nnet", quietly = TRUE))
        stop("family='multinomial' requires {nnet}. Install it or convert x to binary.")
      x_fac <- if (!is.factor(x)) factor(x) else x
      fit   <- nnet::multinom(x_fac ~ ., data = data.frame(z), trace = FALSE)
      P     <- stats::predict(fit, type = "probs")
      if (is.null(dim(P))) P <- cbind(P)
      levs  <- levels(x_fac)
      pi_i  <- P[cbind(seq_len(nrow(P)), match(x_fac, levs))]
      p_m   <- prop.table(table(x_fac))[as.character(x_fac)]
      w     <- as.numeric(p_m / pmax(pi_i, .Machine$double.eps))
    } else { # gaussian working model
      fit <- stats::lm(as.numeric(x) ~ ., data = data.frame(z))
      mu  <- stats::fitted(fit)
      sig <- sqrt(mean(stats::residuals(fit)^2))
      fyz <- stats::dnorm(as.numeric(x), mean = mu, sd = max(sig, .Machine$double.eps))
      fy  <- stats::dnorm(as.numeric(x),
                          mean = mean(as.numeric(x)),
                          sd   = stats::sd(as.numeric(x)) + .Machine$double.eps)
      w   <- if (stabilise) fy / pmax(fyz, .Machine$double.eps) else 1 / pmax(fyz, .Machine$double.eps)
    }
    if (is.numeric(trim) && trim > 0) {
      lo <- stats::quantile(w, trim / 2)
      hi <- stats::quantile(w, 1 - trim / 2)
      w  <- pmin(pmax(w, lo), hi)
    }
    w
  }
  ## ---- set up ---------------------------------------------------------------
  D  <- as_matrix(D)
  n  <- nrow(D)
  if (!is.null(x) && length(x) != n) stop("`x` must have length nrow(D).")
  if (!is.null(z) && nrow(as.data.frame(z)) != n) stop("`z` must have n rows.")
  # Gower double-centering (unweighted) to obtain Gram matrix G
  D2 <- D^2
  J  <- diag(n) - matrix(1/n, n, n)
  G  <- -0.5 * J %*% D2 %*% J
  # Build / get weights
  fam <- match.arg(family)
  if (is.null(weights)) {
    if (!is.null(x) && !is.null(z)) {
      w_raw <- est_weights(x, z, fam, trim, stabilise)
    } else {
      w_raw <- rep(1, n)
    }
  } else {
    if (length(weights) != n) stop("`weights` must have length n.")
    if (any(weights <= 0))   stop("All weights must be positive.")
    w_raw <- as.numeric(weights)
  }
  if (max(w_raw) > weight_warn_cutoff) {
    warning("Extreme weights detected (max = ", round(max(w_raw), 2),
            "). Check overlap / consider stronger trimming.")
  }
  v <- as.numeric(w_raw / sum(w_raw))  # normalized to sum 1 (needed by weighted centering)
  ## ---- weighted centering of Gram ------------------------------------------
  # G_w = (I - v 1^T) G (I - 1 v^T)
  One <- matrix(1, n, 1)
  H_L <- diag(n) - v %*% t(One)
  H_R <- diag(n) - One %*% t(v)
  Gw  <- H_L %*% G %*% H_R
  # Numerical symmetrization (very mild) before eigen
  Gw  <- (Gw + t(Gw)) / 2
  ## ---- eigendecomposition ---------------------------------------------------
  ei  <- eigen(Gw, symmetric = TRUE)
  vals <- ei$values
  vecs <- ei$vectors
  # keep nonnegative (tolerance) eigenvalues to avoid spurious negatives
  tol  <- max(1e-12, 1e-10 * max(abs(vals)))
  keep <- which(vals > tol)
  if (length(keep) == 0L) stop("No positive eigenvalues after IPW-centering; check distances/weights.")
  vals_pos <- vals[keep]
  vecs_pos <- vecs[, keep, drop = FALSE]
  # principal coordinates: U * sqrt(Lambda)
  scores <- vecs_pos %*% diag(sqrt(vals_pos), nrow = length(vals_pos))
  # explained variance proportions among positive eigenvalues
  prop  <- vals_pos / sum(vals_pos)
  cum   <- cumsum(prop)
  # optionally truncate to k axes
  if (!is.null(k)) {
    k <- min(k, ncol(scores))
    scores   <- scores[, seq_len(k), drop = FALSE]
    vals_pos <- vals_pos[seq_len(k)]
    prop     <- prop[seq_len(k)]
    cum      <- cum[seq_len(k)]
  }
  colnames(scores) <- paste0("PCo", seq_len(ncol(scores)))
  list(
    scores         = scores,
    eigenvalues    = vals_pos,
    prop_explained = prop,
    cumulative     = cum,
    weights        = as.numeric(v),     # the normalized weights used in centering
    gram           = G,
    gram_weighted  = Gw
  )
}

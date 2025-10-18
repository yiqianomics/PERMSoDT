#' IPW-adjusted & residualized PCoA (ordination with inverse-probability weights)
#'
#' Computes an inverse-probability-weighted principal coordinates analysis with
#' optional confounder residualization in the IPW-induced geometry.
#'
#' Starting from a squared distance matrix D^2, we apply *weighted* Gower
#' double-centering to obtain a Gram matrix:
#'   B_w = -1/2 * J_w D^2 J_w^T,  with  J_w = I - 1 v^T,  v = w / sum(w).
#' If confounders z are provided and residualize=TRUE, we remove the linear
#' component explained by z under the IPW geometry:
#'   B_res = W^{-1/2} M_Z~ ( W^{1/2} B_w W^{1/2} ) M_Z~ W^{-1/2},
#' where  Z~ = W^{1/2} Z and M_Z~ is the orthogonal projector onto the residual
#' space w.r.t. Z~. An eigendecomposition of B_res (or B_w if no residualization)
#' yields the principal coordinates.
#'
#' If z=NULL and weights=NULL, the method reduces to classical PCoA.
#' If weights are supplied, they are used directly (rescaled to sum to 1).
#'
#' @param D A square distance object (`dist`) or numeric matrix (n x n).
#' @param x Optional target/label (length n) used only to estimate IPW from `z`
#'   when `weights` are not supplied. Numeric or factor. If `x=NULL` and
#'   `weights=NULL`, all weights are 1.
#' @param z Optional confounder matrix or data.frame with n rows for residualization.
#'   If provided and `residualize=TRUE`, confounder effects are projected out in
#'   the IPW-induced geometry.
#' @param weights Optional positive inverse-probability weights (length n).
#'   If supplied, `x`, `z` and `family` are ignored for weighting. Internally
#'   rescaled to sum 1.
#' @param family Character; one of "auto","binomial","gaussian".
#' @param trim Numeric in [0,1); symmetric trim proportion for estimated weights.
#' @param stabilise Logical; return stabilised IPW (marginal over conditional).
#' @param residualize Logical; if TRUE and z is provided, perform weighted residualization.
#' @param add_intercept Logical; include intercept in Z when residualizing (default TRUE).
#' @param k Optional integer; number of axes to return (default: all positive-eig axes).
#' @param weight_warn_cutoff Numeric; warn if max(weights) exceeds this threshold.
#'
#' @return A list with components:
#'   - scores: n x k matrix of principal coordinates.
#'   - eigenvalues: retained nonnegative eigenvalues.
#'   - prop_explained: fraction of sum of positive eigenvalues per axis.
#'   - cumulative: cumulative explained fractions.
#'   - weights: normalized weights used for centering (sum to 1).
#'   - gram: weighted Gram before residualization (B_w).
#'   - gram_resid: residualized Gram (if residualize), else NULL.
#'
#' @importFrom stats dist glm predict lm fitted residuals dnorm var quantile model.matrix
#' @export
ipw_pcoa <- function(
    D, x = NULL, z = NULL, weights = NULL,
    family = c("multinomial","binomial","gaussian"),
    trim = 0.01, stabilise = TRUE,
    residualize = TRUE, add_intercept = TRUE,
    k = NULL, weight_warn_cutoff = 10
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
  mmatrix_Z <- function(z, add_intercept = TRUE) {
    z <- as.data.frame(z)
    if (add_intercept) {
      stats::model.matrix(~ ., data = z)  # includes intercept
    } else {
      M <- stats::model.matrix(~ . - 1, data = z)
      # if all numeric and no column names, ensure matrix
      as.matrix(M)
    }
  }
  proj_residualizer <- function(Ztil) {
    # returns MZtil = I - P_Ztil, with numerical stability via qr.solve
    if (is.null(Ztil) || ncol(Ztil) == 0) return(diag(nrow(Ztil)))
    R <- qr(t(Ztil) %*% Ztil)
    P <- Ztil %*% qr.solve(R, t(Ztil))  # (Z~ (Z~'Z~)^-1 Z~')
    diag(nrow(Ztil)) - P
  }
  
  ## ---- set up ---------------------------------------------------------------
  D  <- as_matrix(D)
  n  <- nrow(D)
  if (!is.null(x) && length(x) != n) stop("`x` must have length nrow(D).")
  if (!is.null(z) && nrow(as.data.frame(z)) != n) stop("`z` must have n rows.")
  
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
  v  <- as.numeric(w_raw / sum(w_raw))  # normalized to sum 1
  W  <- diag(as.numeric(w_raw), n, n)
  Wt <- diag(sqrt(as.numeric(w_raw)), n, n)
  
  ## ---- weighted double-centering (Gower) ------------------------------------
  D2 <- D^2
  Jw <- diag(n) - matrix(1, n, 1) %*% t(matrix(v, n, 1))  # I - 1 v^T
  Bw <- -0.5 * Jw %*% D2 %*% t(Jw)                        # weighted Gram
  
  ## ---- weighted residualization in IPW geometry -----------------------------
  Bres <- NULL
  if (!is.null(z) && isTRUE(residualize)) {
    Zmm  <- mmatrix_Z(z, add_intercept = add_intercept)
    Ztil <- Wt %*% Zmm
    Btil <- Wt %*% Bw %*% Wt
    MZt  <- proj_residualizer(Ztil)
    Bres <- solve(Wt) %*% MZt %*% Btil %*% MZt %*% solve(Wt)
    # mild symmetrization
    Bres <- (Bres + t(Bres)) / 2
    Guse <- Bres
  } else {
    Bw   <- (Bw + t(Bw)) / 2
    Guse <- Bw
  }
  
  ## ---- eigendecomposition ---------------------------------------------------
  ei   <- eigen(Guse, symmetric = TRUE)
  vals <- ei$values
  vecs <- ei$vectors
  tol  <- max(1e-12, 1e-10 * max(abs(vals)))
  keep <- which(vals > tol)
  if (length(keep) == 0L)
    stop("No positive eigenvalues after weighting/residualization; check distances/weights/z.")
  vals_pos <- vals[keep]
  vecs_pos <- vecs[, keep, drop = FALSE]
  scores   <- vecs_pos %*% diag(sqrt(vals_pos), nrow = length(vals_pos))
  prop     <- vals_pos / sum(vals_pos)
  cum      <- cumsum(prop)
  
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
    weights        = as.numeric(v),
    gram           = Bw,
    gram_resid     = Bres
  )
}

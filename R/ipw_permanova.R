#' IPW-adjusted PERMANOVA for a single trait (binary, continuous, or multiclass)
#'
#' Computes the inverse-probability-weighted PERMANOVA statistic
#' \eqn{F_{\mathrm{IPW}}=\dfrac{\mathrm{tr}(H_W G H_W)}
#' {\mathrm{tr}\{(I-H_W)\,G\,(I-H_W)\}}}, where
#' \eqn{H_W = Y (Y^\top W Y)^{-1} Y^\top W} is the weighted projector built from
#' inverse-probability weights that balance the trait \code{x} with respect to
#' confounders \code{z}. The distance Gram matrix \eqn{G} is computed by
#' Gower double-centering of the squared distances and is \emph{not} weighted.
#'
#' If \code{z = NULL}, the weights default to ones, reproducing classical one-way
#' PERMANOVA. A permutation p-value is available; at each shuffle of \code{x},
#' the propensity model is re-fit and weights are re-estimated.
#'
#' @param D A square distance object (\code{dist}) or numeric matrix of size \eqn{n\times n}.
#' @param x Trait of interest (length \eqn{n}): numeric (continuous) or factor
#'   (binary or multiclass). Factors are encoded via \code{model.matrix(~ x)}
#'   so the design \eqn{Y} includes an intercept.
#' @param z Optional confounder matrix or \code{data.frame} with \eqn{n} rows.
#'   If provided, inverse-probability weights are estimated from a propensity model.
#' @param weights Optional vector of positive inverse-probability weights of length \eqn{n}.
#'   If supplied, \code{z} and \code{family} are ignored for the observed statistic
#'   (but are still used during permutations if \code{permutation = TRUE} and \code{z} is not \code{NULL}).
#'   Internally, weights are rescaled to have mean 1 (scale cancels in \eqn{F_{\mathrm{IPW}}}).
#' @param family Character; one of \code{"auto"}, \code{"binomial"}, or \code{"gaussian"}.
#'   With \code{"auto"}, a logistic model is used when \code{x} is binary, a Gaussian
#'   working model when \code{x} is continuous, and (automatically) a multinomial
#'   model via \pkg{nnet} when \code{x} is a factor with \eqn{>2} levels.
#' @param permutation Logical; if \code{TRUE}, compute a permutation p-value with weight
#'   re-estimation at each shuffle of \code{x}.
#' @param B Integer; number of permutations (default \code{999}).
#' @param seed Integer seed for reproducibility.
#' @param trim Numeric in \([0, 1)\); symmetric trim proportion for weights
#'   (e.g., \code{0.01} trims the lower/upper 0.5\%).
#' @param stabilise Logical; if \code{TRUE}, return stabilised IPW
#'   (marginal \eqn{f_Y(y)} over conditional \eqn{f_{Y|Z}(y|z)}). See Details.
#' @param weight_warn_cutoff Numeric; if \code{max(weights)} exceeds this threshold,
#'   a warning is issued to flag poor overlap / extreme weights.
#'
#' @details
#' \strong{Propensity and weights.} For binary \code{x}, a logistic regression yields
#' \eqn{w_i \propto \Pr(Y=y_i)/\Pr(Y=y_i\mid Z=z_i)} with optional stabilisation;
#' for continuous \code{x}, a Gaussian working model provides
#' \eqn{w_i \propto f_Y(y_i)/f_{Y|Z}(y_i\mid z_i)}. When \code{x} is a factor with more
#' than two levels, a multinomial model (via \code{nnet::multinom}) is used under
#' \code{family = "auto"}. Weights are trimmed to improve numerical stability and
#' rescaled to mean 1.
#'
#' \strong{Design and centering.} The design matrix \eqn{Y} always includes an
#' intercept (or, equivalently, satisfies weighted centering), as required by the
#' large-sample arguments for IPW-PERMANOVA. The Gram matrix \eqn{G} is computed
#' from \code{D} by Gower double-centering and is not altered by the weights.
#'
#' \strong{Statistic.} The numerator and denominator are evaluated via trace identities:
#' \eqn{\mathrm{num}=\mathrm{tr}\{(Y^\top W Y)^{-1}(Y^\top W G Y)\}}, and
#' \eqn{\mathrm{den}=\mathrm{tr}(G)-\mathrm{num}}. This avoids forming \eqn{H_W}
#' explicitly and is efficient and numerically stable.
#'
#' \strong{Permutation.} Under \code{permutation = TRUE}, the trait \code{x} is
#' permuted without replacement; for each permutation, a new \eqn{Y} and a new set
#' of weights are fit from the shuffled data (if \code{z} is provided), and the
#' statistic is recomputed. The reported p-value uses the standard \code{(1 + #\{F_b \ge F_{\mathrm{obs}}\})/(B+1)}
#' correction.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{F_statistic}: observed \eqn{F_{\mathrm{IPW}}}.
#'   \item \code{p_value}: permutation p-value (or \code{NA} if \code{permutation = FALSE}).
#'   \item \code{R2}: \eqn{F_{\mathrm{IPW}}/(1+F_{\mathrm{IPW}})}.
#'   \item \code{weights}: the (rescaled) weights used for the observed statistic.
#'   \item \code{perm_stats}: vector of permutation statistics (if computed).
#' }
#'
#' @examples
#' ## Binary trait with confounding
#' set.seed(1)
#' n <- 60
#' z <- data.frame(batch = sample(1:3, n, TRUE))
#' x <- rbinom(n, 1, c(0.2, 0.5, 0.8)[z$batch])  # confounded with batch
#' Y <- matrix(rnorm(n * 30, mean = 0.6 * x), n)
#' fit <- ipw_permanova(dist(Y), x, z, B = 199, seed = 2025)
#' fit$F_statistic; fit$p_value; fit$R2
#'
#' ## Continuous trait
#' set.seed(2)
#' x2  <- scale(0.3 * z$batch + rnorm(n))
#' Y2  <- matrix(rnorm(n * 25, mean = 0.5 * x2), n)
#' fit2 <- ipw_permanova(dist(Y2), as.numeric(x2), z, B = 199, seed = 2025)
#' fit2$F_statistic; fit2$p_value
#'
#' ## Multiclass trait (auto-selects multinomial when nnet is available)
#' if (requireNamespace("nnet", quietly = TRUE)) {
#'   set.seed(3)
#'   gx <- factor(sample(letters[1:3], n, TRUE, prob = c(0.2,0.5,0.3)))
#'   Y3 <- matrix(rnorm(n * 20, mean = c(0.8, 0.4, 0)[gx]), n)
#'   fit3 <- ipw_permanova(dist(Y3), gx, z, B = 99, seed = 2025)
#'   fit3$R2
#' }
#'
#' @importFrom stats dist glm predict lm fitted residuals dnorm var quantile model.matrix
#' @importFrom utils head
#' @seealso \code{\link[vegan]{adonis2}} for classical PERMANOVA (unweighted).
#' @export
ipw_permanova <- function(
    D, x, z = NULL, weights = NULL,
    family = c("auto", "binomial", "gaussian"),
    permutation = TRUE, B = 999, seed = 2025,
    trim = 0.01, stabilise = TRUE,
    weight_warn_cutoff = 10
) {
  ## ---- helpers --------------------------------------------------------------
  as_matrix <- function(D) {
    if (inherits(D, "dist")) return(as.matrix(D))
    D <- as.matrix(D)
    if (nrow(D) != ncol(D)) stop("Distance matrix must be square.")
    if (!isTRUE(all.equal(D, t(D)))) warning("D is not exactly symmetric; symmetrising.")
    (D + t(D)) / 2
  }
  
  mm_Y <- function(x) {
    if (is.factor(x)) {
      # Intercept + (K-1) dummies
      model.matrix(~ x)
    } else {
      # Intercept + numeric covariate
      cbind(Intercept = 1, x = as.numeric(x))
    }
  }
  
  est_weights <- function(x, z, fam, trim, stabilise) {
    z <- as.data.frame(z)
    if (identical(fam, "auto")) {
      if (is.factor(x) && nlevels(x) > 2) fam <- "multinomial"
      else if (is.factor(x) || all(x %in% c(0, 1))) fam <- "binomial"
      else fam <- "gaussian"
    }
    
    if (fam == "binomial") {
      # Logistic PS; stabilized by marginal Pr(Y=1)
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
      if (is.null(dim(P))) P <- cbind(P)  # 2-level edge case
      levs  <- levels(x_fac)
      pi_i  <- P[cbind(seq_len(nrow(P)), match(x_fac, levs))]
      p_m   <- prop.table(table(x_fac))[as.character(x_fac)]
      w     <- as.numeric(p_m / pmax(pi_i, .Machine$double.eps))
      # (Always stabilized by marginal class probabilities)
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
    
    # Trim extremes and rescale to mean 1
    if (is.numeric(trim) && trim > 0) {
      lo <- stats::quantile(w, trim / 2)
      hi <- stats::quantile(w, 1 - trim / 2)
      w  <- pmin(pmax(w, lo), hi)
    }
    w / mean(w)
  }
  
  # Efficient traces:
  # num = tr( (Y'WY)^(-1) (Y'WG Y) ), den = tr(G) - num
  tr_num <- function(G, Y, w) {
    WY  <- Y * w                       # each column scaled by w
    B   <- crossprod(Y, WY)            # Y' W Y
    A   <- crossprod(WY, G %*% Y)      # Y' W G Y
    sum(diag(solve(B, A)))             # tr(B^{-1} A)
  }
  
  ## ---- checks & setup -------------------------------------------------------
  set.seed(seed)
  
  D <- as_matrix(D)
  n <- nrow(D)
  
  if (length(x) != n) stop("`x` must have length nrow(D).")
  if (!is.null(z) && nrow(as.data.frame(z)) != n) stop("`z` must have n rows.")
  if (is.numeric(x) && var(x) == 0) stop("Trait `x` is constant; F statistic undefined.")
  
  # Gower double-centering (unweighted, as in manuscript)
  D2 <- D^2
  J  <- diag(n) - matrix(1 / n, n, n)
  G  <- -0.5 * J %*% D2 %*% J
  
  # Design matrix Y (include intercept per Assumption 5)
  Y <- mm_Y(x)
  
  # Weights
  fam <- match.arg(family)
  if (is.null(weights)) {
    weights <- if (is.null(z)) rep(1, n) else est_weights(x, z, fam, trim, stabilise)
  } else {
    if (length(weights) != n) stop("`weights` must have length n.")
    if (any(weights <= 0))   stop("All weights must be positive.")
    weights <- weights / mean(weights)
  }
  if (max(weights) > weight_warn_cutoff) {
    warning("Extreme weights detected (max = ", round(max(weights), 2),
            "). Check overlap / consider stronger trimming.")
  }
  
  # Statistic
  num   <- tr_num(G, Y, weights)
  den   <- sum(diag(G)) - num
  F_obs <- num / den
  R2    <- F_obs / (1 + F_obs)
  
  if (!permutation) {
    return(list(F_statistic = F_obs, p_value = NA_real_, R2 = R2,
                weights = as.numeric(weights), perm_stats = NULL))
  }
  
  # Permutations: shuffle x, recompute weights & Y each time
  perm_F <- numeric(B)
  for (b in seq_len(B)) {
    xb <- sample(x, replace = FALSE)
    Yb <- mm_Y(xb)
    wb <- if (is.null(z)) rep(1, n) else est_weights(xb, z, fam, trim, stabilise)
    perm_F[b] <- tr_num(G, Yb, wb) / (sum(diag(G)) - tr_num(G, Yb, wb))
  }
  
  list(
    F_statistic = F_obs,
    p_value     = (1 + sum(perm_F >= F_obs)) / (B + 1),
    R2          = R2,
    weights     = as.numeric(weights),
    perm_stats  = perm_F
  )
}

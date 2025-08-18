#' Type III Sum of Dissimilarity Table (SoDT) for PERMANOVA (multi-term)
#'
#' Computes the \strong{Type III (marginal) contribution} of a target term to the
#' PERMANOVA sum of dissimilarities, adjusting for all other terms in a multi-term
#' design. Let \eqn{G} be the Gower double-centered Gram matrix from the distance
#' matrix \eqn{D}. For a full design \eqn{F=[F_1,\ldots,F_K]}, the full-model
#' projector is \eqn{H = F(F^\top F)^{-1}F^\top}, the residual projector is
#' \eqn{R = I - H}, and the within (residual) sum of dissimilarities is
#' \eqn{\mathrm{SSW} = \mathrm{tr}(R G R)}.
#'
#' For a specific term \eqn{k}, the Type III projector given all other terms is
#' \deqn{
#'   H_{k\mid -k} \;=\; R_{-k}\,F_k\,\big(F_k^\top R_{-k} F_k\big)^{-1} F_k^\top R_{-k},
#'   \quad R_{-k} = I - H_{-k},\; H_{-k} = F_{-k}(F_{-k}^\top F_{-k})^{-1}F_{-k}^\top,
#' }
#' and the term’s Type III sum of dissimilarities is
#' \eqn{\mathrm{SS}_k = \mathrm{tr}(H_{k\mid -k}\, G)}, with degrees of freedom
#' \eqn{df_k = \mathrm{rank}(R_{-k} F_k)}. The identity
#' \eqn{\mathrm{SST}=\mathrm{tr}(G)=\sum_k \mathrm{SS}_k + \mathrm{SSW}} holds.
#'
#' Optionally, group-wise residual shares \eqn{\mathrm{SSW}_g} are returned via
#' \eqn{A = R G R} and \eqn{\mathrm{SSW}_g = \mathrm{tr}(P_g A P_g)
#'  = \sum_{i \in g} A_{ii}}, where \eqn{P_g} selects samples in group \eqn{g}.
#'
#' A permutation \eqn{F}-test for the chosen term is provided:
#' \eqn{F_k = (\mathrm{SS}_k/df_k)\,/\,(\mathrm{SSW}/df_W)}, with
#' \eqn{df_W = n - \mathrm{rank}(F)}. Permutations shuffle sample rows jointly across
#' all term matrices in \code{terms}.
#'
#' @param D A square distance object (\code{dist}) or numeric \eqn{n \times n} matrix
#'   with zero diagonal and symmetry.
#' @param terms A \emph{named} list of model matrices \eqn{F_k} (each \eqn{n \times p_k}),
#'   one per term in the design. The column encoding (e.g., treatment or sum contrasts)
#'   is up to the caller.
#' @param term Character scalar giving the name of the term in \code{terms} to test
#'   (the “main variable”). All other entries in \code{terms} are treated as covariates
#'   to adjust for (confounders).
#' @param group Optional vector of group labels (length \eqn{n}) to report group-wise
#'   residual shares \eqn{\mathrm{SSW}_g}.
#' @param nperm Integer; number of permutations for the \eqn{F}-test (default \code{999}).
#' @param seed Integer random seed for reproducibility.
#'
#' @details
#' \strong{Design and ranks.} The function is agnostic to how each \eqn{F_k} is built.
#' The degrees of freedom are computed as \eqn{df_k=\mathrm{rank}(R_{-k}F_k)} and
#' \eqn{df_W = n - \mathrm{rank}(F)} via QR ranks. If \eqn{df_k=0} or \eqn{df_W=0},
#' the \eqn{F}-statistic is reported as \code{NA}.
#'
#' \strong{Permutation scheme.} Each permutation shuffles sample indices once and applies
#' the same permutation to the rows of every term matrix in \code{terms}, then recomputes
#' \eqn{H}, \eqn{R}, \eqn{H_{k\mid -k}}, \eqn{\mathrm{SSW}}, and \eqn{\mathrm{SS}_k}.
#' The p-value uses the standard \code{(1 + #\{F_b \ge F_{\mathrm{obs}}\})/(1 + B)} correction.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{table}: data frame with \code{Term}, \code{SumOfDissimilarity}
#'         (each \eqn{\mathrm{SS}_k}, plus \code{Residual} = \eqn{\mathrm{SSW}} and
#'         \code{Total} = \eqn{\mathrm{SST}}), and \code{df} (each \eqn{df_k}, plus
#'         \eqn{df_W} and \eqn{n-1}).
#'   \item \code{F_stat}: observed \eqn{F_k} for \code{term}.
#'   \item \code{p_value}: permutation p-value for \code{term}.
#'   \item \code{df_term}: \eqn{df_k} for \code{term}.
#'   \item \code{df_within}: \eqn{df_W}.
#'   \item \code{group_within}: optional named vector of \eqn{\mathrm{SSW}_g} (if \code{group} supplied).
#'   \item \code{perm_stats}: vector of permutation \eqn{F_k} values.
#' }
#'
#' @examples
#' ## Synthetic example with a grouping factor (to test), plus batch + age as covariates
#' set.seed(42)
#' n <- 60
#' batch <- factor(sample(1:3, n, TRUE))
#' age   <- scale(rnorm(n))
#' group <- factor(sample(letters[1:2], n, TRUE))
#'
#' ## Build term matrices (user controls contrasts and columns)
#' F_group <- model.matrix(~ 0 + group)     # main effect to test
#' F_batch <- model.matrix(~ 0 + batch)     # confounder
#' F_age   <- model.matrix(~ age)           # confounder
#'
#' ## Simulate an embedding and distances with a group effect
#' X  <- matrix(rnorm(n * 20), n, 20) + 0.6 * F_group %*% c(1, -1)
#' D  <- dist(X)
#'
#' terms <- list(group = F_group, batch = F_batch, age = F_age)
#' out   <- compute_sodt_type3(D, terms, term = "group", group = group, nperm = 199)
#' out$table
#' out$F_stat; out$p_value
#'
#' @importFrom stats qr
#' @export
compute_sodt_type3 <- function(
    D,
    terms,                 # named list of n x p_k model matrices for each term
    term,                  # name of the term to test (must be in names(terms))
    group = NULL,          # optional vector of group labels for SSW_g
    nperm = 999,
    seed = 2025
) {
  set.seed(seed)
  
  # --- Input checks ----------------------------------------------------------
  if (inherits(D, "dist")) D <- as.matrix(D) else D <- as.matrix(D)
  if (!all(D == t(D))) stop("D must be symmetric.")
  if (any(abs(diag(D)) > .Machine$double.eps)) stop("D must have 0 diagonal.")
  if (!is.list(terms) || is.null(names(terms)))
    stop("`terms` must be a named list of model matrices.")
  n <- nrow(D)
  if (any(vapply(terms, function(M) nrow(M) == n, logical(1)) == FALSE))
    stop("Each term matrix in `terms` must have n rows.")
  if (!(term %in% names(terms)))
    stop("`term` must be one of names(terms).")
  
  # --- Gower double-centering: G = -1/2 J D^2 J -----------------------------
  D2 <- D^2
  J  <- diag(n) - matrix(1 / n, n, n)
  G  <- -0.5 * J %*% D2 %*% J
  trG <- sum(diag(G))
  
  # --- Build full design F and projectors -----------------------------------
  F_list <- terms
  F_full <- do.call(cbind, F_list)
  
  # rank-safe projector H = F (F'F)^-1 F'
  XTX   <- crossprod(F_full)
  H     <- F_full %*% solve(XTX, t(F_full))
  R     <- diag(n) - H
  
  # SSW and dfW
  SST   <- trG
  SSW   <- sum(diag(R %*% G %*% R))     # equivalently tr(G) - tr(H G)
  dfW   <- n - qr(F_full)$rank
  
  # --- Type III for each term: SS_k, df_k -----------------------------------
  term_names <- names(F_list)
  SS_k  <- setNames(numeric(length(F_list)), term_names)
  df_k  <- setNames(integer(length(F_list)), term_names)
  
  for (nm in term_names) {
    # F_-k and its projector
    F_minus <- do.call(cbind, F_list[setdiff(term_names, nm)])
    H_minus <- if (is.null(F_minus) || length(F_minus) == 0) {
      matrix(0, n, n)
    } else {
      F_minus %*% solve(crossprod(F_minus), t(F_minus))
    }
    R_minus <- diag(n) - H_minus
    
    # Type III projector for term k given others:
    Fk      <- F_list[[nm]]
    Zk      <- R_minus %*% Fk                    # residualized columns
    B       <- crossprod(Zk)                     # Fk' R_-k Fk
    # rank of H_{k|-k} equals rank(Zk)
    rk      <- qr(Zk)$rank
    if (rk == 0) {
      SS_k[nm] <- 0
      df_k[nm] <- 0
    } else {
      Hk <- if (ncol(Fk) == 0) matrix(0, n, n) else Zk %*% solve(B, t(Zk))
      # SS_k = tr(H_{k|-k} G) (trace identity for symmetric idempotent Hk)
      SS_k[nm] <- sum(diag(Hk %*% G %*% Hk))
      df_k[nm] <- rk
    }
  }
  
  # --- Check identity (optional small tolerance) -----------------------------
  # identity: SST == sum_k SS_k + SSW
  if (abs(SST - (sum(SS_k) + SSW)) > 1e-6 * max(1, abs(SST))) warning("SoDT identity off beyond tolerance.")
  
  # --- F-stat for requested term --------------------------------------------
  if (df_k[term] == 0 || dfW == 0) {
    F_obs <- NA_real_
  } else {
    F_obs <- (SS_k[term] / df_k[term]) / (SSW / dfW)
  }
  
  # --- Group-wise residual shares (optional) ---------------------------------
  SSW_g <- NULL
  if (!is.null(group)) {
    fgrp <- factor(group)
    A    <- R %*% G %*% R
    SSW_g <- setNames(vapply(levels(fgrp), function(g) {
      idx <- which(fgrp == g)
      sum(diag(A[idx, idx, drop = FALSE]))
    }, numeric(1)), levels(fgrp))
  }
  
  # --- Permutations: permute rows of all term matrices together --------------
  F_perm <- rep(NA_real_, nperm)
  if (!is.na(F_obs)) {
    for (b in seq_len(nperm)) {
      perm <- sample.int(n)
      # permute rows of each term matrix
      F_list_b <- lapply(F_list, function(M) M[perm, , drop = FALSE])
      F_full_b <- do.call(cbind, F_list_b)
      XTX_b    <- crossprod(F_full_b)
      H_b      <- F_full_b %*% solve(XTX_b, t(F_full_b))
      R_b      <- diag(n) - H_b
      
      # SSW_b
      SSW_b    <- sum(diag(R_b %*% G %*% R_b))
      dfW_b    <- n - qr(F_full_b)$rank
      
      # Type III for the target term only (faster)
      F_minus_b <- do.call(cbind, F_list_b[setdiff(term_names, term)])
      H_minus_b <- if (is.null(F_minus_b) || length(F_minus_b) == 0) {
        matrix(0, n, n)
      } else {
        F_minus_b %*% solve(crossprod(F_minus_b), t(F_minus_b))
      }
      R_minus_b <- diag(n) - H_minus_b
      
      Fk_b   <- F_list_b[[term]]
      Zk_b   <- R_minus_b %*% Fk_b
      rk_b   <- qr(Zk_b)$rank
      if (rk_b == 0 || dfW_b == 0) {
        F_perm[b] <- NA_real_
      } else {
        Hk_b   <- Zk_b %*% solve(crossprod(Zk_b), t(Zk_b))
        SS_k_b <- sum(Hk_b * G)
        F_perm[b] <- (SS_k_b / rk_b) / (SSW_b / dfW_b)
      }
    }
  }
  
  # p-value with +1 correction, dropping NAs safely
  valid <- is.finite(F_perm)
  p_val <- if (is.na(F_obs) || !any(valid)) NA_real_ else {
    (1 + sum(F_perm[valid] >= F_obs)) / (1 + sum(valid))
  }
  
  # --- Assemble output tables ------------------------------------------------
  term_table <- data.frame(
    Term = c(term_names, "Residual", "Total"),
    SumOfDissimilarity = c(unname(SS_k), SSW, SST),
    df = c(unname(df_k), dfW, n - 1L),
    stringsAsFactors = FALSE
  )
  
  list(
    table = term_table,
    F_stat = F_obs,
    p_value = p_val,
    df_term = df_k[term],
    df_within = dfW,
    group_within = SSW_g,
    perm_stats = F_perm
  )
}

################################ Replicate #########################################
args    <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.numeric(args[1]); if (is.na(rep_id)) rep_id <- 1L
setwd("/home/yxz3116/IPWPERMANOVA/Ordinal")
################################ Simulation ########################################
user_lib <- "/home/yxz3116/R_Package"

library(phyloseq, lib.loc = user_lib)
library(vegan, lib.loc = user_lib)
# install.packages("hexbin", lib = user_lib)
library(hexbin, lib.loc = user_lib)
library(GUniFrac, lib.loc = user_lib)
# install.packages("RSpectra", lib = user_lib)
library(RSpectra, lib.loc = user_lib)
library(Matrix, lib.loc = user_lib)
library(ape, lib.loc = user_lib)
# install.packages("phangorn", lib = user_lib)
library(phangorn, lib.loc = user_lib)
library(nnet, lib.loc = user_lib)
# Tree only if build_tree=TRUE:
#   install.packages(c("ape","phangorn"))
#   library(ape); library(phangorn)

## -- Build a global hex grid once from a bounding box (no per-sample binning) --
# We generate a dense uniform cloud inside the bbox only to let hexbin produce a grid,
# then take hex cell centers as fixed OTUs.
.hex_grid_centers <- function(xmin, xmax, ymin, ymax, bin_size = 0.02){
  dx <- bin_size; dy <- sqrt(3)/2 * dx
  xs <- seq(xmin, xmax, by = dx); ys <- seq(ymin, ymax, by = dy)
  centers <- vector("list", length(ys))
  for (j in seq_along(ys)){
    yj <- ys[j]; xoff <- if ((j %% 2) == 0) dx/2 else 0
    xrow <- xs + xoff; keep <- xrow >= xmin & xrow <= xmax
    centers[[j]] <- cbind(cx = xrow[keep], cy = rep(yj, sum(keep)))
  }
  centers <- do.call(rbind, centers)
  keep <- centers[,1] >= xmin & centers[,1] <= xmax & centers[,2] >= ymin & centers[,2] <= ymax
  centers[keep, , drop=FALSE]
}

.solve_bin_size_for_K <- function(xmin, xmax, ymin, ymax, target_K){
  area_bbox <- (xmax-xmin)*(ymax-ymin)
  sqrt((2*area_bbox)/(sqrt(3)*max(target_K,1)))
}

.center_probs <- function(centers, o, r, alpha, k_abund=1, k_mod=1, k_rare=1, eps=1e-15){
  dx <- centers[,1]-o[1]; dy <- centers[,2]-o[2]
  d  <- sqrt(dx*dx + dy*dy)
  in_ball <- d <= r
  if(!any(in_ball)){ j <- which.min(d); p <- rep(0, nrow(centers)); p[j] <- 1; return(p) }
  x <- pmin(1, pmax(0, d / r))
  base <- x^alpha
  band <- rep(1, length(d))
  band[x <= 0.20]             <- k_abund
  band[x >= 0.40 & x <= 0.80] <- k_mod
  band[x >= 0.80 & x <= 1.00] <- k_rare
  prob <- base * band
  prob[!in_ball] <- 0
  prob <- prob + eps
  prob / sum(prob)
}

.build_midpoint_tree <- function(OTU){
  otu_rel <- sweep(OTU, 1, pmax(rowSums(OTU),1), "/")
  D_otu   <- dist(t(otu_rel))
  tr      <- ape::nj(D_otu)
  tr$tip.label <- colnames(otu_rel)
  phangorn::midpoint(tr)
}

## ---------- MAIN: 3-level ordered factor X with ordinal GPS ----------
ConMiseq_fixedGrid_Ordinal <- function(
    n = 60, N_pts = 300,
    # Confounder
    Z_sd = 1.0,                     # Z ~ N(0, Z_sd^2)
    # Ordinal GPS model: proportional-odds  (three levels: L<M<H)
    # P(X <= k | Z) = plogis(theta_k - gamma * Z), k=1 (L), 2 (M). H is remainder.
    theta = c(L=-0.5, M=0.5),       # cutpoints (ordered: theta_L < theta_M)
    gamma = 1.5,                    # slope (larger => stronger X~Z association)
    # Marginal stabilization: use empirical pi_k or pass pi_target; default = empirical
    pi_target = NULL,               # optional c(pi_L, pi_M, pi_H); else estimated from sample
    # Geometry channels driven by ordered level score s in {-1,0,1}
    level_scores = c(L=-1, M=0, H=1), # numeric scores for monotone mapping
    r0 = 0.15,   r_slope = 0.00,    # r = r0 + r_slope * s
    a0 = 0.5,    a_slope = 0.00,    # a = a0 + a_slope * s
    k0 = 0.0,    k_slope = 0.00,    # log kR = k0 + k_slope * s
    # Optional centroid separation by levels (set 0 for none)
    Delta_levels = 0.00,            # baseline separation across levels along x (monotone in s)
    bo_z = 0.00,                    # extra centroid shift by Z along x (membership confounding)
    # Grid controls
    centers = NULL, bin_size = 0.02, bbox_margin = 0.10, target_K = NULL,
    min_prev = 0.05, min_var = 1e-8,
    # IPTW stabilization & trimming
    w_trim_q = 0.99,
    # Tree
    build_tree = TRUE,
    # Reproducibility
    seed = 1,
    return_centers = FALSE
){
  set.seed(seed)
  
  ## 1) Confounder
  Z <- rnorm(n, sd = Z_sd)
  
  ## 2) Ordinal exposure via proportional-odds model
  # cumulative probs
  eta_L <- theta["L"] - gamma * Z
  eta_M <- theta["M"] - gamma * Z
  cdf_L <- plogis(eta_L)
  cdf_M <- plogis(eta_M)
  # ensure monotone ordering
  cdf_L <- pmin(pmax(cdf_L, 1e-12), 1-1e-12)
  cdf_M <- pmin(pmax(pmax(cdf_M, cdf_L + 1e-9), 1e-12), 1-1e-12)
  # category probabilities
  pL <- cdf_L
  pM <- pmax(0, cdf_M - cdf_L)
  pH <- pmax(0, 1 - cdf_M)
  P  <- cbind(L=pL, M=pM, H=pH)
  
  # sample X
  u <- runif(n)
  X_idx <- ifelse(u < pL, 1L, ifelse(u < pL+pM, 2L, 3L))
  X_ord <- factor(c("L","M","H")[X_idx], levels = c("L","M","H"), ordered = TRUE)
  
  # GPS for realized category; stabilized IPTW
  gps <- P[cbind(seq_len(n), X_idx)]
  if (is.null(pi_target)){
    pi_hat <- colMeans(P)  # marginal via model (or use table(X_ord)/n)
  } else {
    stopifnot(length(pi_target) == 3)
    pi_hat <- pi_target / sum(pi_target)
  }
  w_stab <- pi_hat[X_idx] / pmax(1e-12, gps)
  w_trim <- pmin(w_stab, quantile(w_stab, w_trim_q, names = FALSE))
  
  ## 3) Map ordered level to numeric score s and build geometry
  s <- as.numeric(level_scores[X_ord])  # -1,0,1 by default
  # Centroids: place levels along x if Delta_levels>0 (monotone in s), plus Z-shift
  o_list <- lapply(seq_len(n), function(i){
    c( (Delta_levels/2) * s[i], 0 ) + c(bo_z * Z[i], 0)
  })
  # Radius / evenness / rare band
  r  <- pmax(0.05, r0 + r_slope * s)
  a  <- a0 + a_slope * s
  kR <- exp(k0 + k_slope * s)
  kM <- 1.0; kA <- 1.0
  
  ## 4) Build/use grid
  xs_min <- min(vapply(seq_len(n), function(i) o_list[[i]][1] - r[i], 0.0)) - bbox_margin
  xs_max <- max(vapply(seq_len(n), function(i) o_list[[i]][1] + r[i], 0.0)) + bbox_margin
  ys_min <- min(vapply(seq_len(n), function(i) o_list[[i]][2] - r[i], 0.0)) - bbox_margin
  ys_max <- max(vapply(seq_len(n), function(i) o_list[[i]][2] + r[i], 0.0)) + bbox_margin
  if (is.null(centers)){
    if(!is.null(target_K)){
      bin_size <- .solve_bin_size_for_K(xs_min, xs_max, ys_min, ys_max, target_K)
    }
    centers <- .hex_grid_centers(xs_min, xs_max, ys_min, ys_max, bin_size)
  } else {
    centers <- as.matrix(centers); stopifnot(ncol(centers)==2); colnames(centers) <- c("cx","cy")
  }
  
  ## 5) Counts via multinomial per sample
  K <- nrow(centers)
  OTU <- matrix(0L, nrow = n, ncol = K)
  for(i in seq_len(n)){
    p <- .center_probs(centers, o=o_list[[i]], r=r[i], alpha=a[i],
                       k_abund=kA, k_mod=kM, k_rare=kR[i])
    OTU[i, ] <- as.vector(rmultinom(1, size = N_pts, prob = p))
  }
  colnames(OTU) <- paste0("OTU_", seq_len(K))
  rownames(OTU) <- paste0("S", seq_len(n))
  
  ## 6) Filter cols to avoid zeros-only taxa
  keep <- colSums(OTU) > 0
  OTU  <- OTU[, keep, drop=FALSE]
  centers <- centers[keep, , drop=FALSE]
  if(!is.null(min_prev)){
    thr <- if(min_prev < 1) ceiling(min_prev * n) else as.integer(min_prev); thr <- max(thr, 2L)
    keep2 <- colSums(OTU > 0) >= thr
    OTU <- OTU[, keep2, drop=FALSE]; centers <- centers[keep2,,drop=FALSE]
  }
  if(!is.null(min_var) && min_var > 0){
    keep3 <- apply(OTU, 2, function(z) var(as.numeric(z)) > min_var)
    OTU <- OTU[, keep3, drop=FALSE]; centers <- centers[keep3,,drop=FALSE]
  }
  
  ## 7) Build phyloseq (optionally with tree)
  ps_otu <- otu_table(OTU, taxa_are_rows = FALSE)
  
  # add polynomial design for ordered factor
  Xmm <- model.matrix(~ X_ord)  # (Intercept), X_ord.L, X_ord.Q
  colnames(Xmm) <- make.names(colnames(Xmm))
  
  ps_sd <- sample_data(data.frame(
    SampleID = rownames(OTU),
    X_ord = X_ord,
    X_ord_L = Xmm[,"X_ord.L", drop=TRUE],
    X_ord_Q = if("X_ord.Q" %in% colnames(Xmm)) Xmm[,"X_ord.Q"] else 0,
    Z = Z,
    gps_L = pL, gps_M = pM, gps_H = pH,
    gps_realized = pmax(1e-12, gps),
    w_stab = w_stab,
    w_trim = w_trim,
    s = s,
    r = r, alpha = a, kRare = kR,
    row.names = rownames(OTU),
    check.names = FALSE
  ))
  
  if(!build_tree){
    out <- phyloseq(ps_otu, ps_sd)
    if(return_centers) return(list(ps = out, centers = centers))
    return(out)
  }
  tr <- .build_midpoint_tree(OTU)
  out <- phyloseq(ps_otu, ps_sd, phy_tree(tr))
  if(return_centers) return(list(ps = out, centers = centers))
  return(out)
}

## =================== END FIXED-GRID FAST SIMULATOR (TREE OPTIONAL) =================== ##

################################## IPW- Permanova ####################################

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



################################## GLaD ##############################################

GLaD_eigen <- function(physeq, rho = 1, weighted = TRUE) {
  tree <- phy_tree(physeq)
  tip_labels <- tree$tip.label
  
  if (any(duplicated(tip_labels))) {
    tip_labels_fixed <- make.unique(tip_labels)
    tree$tip.label <- tip_labels_fixed
    phy_tree(physeq) <- tree
    if (taxa_are_rows(physeq)) {
      rownames(otu_table(physeq)) <- tip_labels_fixed
    } else {
      colnames(otu_table(physeq)) <- tip_labels_fixed
    }
  }
  
  if (!is.rooted(tree)) stop("The phylogenetic tree must be rooted.")
  
  A <- build_adjacency_matrix(tree)
  D <- Matrix::Diagonal(x = rowSums(A))
  L <- D - rho * A
  L <- L + Matrix::Diagonal(n = nrow(L), x = 1e-8)
  L_sparse <- as(L, "dgCMatrix")
  
  k <- min(100, nrow(L) - 1)
  eigs_res <- RSpectra::eigs_sym(L_sparse, k = k, which = "SM")
  keep <- which(abs(eigs_res$values) > 1e-8)
  lambda <- eigs_res$values[keep]
  U <- eigs_res$vectors[, keep]
  
  ntips <- Ntip(tree)
  nnodes <- Nnode(tree)
  tip_names <- tree$tip.label
  internal_ids <- (ntips + 1):(ntips + nnodes)
  internal_names <- paste0("Node_", internal_ids)
  all_node_names <- c(tip_names, internal_names)
  
  otu_mat <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  full_mat <- matrix(0, nrow = nsamples(physeq), ncol = length(all_node_names))
  rownames(full_mat) <- sample_names(physeq)
  colnames(full_mat) <- all_node_names
  full_mat[, tip_names] <- otu_mat[, tip_names, drop = FALSE]
  
  for (i in seq_along(internal_ids)) {
    node_id <- internal_ids[i]
    node_name <- internal_names[i]
    desc <- Descendants(tree, node_id, type = "tips")[[1]]
    desc_names <- tip_names[desc]
    full_mat[, node_name] <- rowSums(otu_mat[, desc_names, drop = FALSE])
  }
  
  rel_abund <- if (weighted) {
    sweep(full_mat, 1, rowSums(full_mat), FUN = "/")
  } else {
    (full_mat > 0) * 1
  }
  
  proj <- rel_abund %*% U
  scaled_proj <- sweep(proj, 2, sqrt(lambda), FUN = "/")
  
  dist(scaled_proj, method = "euclidean")
}

# Internal utility: build adjacency matrix from phylo tree
build_adjacency_matrix <- function(tree) {
  N <- max(tree$edge)
  A <- matrix(0, nrow = N, ncol = N)
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child  <- tree$edge[i, 2]
    A[parent, child] <- 1
    A[child, parent] <- 1
  }
  A
}


################################## Main ##############################################
set.seed(rep_id)

## -------------------- grids --------------------
Ns        <- c(40L, 60L, 80L)
Deltas    <- seq(0, 0.5, by = 0.05)
Ps_levels <- c(0.2, 1.2)  # weak/strong confounding
Dissims   <- c("UniFrac", "GUniFrac", "weighted UniFrac", "GLaD", "Bray-Curtis", "Aitchison", "Jaccard")

n_total   <- length(Ns) * length(Deltas) * length(Ps_levels) * length(Dissims)

## -------------------- helpers --------------------
as_square <- function(D) {
  if (inherits(D, "dist")) return(as.matrix(D))
  D <- as.matrix(D)
  if (nrow(D) != ncol(D)) stop("Distance matrix must be square.")
  (D + t(D)) / 2
}

to_dist <- function(D) {
  if (inherits(D, "dist")) return(D)
  as.dist(as_square(D))
}

## Aitchison distance: Euclidean distance on CLR of row-wise proportions (with small pseudocount)
aitchison_dist <- function(otu, pseudo = 1e-6) {
  # proportions
  rs <- rowSums(otu)
  P  <- sweep(otu, 1, pmax(rs, 1), "/")
  # small stabilization & renormalize
  P  <- P + pseudo
  P  <- sweep(P, 1, rowSums(P), "/")
  # CLR
  logP <- log(P)
  clrP <- logP - rowMeans(logP)
  dist(clrP)  # Euclidean on CLR
}

## GUniFrac extractor with robust name matching
pick_unifrac <- function(unifracs_array, name_candidates) {
  have <- dimnames(unifracs_array)[[3]]
  for (nm in name_candidates) {
    if (nm %in% have) return(unifracs_array[, , nm, drop = TRUE])
  }
  stop(sprintf("None of the UniFrac names found: %s; available: %s",
               paste(name_candidates, collapse = ", "),
               paste(have, collapse = ", ")))
}

## Dissimilarity switch (assumes GLaD() exists in your env when chosen)
make_distance <- function(otu, tr, psA, which) {
  switch(which,
         "UniFrac" = {
           uf <- GUniFrac::GUniFrac(otu, tr, alpha = c(0, 0.5, 1))$unifracs
           D  <- uf[,,"d_0"]
           to_dist(D)
         },
         "GUniFrac" = {
           uf <- GUniFrac::GUniFrac(otu, tr, alpha = c(0, 0.5, 1))$unifracs
           D  <- uf[,,"d_0.5"]
           to_dist(D)
         },
         "weighted UniFrac" = {
           uf <- GUniFrac::GUniFrac(otu, tr, alpha = c(0, 0.5, 1))$unifracs
           # Common names seen in GUniFrac outputs; this list covers typical builds
           D  <- uf[,,"d_1"]
           to_dist(D)
         },
         "GLaD" = {
           D = GLaD_eigen(psA)
           to_dist(D)
         },
         "Bray-Curtis" = vegan::vegdist(otu, method = "bray"),
         "Aitchison"   = vegan::vegdist(otu, method = "aitchison", pseudocount = 0.5),
         "Jaccard"     = vegan::vegdist(otu, method = "jaccard"),
         stop(sprintf("Unknown dissimilarity: %s", which))
  )
}

## Single analysis run (returns a named list of stats)
analyze_one <- function(D, meta) {
  # Ensure alignment
  meta <- meta[labels(D), , drop = FALSE]
  
  # Unadjusted (confounded)
  t1 <- try(vegan::adonis2(D ~ X_ord, data = meta, permutations = 999), silent = TRUE)
  F1 <- P1 <- NA_real_
  if (!inherits(t1, "try-error")) {
    F1 <- suppressWarnings(t1$F[1])
    P1 <- suppressWarnings(t1$`Pr(>F)`[1])
  }
  
  # Pooled (include Z)
  t2 <- try(vegan::adonis2(D ~ X_ord + Z, data = meta, permutations = 999), silent = TRUE)
  F2 <- P2 <- NA_real_
  if (!inherits(t2, "try-error")) {
    F2 <- suppressWarnings(t2$F[1])
    P2 <- suppressWarnings(t2$`Pr(>F)`[1])
  }
  
  # IPW (deconfounded)
  t3 <- try(ipw_permanova(as_square(D), meta$X_ord, meta$Z), silent = TRUE)
  F3 <- P3 <- NA_real_
  if (!inherits(t3, "try-error")) {
    F3 <- suppressWarnings(unname(t3$F_statistic))
    P3 <- suppressWarnings(unname(t3$p_value))
  }
  
  list(F_unadj = F1, p_unadj = P1,
       F_pool  = F2, p_pool  = P2,
       F_ipw   = F3, p_ipw   = P3)
}

## -------------------- main loop --------------------
results <- vector("list", n_total)
row_id  <- 0L

pb <- txtProgressBar(min = 0, max = n_total, style = 3)

## We'll vary the simulation seed per dataset for stability: base on rep_id and a running index
dataset_counter <- 0L

for (n in Ns) {
  for (Delta in Deltas) {
    for (ps in Ps_levels) {
      
      dataset_counter <- dataset_counter + 1L
      sim_seed <- as.integer(1e6 + rep_id * 1e4 + dataset_counter)
      
      ## Simulate one dataset for this (n, Delta, ps)
      psA <- ConMiseq_fixedGrid_Ordinal(
        n=n, N_pts=400,
        theta=c(L=-0.5, M=0.5), gamma=ps,          # strong X~Z
        r_slope=0, a_slope=0, k_slope=0,             # X has NO direct effect
        Delta_levels=Delta, bo_z=0.5,                      # no centroid effect
        build_tree=T,
        seed = sim_seed
      )
      ###### filter all 0 taxa and sample
      psA <- phyloseq::prune_taxa(phyloseq::taxa_sums(psA) > 0, psA)
      psA <- phyloseq::prune_samples(phyloseq::sample_sums(psA) > 0, psA)
      ####################################
      otu <- as(otu_table(psA), "matrix")
      if (taxa_are_rows(otu_table(psA))) otu <- t(otu)
      tr   <- phy_tree(psA)
      meta <- as(sample_data(psA), "data.frame")
      
      ## Make sure meta rows align to OTU rows for distances
      meta <- meta[rownames(otu), , drop = FALSE]
      
      ## Pre-compute UniFrac cube once and pass into accessor to avoid recompute per metric if preferred
      ## (Here we compute inside make_distance for clarity.)
      
      for (dname in Dissims) {
        row_id <- row_id + 1L
        
        stats <- try({
          D <- make_distance(otu, tr, psA, dname)
          analyze_one(D, meta)
        }, silent = TRUE)
        
        if (inherits(stats, "try-error")) {
          stats <- list(F_unadj = NA_real_, p_unadj = NA_real_,
                        F_pool  = NA_real_, p_pool  = NA_real_,
                        F_ipw   = NA_real_, p_ipw   = NA_real_)
        }
        
        results[[row_id]] <- data.frame(
          rep_id = rep_id,
          n = n,
          Delta = Delta,
          ps = ps,
          dissimilarity = dname,
          F_unadj = stats$F_unadj, p_unadj = stats$p_unadj,
          F_pool  = stats$F_pool,  p_pool  = stats$p_pool,
          F_ipw   = stats$F_ipw,   p_ipw   = stats$p_ipw,
          stringsAsFactors = FALSE
        )
        
        setTxtProgressBar(pb, row_id)
      }
    }
  }
}
close(pb)

## -------------------- save --------------------
out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outfile <- file.path(out_dir, sprintf("sim_results_rep_%03d.csv", rep_id))
outdf   <- do.call(rbind, results)
row.names(outdf) <- NULL
write.csv(outdf, outfile, row.names = FALSE)

cat("\nSaved:", outfile, "\n")
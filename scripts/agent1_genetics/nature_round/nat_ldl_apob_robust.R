#!/usr/bin/env Rscript
# ============================================================
# Weak-instrument-robust inference on the LDL-C minus ApoB contrast.
# Section C of nature_round/PREREGISTERED_READOUT.md.
#
# WHAT THE CONTRAST IS
#   delta = theta_LDL - theta_ApoB, in log OR per SD. The exposures are rotated so that
#   delta is a single coefficient: with beta_D = (b_LDL - b_ApoB)/2 and beta_S = b_LDL + b_ApoB,
#       b_LDL*theta_LDL + b_ApoB*theta_ApoB = beta_D*delta + beta_S*sigma,
#   where sigma = (theta_LDL + theta_ApoB)/2. The SNP-exposure covariance is carried through
#   the rotation: T' [Pcor * outer(se_j, se_j)] T, the phenocov_mvmr form. The rotation does
#   not change the model. It only makes the parameter of interest one coordinate, so that
#   inference on it is inference on a scalar.
#
# WHERE THE STATISTICS COME FROM (read, not recalled)
#   Patel, Lane & Burgess, arXiv:2408.09868v1, and their R code, github.com/ash-res/mvmr-weakiv.
#     S(theta)  = g'Omega^-1 g, chi2_J                                 (their section 3.2)
#     K(theta)  = g'Omega^-1 D (D'Omega^-1 D)^-1 D'Omega^-1 g, chi2_K  (section 3.3)
#     K + a*S vs (1+a)chi2_K + a*chi2_(J-K), Andrews (2018)            (section 3.4)
#     Distortion cutoff, a_tilde = sup_{W>c} (c-K)/S                   (their code)
#     Kleibergen-OH: Omega + kappa2, kappa2 from S(theta, kappa2) = J  (section 3.5)
#   Section 3.4 prints K* with (D'Omega^-1 D)^-1 in the middle. That cannot be right
#   dimensionally, and their code uses the inverse of the variance, which is identical to K.
#   We use K.
#   Their recommendation: report robust and non-robust sets with the cutoff; if AR is empty
#   the cause may be heterogeneity, and then read Kleibergen-OH. An AR set much tighter than
#   K and Andrews is to be distrusted (Davidson & MacKinnon 2014).
#
# WHAT IS OURS, NOT THEIRS
#   They give full-vector sets, with projection onto single coefficients. Projection from a
#   3-D set is conservative, so the subvector step here is a profile: at each fixed delta the
#   nuisance (sigma, theta_TG) is set to its continuously-updated value. Kleibergen (2005) and
#   Stock & Wright (2000) show the profiled K is chi2_1 and the profiled S is chi2_(J-K+1) when
#   the nuisance is strongly identified. That condition is checkable and is expected here:
#   with delta fixed, identifying the nuisance needs beta_S and beta_TG not to be collinear with
#   each other, and needs NOTHING from the LDL-C/ApoB difference. F_nuisance is reported, and
#   the projection sets (AR_projection, KOH_projection), which do not need the condition, are
#   reported beside them. The Andrews combination on the profile, (1+a)chi2_1 + a*chi2_(J-K), is
#   likewise our adaptation. None of this is taken on trust: --calibrate tests each one.
#
# PRE-SPECIFIED SET FOR THE VERDICT: KOH_subset, FinnGen, biochemistry arm.
#   Lipid-CAD instruments are heterogeneous (Q_A 506 on 148 df on UK Biobank CAD), and Patel et
#   al. show AR, K and Andrews lose coverage under heterogeneity. KOH reduces to K when kappa2=0.
#   The choice was fixed before any real result was computed.
#
# MVMR 0.4.8 HAS A qhet_mvmr DEFECT, SO IT CANNOT BE THE ESTIMATOR HERE
#   MVMR 0.4.8 (the TACC version) searches tau2 over (-10, 10), which allows negative weights and
#   is the wrong scale by orders of magnitude for summary data, and targets L-2, right only for
#   two exposures. The package's own 0.4.9 NEWS records both. qhet_fit() below follows Sanderson,
#   Spiller & Bowden (2021) as 0.4.9 does (tau2>=0, min Q = L-K) and is checked against 0.4.9 in
#   the selftest when that version is available. MVMR::qhet_mvmr is still called and its output
#   saved, labelled as known-defective, and it feeds nothing.
#
# SAMPLE OVERLAP. FinnGen does not overlap the exposures. UK Biobank CAD does, every statistic
#   here assumes it does not, and the rows carry a flag. They can inform sign, not the verdict.
#
# Modes: --selftest | --calibrate | (no flag) real run, which refuses to start unless
#        calibration_status.txt says PASSED for this exact version of this script.
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})
set.seed(42)

ALPHA    <- 0.05
COV_LO   <- 0.93   # readout section C
COV_HI   <- 0.97
BIG_MULT <- 1e6    # delta = center +- BIG_MULT*s stands in for the point at infinity
MAX_NA   <- 0.005  # coverage among replicates where a statistic computed is conditional on that;
                   # above this fraction of failures, coverage is not a readable number
SELECT_PASS_FRAC <- 0.875  # measured: share of candidate instruments passing z > 5.45 in sim_gamma

# ------------------------------------------------------------
# Problem setup
# ------------------------------------------------------------
# bx, bxse are J x K with the contrast pair in columns 1 and 2. rotate=FALSE keeps the original
# exposures, which the selftest uses to compare against the packages.
mk_prob <- function(bx, bxse, by, bys, pcor, rotate = TRUE) {
  bx <- as.matrix(bx); bxse <- as.matrix(bxse)
  J <- nrow(bx); K <- ncol(bx)
  bad <- function(why) list(degenerate = TRUE, reason = why, J = J, K = K)
  if (J <= K) return(bad(sprintf("J=%d not greater than K=%d", J, K)))
  if (anyNA(bx) || anyNA(bxse) || anyNA(by) || anyNA(bys)) return(bad("missing values"))
  if (any(bxse <= 0) || any(bys <= 0)) return(bad("non-positive standard error"))
  ev <- eigen(pcor, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-6) return(bad(sprintf("phenotypic correlation not positive definite (min eig %.2e)", min(ev))))
  Tm <- diag(K)
  if (rotate) {
    if (K < 2) return(bad("rotation needs two exposures"))
    Tm[, 1] <- 0; Tm[1, 1] <- 0.5; Tm[2, 1] <- -0.5
    Tm[, 2] <- 0; Tm[1:2, 2] <- 1
  }
  # SG holds Sigma*_j row-wise: J x K^2, entry (a,b) in position a + (b-1)K. Stored flat because
  # theta'Sigma_j theta is then a single matrix product, which is what makes the grid affordable.
  SG <- matrix(0, J, K^2)
  for (a in 1:K) for (b in 1:K) {
    v <- numeric(J)
    for (k in 1:K) for (l in 1:K) {
      w <- Tm[k, a] * Tm[l, b] * pcor[k, l]
      if (w != 0) v <- v + w * bxse[, k] * bxse[, l]
    }
    SG[, a + (b - 1L) * K] <- v
  }
  SGa <- lapply(1:K, function(a) SG[, a + (0:(K - 1L)) * K, drop = FALSE])
  list(degenerate = FALSE, J = J, K = K, bs = bx %*% Tm, by = as.numeric(by),
       s2y = as.numeric(bys)^2, SG = SG, SGa = SGa, Tm = Tm)
}

sub_prob <- function(p, i) {
  p$bs <- p$bs[i, , drop = FALSE]; p$by <- p$by[i]; p$s2y <- p$s2y[i]
  p$SG <- p$SG[i, , drop = FALSE]
  p$SGa <- lapply(p$SGa, function(M) M[i, , drop = FALSE])
  p$J <- length(i)
  p
}

om_fn  <- function(p, th, tau2 = 0) p$s2y + tau2 + as.numeric(p$SG %*% as.numeric(tcrossprod(th)))
g_fn   <- function(p, th) p$by - as.numeric(p$bs %*% th)
S_fn   <- function(p, th, tau2 = 0) { g <- g_fn(p, th); sum(g^2 / om_fn(p, th, tau2)) }
Sth_fn <- function(p, th) do.call(cbind, lapply(p$SGa, function(M) M %*% th))

# Solves A x = b after scaling A to a correlation matrix: bs columns differ in scale by orders
# of magnitude once delta is large. Returns NULL when singular so the caller emits NA.
solve_scaled <- function(A, b) {
  d <- sqrt(diag(A))
  if (any(!is.finite(d)) || any(d <= 0)) return(NULL)
  Ai <- A / tcrossprod(d)
  ch <- tryCatch(chol(Ai), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  if (min(diag(ch)) < 1e-7) return(NULL)
  backsolve(ch, backsolve(ch, b / d, transpose = TRUE)) / d
}

# The Kleibergen statistic and its parts. D is g-orthogonalised gamma-hat; its first-order role
# is that -2 D'Omega^-1 g is the exact gradient of S (checked numerically in the selftest).
kstat <- function(p, th, tau2 = 0) {
  g <- g_fn(p, th); Om <- om_fn(p, th, tau2); r <- g / Om
  D <- p$bs + Sth_fn(p, th) * r
  sc <- crossprod(D, r)
  info <- crossprod(D / sqrt(Om))
  x <- solve_scaled(info, sc)
  list(K = if (is.null(x)) NA_real_ else sum(sc * x), S = sum(g * r), D = D, g = g, Om = Om,
       sc = sc, info = info)
}

# Gauss-Newton on S with step (D'Omega^-1 D)^-1 D'Omega^-1 g, the theta* of their section 3.4,
# with Armijo backtracking. Convergence is in chi2 units (the score statistic for the free
# coordinates), which is scale free, so the same tolerance holds at delta = 1 and delta = 1e6.
gn_min <- function(p, th0, free, tau2 = 0, maxit = 60L, tol = 1e-9) {
  th <- th0; S <- S_fn(p, th, tau2); conv <- FALSE; dec <- NA_real_
  for (it in seq_len(maxit)) {
    g <- g_fn(p, th); Om <- om_fn(p, th, tau2); r <- g / Om
    D <- (p$bs + Sth_fn(p, th) * r)[, free, drop = FALSE]
    sc <- crossprod(D, r)
    step <- solve_scaled(crossprod(D / sqrt(Om)), sc)
    if (is.null(step)) break
    dec <- sum(sc * step)
    if (!is.finite(dec)) break
    if (dec < tol) { conv <- TRUE; break }
    lam <- 1; ok <- FALSE
    for (ls in 1:30) {
      tn <- th; tn[free] <- th[free] + lam * step
      Sn <- S_fn(p, tn, tau2)
      if (is.finite(Sn) && Sn <= S - 1e-4 * lam * 2 * dec) { ok <- TRUE; break }
      lam <- lam / 2
    }
    if (!ok) break
    th <- tn; S <- Sn
  }
  if (!conv) {
    # Fallback: when GN cannot make progress, typically a nearly flat direction. Only a score
    # statistic close to zero counts as converged; otherwise the caller gets conv=FALSE and NA.
    sc0 <- pmax(abs(th[free]), 1e-3)
    ob <- function(u) { tn <- th; tn[free] <- th[free] + sc0 * u; S_fn(p, tn, tau2) }
    gr <- function(u) {
      tn <- th; tn[free] <- th[free] + sc0 * u
      g <- g_fn(p, tn); Om <- om_fn(p, tn, tau2); r <- g / Om
      D <- (p$bs + Sth_fn(p, tn) * r)[, free, drop = FALSE]
      as.numeric(-2 * crossprod(D, r)) * sc0
    }
    fit <- tryCatch(nlminb(rep(0, length(free)), ob, gr,
                           control = list(iter.max = 300, eval.max = 600, rel.tol = 1e-13)),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      tn <- th; tn[free] <- th[free] + sc0 * fit$par
      if (is.finite(fit$objective) && fit$objective <= S + 1e-10) { th <- tn; S <- fit$objective }
    }
    g <- g_fn(p, th); Om <- om_fn(p, th, tau2); r <- g / Om
    D <- (p$bs + Sth_fn(p, th) * r)[, free, drop = FALSE]
    sc <- crossprod(D, r); x <- solve_scaled(crossprod(D / sqrt(Om)), sc)
    dec <- if (is.null(x)) NA_real_ else sum(sc * x)
    conv <- is.finite(dec) && dec < 1e-5
  }
  list(th = th, S = S, conv = conv, dec = dec)
}

wls_init <- function(p, delta, tau2 = 0) {
  free <- 2:p$K
  X <- p$bs[, free, drop = FALSE]; y <- p$by - delta * p$bs[, 1]
  nu <- rep(0, length(free))
  for (it in 1:4) {
    w <- 1 / om_fn(p, c(delta, nu), tau2)
    x <- solve_scaled(crossprod(X * sqrt(w)), crossprod(X, w * y))
    if (is.null(x)) return(NULL)
    nu <- as.numeric(x)
  }
  nu
}

# S profiled over the nuisance at fixed delta. Two starts: the continuation from the neighbouring
# grid point and a fresh weighted regression, because a warm start can stay on a wrong branch.
prof_best <- function(p, delta, tau2 = 0, warm = NULL) {
  free <- 2:p$K
  best <- NULL
  st <- wls_init(p, delta, tau2)
  for (nu0 in list(warm, st)) {
    if (is.null(nu0)) next
    f <- gn_min(p, c(delta, nu0), free, tau2)
    if (f$conv && (is.null(best) || f$S < best$S)) best <- f
  }
  if (is.null(best)) return(list(conv = FALSE, S = NA_real_, nu = NULL))
  list(conv = TRUE, S = best$S, nu = best$th[free])
}

# kappa2 at fixed delta. Target is J-K+1: the nuisance has been fitted, so that is the mean of
# the profiled S at the truth (the J of their section 3.5 applies to the full-vector statistic).
solve_tau <- function(p, delta, pr0, warm = NULL) {
  target <- p$J - p$K + 1
  if (!pr0$conv) return(list(conv = FALSE))
  if (pr0$S <= target) return(list(conv = TRUE, tau2 = 0, nu = pr0$nu, S = pr0$S))
  cur <- if (is.null(warm)) pr0$nu else warm
  fail <- FALSE
  h <- function(t) {
    pr <- prof_best(p, delta, t, cur)
    if (!pr$conv) { fail <<- TRUE; return(NA_real_) }
    cur <<- pr$nu
    pr$S - target
  }
  lo <- 0
  hi <- max(median(om_fn(p, c(delta, pr0$nu))), 1e-300)
  hh <- h(hi); n <- 0L
  while (!fail && hh > 0 && n < 200L) { lo <- hi; hi <- hi * 4; hh <- h(hi); n <- n + 1L }
  if (fail || hh > 0) return(list(conv = FALSE))
  rt <- tryCatch(uniroot(h, c(lo, hi), f.upper = hh, tol = hi * 1e-7, maxiter = 100L),
                 error = function(e) NULL)
  if (is.null(rt) || fail) return(list(conv = FALSE))
  pr <- prof_best(p, delta, rt$root, cur)
  if (!pr$conv) return(list(conv = FALSE))
  list(conv = TRUE, tau2 = rt$root, nu = pr$nu, S = pr$S)
}

# kappa2 at a full theta (their section 3.5). S is convex decreasing in kappa2, so Newton from
# zero is monotone and needs no bracket.
tau_newton <- function(p, th, target = p$J) {
  g <- g_fn(p, th); Om <- om_fn(p, th); g2 <- g^2; t <- 0
  for (i in 1:100) {
    v <- Om + t; S <- sum(g2 / v)
    if (S <= target * (1 + 1e-10)) break
    t <- t + (S - target) / sum(g2 / v^2)
  }
  t
}

# Everything measured at one delta.
eval_delta <- function(p, delta, warm = list(), tau_global = NULL, Vdd = NULL, do_oh = TRUE) {
  out <- list(delta = delta, S0 = NA_real_, K0 = NA_real_, Koh = NA_real_, tauoh = NA_real_,
              W = NA_real_, Stg = NA_real_, nu0 = NULL, nuoh = NULL, nutg = NULL)
  pr <- prof_best(p, delta, 0, warm$nu0)
  if (!pr$conv) return(out)
  out$S0 <- pr$S; out$nu0 <- pr$nu
  th <- c(delta, pr$nu)
  out$K0 <- kstat(p, th, 0)$K
  if (!is.null(Vdd)) {
    Om <- om_fn(p, th); G <- p$bs
    x <- solve_scaled(crossprod(G / sqrt(Om)), crossprod(G, g_fn(p, th) / Om))
    if (!is.null(x)) out$W <- x[1]^2 / Vdd
  }
  if (do_oh) {
    so <- solve_tau(p, delta, pr, warm$nuoh)
    if (isTRUE(so$conv)) {
      out$tauoh <- so$tau2; out$nuoh <- so$nu
      out$Koh <- kstat(p, c(delta, so$nu), so$tau2)$K
    }
  }
  if (!is.null(tau_global)) {
    pt <- prof_best(p, delta, tau_global, warm$nutg)
    if (pt$conv) { out$Stg <- pt$S; out$nutg <- pt$nu }
  }
  out
}

# Kleibergen-OH projection: min over the nuisance of K(theta) with kappa2 refitted at each theta,
# against chi2_K. Valid without strong nuisance, conservative because of the K degrees of freedom.
# K is not convex in the nuisance, so several starts; a missed lower value would make the reported
# projection set too small, which is why it is a sensitivity row and not the primary set.
koh_proj_stat <- function(p, delta, nu_start, crit) {
  fn <- function(nu) {
    th <- c(delta, nu); kstat(p, th, tau_newton(p, th))$K
  }
  v0 <- fn(nu_start)
  if (!is.na(v0) && v0 <= crit) return(v0)
  ks <- kstat(p, c(delta, nu_start), tau_newton(p, c(delta, nu_start)))
  free <- 2:p$K
  x <- tryCatch(solve(ks$info[free, free, drop = FALSE]), error = function(e) NULL)
  se <- if (is.null(x)) pmax(abs(nu_start), 1e-3) else sqrt(pmax(diag(x), 1e-30))
  starts <- list(rep(0, length(free)))
  for (i in seq_along(free)) for (s in c(-2, 2)) {
    u <- rep(0, length(free)); u[i] <- s; starts <- c(starts, list(u))
  }
  best <- if (is.na(v0)) Inf else v0
  for (u0 in starts) {
    fit <- tryCatch(nlminb(u0, function(u) {
      v <- fn(nu_start + se * u); if (is.na(v)) 1e10 else v
    }, control = list(iter.max = 40, eval.max = 80)), error = function(e) NULL)
    if (!is.null(fit) && fit$objective < best) best <- fit$objective
    if (best <= crit) break
  }
  if (is.finite(best)) best else NA_real_
}

# ------------------------------------------------------------
# Andrews (2018) linear combination, with df1 = 1 after profiling
# ------------------------------------------------------------
# P((1+a) X + a Y <= q) for X ~ chi2_1, Y ~ chi2_df2, by quadrature: exact and reproducible where
# their code simulates. x = z^2 removes the chi2_1 singularity at zero.
plc <- function(q, a, df2) {
  if (a <= 0) return(pchisq(q, 1))
  zmax <- sqrt(q / (1 + a))
  integrate(function(z) 2 * dnorm(z) * pchisq((q - (1 + a) * z^2) / a, df2),
            0, zmax, rel.tol = 1e-9, subdivisions = 400L, stop.on.error = FALSE)$value
}
qlc <- function(prob, a, df2) {
  if (a <= 0) return(qchisq(prob, 1))
  m <- (1 + a) + a * df2
  uniroot(function(q) plc(q, a, df2) - prob, c(1e-8, m + 40 * (1 + a) + 40 * a * sqrt(2 * df2) + 10),
          tol = 1e-8, maxiter = 200L)$root
}
a_for_gamma <- function(gamma, df2) {
  c1 <- qchisq(1 - ALPHA, 1)
  uniroot(function(a) plc(c1, a, df2) - (1 - ALPHA - gamma), c(1e-10, 50), tol = 1e-10,
          maxiter = 200L)$root
}
distortion_cutoff <- function(K0, S0, W, df2, gamma_min) {
  c1 <- qchisq(1 - ALPHA, 1)
  out <- which(!is.na(W) & !is.na(K0) & !is.na(S0) & W > c1)
  a <- if (length(out)) max(0, max((c1 - K0[out]) / S0[out])) else 0
  gt <- (1 - ALPHA) - plc(c1, a, df2)
  max(gt, gamma_min)
}

# Critical values for one problem; a single place so that simulation and real run cannot drift.
crit_table <- function(p, gamma_min = 0.05) {
  df2 <- p$J - p$K
  a <- a_for_gamma(gamma_min, df2)
  list(a = a,
       AR_subset = qchisq(1 - ALPHA, p$J - p$K + 1),
       AR_projection = qchisq(1 - ALPHA, p$J),
       K_subset = qchisq(1 - ALPHA, 1),
       LC_subset = qlc(1 - ALPHA, a, df2),
       KOH_subset = qchisq(1 - ALPHA, 1),
       KOH_projection = qchisq(1 - ALPHA, p$K),
       WALD_GMM = qchisq(1 - ALPHA, 1))
}

stat_vec <- function(e, a) c(AR_subset = e$S0, AR_projection = e$S0, K_subset = e$K0,
                             LC_subset = e$K0 + a * e$S0, KOH_subset = e$Koh, WALD_GMM = e$W)

# ------------------------------------------------------------
# Non-robust and heterogeneity-aware point estimators
# ------------------------------------------------------------
# Weighted regression with weights 1/se_y^2: what mr_mvivw computes in rotated coordinates. The
# random-effects standard error uses max(1, residual SE), the package default for more than 3 SNPs.
ivw_fit <- function(p) {
  w <- 1 / p$s2y
  A <- crossprod(p$bs * sqrt(w)); x <- solve_scaled(A, crossprod(p$bs, w * p$by))
  if (is.null(x)) return(NULL)
  th <- as.numeric(x)
  res <- p$by - as.numeric(p$bs %*% th)
  Q <- sum(w * res^2)
  V <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  list(th = th, V = V * max(1, Q / (p$J - p$K)), Q = Q)
}

fit_cue <- function(p, starts) {
  best <- NULL
  for (s in starts) {
    f <- gn_min(p, s, 1:p$K, 0)
    if (f$conv && (is.null(best) || f$S < best$S)) best <- f
  }
  if (is.null(best)) return(NULL)
  Om <- om_fn(p, best$th)
  V <- tryCatch(solve(crossprod(p$bs / sqrt(Om))), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  list(th = best$th, S = best$S, V = V)
}

# Sanderson, Spiller & Bowden (2021): min_theta Q(theta, tau2) = L - K, tau2 >= 0.
qhet_fit <- function(p, starts) {
  target <- p$J - p$K
  mn <- function(t, st) {
    best <- NULL
    for (s in st) {
      f <- gn_min(p, s, 1:p$K, t)
      if (f$conv && (is.null(best) || f$S < best$S)) best <- f
    }
    best
  }
  f0 <- mn(0, starts)
  if (is.null(f0)) return(NULL)
  if (f0$S <= target) return(list(th = f0$th, tau2 = 0, S = f0$S))
  cur <- f0$th; fail <- FALSE
  h <- function(t) {
    f <- mn(t, list(cur))
    if (is.null(f)) { fail <<- TRUE; return(NA_real_) }
    cur <<- f$th
    f$S - target
  }
  lo <- 0; hi <- max(median(om_fn(p, f0$th)), 1e-300)
  hh <- h(hi); n <- 0L
  while (!fail && hh > 0 && n < 200L) { lo <- hi; hi <- hi * 4; hh <- h(hi); n <- n + 1L }
  if (fail || hh > 0) return(NULL)
  rt <- tryCatch(uniroot(h, c(lo, hi), f.upper = hh, tol = hi * 1e-7, maxiter = 100L),
                 error = function(e) NULL)
  if (is.null(rt) || fail) return(NULL)
  f <- mn(rt$root, list(cur, f0$th))
  if (is.null(f)) return(NULL)
  list(th = f$th, tau2 = rt$root, S = f$S)
}

# Sanderson-Windmeijer conditional F for one column on others (their section 2.4), by the same
# machinery. With the outcome set to zero, S at theta=(1,-c) is exactly their quadratic form.
cond_F <- function(p, target, others) {
  keep <- c(target, others); K <- length(keep)
  q <- p
  q$K <- K
  q$bs <- p$bs[, keep, drop = FALSE]
  q$by <- rep(0, p$J)
  q$s2y <- rep(1e-300, p$J)
  ii <- as.vector(outer(keep, keep, function(a, b) a + (b - 1L) * p$K))
  q$SG <- p$SG[, ii, drop = FALSE]
  q$SGa <- lapply(1:K, function(a) q$SG[, a + (0:(K - 1L)) * K, drop = FALSE])
  if (K == 1) return(S_fn(q, 1) / p$J)
  # With delta fixed at 1 the free coefficients are the -c of their section 2.4.
  pr <- prof_best(q, 1, 0, NULL)
  if (!pr$conv) return(NA_real_)
  pr$S / (p$J - K + 1)
}

# ------------------------------------------------------------
# Grid inversion. The grid is uniform in an angle, so it covers the whole real line.
# ------------------------------------------------------------
# delta = center + s*tan(phi). A set cannot be clipped by an edge, because there is none: the
# endpoints stand for infinity, and acceptance there is what "unbounded" means. On the projective
# line +infinity and -infinity are the same point, since S depends on delta only through a ratio
# at that limit, so an unbounded set is normally two-sided: (-inf, a] U [b, inf).
make_grid <- function(center, s, n) {
  eta <- pi / (2 * n)
  phi <- seq(-pi / 2 + eta, pi / 2 - eta, length.out = n)
  far <- s * 10^seq(log10(2 * n / pi) + 0.2, log10(BIG_MULT) - 0.2, length.out = 8)
  sort(unique(c(center + s * tan(phi), center + far, center - far,
                center - s * BIG_MULT, center + s * BIG_MULT, 0)))
}

# Interval end at the sign change in (a, b).
find_root <- function(fn, crit, a, b, va, vb, s) {
  f <- function(x) { v <- fn(x); if (is.na(v)) NA_real_ else v - crit }
  rt <- tryCatch(uniroot(f, c(a, b), f.lower = va - crit, f.upper = vb - crit,
                         tol = 1e-6 * s, maxiter = 100L), error = function(e) NULL)
  if (is.null(rt)) return(if (va <= crit) a else b)
  rt$root
}

invert_stat <- function(deltas, vals, crit, fn, s, zero_val, refine = TRUE, max_na = 0.005) {
  res <- list(status = "OK", class = NA_character_, intervals = matrix(numeric(0), 0, 2),
              includes_zero = NA, frac_na = mean(is.na(vals)))
  if (res$frac_na > max_na || is.na(vals[1]) || is.na(vals[length(vals)]) || is.na(zero_val)) {
    res$status <- "NUMERIC_FAILURE"
    return(res)
  }
  keep <- !is.na(vals); d <- deltas[keep]; v <- vals[keep]
  # A narrow acceptance region can fall between grid points, and K is prone to them (it is zero
  # at every stationary point of S, maxima included), so each local minimum above crit is polished.
  if (refine && length(v) >= 3) {
    im <- which(v[-c(1, length(v))] < v[-c(length(v) - 1, length(v))] &
                  v[-c(1, length(v))] < v[-c(1, 2)] & v[-c(1, length(v))] > crit) + 1L
    if (length(im) > 50L) im <- im[order(v[im])[1:50]]
    add_d <- numeric(0); add_v <- numeric(0)
    for (i in im) {
      o <- tryCatch(optimize(function(x) { y <- fn(x); if (is.na(y)) 1e12 else y },
                             c(d[i - 1], d[i + 1]), tol = 1e-6 * s),
                    error = function(e) NULL)
      if (!is.null(o) && o$objective <= crit) { add_d <- c(add_d, o$minimum); add_v <- c(add_v, o$objective) }
    }
    if (length(add_d)) { o <- order(c(d, add_d)); d <- c(d, add_d)[o]; v <- c(v, add_v)[o] }
  }
  acc <- v <= crit; n <- length(v)
  if (!any(acc)) { res$class <- "empty"; res$includes_zero <- FALSE; return(res) }
  if (all(acc)) {
    res$class <- "whole_line"; res$includes_zero <- TRUE
    res$intervals <- matrix(c(-Inf, Inf), 1, 2)
    return(res)
  }
  rr <- rle(acc); ends <- cumsum(rr$lengths); starts <- ends - rr$lengths + 1L
  iv <- matrix(NA_real_, sum(rr$values), 2); k <- 0L
  for (r in which(rr$values)) {
    k <- k + 1L; i <- starts[r]; j <- ends[r]
    iv[k, 1] <- if (i == 1L) -Inf else find_root(fn, crit, d[i - 1], d[i], v[i - 1], v[i], s)
    iv[k, 2] <- if (j == n) Inf else find_root(fn, crit, d[j], d[j + 1], v[j], v[j + 1], s)
  }
  res$intervals <- iv
  res$class <- if (any(!is.finite(iv))) "unbounded" else "bounded"
  res$includes_zero <- zero_val <= crit
  res
}

fmt_intervals <- function(iv) {
  if (!nrow(iv)) return("empty")
  paste(apply(iv, 1, function(r) sprintf("[%s, %s]",
      if (is.infinite(r[1])) "-Inf" else sprintf("%.4f", r[1]),
      if (is.infinite(r[2])) "Inf" else sprintf("%.4f", r[2]))), collapse = " U ")
}

# ------------------------------------------------------------
# One full analysis of one arm and one outcome
# ------------------------------------------------------------
robust_analysis <- function(p, grid_n = 2001L, gamma_min = 0.05, proj_every = 10L, label = "") {
  iv <- ivw_fit(p)
  if (is.null(iv)) return(list(degenerate = TRUE, reason = "IVW design singular"))
  dfree <- p$J - p$K
  ct <- crit_table(p, gamma_min)

  cue <- fit_cue(p, list(iv$th))
  # A second CUE start from the delta=0 profile: S can have more than one mode under collinearity.
  z0 <- prof_best(p, 0, 0, NULL)
  if (z0$conv) {
    cue2 <- fit_cue(p, list(c(0, z0$nu)))
    if (!is.null(cue2) && (is.null(cue) || cue2$S < cue$S)) cue <- cue2
  }
  # S can decrease all the way to infinity when delta is unidentified, and then no finite CUE or
  # qhet exists. That is a result and not a failure: the robust sets still stand and those
  # estimators read NA.
  qh <- qhet_fit(p, if (!is.null(cue)) list(cue$th, iv$th) else list(iv$th))
  Vdd <- if (!is.null(cue)) cue$V[1, 1]

  s <- sqrt(iv$V[1, 1])
  if (!is.finite(s) || s <= 0) return(list(degenerate = TRUE, reason = "contrast SE undefined"))
  center <- if (!is.null(cue)) cue$th[1] else iv$th[1]
  grid <- make_grid(center, s, grid_n)

  # Sweep outward from the centre; the nuisance path is nearly linear in delta (it tends to a
  # straight line at infinity), so a secant predictor keeps the warm start good across big gaps.
  ng <- length(grid)
  ev <- vector("list", ng)
  ic <- which.min(abs(grid - center))
  sweep <- function(idx) {
    prev <- list(); pd <- NA_real_; prev2 <- list(); pd2 <- NA_real_
    for (i in idx) {
      d <- grid[i]; warm <- list()
      for (nm in c("nu0", "nuoh", "nutg")) {
        a <- prev[[nm]]; b <- prev2[[nm]]
        warm[[nm]] <- if (!is.null(a) && !is.null(b) && is.finite(pd2) && pd != pd2)
          a + (a - b) * (d - pd) / (pd - pd2) else a
      }
      e <- eval_delta(p, d, warm, tau_global = if (!is.null(qh)) qh$tau2 else 0, Vdd = Vdd, do_oh = TRUE)
      ev[[i]] <<- e
      if (!is.null(e$nu0)) { prev2 <- prev; pd2 <- pd; prev <- e; pd <- d }
    }
  }
  sweep(ic:ng)
  if (ic > 1L) sweep((ic - 1L):1L)

  col <- function(nm) vapply(ev, function(e) as.numeric(e[[nm]])[1], numeric(1))
  S0 <- col("S0"); K0 <- col("K0"); Koh <- col("Koh"); W <- col("W"); Stg <- col("Stg"); tauoh <- col("tauoh")
  nearest <- function(x) ev[[which.min(abs(grid - x))]]
  warm_of <- function(x) { e <- nearest(x); list(nu0 = e$nu0, nuoh = e$nuoh, nutg = e$nutg) }
  fns <- list(
    AR_subset = function(x) eval_delta(p, x, warm_of(x), do_oh = FALSE)$S0,
    AR_projection = function(x) eval_delta(p, x, warm_of(x), do_oh = FALSE)$S0,
    K_subset = function(x) eval_delta(p, x, warm_of(x), do_oh = FALSE)$K0,
    LC_subset = function(x) { e <- eval_delta(p, x, warm_of(x), do_oh = FALSE); e$K0 + ct$a * e$S0 },
    KOH_subset = function(x) eval_delta(p, x, warm_of(x), do_oh = TRUE)$Koh,
    WALD_GMM = function(x) eval_delta(p, x, warm_of(x), Vdd = Vdd, do_oh = FALSE)$W)
  vals <- list(AR_subset = S0, AR_projection = S0, K_subset = K0, LC_subset = K0 + ct$a * S0,
               KOH_subset = Koh, WALD_GMM = W)
  iz <- which(grid == 0)
  sets <- lapply(names(fns), function(nm)
    invert_stat(grid, vals[[nm]], ct[[nm]], fns[[nm]], s, vals[[nm]][iz]))
  names(sets) <- names(fns)

  # Projection on a thinned grid; it is costly and only a sensitivity row.
  ip <- sort(unique(c(seq(1L, ng, by = proj_every), iz, 1L, ng)))
  pfn <- function(x) {
    e <- eval_delta(p, x, warm_of(x), do_oh = FALSE)
    if (is.null(e$nu0)) return(NA_real_)
    koh_proj_stat(p, x, e$nu0, ct$KOH_projection)
  }
  pv <- vapply(ip, function(i) if (is.null(ev[[i]]$nu0)) NA_real_ else
    koh_proj_stat(p, grid[i], ev[[i]]$nu0, ct$KOH_projection), numeric(1))
  sets$KOH_projection <- invert_stat(grid[ip], pv, ct$KOH_projection, pfn, s,
                                     pv[which(ip == iz)], refine = FALSE)

  # Does a component contain a minimum of the profiled S? K cannot tell a maximum of S from a
  # minimum, so a K component with no S minimum in it is a known artefact of the statistic.
  lm <- which(diff(sign(diff(Stg))) > 0) + 1L
  smins <- grid[lm]
  in_comp <- function(iv, x) any(iv[, 1] <= x & x <= iv[, 2])
  zero_comp_has_min <- function(st) {
    if (st$status != "OK" || !isTRUE(st$includes_zero) || !nrow(st$intervals)) return(NA)
    z <- st$intervals[st$intervals[, 1] <= 0 & 0 <= st$intervals[, 2], , drop = FALSE]
    if (!nrow(z)) return(NA)
    any(vapply(smins, function(x) in_comp(z, x), logical(1)))
  }

  rows <- rbindlist(lapply(names(sets), function(nm) {
    st <- sets[[nm]]
    role <- if (nm == "KOH_subset") "PRE_SPECIFIED_FOR_VERDICT" else if (nm == "WALD_GMM") "NON_ROBUST" else "REPORTED"
    fin <- st$intervals[is.finite(st$intervals)]
    data.table(statistic = nm, role = role, status = st$status, set_class = st$class,
               n_intervals = nrow(st$intervals),
               intervals = if (st$status == "OK") fmt_intervals(st$intervals) else NA_character_,
               lower = if (nrow(st$intervals)) min(st$intervals) else NA_real_,
               upper = if (nrow(st$intervals)) max(st$intervals) else NA_real_,
               includes_zero = st$includes_zero, frac_grid_failed = st$frac_na, critical_value = ct[[nm]],
               zero_component_has_S_minimum = if (grepl("^K", nm)) zero_comp_has_min(st) else NA,
               set_sign = if (st$status == "OK" && st$class == "bounded" && !isTRUE(st$includes_zero)) {
                 if (all(fin > 0)) "+" else if (all(fin < 0)) "-" else "mixed" } else NA_character_)
  }))

  gh <- distortion_cutoff(K0, S0, W, dfree, gamma_min)
  ests <- list(
    J = p$J, K = p$K,
    ivw_contrast = iv$th[1], ivw_contrast_se_re = s, Q_ivw = iv$Q,
    cue_contrast = if (!is.null(cue)) cue$th[1] else NA_real_,
    cue_contrast_se = if (!is.null(cue)) sqrt(Vdd) else NA_real_,
    S_min = if (!is.null(cue)) cue$S else NA_real_,
    S_min_p_chi2_J_minus_K = if (!is.null(cue)) pchisq(cue$S, dfree, lower.tail = FALSE) else NA_real_,
    qhet_contrast = if (!is.null(qh)) qh$th[1] else NA_real_,
    qhet_tau2 = if (!is.null(qh)) qh$tau2 else NA_real_,
    qhet_tau_over_median_sey = if (!is.null(qh)) sqrt(qh$tau2 / median(p$s2y)) else NA_real_,
    # S_prof at infinity divided by its df is the Sanderson-Windmeijer conditional F of D given
    # the nuisance, which is why AR unboundedness and F are the same event.
    F_contrast_given_nuisance = cond_F(p, 1L, 2:p$K),
    F_nuisance_min = min(vapply(2:p$K, function(k)
      cond_F(p, k, setdiff(2:p$K, k)), numeric(1))),
    S_prof_at_inf = S0[ng],
    koh_tau2_at_zero = tauoh[iz], lc_a = ct$a, lc_gamma_min = gamma_min,
    distortion_cutoff = gh,
    # The cutoff shrinks with S, so heterogeneity can make it look reassuring for the wrong reason.
    distortion_cutoff_interpretable = !is.null(cue) && cue$S <= qchisq(1 - ALPHA, dfree))
  list(degenerate = FALSE, sets = rows, ests = ests, th_ivw = iv$th,
       grid = data.table(delta = grid, S0 = S0, K0 = K0, Koh = Koh, W = W, Stg = Stg, tauoh = tauoh))
}

# ------------------------------------------------------------
# Simulation
# ------------------------------------------------------------
# A scenario gives the generating parameters; sim_gamma/sim_noise are separate so the plasmode
# can reuse the noise part with the observed SNP effects as truth.
sim_gamma <- function(sc) {
  J0 <- if (isTRUE(sc$select)) round(sc$J / SELECT_PASS_FRAC) else sc$J
  f <- runif(J0, 0.05, 0.5)
  base <- 1 / sqrt(2 * f * (1 - f))
  sx <- base / sqrt(sc$nx)
  # Lead-variant z-scores for lipids: a floor near genome-wide significance and a long right tail.
  # In the selection scenario there is no floor and a smaller mean, so that many instruments sit
  # near the threshold and it bites (winner's curse); selftest section 6 checks J stays near sc$J.
  # zscale sets instrument strength without touching the 0.97 correlation of the SNP effects.
  zs <- if (is.null(sc$zscale)) 1 else sc$zscale
  z0 <- if (isTRUE(sc$select)) rexp(J0, 1 / 5) + rbinom(J0, 1, 0.08) * rexp(J0, 1 / 40) else
    zs * (6 + rexp(J0, 1 / 9) + rbinom(J0, 1, 0.08) * rexp(J0, 1 / 40))
  z0 <- z0 * sample(c(-1, 1), J0, TRUE)
  g1 <- z0 * sx
  # Exact sample correlation sc$corg between the first two columns.
  if (sc$corg >= 1) g2 <- g1 else {
    e <- rnorm(J0); e <- resid(lm(e ~ g1)); e <- e / sd(e) * sd(g1)
    g2 <- sc$corg * g1 + sqrt(1 - sc$corg^2) * e
  }
  G <- cbind(g1, g2)
  if (sc$K == 3L) {
    g3 <- 0.3 * g1 + sqrt(1 - 0.3^2) * sd(g1) * rnorm(J0)
    G <- cbind(G, g3)
  }
  list(G = unname(G), sx = sx, base = base)
}

scenario_rho <- function(sc, rho = sc$rho) {
  R <- diag(sc$K); R[1, 2] <- R[2, 1] <- rho
  if (sc$K == 3L) { R[1, 3] <- R[3, 1] <- 0.25; R[2, 3] <- R[3, 2] <- 0.35 }
  R
}

sim_rep <- function(sc, seed, gamma_min = 0.05) {
  set.seed(seed)
  gm <- sim_gamma(sc)
  G <- gm$G; K <- sc$K
  sxm <- gm$sx %o% rep(1, K)
  sxm[, 2] <- sxm[, 2] * 1.01
  L <- chol(scenario_rho(sc))
  bx <- G + (matrix(rnorm(length(G)), nrow(G), K) %*% L) * sxm
  sy <- gm$base / sqrt((sc$ncase + sc$nctrl) * (sc$ncase / (sc$ncase + sc$nctrl)) * (sc$nctrl / (sc$ncase + sc$nctrl)))
  tau2 <- sc$tau_mult * median(sy^2)
  nu <- if (sc$tau_mult == 0) 0 else if (identical(sc$het, "t3"))
    sqrt(tau2) * rt(nrow(G), 3) / sqrt(3) else rnorm(nrow(G), 0, sqrt(tau2))
  by <- as.numeric(G %*% unlist(sc$theta)) + nu + rnorm(nrow(G), 0, sy)
  if (isTRUE(sc$select)) {
    keep <- apply(abs(bx / sxm), 1, max) > 5.45
    bx <- bx[keep, , drop = FALSE]; sxm <- sxm[keep, , drop = FALSE]; by <- by[keep]; sy <- sy[keep]
  }
  pc <- if (is.null(sc$rho_analysis) || is.na(sc$rho_analysis)) scenario_rho(sc) else scenario_rho(sc, sc$rho_analysis)
  p <- mk_prob(bx, sxm, by, sy, pc)
  test_prob(p, sc$delta, sc$rep_proj, gamma_min)
}

# Fill in test results at the true value, zero, and infinity. Coverage is membership of the truth,
# so it needs the statistic only there; full inversion is exercised separately (check_inversion).
test_prob <- function(p, delta_true, do_proj = FALSE, gamma_min = 0.05) {
  stats <- c("AR_subset", "AR_projection", "K_subset", "LC_subset", "KOH_subset", "WALD_GMM")
  nm <- c(paste0("cover_", c(stats, "KOH_projection", "IVW_RE")),
          paste0("accept0_", stats), paste0("acceptinf_", stats), "F_D", "S_min_over_df", "tau2hat", "J", "ok")
  out <- setNames(rep(NA_real_, length(nm)), nm)
  if (isTRUE(p$degenerate)) { out["ok"] <- 0; return(out) }
  iv <- ivw_fit(p)
  if (is.null(iv)) { out["ok"] <- 0; return(out) }
  ct <- crit_table(p, gamma_min)
  # No CUE only rules out the Wald row: dropping the replicate would condition coverage on the
  # very event that identification has failed, the case a robust set exists for.
  cue <- fit_cue(p, list(iv$th))
  s <- sqrt(iv$V[1, 1])
  ev <- function(d) stat_vec(eval_delta(p, d, list(), Vdd = if (!is.null(cue)) cue$V[1, 1], do_oh = TRUE), ct$a)
  v1 <- ev(delta_true); v0 <- ev(0); vi <- ev(iv$th[1] + s * BIG_MULT)
  cr <- unlist(ct[stats])
  out[paste0("cover_", stats)] <- as.numeric(v1 <= cr)
  out[paste0("accept0_", stats)] <- as.numeric(v0 <= cr)
  out[paste0("acceptinf_", stats)] <- as.numeric(vi <= cr)
  out["cover_IVW_RE"] <- as.numeric(abs(iv$th[1] - delta_true) <= qnorm(0.975) * s)
  if (isTRUE(do_proj)) {
    e <- eval_delta(p, delta_true, list(), do_oh = FALSE)
    if (!is.null(e$nu0)) out["cover_KOH_projection"] <- as.numeric(koh_proj_stat(p, delta_true, e$nu0, ct$KOH_projection) <= ct$KOH_projection)
  }
  out["F_D"] <- cond_F(p, 1L, 2:p$K)
  if (!is.null(cue)) out["S_min_over_df"] <- cue$S / (p$J - p$K)
  qh <- qhet_fit(p, if (!is.null(cue)) list(cue$th, iv$th) else list(iv$th))
  if (!is.null(qh)) out["tau2hat"] <- qh$tau2
  out["J"] <- p$J
  out["ok"] <- 1
  out
}

# Scenario table. gated=TRUE ones decide the calibration.
# rho=0.851 and 0.958 are the phenotypic correlations of the biochemistry and statin-naive arms;
# tau_mult=2.5 makes Q/df about 3.4, the value on UK Biobank CAD (506 on 148 df). The FinnGen case
# and control counts are approximate; the plasmode on the real standard errors does not depend on them.
scenarios <- function() {
  base <- list(J = 150L, K = 3L, nx = 400000, ncase = 47000, nctrl = 370000, corg = 0.97, rho = 0.851, zscale = 1,
               tau_mult = 0, het = "normal", select = FALSE, rho_analysis = NA_real_, gated = TRUE)
  th0 <- c(0.13, 0.13, 0.08); th1 <- c(0.593, -0.329, 0.079)
  mk <- function(id, ..., theta = th0) {
    s <- base; a <- list(...); for (n in names(a)) s[[n]] <- a[[n]]
    s$id <- id; s$theta <- theta; s$delta <- theta[1] - theta[2]
    s
  }
  list(
    mk("A1_cor97_rho851_null"),
    mk("A2_cor97_rho851_alt", theta = th1),
    mk("A3_cor97_rho958_null", rho = 0.958, nx = 330000),
    mk("A4_cor97_rho958_alt", rho = 0.958, nx = 330000, theta = th1),
    mk("B1_het_cor97_rho851_null", tau_mult = 2.5),
    mk("B2_het_cor97_rho851_alt", tau_mult = 2.5, theta = th1),
    mk("B3_het_cor97_rho958_null", tau_mult = 2.5, rho = 0.958, nx = 330000),
    mk("B4_het_cor97_rho958_alt", tau_mult = 2.5, rho = 0.958, nx = 330000, theta = th1),
    # No conditional strength at all: the difference of the true effects is exactly zero. Delta
    # is unidentified, so a valid set must still cover, and should usually be unbounded.
    mk("C1_cor100_rho958_null", corg = 1, rho = 0.958, nx = 330000),
    mk("C2_cor100_rho958_alt", corg = 1, rho = 0.958, nx = 330000, theta = th1),
    mk("C3_het_cor100_rho958_alt", corg = 1, rho = 0.958, nx = 330000, theta = th1, tau_mult = 2.5),
    # NMR arm: two exposures, 62 instruments, N = 212,797. zscale puts the conditional F near the
    # published 76, not the ~390 that the lipid z-score distribution would give at this J.
    mk("D1_nmr_J62_null", J = 62L, K = 2L, nx = 212797, rho = 0.955, zscale = 0.44, theta = c(0.3, 0.3)),
    mk("D2_nmr_J62_het_alt", J = 62L, K = 2L, nx = 212797, rho = 0.955, zscale = 0.44, theta = c(-0.37, 0.91), tau_mult = 2.5),
    # The A and B scenarios keep the 0.97 SNP-effect correlation but have conditional F near 146 and
    # 430, about three to four times the real arms (48 biochemistry, 104-117 statin-naive), and
    # coverage where the instruments are strong says little about where they are not. F and G
    # match the real strength, still at correlation 0.97, with the same exposure sample sizes.
    mk("F1_F48_rho851_null", zscale = 0.57),
    mk("F2_F48_rho851_alt", zscale = 0.57, theta = th1),
    mk("F3_het_F48_rho851_null", zscale = 0.57, tau_mult = 2.5),
    mk("F4_het_F48_rho851_alt", zscale = 0.57, tau_mult = 2.5, theta = th1),
    mk("G1_F105_rho958_null", zscale = 0.49, rho = 0.958, nx = 330000),
    mk("G2_het_F105_rho958_alt", zscale = 0.49, rho = 0.958, nx = 330000, tau_mult = 2.5, theta = th1),
    # Weaker than any real arm: conditional F near 10 and near 3, reached by raising the SNP-effect
    # correlation so that the nuisance (sum of LDL-C and ApoB, TG) stays as strong as it really is.
    # This is the regime the robust statistics exist for and where the profile could break.
    mk("H1_F10_rho851_null", corg = 0.998),
    mk("H2_het_F10_rho851_alt", corg = 0.998, tau_mult = 2.5, theta = th1),
    mk("H3_F3_rho851_null", corg = 0.9996),
    mk("H4_het_F3_rho851_alt", corg = 0.9996, tau_mult = 2.5, theta = th1),
    # Instruments chosen on the same data, as in the real analysis (winner's curse in the exposures).
    mk("E1_selected_rho851_het_alt", select = TRUE, tau_mult = 2.5, theta = th1),
    mk("E2_selected_rho958_null", select = TRUE, rho = 0.958, nx = 330000),
    # Diagnostics, not gated: departures that no statistic here claims to handle.
    mk("Z1_heavytail_het_alt", het = "t3", tau_mult = 2.5, theta = th1, gated = FALSE),
    mk("Z2_rho_true851_used871", rho_analysis = 0.871, gated = FALSE, theta = th1),
    mk("Z3_rho_true851_used831", rho_analysis = 0.831, gated = FALSE, theta = th1),
    mk("Z4_rho_true958_used968", rho = 0.958, nx = 330000, rho_analysis = 0.968, gated = FALSE, theta = th1),
    mk("Z5_rho_true958_used948", rho = 0.958, nx = 330000, rho_analysis = 0.948, gated = FALSE, theta = th1))
}

# Gate: heterogeneity-free scenarios gate all four profile statistics; under heterogeneity only KOH
# is claimed to hold, and the others are expected to fail (that is the point of KOH).
gated_stats <- function(sc) if (sc$tau_mult > 0) "KOH_subset" else c("AR_subset", "K_subset", "LC_subset", "KOH_subset")

rep_ok <- function(r) is.numeric(r) && "ok" %in% names(r) && isTRUE(r[["ok"]] == 1)

run_scenario <- function(sc, n_reps, ncores, seed_base, n_proj) {
  res <- parallel::mclapply(seq_len(n_reps), function(i) {
    s <- sc; s$rep_proj <- i <= n_proj
    tryCatch(sim_rep(s, seed_base + i), error = function(e) {
      m <- rep(NA_real_, 1); names(m) <- "ok"; m["ok"] <- 0; m })
  }, mc.cores = ncores, mc.preschedule = TRUE)
  # A rep that fell over and a rep that returned the short vector are both counted as failed.
  # A dead fork comes back as NULL or a try-error, so the check has to tolerate both.
  full <- res[vapply(res, rep_ok, logical(1))]
  list(mat = if (length(full)) do.call(rbind, full) else NULL, n_failed = n_reps - length(full))
}

summarise_scenario <- function(sc, run, n_reps) {
  m <- run$mat
  stats <- c("AR_subset", "AR_projection", "K_subset", "LC_subset", "KOH_subset", "KOH_projection", "WALD_GMM", "IVW_RE")
  rbindlist(lapply(stats, function(st) {
    cv <- m[, paste0("cover_", st)]; n <- sum(!is.na(cv))
    a0 <- paste0("accept0_", st); ai <- paste0("acceptinf_", st)
    gate <- isTRUE(sc$gated) && st %in% gated_stats(sc)
    cov <- if (n) mean(cv, na.rm = TRUE) else NA_real_
    # Replicates where a statistic failed to compute cannot count as covered or as not covered, so
    # a statistic that often fails has no readable coverage: failure may track the hard data sets.
    frac_na <- 1 - n / n_reps
    data.table(scenario = sc$id, gated = gate, statistic = st, n_rep = n, n_failed = run$n_failed, frac_na = frac_na,
               coverage = cov, mc_se = if (n) sqrt(cov * (1 - cov) / n) else NA_real_,
               rate_excludes_zero = if (a0 %in% colnames(m)) 1 - mean(m[, a0], na.rm = TRUE) else NA_real_,
               rate_unbounded = if (ai %in% colnames(m)) mean(m[, ai], na.rm = TRUE) else NA_real_,
               mean_F_contrast = mean(m[, "F_D"], na.rm = TRUE), mean_J = mean(m[, "J"], na.rm = TRUE),
               mean_Smin_over_df = mean(m[, "S_min_over_df"], na.rm = TRUE),
               pass = if (gate) !is.na(cov) && cov >= COV_LO && cov <= COV_HI && frac_na <= MAX_NA else NA)
  }))
}

script_path <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f[1]) else NA_character_
}
script_md5 <- function() { f <- script_path(); if (is.na(f)) "unknown" else unname(tools::md5sum(f)) }

run_calibration <- function(out_dir, n_reps, ncores, scens = scenarios(), n_proj = 200L, min_reps_for_pass = 1000L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  all <- list()
  for (k in seq_along(scens)) {
    sc <- scens[[k]]
    t0 <- Sys.time()
    run <- run_scenario(sc, n_reps, ncores, 42L + k * 100000L, n_proj)
    if (is.null(run$mat)) stop("calibration scenario produced nothing: ", sc$id)
    all[[k]] <- summarise_scenario(sc, run, n_reps)
    message(sprintf("[calib] %-30s %5.0fs failed=%d KOH=%.3f K=%.3f AR=%.3f LC=%.3f Wald=%.3f IVW=%.3f unb(KOH)=%.2f F=%.1f",
      sc$id, as.numeric(difftime(Sys.time(), t0, units = "secs")), run$n_failed,
      all[[k]][statistic == "KOH_subset", coverage], all[[k]][statistic == "K_subset", coverage],
      all[[k]][statistic == "AR_subset", coverage], all[[k]][statistic == "LC_subset", coverage],
      all[[k]][statistic == "WALD_GMM", coverage], all[[k]][statistic == "IVW_RE", coverage],
      all[[k]][statistic == "KOH_subset", rate_unbounded], all[[k]][statistic == "KOH_subset", mean_F_contrast]))
  }
  tab <- rbindlist(all)
  fwrite(tab, file.path(out_dir, "calibration_coverage.csv"))
  g <- tab[gated == TRUE]
  ok <- nrow(g) > 0 && !anyNA(g$pass) && all(g$pass) && all(g$n_failed <= MAX_NA * n_reps)
  passed <- ok && n_reps >= min_reps_for_pass
  writeLines(c(sprintf("PASSED=%s", if (passed) "TRUE" else "FALSE"),
               sprintf("script_md5=%s", script_md5()),
               sprintf("n_reps=%d", n_reps), sprintf("n_scenarios=%d", length(scens)), sprintf("date=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               sprintf("window=%.2f-%.2f", COV_LO, COV_HI)),
             file.path(out_dir, "calibration_status.txt"))
  if (nrow(g[pass == FALSE | is.na(pass)])) {
    message("[calib] gated failures:")
    print(g[pass == FALSE | is.na(pass), .(scenario, statistic, coverage, mc_se, n_failed, frac_na)])
  }
  message(sprintf("[calib] overall %s (n_reps=%d)", if (passed) "PASSED" else "NOT PASSED", n_reps))
  invisible(list(tab = tab, passed = passed))
}

# The gate read by the real run. Anything but an exact match on this script's md5 refuses: this
# project has been bitten repeatedly by results from stale artifacts.
read_gate <- function(path, md5 = script_md5()) {
  if (!file.exists(path)) return(list(ok = FALSE, why = sprintf("no calibration file at %s; run --calibrate first", path)))
  x <- readLines(path)
  kv <- setNames(sub("^[^=]*=", "", x), sub("=.*$", "", x))
  if (!identical(unname(kv["PASSED"]), "TRUE")) return(list(ok = FALSE, why = sprintf("calibration not passed (%s)", path)))
  if (!identical(unname(kv["script_md5"]), unname(md5))) return(list(ok = FALSE, why = "calibration was run with a different version of this script; rerun --calibrate"))
  if (as.integer(kv["n_reps"]) < 1000L) return(list(ok = FALSE, why = "calibration used fewer than 1000 replicates"))
  if (!identical(as.integer(kv["n_scenarios"]), length(scenarios()))) return(list(ok = FALSE, why = "calibration did not cover all scenarios"))
  list(ok = TRUE, why = "")
}

# Plasmode: the real standard errors and correlation, with the observed SNP effects taken as truth
# and an outcome simulated from them. Reads the real outcome only to set tau2 and the nuisance.
plasmode <- function(bx, bxse, bys, pcor, p_real, th_qhet, tau2, delta_alt, n_reps, ncores, window = c(COV_LO, COV_HI)) {
  L <- chol(pcor); J <- nrow(bx); K <- ncol(bx)
  # Null truth keeps the qhet nuisance values and sets delta to zero; alternative uses delta_alt.
  run <- function(delta, seed_base) {
    th <- th_qhet; th[1] <- delta
    orig <- p_real$Tm %*% th
    Yt <- as.numeric(bx %*% orig)
    res <- parallel::mclapply(seq_len(n_reps), function(i) {
      set.seed(seed_base + i)
      bxs <- bx + (matrix(rnorm(J * K), J, K) %*% L) * bxse
      bys2 <- Yt + rnorm(J, 0, sqrt(tau2)) + rnorm(J, 0, bys)
      tryCatch(test_prob(mk_prob(bxs, bxse, bys2, bys, pcor), delta), error = function(e) c(ok = 0))
    }, mc.cores = ncores)
    res <- res[vapply(res, rep_ok, logical(1))]
    if (length(res) < (1 - MAX_NA) * n_reps) return(list(cov = NA_real_, row = NULL))
    do.call(rbind, res)
  }
  tabs <- list(); pass <- TRUE
  for (w in list(list("null", 0, 910000L), list("alt", delta_alt, 920000L))) {
    m <- run(w[[2]], w[[3]])
    if (!is.matrix(m)) { pass <- FALSE; tabs[[w[[1]]]] <- data.table(truth = w[[1]], statistic = "KOH_subset", coverage = NA_real_, n_rep = 0L); next }
    cv <- mean(m[, "cover_KOH_subset"], na.rm = TRUE)
    # Same rule as the calibration gate: failures of the statistic may not be averaged away.
    fna <- mean(is.na(m[, "cover_KOH_subset"]))
    if (is.na(cv) || cv < window[1] || cv > window[2] || fna > MAX_NA) pass <- FALSE
    tabs[[w[[1]]]] <- data.table(truth = w[[1]], delta_true = w[[2]],
      statistic = c("KOH_subset", "K_subset", "AR_subset", "LC_subset", "WALD_GMM", "IVW_RE"),
      coverage = c(cv, mean(m[, "cover_K_subset"], na.rm = TRUE), mean(m[, "cover_AR_subset"], na.rm = TRUE),
                   mean(m[, "cover_LC_subset"], na.rm = TRUE), mean(m[, "cover_WALD_GMM"], na.rm = TRUE),
                   mean(m[, "cover_IVW_RE"], na.rm = TRUE)),
      rate_unbounded_KOH = mean(m[, "acceptinf_KOH_subset"], na.rm = TRUE),
      rate_excludes_zero_KOH = 1 - mean(m[, "accept0_KOH_subset"], na.rm = TRUE),
      frac_na_KOH = fna, n_rep = nrow(m))
  }
  list(pass = pass, tab = rbindlist(tabs, fill = TRUE))
}

# ------------------------------------------------------------
# Verdict under the readout rules
# ------------------------------------------------------------
PRIMARY_ARM <- "biochem"; PRIMARY_OUT <- "FINNGEN"
# GRCh37 regions, chromosome / start / end. 157-163 Mb is the cis region of run_mr.R.
LOCI_OUT <- list(LPA = list(c(6, 157e6, 163e6)), APOE = list(c(19, 44e6, 47e6)),
                 LPA_APOE = list(c(6, 157e6, 163e6), c(19, 44e6, 47e6)))
locus_mask <- function(chrom, pos, regs) Reduce(`|`, lapply(regs, function(r) as.character(chrom) == as.character(r[1]) & pos >= r[2] & pos <= r[3]))

readout_verdict <- function(sets, ests, invalid_arms = character(0)) {
  key <- sets[statistic == "KOH_subset"]
  pr <- key[arm == PRIMARY_ARM & outcome == PRIMARY_OUT]
  lines <- character(0)
  if (nrow(pr) != 1L) return(list(label = "INVALID", lines = "No primary KOH row: nothing may be read."))
  if (pr$status != "OK") return(list(label = "INVALID", lines = sprintf("Primary set status %s: nothing may be read.", pr$status)))
  label <- NA_character_
  if (pr$set_class == "empty") {
    label <- "INVALID"
    lines <- c(lines, "Primary KOH set is empty, which KOH should not be able to produce. Investigate; do not read as precision.")
  } else if (pr$set_class %in% c("unbounded", "whole_line")) {
    label <- "SUPPORTS_NOT_IDENTIFIABLE"
    lines <- c(lines, sprintf("Primary robust set is %s: %s. Readout: SUPPORTS 'not identifiable'.", pr$set_class, pr$intervals))
  } else if (isTRUE(pr$includes_zero)) {
    label <- "SUPPORTS_NOT_IDENTIFIABLE"
    lines <- c(lines, sprintf("Primary robust set %s includes 0. Readout: SUPPORTS 'not identifiable'.", pr$intervals))
    if (identical(pr$zero_component_has_S_minimum, FALSE))
      lines <- c(lines, "CAUTION: the component containing 0 has no minimum of S. K is also zero at maxima of S, so 0 may be in the set through that artefact. The rule is applied as written; do not lean on this inclusion alone.")
  } else {
    sg <- pr$set_sign
    if (!sg %in% c("+", "-")) {
      label <- "UNRESOLVED"
      lines <- c(lines, sprintf("Primary set %s excludes 0 but lies on both sides of it, so it has no sign. The readout has no rule for this.", pr$intervals))
    } else {
      s <- if (sg == "+") 1 else -1
      others <- ests[!(arm == PRIMARY_ARM & outcome == PRIMARY_OUT)]
      okk <- others[status == "OK"]
      conf_point <- okk[sign(qhet_contrast) != s]
      oppo <- key[!(arm == PRIMARY_ARM & outcome == PRIMARY_OUT) & status == "OK" & !is.na(set_sign) & set_sign != sg]
      # The primary arm is present by construction here; only the secondary arms can be missing.
      missing_arms <- setdiff(c("statinnaive", "nmr"), unique(okk$arm))
      if (nrow(conf_point) || nrow(oppo)) {
        label <- "UNRESOLVED"
        lines <- c(lines, sprintf("Primary set %s excludes 0 (sign %s), but sign disagrees in: %s. CONTRADICTS requires agreement across arms; the readout has no rule for this case.",
                                  pr$intervals, sg, paste(unique(c(paste(conf_point$arm, conf_point$outcome), paste(oppo$arm, oppo$outcome))), collapse = "; ")))
      } else if (length(missing_arms) || length(invalid_arms)) {
        label <- "UNRESOLVED"
        lines <- c(lines, sprintf("Primary set %s excludes 0 (sign %s) and available arms agree, but arms missing or invalid: %s. Agreement across arms is not established.",
                                  pr$intervals, sg, paste(unique(c(missing_arms, invalid_arms)), collapse = ", ")))
      } else {
        label <- "CONTRADICTS_NOT_IDENTIFIABLE"
        lines <- c(lines, sprintf("Primary set %s is bounded, excludes 0 and sign (%s) agrees across arms. Readout: CONTRADICTS; the 'not identifiable' claim must be rewritten.", pr$intervals, sg))
      }
    }
  }
  # Reading aids the readout requires, on every row they apply to.
  ar <- sets[statistic == "AR_subset" & status == "OK" & set_class == "empty"]
  for (i in seq_len(nrow(ar)))
    lines <- c(lines, sprintf("AR set empty in %s / %s: evidence of heterogeneity or misspecification, never of a precise estimate.", ar$arm[i], ar$outcome[i]))
  if (length(invalid_arms))
    lines <- c(lines, sprintf("Arms invalidated by a failed control, their numbers are NA: %s.", paste(invalid_arms, collapse = ", ")))
  list(label = label, lines = lines)
}

# The label rests on one set, chosen before any data were seen. These lines put the rest of the
# primary analysis beside it, so a reader can see when the choice of statistic is doing the work.
primary_context <- function(sets, ests) {
  pr <- sets[arm == PRIMARY_ARM & outcome == PRIMARY_OUT]
  e <- ests[arm == PRIMARY_ARM & outcome == PRIMARY_OUT]
  if (!nrow(pr) || !"intervals" %in% names(pr)) return(character(0))
  ln <- sprintf("Primary, %s (%s): %s; includes 0: %s", pr$statistic, pr$role, pr$intervals, pr$includes_zero)
  key <- pr[statistic == "KOH_subset"]
  dis <- pr[statistic != "KOH_subset" & role != "NON_ROBUST" & status == "OK" & !is.na(includes_zero) & includes_zero != key$includes_zero]
  if (nrow(key) == 1L && nrow(dis))
    ln <- c(ln, sprintf("DISAGREEMENT on whether 0 is in the set: KOH_subset says %s; %s say otherwise. AR, K and Andrews lose coverage under heterogeneity (S_min/df %.2f here) and KOH_projection is conservative, so read the reason before the conclusion.",
                        key$includes_zero, paste(dis$statistic, collapse = ", "), if (nrow(e) && "S_min" %in% names(e)) e$S_min / (e$J - e$K) else NA_real_))
  if (nrow(e) && "F_nuisance_min" %in% names(e)) {
    ln <- c(ln, sprintf("Primary, conditional F of the contrast given the nuisance %.1f; minimum F of the nuisance %.1f; Andrews distortion cutoff %.3f (interpretable: %s).",
                        e$F_contrast_given_nuisance, e$F_nuisance_min, e$distortion_cutoff, e$distortion_cutoff_interpretable))
    # The profile statistics are justified only when the nuisance is strong; projection is not.
    if (!is.na(e$F_nuisance_min) && e$F_nuisance_min < 10)
      ln <- c(ln, "WARNING: the nuisance is not strongly identified (F < 10), which the profile statistics assume. Read the projection rows.")
  }
  if (nrow(e) && "qhet_contrast" %in% names(e) && "qhet_boot_perc_lo" %in% names(e))
    ln <- c(ln, sprintf("Primary, qhet contrast %.3f, tau2 %.3g, bootstrap percentile CI %.3f to %.3f, jackknife SE %.3f.",
                        e$qhet_contrast, e$qhet_tau2, e$qhet_boot_perc_lo, e$qhet_boot_perc_hi, e$qhet_jack_se))
  ln
}

# ------------------------------------------------------------
# Reading summary statistics without loading whole GWAS files
# ------------------------------------------------------------
# bash, not sh, so pipefail catches a truncated gzip: a silently short extract would otherwise
# look like a clean file with fewer instruments.
run_bash <- function(cmd) {
  st <- system2("bash", c("-c", shQuote(paste("set -o pipefail;", cmd))))
  if (st != 0) stop("shell step failed (status ", st, "): ", substr(cmd, 1, 300))
}
reader <- function(f) if (grepl("\\.gz$", f)) paste("gunzip -c", shQuote(f)) else paste("cat", shQuote(f))

# Cache files carry a hash of what produced them: the input paths with size and mtime, and the rule
# or ID list. A name built from the arm alone would serve the earlier extract when --sumstat_dir
# changes (the censored Lp(a) route) or when the instrument list moves, and stale artifacts are
# this project's recurring failure.
cache_file <- function(base, files, ...) {
  tf <- tempfile(); on.exit(unlink(tf))
  writeLines(c(normalizePath(files), as.character(file.size(files)), format(file.mtime(files), "%s"), as.character(unlist(list(...)))), tf)
  sprintf("%s.%s.tsv", base, substr(unname(tools::md5sum(tf)), 1, 12))
}

# Rows passing LOG10P, from one or many files (the NMR GWAS is per chromosome). The threshold is
# on LOG10P, never on a p-value, which underflows to zero above LOG10P ~ 323.
extract_gws <- function(files, base, thr = -log10(5e-8)) {
  out <- cache_file(base, files, "gws", sprintf("%.10f", thr))
  if (file.exists(out) && file.size(out) > 0) return(fread(out, na.strings = "NA"))
  tmp <- paste0(out, ".tmp"); unlink(tmp)
  for (k in seq_along(files)) {
    run_bash(sprintf("%s | awk -v h=%d -v t=%.10f 'FNR==1{for(i=1;i<=NF;i++) if($i==\"LOG10P\") c=i; if(!c){print \"no LOG10P column\" > \"/dev/stderr\"; exit 3} if(h)print; next} ($c+0)>t' >> %s",
                     reader(files[k]), as.integer(k == 1L), thr, shQuote(tmp)))
  }
  file.rename(tmp, out)
  fread(out, na.strings = "NA")
}

# Rows whose ID is in a list. For FinnGen the ID is the FIRST rsID in a comma list, which is what
# run_external_mr.R uses, so that the same rows are matched as in the published replication.
extract_ids <- function(files, ids, base, idcol = "ID", sep = NULL, first_of_list = FALSE) {
  out <- cache_file(base, files, "ids", sort(ids), idcol, first_of_list)
  if (file.exists(out) && file.size(out) > 0) return(fread(out, na.strings = c("NA", ".")))
  idf <- paste0(out, ".ids"); writeLines(ids, idf)
  tmp <- paste0(out, ".tmp"); unlink(tmp)
  fs <- if (is.null(sep)) "" else sprintf("-F'%s'", sep)
  val <- if (first_of_list) "split($c,r,\",\"); v=r[1]" else "v=$c"
  for (k in seq_along(files)) {
    run_bash(sprintf("%s | awk %s -v h=%d 'NR==FNR{a[$1];next} FNR==1{for(i=1;i<=NF;i++) if($i==\"%s\") c=i; if(!c){print \"no %s column\" > \"/dev/stderr\"; exit 3} if(h)print; next} {%s; if(v in a) print}' %s - >> %s",
                     reader(files[k]), fs, as.integer(k == 1L), idcol, idcol, val, shQuote(idf), shQuote(tmp)))
  }
  file.rename(tmp, out)
  fread(out, na.strings = c("NA", "."))
}

# Exact copy of the greedy 10 Mb distance clump in run_mvmr_apob.R, including its two bug fixes:
# rank on LOG10P, and subset to the chromosome before comparing positions.
greedy_clump <- function(dt, kb_window = 10000L) {
  gws <- dt[!is.na(BETA) & !is.na(SE)]
  if (nrow(gws) == 0) return(character(0))
  gws <- gws[order(-LOG10P)]
  keep_snps <- character(0); keep_pos <- integer(0); keep_chr <- character(0)
  for (i in seq_len(nrow(gws))) {
    chr <- as.character(gws$CHROM[i]); pos <- as.integer(gws$GENPOS[i])
    same_chr <- keep_chr == chr
    if (length(keep_pos) == 0 || !any(abs(keep_pos[same_chr] - pos) < kb_window * 1000L)) {
      keep_snps <- c(keep_snps, gws$ID[i]); keep_pos <- c(keep_pos, pos); keep_chr <- c(keep_chr, chr)
    }
  }
  keep_snps
}

select_instruments <- function(gws, kb = 10000L) {
  trait_ivs <- lapply(gws, greedy_clump, kb_window = kb)
  union_ivs <- unique(unlist(trait_ivs))
  ref <- rbindlist(lapply(gws, function(d) d[ID %in% union_ivs, .(ID, CHROM, GENPOS, LOG10P)]))
  ref <- ref[ref[, .I[which.max(LOG10P)], by = ID]$V1][order(-LOG10P)]
  ids <- character(0); pos <- integer(0); chrs <- character(0)
  for (i in seq_len(nrow(ref))) {
    chr <- as.character(ref$CHROM[i]); po <- as.integer(ref$GENPOS[i]); sm <- chrs == chr
    if (length(pos) == 0 || !any(abs(pos[sm] - po) < kb * 1000L)) {
      ids <- c(ids, ref$ID[i]); pos <- c(pos, po); chrs <- c(chrs, chr)
    }
  }
  ids
}

# Same alignment as run_mvmr_apob.R: first trait is the allele reference, complete cases only.
align_exposures <- function(rows, ids) {
  tr <- names(rows)
  ref <- rows[[1]][ID %in% ids, .(ID, ALLELE1, ALLELE0, A1FREQ, CHROM, GENPOS)]
  if (anyDuplicated(ref$ID)) message(sprintf("[arm] note: %d rsIDs on more than one exposure row; last row wins, as in run_mvmr_apob.R", sum(duplicated(ref$ID))))
  bx <- se <- matrix(NA_real_, length(ids), length(tr), dimnames = list(ids, tr))
  for (t in tr) {
    d <- merge(rows[[t]][ID %in% ids], ref[, .(ID, ref_A1 = ALLELE1, ref_A0 = ALLELE0)], by = "ID", all.x = TRUE, allow.cartesian = TRUE)
    d[, beta_aligned := fifelse(ALLELE1 == ref_A1, BETA, fifelse(ALLELE1 == ref_A0, -BETA, NA_real_))]
    i <- match(d$ID, ids)
    bx[i, t] <- d$beta_aligned; se[i, t] <- d$SE
  }
  ok <- complete.cases(bx) & complete.cases(se)
  list(bx = bx[ok, , drop = FALSE], se = se[ok, , drop = FALSE], ref = ref[ID %in% rownames(bx)[ok]])
}

outcome_regenie <- function(o, ref, ids) {
  o <- merge(o, ref[, .(ID, ref_A1 = ALLELE1, ref_A0 = ALLELE0)], by = "ID", all.x = TRUE, allow.cartesian = TRUE)
  o[, beta_aligned := fifelse(ALLELE1 == ref_A1, BETA, fifelse(ALLELE1 == ref_A0, -BETA, NA_real_))]
  o <- o[match(ids, ID)]
  list(b = o$beta_aligned, s = o$SE)
}

# Harmonisation from run_external_mr.R: exact and strand flips, palindromes resolved by allele
# frequency and dropped when ambiguous (MAF 0.42-0.58), one row per rsID. palin="drop" removes all
# palindromes, which is what nmr_mvmr.R did for the published 62-instrument NMR analysis.
outcome_external <- function(o, ref, ids, FM, palin = "freq", first_of_list = FALSE) {
  setnames(o, c(FM$snp, FM$ea, FM$oa, FM$beta, FM$se, FM$eaf), c("SNP", "oEA", "oOA", "oBETA", "oSE", "oEAF"), skip_absent = FALSE)
  if (first_of_list) o[, SNP := tstrsplit(SNP, ",", fixed = TRUE)[[1]]]
  o <- o[!is.na(SNP) & SNP != "" & !is.na(oBETA) & !is.na(oSE)]
  m <- merge(ref[, .(SNP = ID, eEA = toupper(ALLELE1), eOA = toupper(ALLELE0), eEAF = A1FREQ)],
             o[, .(SNP, oEA = toupper(oEA), oOA = toupper(oOA), oBETA, oSE, oEAF)], by = "SNP", allow.cartesian = TRUE)
  flip_base <- function(x) chartr("ACGT", "TGCA", x)
  m[, palin_ := paste0(eEA, eOA) %chin% c("AT", "TA", "CG", "GC")]
  m[, aligned := fifelse(eEA == oEA & eOA == oOA, 1,
               fifelse(eEA == oOA & eOA == oEA, -1,
               fifelse(eEA == flip_base(oEA) & eOA == flip_base(oOA), 1,
               fifelse(eEA == flip_base(oOA) & eOA == flip_base(oEA), -1, NA_real_))))]
  # Ambiguous in EITHER source is dropped: UKB 0.40 with FinnGen 0.53 would otherwise be flipped
  # although its alleles match. A row that failed the allele match stays dropped. run_external_mr.R
  # has neither guard, so the count below is how many instruments its rule treats differently.
  if (palin == "drop") m[palin_ == TRUE, aligned := NA_real_] else {
    amb_o <- m[palin_ == TRUE & !is.na(aligned) & !is.na(oEAF) & oEAF > 0.42 & oEAF < 0.58 & !is.na(eEAF) & !(eEAF > 0.42 & eEAF < 0.58), uniqueN(SNP)]
    bad_o <- m[palin_ == TRUE & is.na(aligned) & !is.na(eEAF) & !is.na(oEAF) & !(eEAF > 0.42 & eEAF < 0.58), uniqueN(SNP)]
    if (amb_o + bad_o > 0L) message(sprintf("[harmonise] palindromes treated differently from run_external_mr.R: %d ambiguous only in the outcome, %d failing the allele match", amb_o, bad_o))
    m[palin_ == TRUE & !is.na(aligned), aligned := fifelse(is.na(eEAF) | is.na(oEAF) | (eEAF > 0.42 & eEAF < 0.58) | (oEAF > 0.42 & oEAF < 0.58), NA_real_,
                                                           fifelse((eEAF - 0.5) * (oEAF - 0.5) > 0, 1, -1))]
  }
  m <- m[!is.na(aligned)]
  m <- m[!duplicated(SNP)]
  i <- match(ids, m$SNP)
  list(b = m$oBETA[i] * m$aligned[i], s = m$oSE[i])
}

FORMATS <- list(
  cardiogram = list(snp = "markername", ea = "effect_allele", oa = "noneffect_allele", beta = "beta", se = "se_dgc", eaf = "effect_allele_freq"),
  finngen = list(snp = "rsids", ea = "alt", oa = "ref", beta = "beta", se = "sebeta", eaf = "af_alt"))

phen_cor <- function(pheno, traits, keep = NULL) {
  ph <- fread(pheno, select = c("FID", "IID", traits), na.strings = c("NA", ""))
  if (!is.null(keep)) {
    k <- fread(keep, header = FALSE, col.names = c("FID", "IID"))
    ph <- ph[as.character(IID) %chin% as.character(k$IID)]
  }
  m <- as.matrix(ph[, ..traits])
  # Spearman: the GWAS rank-normalised the traits, and Spearman is invariant to that transform.
  # A CORRELATION matrix, never a covariance: see the bug comment in R/mvmr_covariance.R.
  R <- cor(m, use = "pairwise.complete.obs", method = "spearman")
  n1 <- colSums(!is.na(m))
  n12 <- crossprod(!is.na(m))
  # What the GWAS actually correlated: Pearson on rank-normal scores. It is not quite Spearman (for
  # bivariate normal data at r = 0.85, Spearman is 0.84), so it is reported beside it and the rho
  # sensitivity rows cover the gap. Spearman stays primary to reproduce the published phenocov values.
  rn <- apply(m, 2, function(x) { o <- rank(x, na.last = "keep"); qnorm((o - 0.5) / sum(!is.na(x))) })
  list(R = R, rho_rint = cor(rn[, 1], rn[, 2], use = "pairwise.complete.obs"),
       overlap_factor = min((n12 / sqrt(outer(n1, n1)))[upper.tri(n12)]), n = nrow(m))
}

# ------------------------------------------------------------
# Real run
# ------------------------------------------------------------
# An arm whose instruments and harmonised effects are fixed in a file rather than selected here.
# The NMR arm uses this. Supplementary Table oS3b came from exports/agent1/nmr_gencov/nmr_mvmr_snps.csv
# (62 instruments; columns SNP bxp sxp bxc sxc by sy; outcome CARDIoGRAM), and reselecting from the
# archived per-chromosome GWAS with this script's distance-only clump does not reproduce it (job 3549234:
# LDL-P 2.94 and LDL-C(NMR) 0.536 against the published 2.48 and 0.69). The file carries no alleles, so
# no other outcome can be harmonised to it and such an arm has its one outcome only.
load_precomputed <- function(pc, traits) {
  x <- fread(pc$file)
  need <- c("SNP", unname(pc$cols$bx[traits]), unname(pc$cols$se[traits]), pc$cols$by, pc$cols$sy)
  gone <- setdiff(need, names(x))
  if (length(gone)) stop("precomputed instrument file lacks columns: ", paste(gone, collapse = ", "))
  if (anyDuplicated(x$SNP)) stop("precomputed instrument file has duplicated SNPs: ", pc$file)
  x <- x[complete.cases(x[, need, with = FALSE])]
  bx <- as.matrix(x[, unname(pc$cols$bx[traits]), with = FALSE]); se <- as.matrix(x[, unname(pc$cols$se[traits]), with = FALSE])
  dimnames(bx) <- dimnames(se) <- list(x$SNP, traits)
  list(O = pc$outcome, bx = bx, se = se, by = x[[pc$cols$by]], bys = x[[pc$cols$sy]], snps = x$SNP)
}

arm_config <- function(o) {
  R <- file.path(o$repo, "results/agent1_genetics")
  sumstat <- function(dir) function(t) file.path(dir, sprintf("lipids_%s.txt.gz", t))
  ukb <- function(dir) list(id = "UKB_CAD", kind = "regenie", files = file.path(dir, "outcomes_CAD.txt.gz"),
                            overlap = "UKB_ON_UKB_UNCORRECTED")
  fg <- list(id = "FINNGEN", kind = "finngen", files = o$finngen, overlap = "none", palin = "freq")
  cg <- function(p) list(id = "CARDIOGRAM", kind = "cardiogram", files = o$cardiogram, overlap = "none", palin = p)
  list(
    biochem = list(name = "biochem", traits = c("LDL_C", "APOB", "TRIGLY"), files = sumstat(o$sumstat_dir),
                   required = TRUE, pheno = o$pheno, keep = NULL, outcomes = list(fg, ukb(o$sumstat_dir)),
                   # From exports/agent1/mvmr_apob_prunefix/mvmr_apob_CAD.csv.
                   control = list(outcome = "UKB_CAD", n = 151L, b = c(LDL_C = 0.59246, APOB = -0.32881),
                                  rho = 0.851, condF = c(LDL_C = 48.2, APOB = 49.0))),
    statinnaive = list(name = "statinnaive", traits = c("LDL_C", "APOB", "TRIGLY"), files = sumstat(o$statinnaive_dir),
                       required = FALSE, pheno = o$pheno_statinnaive, keep = o$keep_statinnaive,
                       outcomes = list(fg, ukb(o$statinnaive_dir)),
                       control = list(outcome = "UKB_CAD", n = 148L, b = NULL, rho = 0.958, condF = c(LDL_C = 103.9, APOB = 117.4))),
    # Contrast oriented LDL-C minus particle number, to match LDL-C minus ApoB.
    nmr = list(name = "nmr", traits = c("LDL_C_NMR", "LDL_P"),
               # The reselection route below is kept for the selftest; precomputed takes precedence in a run
               # (PREREGISTERED_READOUT.md, amendment of 2026-11-05).
               files = function(t) file.path(o$nmr_dir, sprintf("nmr_chr%d_%s.regenie.gz", 1:22, t)),
               required = FALSE, pheno = o$pheno_nmr, keep = NULL,
               outcomes = list(cg("drop"), fg),
               precomputed = list(file = o$nmr_instruments,
                                  outcome = list(id = "CARDIOGRAM", kind = "precomputed", overlap = "none"),
                                  cols = list(bx = c(LDL_C_NMR = "bxc", LDL_P = "bxp"), se = c(LDL_C_NMR = "sxc", LDL_P = "sxp"),
                                              by = "by", sy = "sy")),
               # Published on CARDIoGRAM (Table oS3b): LDL_P OR 2.48, LDL_C_NMR OR 0.69, 62 instruments. Full
               # precision from exports/agent1/nmr_gencov/nmr_mvmr_results.csv, the MVMR-IVW row.
               control = list(outcome = "CARDIOGRAM", n = 62L, b = c(LDL_C_NMR = log(0.686484624260781), LDL_P = log(2.47578092791507)),
                              rho = 0.955, condF = c(LDL_P = 75.6, LDL_C_NMR = 72.1)))
  )
}

# The positive control. Returns problems as strings; the caller decides whether to stop or to
# invalidate. Tolerance on b is about half a rounding unit of the published OR.
check_control <- function(ctl, n, fit_b, rho, condF, b_tol = 0.01) {
  hard <- character(0); soft <- character(0)
  if (n != ctl$n) hard <- c(hard, sprintf("instrument count %d, expected %d", n, ctl$n))
  if (!is.null(ctl$b)) for (nm in names(ctl$b)) {
    # A missing estimate fails the control; it must not surface as an R error in the if().
    if (!isTRUE(abs(fit_b[nm] - ctl$b[nm]) <= b_tol))
      hard <- c(hard, sprintf("%s IVW log OR %.5f, expected %.5f", nm, fit_b[nm], ctl$b[nm]))
  }
  # These two come from the project's notes rather than an export, so they warn instead of stopping.
  if (!is.na(ctl$rho) && abs(rho - ctl$rho) > 0.01) soft <- c(soft, sprintf("rho %.3f, notes say %.3f", rho, ctl$rho))
  for (nm in names(ctl$condF)) if (abs(condF[nm] / ctl$condF[nm] - 1) > 0.10)
    soft <- c(soft, sprintf("%s conditional F %.1f, notes say %.1f", nm, condF[nm], ctl$condF[nm]))
  list(hard = hard, soft = soft)
}

run_real <- function(o, gate_path, plasmode_window = c(COV_LO, COV_HI), arms_override = NULL, md5 = script_md5()) {
  # Before any data is read: the readout says calibration passes first.
  gate <- read_gate(gate_path, md5)
  if (!gate$ok) stop("[gate] refusing to run: ", gate$why)
  for (pk in c("MendelianRandomization", "MVMR", "boot"))
    if (!requireNamespace(pk, quietly = TRUE)) stop("required package missing: ", pk)
  dir.create(o$out_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- if (nzchar(o$cache_dir)) o$cache_dir else file.path(o$out_dir, "cache")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  arms <- if (is.null(arms_override)) arm_config(o)[strsplit(o$arms, ",")[[1]]] else arms_override
  if (is.null(arms$biochem)) stop("the biochemistry arm is the primary analysis and cannot be left out")

  # After the gate, before any data: a rerun that stops early must not leave an older verdict behind.
  unlink(file.path(o$out_dir, c("robust_contrast_verdict.txt", "robust_contrast_sets.csv", "robust_contrast_estimates.csv",
                                "plasmode_calibration.csv", "rho_sensitivity.csv", "robust_contrast_leave_locus.csv")))
  all_sets <- list(); all_ests <- list(); all_pl <- list(); all_rho <- list(); all_locus <- list(); invalid <- character(0); notes <- character(0)
  flush <- function() {
    if (length(all_sets)) fwrite(rbindlist(all_sets, fill = TRUE), file.path(o$out_dir, "robust_contrast_sets.csv"))
    if (length(all_ests)) fwrite(rbindlist(all_ests, fill = TRUE), file.path(o$out_dir, "robust_contrast_estimates.csv"))
    if (length(all_pl)) fwrite(rbindlist(all_pl, fill = TRUE), file.path(o$out_dir, "plasmode_calibration.csv"))
    if (length(all_rho)) fwrite(rbindlist(all_rho, fill = TRUE), file.path(o$out_dir, "rho_sensitivity.csv"))
    if (length(all_locus)) fwrite(rbindlist(all_locus, fill = TRUE), file.path(o$out_dir, "robust_contrast_leave_locus.csv"))
  }

  for (A in arms) {
    message(sprintf("\n===== arm %s =====", A$name))
    pc <- A$precomputed
    trfiles <- if (is.null(pc)) lapply(setNames(A$traits, A$traits), A$files) else list()
    miss <- unique(c(unlist(lapply(trfiles, function(f) f[!file.exists(f)])),
                     if (!is.null(pc) && !isTRUE(file.exists(pc$file))) (if (is.null(pc$file)) "precomputed instrument file not set" else pc$file),
                     if (!file.exists(A$pheno)) A$pheno,
                     if (!is.null(A$keep) && !file.exists(A$keep)) A$keep))
    if (length(miss)) {
      if (A$required) stop(sprintf("primary arm is missing files: %s", paste(head(miss, 5), collapse = ", ")))
      message(sprintf("[%s] SKIPPED, files absent (%d), e.g. %s", A$name, length(miss), miss[1]))
      notes <- c(notes, sprintf("Arm %s skipped: files absent, e.g. %s", A$name, miss[1]))
      invalid <- c(invalid, A$name)
      next
    }
    if (is.null(pc)) {
      gws <- lapply(A$traits, function(t) extract_gws(trfiles[[t]], file.path(cache, sprintf("%s_%s_gws", A$name, t))))
      names(gws) <- A$traits
      ids <- select_instruments(gws)
      message(sprintf("[%s] %d instruments after second-pass clump", A$name, length(ids)))
      rows <- lapply(A$traits, function(t) extract_ids(trfiles[[t]], ids, file.path(cache, sprintf("%s_%s_ivrows", A$name, t))))
      names(rows) <- A$traits
      al <- align_exposures(rows, ids)
      ph <- phen_cor(A$pheno, A$traits, A$keep)
      message(sprintf("[%s] phenotypic Spearman correlation of the pair %.3f; missingness factor %.3f", A$name, ph$R[1, 2], ph$overlap_factor))
      if (min(eigen(ph$R, only.values = TRUE)$values) < 1e-6) stop("phenotypic correlation is singular for ", A$name)

      outs <- list()
      for (O in A$outcomes) {
        if (!file.exists(O$files)) {
          if (A$required && O$id == PRIMARY_OUT) stop("primary outcome file missing: ", O$files)
          message(sprintf("[%s] outcome %s absent, skipped: %s", A$name, O$id, O$files)); next
        }
        tag <- file.path(cache, sprintf("%s_%s_ivrows", A$name, O$id))
        ivs <- rownames(al$bx)
        if (O$kind == "regenie") {
          raw <- extract_ids(O$files, ivs, tag)
          oo <- outcome_regenie(raw, al$ref, ivs)
        } else if (O$kind == "finngen") {
          raw <- extract_ids(O$files, ivs, tag, idcol = "rsids", sep = "\\t", first_of_list = TRUE)
          oo <- outcome_external(raw, al$ref, ivs, FORMATS$finngen, O$palin, TRUE)
        } else {
          raw <- extract_ids(O$files, ivs, tag, idcol = "markername")
          oo <- outcome_external(raw, al$ref, ivs, FORMATS$cardiogram, O$palin)
        }
        keep <- !is.na(oo$b) & !is.na(oo$s)
        outs[[O$id]] <- list(O = O, bx = al$bx[keep, , drop = FALSE], se = al$se[keep, , drop = FALSE],
                             by = oo$b[keep], bys = oo$s[keep], snps = ivs[keep])
        message(sprintf("[%s] %-10s %d of %d instruments harmonised", A$name, O$id, sum(keep), length(ivs)))
      }
    } else {
      ph <- phen_cor(A$pheno, A$traits, A$keep)
      message(sprintf("[%s] phenotypic Spearman correlation of the pair %.3f; missingness factor %.3f", A$name, ph$R[1, 2], ph$overlap_factor))
      if (min(eigen(ph$R, only.values = TRUE)$values) < 1e-6) stop("phenotypic correlation is singular for ", A$name)
      al <- NULL  # no selected instruments here, and a previous arm's must not leak into this one
      cur <- load_precomputed(pc, A$traits)
      outs <- setNames(list(cur), pc$outcome$id)
      message(sprintf("[%s] %d instruments from the published file %s; outcome %s only, as published", A$name, nrow(cur$bx), basename(pc$file), pc$outcome$id))
    }

    # Positive control first, so that a wrong arm costs minutes and not the whole job.
    ctl <- A$control; cd <- outs[[ctl$outcome]]
    arm_ok <- TRUE
    if (is.null(cd)) {
      probs <- list(hard = sprintf("control outcome %s not available", ctl$outcome), soft = character(0))
    } else {
      fit0 <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = cd$bx, bxse = cd$se, by = cd$by, byse = cd$bys))
      fb <- setNames(fit0@Estimate, A$traits)
      fmt <- MVMR::format_mvmr(BXGs = cd$bx, BYG = cd$by, seBXGs = cd$se, seBYG = cd$bys, RSID = rownames(cd$bx))
      sf <- suppressWarnings(MVMR::strength_mvmr(fmt, gencov = MVMR::phenocov_mvmr(ph$R, cd$se)))
      cf <- setNames(as.numeric(sf[1, seq_along(A$traits)]), A$traits)
      probs <- check_control(ctl, nrow(cd$bx), fb, ph$R[1, 2], cf)
      message(sprintf("[%s] control: n=%d, IVW log OR %s, conditional F %s", A$name, nrow(cd$bx),
                      paste(sprintf("%s %.4f", names(fb), fb), collapse = ", "),
                      paste(sprintf("%s %.1f", names(cf), cf), collapse = ", ")))
    }
    for (s in probs$soft) { message(sprintf("[%s] CONTROL WARNING: %s", A$name, s)); notes <- c(notes, sprintf("Arm %s control warning: %s", A$name, s)) }
    if (length(probs$hard)) {
      msg <- sprintf("POSITIVE CONTROL FAILED in arm %s: %s", A$name, paste(probs$hard, collapse = "; "))
      if (A$required) stop(msg)
      message(msg); notes <- c(notes, msg); invalid <- c(invalid, A$name); arm_ok <- FALSE
    }

    for (id in names(outs)) {
      cur <- outs[[id]]; O <- cur$O
      lab <- sprintf("%s/%s", A$name, id)
      base_row <- data.table(arm = A$name, outcome = id, overlap_flag = O$overlap,
                             rho_pair = ph$R[1, 2], rho_pair_rint_pearson = ph$rho_rint, rho_missingness_factor = ph$overlap_factor,
                             n_instruments = nrow(cur$bx))
      if (!arm_ok) {
        all_ests[[lab]] <- cbind(base_row, data.table(status = "ARM_CONTROL_FAILED", qhet_contrast = NA_real_))
        all_sets[[lab]] <- cbind(base_row, data.table(statistic = "KOH_subset", role = "PRE_SPECIFIED_FOR_VERDICT", status = "ARM_CONTROL_FAILED"))
        flush(); next
      }
      p <- mk_prob(cur$bx, cur$se, cur$by, cur$bys, ph$R)
      if (isTRUE(p$degenerate)) {
        message(sprintf("[%s] degenerate design: %s", lab, p$reason))
        all_ests[[lab]] <- cbind(base_row, data.table(status = paste("DEGENERATE:", p$reason), qhet_contrast = NA_real_))
        all_sets[[lab]] <- cbind(base_row, data.table(statistic = "KOH_subset", role = "PRE_SPECIFIED_FOR_VERDICT", status = "DEGENERATE"))
        flush(); next
      }
      # Harmonisation sign check on the external outcomes: the sum of the LDL-C and ApoB effects is
      # the total LDL effect, and it must be clearly causal. A flipped allele fails this.
      iv <- ivw_fit(p)
      total <- 2 * iv$th[2]
      if (O$kind != "regenie" && !(exp(total) > 1.1)) {
        msg <- sprintf("HARMONISATION CONTROL FAILED in %s: total LDL OR %.3f; expected clearly above 1", lab, exp(total))
        if (A$required) stop(msg)
        message(msg); notes <- c(notes, msg); invalid <- c(invalid, A$name)
        all_ests[[lab]] <- cbind(base_row, data.table(status = "HARMONISATION_FAILED", qhet_contrast = NA_real_))
        all_sets[[lab]] <- cbind(base_row, data.table(statistic = "KOH_subset", role = "PRE_SPECIFIED_FOR_VERDICT", status = "HARMONISATION_FAILED"))
        flush(); next
      }

      t0 <- Sys.time()
      qh0 <- qhet_fit(p, list(iv$th))
      if (is.null(qh0)) {
        all_ests[[lab]] <- cbind(base_row, data.table(status = "QHET_NO_CONVERGENCE", qhet_contrast = NA_real_))
        all_sets[[lab]] <- cbind(base_row, data.table(statistic = "KOH_subset", role = "PRE_SPECIFIED_FOR_VERDICT", status = "QHET_NO_CONVERGENCE"))
        flush(); next
      }
      pl <- plasmode(cur$bx, cur$se, cur$bys, ph$R, p, qh0$th, qh0$tau2, iv$th[1], o$plasmode_reps, o$ncores, plasmode_window)
      all_pl[[lab]] <- cbind(data.table(arm = A$name, outcome = id, plasmode_pass = pl$pass), pl$tab)
      message(sprintf("[%s] plasmode %s (%.0fs)", lab, if (pl$pass) "passed" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      print(pl$tab[statistic == "KOH_subset"])
      if (!pl$pass) {
        all_ests[[lab]] <- cbind(base_row, data.table(status = "PLASMODE_FAILED", qhet_contrast = NA_real_))
        all_sets[[lab]] <- cbind(base_row, data.table(statistic = "KOH_subset", role = "PRE_SPECIFIED_FOR_VERDICT", status = "PLASMODE_FAILED"))
        notes <- c(notes, sprintf("%s: plasmode calibration failed; no robust set is reported for it.", lab))
        flush(); next
      }

      t0 <- Sys.time()
      ra <- robust_analysis(p, grid_n = o$grid, gamma_min = o$gamma_min)
      if (isTRUE(ra$degenerate)) {
        all_ests[[lab]] <- cbind(base_row, data.table(status = paste("DEGENERATE:", ra$reason), qhet_contrast = NA_real_))
        all_sets[[lab]] <- cbind(base_row, data.table(statistic = "KOH_subset", role = "PRE_SPECIFIED_FOR_VERDICT", status = "DEGENERATE"))
        flush(); next
      }
      message(sprintf("[%s] sets done in %.0fs", lab, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      print(ra$sets[, .(statistic, status, set_class, intervals, includes_zero)])

      # qhet with the bootstrap and jackknife. The jackknife is what Sanderson et al. recommend,
      # the bootstrap is what was asked for, and a large gap between them is itself informative.
      RNGkind("L'Ecuyer-CMRG"); set.seed(42)
      bt <- boot::boot(seq_len(p$J), function(d, i) {
        f <- tryCatch(qhet_fit(sub_prob(p, i), list(qh0$th, iv$th)), error = function(e) NULL)
        if (is.null(f)) NA_real_ else f$th[1]
      }, R = o$nboot, parallel = "multicore", ncpus = o$ncores)
      RNGkind("Mersenne-Twister"); set.seed(42)
      bvals <- bt$t[, 1]; bfail <- mean(is.na(bvals))
      perc <- if (bfail <= 0.01) quantile(bvals, c(0.025, 0.975), na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)
      bca <- c(NA_real_, NA_real_)
      if (bfail == 0) {
        b <- tryCatch(boot::boot.ci(bt, type = "bca")$bca[4:5], error = function(e) NULL)
        if (!is.null(b)) bca <- b
      }
      jk <- unlist(parallel::mclapply(seq_len(p$J), function(k) {
        f <- tryCatch(qhet_fit(sub_prob(p, setdiff(seq_len(p$J), k)), list(qh0$th, iv$th)), error = function(e) NULL)
        if (is.null(f)) NA_real_ else f$th[1]
      }, mc.cores = o$ncores))
      jse <- if (any(is.na(jk))) NA_real_ else sqrt((p$J - 1) / p$J * sum((jk - mean(jk))^2))

      pk048 <- tryCatch({
        fmt <- MVMR::format_mvmr(BXGs = cur$bx, BYG = cur$by, seBXGs = cur$se, seBYG = cur$bys, RSID = rownames(cur$bx))
        x <- suppressWarnings(MVMR::qhet_mvmr(fmt, ph$R, CI = FALSE))
        x[1, 1] - x[2, 1]
      }, error = function(e) NA_real_)

      th_o <- as.numeric(p$Tm %*% ra$th_ivw)
      all_ests[[lab]] <- cbind(base_row, data.table(status = "OK"), as.data.table(ra$ests),
        data.table(
          ivw_or_exposure1 = exp(th_o[1]), ivw_or_exposure2 = exp(th_o[2]),
          exposure1 = A$traits[1], exposure2 = A$traits[2],
          qhet_boot_perc_lo = perc[1], qhet_boot_perc_hi = perc[2],
          qhet_boot_bca_lo = bca[1], qhet_boot_bca_hi = bca[2],
          qhet_boot_frac_failed = bfail, qhet_jack_se = jse,
          mvmr_pkg_version = as.character(packageVersion("MVMR")), mvmr_pkg_qhet_contrast = pk048,
          mvmr_pkg_note = "0.4.8 tau2 search is defective (see header); not used"))
      all_sets[[lab]] <- cbind(base_row, ra$sets)
      flush()

      # The rotation makes contrast precision depend directly on rho (var(D) is proportional to
      # 1-rho), and the two arms have very different rho. If the answer turns on rho, say so.
      if (A$name == PRIMARY_ARM && id == PRIMARY_OUT) {
        for (sh in c(-0.03, -0.01, 0.01, 0.03)) {
          R2 <- ph$R; R2[1, 2] <- R2[2, 1] <- min(ph$R[1, 2] + sh, 0.995)
          p2 <- mk_prob(cur$bx, cur$se, cur$by, cur$bys, R2)
          if (isTRUE(p2$degenerate)) next
          r2 <- robust_analysis(p2, grid_n = max(501L, o$grid %/% 4L), gamma_min = o$gamma_min, proj_every = 50L)
          if (isTRUE(r2$degenerate)) next
          all_rho[[paste(lab, sh)]] <- cbind(data.table(arm = A$name, outcome = id, rho_used = R2[1, 2], rho_shift = sh,
                                                        F_contrast = r2$ests$F_contrast_given_nuisance),
                                             r2$sets[statistic %in% c("KOH_subset", "AR_subset"), .(statistic, status, set_class, intervals, includes_zero)])
        }
        # KOH models balanced heterogeneity only, so one locus with a large CAD effect and a
        # non-average LDL/ApoB ratio could carry a bounded set. Reported beside the label, never in it.
        rf <- al$ref[match(rownames(cur$bx), al$ref$ID)]
        for (lc in names(LOCI_OUT)) {
          kp <- which(!locus_mask(rf$CHROM, rf$GENPOS, LOCI_OUT[[lc]]))
          hd <- data.table(arm = A$name, outcome = id, locus_removed = lc, n_removed = p$J - length(kp),
                           ids_removed = paste(setdiff(rownames(cur$bx), rownames(cur$bx)[kp]), collapse = ","))
          r3 <- if (length(kp) == p$J) NULL else robust_analysis(sub_prob(p, kp), grid_n = max(501L, o$grid %/% 4L), gamma_min = o$gamma_min, proj_every = 50L)
          # Keyed by arm/outcome too, or each outcome overwrites the last. Non-OK rows carry the set
          # columns as NA so that the table always has them.
          na_set <- data.table(set_class = NA_character_, intervals = NA_character_, includes_zero = NA)
          all_locus[[paste(lab, lc)]] <- if (is.null(r3)) cbind(hd, data.table(statistic = "KOH_subset", status = "NO_VARIANT_IN_REGION"), na_set)
            else if (isTRUE(r3$degenerate)) cbind(hd, data.table(statistic = "KOH_subset", status = "DEGENERATE"), na_set)
            else cbind(hd, r3$sets[statistic %in% c("KOH_subset", "AR_subset"), .(statistic, status, set_class, intervals, includes_zero)])
        }
        flush()
      }
    }
  }
  if (!length(all_sets)) stop("no arm produced output")

  sets <- rbindlist(all_sets, fill = TRUE); ests <- rbindlist(all_ests, fill = TRUE)
  invalid <- unique(invalid)
  vd <- readout_verdict(sets, ests, invalid)
  rho_note <- character(0)
  if (length(all_rho)) {
    rr <- rbindlist(all_rho, fill = TRUE)[statistic == "KOH_subset"]
    pr <- sets[arm == PRIMARY_ARM & outcome == PRIMARY_OUT & statistic == "KOH_subset"]
    if (nrow(pr) && (any(rr$set_class != pr$set_class, na.rm = TRUE) || any(rr$includes_zero != pr$includes_zero, na.rm = TRUE)))
      rho_note <- "FRAGILE: the primary set changes class or zero-inclusion within +-0.03 of the phenotypic correlation. See rho_sensitivity.csv."
    else if (nrow(pr)) rho_note <- "Primary set class and zero-inclusion are unchanged within +-0.03 of the phenotypic correlation."
  }
  locus_note <- character(0)
  if (length(all_locus)) {
    lr <- rbindlist(all_locus, fill = TRUE)[statistic == "KOH_subset" & arm == PRIMARY_ARM & outcome == PRIMARY_OUT]
    pr <- sets[arm == PRIMARY_ARM & outcome == PRIMARY_OUT & statistic == "KOH_subset"]
    lc <- if (nrow(pr) == 1L) lr[status == "OK" & (set_class != pr$set_class | includes_zero != pr$includes_zero)] else lr[0]
    locus_note <- if (nrow(lc)) sprintf("LOCUS-SENSITIVE: the primary set changes class or zero-inclusion without %s. See robust_contrast_leave_locus.csv.", paste(lc$locus_removed, collapse = ", "))
                  else sprintf("Primary set class and zero-inclusion hold without the LPA and APOE loci (variants removed: %s).", paste(sprintf("%s %d", lr$locus_removed, lr$n_removed), collapse = ", "))
  }
  out <- c(sprintf("VERDICT (readout section C): %s", vd$label), vd$lines, rho_note, locus_note, primary_context(sets, ests), notes,
           sprintf("calibration gate: passed for script md5 %s", md5),
           "UKB rows have uncorrected sample overlap with the exposures and do not enter the verdict except for sign.",
           "KOH and qhet model balanced heterogeneity: random direct effects, mean zero, independent of the instrument strengths. Directional pleiotropy is not handled by any statistic here.")
  writeLines(out, file.path(o$out_dir, "robust_contrast_verdict.txt"))
  message("\n", paste(out, collapse = "\n"))
  invisible(list(sets = sets, ests = ests, verdict = vd, lines = out))
}

# ------------------------------------------------------------
# Selftest
# ------------------------------------------------------------
sim_small <- function(seed, J = 150L, rho = 0.851, corg = 0.97, theta = c(0.593, -0.329, 0.079), tau_mult = 0, K = 3L) {
  set.seed(seed)
  sc <- list(J = J, K = K, nx = 400000, ncase = 47000, nctrl = 370000, corg = corg, rho = rho, tau_mult = tau_mult,
             het = "normal", select = FALSE)
  gm <- sim_gamma(sc); G <- gm$G
  sxm <- gm$sx %o% rep(1, K); sxm[, 2] <- sxm[, 2] * 1.01
  R <- scenario_rho(sc)
  bx <- G + (matrix(rnorm(length(G)), nrow(G), K) %*% chol(R)) * sxm
  sy <- gm$base / sqrt(417000 * (47000 / 417000) * (370000 / 417000))
  by <- as.numeric(G %*% theta) + rnorm(J, 0, sqrt(tau_mult * median(sy^2))) + rnorm(J, 0, sy)
  dimnames(bx) <- dimnames(sxm) <- list(sprintf("rs%d", seq_len(J)), c("LDL_C", "APOB", "TRIGLY")[seq_len(K)])
  list(bx = bx, se = sxm, by = by, sy = sy, R = R, G = G)
}

check_inversion <- function(sc_list, n_each, grid_n = 301L) {
  bad <- 0L; n <- 0L
  for (sc in sc_list) for (i in seq_len(n_each)) {
    set.seed(5000L + i)
    d <- sim_small(5000L + i, rho = sc$rho, corg = sc$corg, theta = sc$theta, tau_mult = sc$tau_mult)
    p <- mk_prob(d$bx, d$se, d$by, d$sy, d$R)
    ra <- robust_analysis(p, grid_n = grid_n, proj_every = 50L)
    if (isTRUE(ra$degenerate)) stop("inversion check: degenerate at i=", i)
    dt <- sc$theta[1] - sc$theta[2]
    ct <- crit_table(p); e <- eval_delta(p, dt, list(), do_oh = TRUE)
    v <- stat_vec(e, ct$a)
    for (st in c("AR_subset", "K_subset", "LC_subset", "KOH_subset")) {
      row <- ra$sets[statistic == st]
      if (row$status != "OK") { bad <- bad + 1L; next }
      m <- regmatches(row$intervals, gregexpr("-?[0-9.]+|-?Inf", row$intervals))[[1]]
      if (row$set_class == "empty") inset <- FALSE else {
        m <- matrix(as.numeric(m), ncol = 2, byrow = TRUE)
        inset <- any(m[, 1] <= dt & dt <= m[, 2])
      }
      direct <- v[st] <= ct[[st]]
      # A truth within rounding of an interval end can differ; everything else is a bug.
      near <- row$set_class != "empty" && min(abs(c(m) - dt)) < 2e-3
      n <- n + 1L
      if (!is.na(direct) && inset != direct && !near) {
        bad <- bad + 1L
        message(sprintf("  inversion mismatch: %s set=%s direct stat %.3f crit %.3f truth %.4f", st, row$intervals, v[st], ct[[st]], dt))
      }
    }
  }
  list(bad = bad, n = n)
}

write_fixture <- function(root, n_iv = 151L, seed = 42L) {
  set.seed(seed)
  dirs <- file.path(root, c("sumstat", "naive", "nmr", "ext", "ph")); for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d <- sim_small(seed, J = n_iv)
  chr <- rep(1:22, length.out = n_iv); pos <- 1e6 + 2e7 * (ave(seq_len(n_iv), chr, FUN = seq_along) - 1)
  # One instrument inside the LPA region and one inside APOE, so the leave-locus refits run.
  # Both stay more than 10 Mb from their chromosome's other instruments.
  if (n_iv >= 138L) { pos[138] <- 160.9e6; pos[63] <- 45.4e6 }
  ids <- sprintf("rs%d", seq_len(n_iv))
  a1 <- rep("A", n_iv); a0 <- rep("G", n_iv)
  # Three palindromic instruments: two resolvable by frequency, one ambiguous and to be dropped.
  a1[c(5, 6, 7)] <- c("A", "C", "A"); a0[c(5, 6, 7)] <- c("T", "G", "T")
  frq <- rep(0.2, n_iv); frq[7] <- 0.48
  # A very strong non-lead a few Mb from each lead (must be clumped away, and its LOG10P is high
  # enough that a p-value would have underflowed); and null rows that must never be selected.
  mkgwas <- function(k, tr, dir, tag = "lipids_%s.txt.gz", extra = TRUE, b = d$bx[, k], s = d$se[, k]) {
    z <- b / s
    lp <- pmin(-(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10), 1000)
    # Every planted instrument is significant in the first trait, so the expected count is exact.
    if (k == 1L) lp <- pmax(lp, 8)
    main <- data.table(CHROM = chr, GENPOS = pos, ID = ids, ALLELE0 = a0, ALLELE1 = a1, A1FREQ = frq, N = 400000L,
                       TEST = "ADD", BETA = b, SE = s, CHISQ = z^2, LOG10P = lp, EXTRA = "NA")
    if (extra) {
      nb <- main[seq(1, n_iv, by = 5)]
      nb[, `:=`(GENPOS = GENPOS + 3e6, ID = paste0("nb", ID), BETA = BETA * 0.2, LOG10P = LOG10P * 0.2)]
      nl <- main[seq(1, n_iv, by = 3)]
      nl[, `:=`(GENPOS = GENPOS + 1.3e7, ID = paste0("nl", ID), BETA = 1e-4, LOG10P = 0.5)]
      main <- rbind(main, nb, nl)
    }
    f <- file.path(dir, sprintf(tag, tr)); fwrite(main, sub("\\.gz$", "", f), sep = " ", na = "NA", quote = FALSE)
    system2("gzip", c("-f", shQuote(sub("\\.gz$", "", f))))
  }
  for (k in 1:3) mkgwas(k, c("LDL_C", "APOB", "TRIGLY")[k], dirs[1])
  mko <- function(dir) {
    b <- d$by; s <- d$sy
    z <- b / s; lp <- -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
    o <- data.table(CHROM = chr, GENPOS = pos, ID = ids, ALLELE0 = a0, ALLELE1 = a1, A1FREQ = frq, N = 400000L, TEST = "ADD",
                    BETA = b, SE = s, CHISQ = z^2, LOG10P = lp, EXTRA = "NA")
    # Half of the outcome rows report the other allele. Harmonisation must flip them.
    fl <- seq(2, n_iv, by = 2)
    o[fl, `:=`(ALLELE0 = a1[fl], ALLELE1 = a0[fl], BETA = -BETA, A1FREQ = 1 - A1FREQ)]
    f <- file.path(dir, "outcomes_CAD.txt"); fwrite(o, f, sep = " "); system2("gzip", c("-f", shQuote(f)))
  }
  mko(dirs[1])
  # FinnGen format: tab-separated, multi-rsID cells, duplicates, an unmatched allele row, and a real
  # "ref/alt" flip on half the rows.
  fl <- seq(2, n_iv, by = 2)
  alt <- a1; ref <- a0; bfg <- d$by; af <- frq
  pal <- 5:7
  fl <- setdiff(fl, pal)
  alt[fl] <- a0[fl]; ref[fl] <- a1[fl]; bfg[fl] <- -bfg[fl]; af[fl] <- 1 - af[fl]
  rsl <- ids; rsl[seq(3, n_iv, by = 7)] <- paste0(rsl[seq(3, n_iv, by = 7)], ",rs9", seq_along(seq(3, n_iv, by = 7)))
  fg <- data.table(`#chrom` = chr, pos = pos, ref = ref, alt = alt, rsids = rsl, nearest_genes = "GENE1,GENE2", pval = 0.5,
                   mlogp = 0.3, beta = bfg, sebeta = d$sy, af_alt = af, af_alt_cases = af, af_alt_controls = af)
  fg <- rbind(fg, fg[10][, `:=`(alt = "T", ref = "C", beta = 9)])  # multi-allelic row that must not match
  ff <- file.path(dirs[4], "finngen.tsv"); fwrite(fg, ff, sep = "\t"); system2("gzip", c("-f", shQuote(ff)))
  # Phenotype file with the right correlation among the traits.
  n <- 20000L
  Z <- matrix(rnorm(n * 3), n, 3) %*% chol(d$R)
  pheno <- data.table(FID = 1:n, IID = 1:n, LDL_C = Z[, 1], APOB = Z[, 2], TRIGLY = Z[, 3], OTHER = rnorm(n))
  fwrite(pheno, file.path(dirs[5], "pheno.txt"), sep = " ")
  list(sumstat = dirs[1], naive = dirs[2], nmr = dirs[3], finngen = paste0(ff, ".gz"), pheno = file.path(dirs[5], "pheno.txt"), data = d)
}

# The NMR arm is the one route with per-chromosome files and the CARDIoGRAM format, so it gets its
# own files: two exposures, three palindromes (dropped under palin="drop"), half the rows flipped.
write_nmr_fixture <- function(root, seed = 42L) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  J <- 62L
  d <- sim_small(seed, J = J, K = 2L, rho = 0.955, theta = c(-0.37, 0.91))
  chr <- rep(1:22, length.out = J); pos <- 1e6 + 2e7 * (ave(seq_len(J), chr, FUN = seq_along) - 1)
  ids <- sprintf("rs%d", seq_len(J))
  a1 <- rep("A", J); a0 <- rep("G", J); a0[c(5, 6, 7)] <- "T"
  for (k in 1:2) {
    z <- d$bx[, k] / d$se[, k]
    lp <- pmin(-(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10), 1000)
    if (k == 1L) lp <- pmax(lp, 8)
    tab <- data.table(CHROM = chr, GENPOS = pos, ID = ids, ALLELE0 = a0, ALLELE1 = a1, A1FREQ = 0.2, N = 212797L, TEST = "ADD",
                      BETA = d$bx[, k], SE = d$se[, k], CHISQ = z^2, LOG10P = lp, EXTRA = "NA")
    for (ch in 1:22) {
      f <- file.path(root, sprintf("nmr_chr%d_%s.regenie", ch, c("LDL_C_NMR", "LDL_P")[k]))
      fwrite(tab[CHROM == ch], f, sep = " "); system2("gzip", c("-f", shQuote(f)))
    }
  }
  fl <- seq(2, J, by = 2)
  ea <- a1; oa <- a0; b <- d$by; frq <- rep(0.2, J)
  fl <- setdiff(fl, 5:7); ea[fl] <- a0[fl]; oa[fl] <- a1[fl]; b[fl] <- -b[fl]; frq[fl] <- 0.8
  cg <- data.table(markername = ids, chr = chr, bp_hg19 = pos, effect_allele = ea, noneffect_allele = oa, effect_allele_freq = frq,
                   median_info = 1, model = "FIXED", beta = b, se_dgc = d$sy, p_dgc = 0.5, het_pvalue = 0.5, n_studies = 48L)
  fwrite(cg, file.path(root, "cardiogram.txt"), sep = "\t")
  n <- 20000L
  Z <- matrix(rnorm(n * 2), n, 2) %*% chol(d$R[1:2, 1:2])
  fwrite(data.table(FID = 1:n, IID = 1:n, LDL_C_NMR = Z[, 1], LDL_P = Z[, 2]), file.path(root, "pheno_nmr.txt"), sep = " ")
  list(data = d, dir = root, cardiogram = file.path(root, "cardiogram.txt"), pheno = file.path(root, "pheno_nmr.txt"))
}

selftest <- function(o) {
  ok <- function(cond, what) {
    if (!isTRUE(cond)) stop("[selftest] FAILED: ", what, call. = FALSE)
    message("[selftest] ok: ", what)
  }
  message("---- 1. algebra ----")
  d <- sim_small(101)
  p <- mk_prob(d$bx, d$se, d$by, d$sy, d$R)
  th <- c(0.7, 0.14, 0.09)
  ks <- kstat(p, th)
  num <- vapply(1:3, function(k) { e <- numeric(3); e[k] <- 1e-6; (S_fn(p, th + e) - S_fn(p, th - e)) / 2e-6 }, numeric(1))
  ok(all(abs(num + 2 * as.numeric(ks$sc)) / (abs(num) + 1) < 1e-4), "-2 D'Omega^-1 g is the gradient of S")
  # Andrews printed K* equals K when the variance is inverted correctly (see header).
  x <- solve_scaled(ks$info, ks$sc)
  ok(abs(sum(x * (ks$info %*% x)) - ks$K) < 1e-8, "K* with inverted middle equals K")
  ok(abs(ks$S - S_fn(p, th)) < 1e-10, "kstat S equals S_fn")
  # Rotation changes nothing: the same S at corresponding parameters; delta recovered as LDL-APOB.
  p0 <- mk_prob(d$bx, d$se, d$by, d$sy, d$R, rotate = FALSE)
  th0 <- as.numeric(p$Tm %*% th)
  ok(abs(S_fn(p, th) - S_fn(p0, th0)) < 1e-8, "rotation leaves S unchanged")
  ok(abs((th0[1] - th0[2]) - th[1]) < 1e-12, "first rotated coordinate is theta1 minus theta2")
  ok(abs(om_fn(p, th)[3] - (d$sy[3]^2 + drop(t(th0) %*% (d$R * outer(d$se[3, ], d$se[3, ])) %*% th0))) < 1e-14, "rotation carries the SNP-exposure covariance (phenocov form)")
  if (requireNamespace("MendelianRandomization", quietly = TRUE)) {
    f <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = d$bx, bxse = d$se, by = d$by, byse = d$sy))
    iv <- ivw_fit(p)
    th_o <- as.numeric(p$Tm %*% iv$th)
    ok(all(abs(th_o - f@Estimate) < 1e-8), "own IVW in rotated coordinates reproduces mr_mvivw estimates")
    ok(abs(sqrt(ivw_fit(p0)$V[1, 1]) - f@StdError[1]) < 1e-8, "own random-effects SE reproduces mr_mvivw")
  } else stop("[selftest] MendelianRandomization needed")
  # The linear combination distribution, against brute force.
  set.seed(1); mc <- (1 + 0.01) * rchisq(4e6, 1) + 0.01 * rchisq(4e6, 148)
  ok(abs(plc(5, 0.01, 148) - mean(mc <= 5)) < 1e-3, "plc matches Monte Carlo")
  ok(abs(qlc(0.95, 0.01, 148) - quantile(mc, 0.95)) < 0.02, "qlc matches Monte Carlo")
  a <- a_for_gamma(0.05, 148)
  ok(abs(plc(qchisq(0.95, 1), a, 148) - 0.90) < 1e-6, "a_for_gamma gives the stated distortion")
  # Stated null distributions of the raw statistics.
  v <- sapply(1:2000, function(i) { d2 <- sim_small(20000L + i, theta = c(0.13, 0.13, 0.08)); p2 <- mk_prob(d2$bx, d2$se, d2$by, d2$sy, d2$R); k <- kstat(p2, c(0, 0.13, 0.08)); c(k$S, k$K) })
  ok(abs(mean(v[1, ]) - 150) < 1.5, "S at the truth has mean J (full vector)")
  ok(abs(mean(v[2, ]) - 3) < 0.3, "K at the truth has mean K (full vector)")
  ok(abs(cor(v[2, ], v[1, ] - v[2, ])) < 0.1, "K and S-K are uncorrelated at the truth")

  message("---- 2. degenerate designs give NA, never a number ----")
  ok(isTRUE(mk_prob(d$bx[1:3, ], d$se[1:3, ], d$by[1:3], d$sy[1:3], d$R)$degenerate), "J<=K is degenerate")
  ok(isTRUE(mk_prob(d$bx, d$se, d$by, d$sy, matrix(c(1, 1, .3, 1, 1, .3, .3, .3, 1), 3))$degenerate), "singular correlation is degenerate")
  dz <- d$se; dz[1, 1] <- 0
  ok(isTRUE(mk_prob(d$bx, dz, d$by, d$sy, d$R)$degenerate), "zero SE is degenerate")
  bcol <- d$bx; bcol[, 2] <- bcol[, 1]
  pc <- mk_prob(bcol, d$se, d$by, d$sy, d$R)
  ra <- robust_analysis(pc, grid_n = 201L, proj_every = 50L)
  ok(isTRUE(ra$degenerate) || all(!is.na(ra$sets$status)), "exact collinearity in the estimates does not crash")
  ok(is.na(kstat(list(by = 1:4 * 0, bs = matrix(0, 4, 2), s2y = rep(1, 4), SG = matrix(0, 4, 4), SGa = list(matrix(0, 4, 2), matrix(0, 4, 2)), K = 2, J = 4), c(0, 0))$K), "K is NA when the Jacobian is zero")

  message("---- 3. known set shapes ----")
  # Exposures with no conditional signal: delta is unidentified. If the set is unbounded it must be
  # reported so, and must not depend on how far the grid looks.
  dw <- sim_small(202, corg = 1, rho = 0.958, theta = c(0.13, 0.13, 0.08))
  pw <- mk_prob(dw$bx, dw$se, dw$by, dw$sy, dw$R)
  r1 <- robust_analysis(pw, grid_n = 801L, proj_every = 40L)
  ok(r1$sets[statistic == "AR_subset", set_class] %in% c("unbounded", "whole_line"), sprintf("no conditional strength: AR is unbounded (F=%.2f, S_inf=%.1f)", r1$ests$F_contrast_given_nuisance, r1$ests$S_prof_at_inf))
  ok(any(grepl("Inf", r1$sets[statistic == "AR_subset", intervals])), "unbounded set is written with Inf and not a grid edge")
  # A bounded set must be stable under the grid; an edge effect would move it with the grid scale.
  ds <- sim_small(203)
  ps <- mk_prob(ds$bx, ds$se, ds$by, ds$sy, ds$R)
  r2 <- robust_analysis(ps, grid_n = 801L, proj_every = 40L)
  r3 <- robust_analysis(ps, grid_n = 1601L, proj_every = 80L)
  a2 <- r2$sets[statistic == "AR_subset"]; a3 <- r3$sets[statistic == "AR_subset"]
  ok(a2$set_class == "bounded", sprintf("strong design: AR bounded (%s, F=%.1f)", a2$intervals, r2$ests$F_contrast_given_nuisance))
  ok(abs(a2$lower - a3$lower) < 1e-3 && abs(a2$upper - a3$upper) < 1e-3, "bounded endpoints do not depend on grid density")
  ok(abs(r2$ests$S_prof_at_inf / (p$J - p$K + 1) - r2$ests$F_contrast_given_nuisance) / r2$ests$F_contrast_given_nuisance < 0.01, "S_prof at infinity equals (J-K+1) times conditional F")
  # Strong design, big contrast: the set excludes zero with the right sign; and non-robust agrees.
  k2 <- r2$sets[statistic == "KOH_subset"]
  ok(k2$status == "OK" && identical(k2$set_sign, "+"), sprintf("true delta 0.92: KOH positive, %s", k2$intervals))
  # Extreme heterogeneity: AR empty must be labelled empty.
  dh <- sim_small(204, tau_mult = 30)
  ph <- mk_prob(dh$bx, dh$se, dh$by, dh$sy, dh$R)
  rh <- robust_analysis(ph, grid_n = 801L, proj_every = 40L)
  ok(rh$sets[statistic == "AR_subset", set_class] == "empty", "extreme heterogeneity: AR is empty")
  ok(is.na(rh$sets[statistic == "AR_subset", lower]), "empty AR gives no endpoints")
  ok(rh$sets[statistic == "KOH_subset", status] == "OK" && rh$sets[statistic == "KOH_subset", set_class] != "empty", "extreme heterogeneity: KOH still gives a set")
  ok(rh$ests$qhet_tau2 > 0, "extreme heterogeneity: qhet tau2 positive")
  # The K=2 path, used by the NMR arm.
  d2e <- sim_small(205, K = 2L, J = 62L, theta = c(-0.37, 0.91), rho = 0.955, tau_mult = 1)
  r2e <- robust_analysis(mk_prob(d2e$bx, d2e$se, d2e$by, d2e$sy, d2e$R[1:2, 1:2]), grid_n = 401L, proj_every = 20L)
  ok(!isTRUE(r2e$degenerate) && r2e$sets[statistic == "KOH_subset", status] == "OK", "two-exposure path runs")

  message("---- 4. set construction agrees with the direct test ----")
  chk <- check_inversion(list(list(rho = 0.851, corg = 0.97, theta = c(0.593, -0.329, 0.079), tau_mult = 0),
                              list(rho = 0.958, corg = 0.97, theta = c(0.13, 0.13, 0.08), tau_mult = 2.5),
                              list(rho = 0.958, corg = 1, theta = c(0.13, 0.13, 0.08), tau_mult = 0)), n_each = 6L)
  ok(chk$bad == 0L, sprintf("grid inversion matches the direct test (%d comparisons)", chk$n))

  message("---- 5. qhet ----")
  dq <- sim_small(301, tau_mult = 2)
  pq <- mk_prob(dq$bx, dq$se, dq$by, dq$sy, dq$R)
  q <- qhet_fit(pq, list(ivw_fit(pq)$th))
  ok(!is.null(q) && q$tau2 > 0 && abs(q$S - (pq$J - pq$K)) < 1e-3, "qhet: minimised Q equals L-K at tau2 > 0")
  # --mvmr049 names a library holding 0.4.9. TACC has only 0.4.8, so the comparison runs where the
  # library exists and is skipped out loud elsewhere; the namespace is unloaded so that everything
  # after this point gets the installed version.
  if (nzchar(o$mvmr049)) {
    if ("MVMR" %in% loadedNamespaces()) unloadNamespace("MVMR")
    ns <- loadNamespace("MVMR", lib.loc = o$mvmr049)
    ok(identical(as.character(getNamespaceVersion(ns)), "0.4.9"), "--mvmr049 holds MVMR 0.4.9")
    fmt <- ns$format_mvmr(BXGs = dq$bx, BYG = dq$by, seBXGs = dq$se, seBYG = dq$sy, RSID = rownames(dq$bx))
    x <- suppressWarnings(ns$qhet_mvmr(fmt, dq$R, CI = FALSE, ncores = 1))
    ok(abs((x[1, 1] - x[2, 1]) - q$th[1]) < 2e-3, sprintf("own qhet matches MVMR 0.4.9 (%.4f vs %.4f)", q$th[1], x[1, 1] - x[2, 1]))
    unloadNamespace("MVMR")
  } else message("[selftest] note: --mvmr049 not given, comparison with MVMR 0.4.9 skipped")
  # A true contrast is recovered better by qhet than a defective tau2 path would: sanity of level.
  rq <- sapply(1:120, function(i) { dd <- sim_small(30000L + i, tau_mult = 2); pp <- mk_prob(dd$bx, dd$se, dd$by, dd$sy, dd$R); qhet_fit(pp, list(ivw_fit(pp)$th))$th[1] })
  ok(abs(mean(rq) - 0.922) < 3 * sd(rq) / sqrt(120) + 0.03, sprintf("qhet unbiased at delta 0.922 (mean %.3f)", mean(rq)))

  message("---- 6. coverage on simulated collinear exposures ----")
  # The strength bands are what make the scenarios mean what their names say.
  bands <- list(F1_F48_rho851_null = c(35, 65), G1_F105_rho958_null = c(75, 140), D1_nmr_J62_null = c(50, 105),
                H1_F10_rho851_null = c(6, 15), H3_F3_rho851_null = c(1.8, 5), C1_cor100_rho958_null = c(0.5, 1.5))
  all_sc <- scenarios()
  ok(!anyDuplicated(vapply(all_sc, function(s) s$id, "")), "scenario ids are unique")
  for (s in Filter(function(s) s$id %in% c(names(bands), "E1_selected_rho851_het_alt"), all_sc)) {
    m <- do.call(rbind, lapply(1:12, function(i) sim_rep(s, 7000L + i)))
    if (!is.null(bands[[s$id]])) ok(mean(m[, "F_D"]) >= bands[[s$id]][1] && mean(m[, "F_D"]) <= bands[[s$id]][2],
                                    sprintf("%s has the conditional F it claims (mean %.1f)", s$id, mean(m[, "F_D"])))
    else ok(abs(mean(m[, "J"]) - 150) < 10, sprintf("selection leaves J near 150 (mean %.0f)", mean(m[, "J"])))
  }
  # Coverage may not be flattered by dropping the replicates where a statistic failed.
  stn <- c("AR_subset", "AR_projection", "K_subset", "LC_subset", "KOH_subset", "WALD_GMM")
  cl <- c(paste0("cover_", c(stn, "KOH_projection", "IVW_RE")), paste0("accept0_", stn), paste0("acceptinf_", stn), "F_D", "S_min_over_df", "tau2hat", "J", "ok")
  fm <- matrix(1, 1000, length(cl), dimnames = list(NULL, cl))
  fm[, grep("^cover_", cl)] <- rep(c(rep(1, 19), 0), 50)
  fm[1:20, "cover_KOH_subset"] <- NA
  fs <- summarise_scenario(all_sc[[1]], list(mat = fm, n_failed = 0L), 1000L)
  ok(isFALSE(fs[statistic == "KOH_subset", pass]), "a gated statistic failing in 2% of replicates does not pass")
  ok(isTRUE(fs[statistic == "K_subset", pass]), "the same coverage without failures passes")
  ok(!rep_ok(NULL) && !rep_ok(structure("x", class = "try-error")) && !rep_ok(c(ok = 0)) && rep_ok(c(a = 1, ok = 1)), "dead forks and failed replicates are not counted")

  sel <- c("A1_cor97_rho851_null", "A2_cor97_rho851_alt", "B2_het_cor97_rho851_alt", "B3_het_cor97_rho958_null",
           "C1_cor100_rho958_null", "D2_nmr_J62_het_alt", "F1_F48_rho851_null", "F4_het_F48_rho851_alt", "H2_het_F10_rho851_alt", "H3_F3_rho851_null")
  sc <- Filter(function(s) s$id %in% sel, all_sc)
  dirc <- file.path(o$selftest_dir, "calib")
  cal <- run_calibration(dirc, o$selftest_reps, o$ncores, sc, n_proj = 100L, min_reps_for_pass = 1000L)
  g <- cal$tab[gated == TRUE]
  print(g[, .(scenario, statistic, coverage, mc_se, n_failed)])
  ok(all(g$coverage >= COV_LO & g$coverage <= COV_HI), sprintf("all gated statistics inside %.2f-%.2f", COV_LO, COV_HI))
  # The point of the exercise: the non-robust interval should break where the robust ones hold.
  cc <- cal$tab[scenario == "C1_cor100_rho958_null"]
  message(sprintf("[selftest] no conditional strength: Wald %.3f, IVW %.3f, KOH %.3f; KOH unbounded %.2f",
                  cc[statistic == "WALD_GMM", coverage], cc[statistic == "IVW_RE", coverage], cc[statistic == "KOH_subset", coverage],
                  cc[statistic == "KOH_subset", rate_unbounded]))
  ok(cc[statistic == "AR_subset", rate_unbounded] > 0.8, "no conditional strength: robust sets are usually unbounded")
  al <- cal$tab[scenario == "A2_cor97_rho851_alt"]
  ok(al[statistic == "KOH_subset", rate_excludes_zero] > 0.5, sprintf("strong design: KOH has power (%.2f)", al[statistic == "KOH_subset", rate_excludes_zero]))

  message("---- 6b. harmonisation cases the fixture does not reach ----")
  rfh <- data.table(ID = c("a", "b", "c", "d"), ALLELE1 = c("A", "A", "A", "A"), ALLELE0 = c("T", "T", "G", "T"), A1FREQ = c(.40, .30, .30, .30))
  oh <- data.table(rsids = c("a", "b", "b", "c", "d"), alt = c("A", "A", "A", "A", "A"), ref = c("T", "C", "T", "G", "T"),
                   beta = c(.1, .2, .3, .4, .5), sebeta = .01, af_alt = c(.53, .3, .3, .3, .7))
  rh <- outcome_external(copy(oh), rfh, c("a", "b", "c", "d"), FORMATS$finngen, "freq", TRUE)
  ok(is.na(rh$b[1]), "palindrome ambiguous only in the outcome is dropped, not flipped")
  ok(isTRUE(all.equal(rh$b[2], .3)), "a row failing the allele match is not rescued by its frequency")
  ok(isTRUE(all.equal(rh$b[3], .4)), "a non-palindromic variant is untouched")
  ok(isTRUE(all.equal(rh$b[4], -.5)), "a palindrome clearly on the other side of 0.5 is flipped")
  ok(identical(locus_mask(c(6, 6, 19, 1), c(160e6, 100e6, 45e6, 160e6), LOCI_OUT$LPA_APOE), c(TRUE, FALSE, TRUE, FALSE)), "locus mask finds LPA and APOE only")
  message("---- 7. verdict rules ----")
  mkset <- function(cl, inc, sg, st = "OK", a = "biochem", o2 = "FINNGEN", iv = "x", zm = NA)
    data.table(arm = a, outcome = o2, statistic = "KOH_subset", status = st, set_class = cl, includes_zero = inc, set_sign = sg, intervals = iv, zero_component_has_S_minimum = zm)
  mkest <- function(q, a = "biochem", o2 = "FINNGEN", st = "OK") data.table(arm = a, outcome = o2, status = st, qhet_contrast = q)
  all3 <- function(pr, q2 = 1, q3 = 1, s2 = "+") {
    s <- rbind(pr, mkset("bounded", FALSE, s2, a = "statinnaive"), mkset("bounded", FALSE, s2, a = "nmr"))
    e <- rbind(mkest(1), mkest(1, o2 = "UKB_CAD"), mkest(q2, "statinnaive"), mkest(q3, "nmr"))
    readout_verdict(s, e)$label
  }
  ok(all3(mkset("unbounded", TRUE, NA_character_)) == "SUPPORTS_NOT_IDENTIFIABLE", "unbounded -> SUPPORTS")
  ok(all3(mkset("whole_line", TRUE, NA_character_)) == "SUPPORTS_NOT_IDENTIFIABLE", "whole line -> SUPPORTS")
  ok(all3(mkset("bounded", TRUE, NA_character_)) == "SUPPORTS_NOT_IDENTIFIABLE", "includes 0 -> SUPPORTS")
  ok(all3(mkset("bounded", FALSE, "+")) == "CONTRADICTS_NOT_IDENTIFIABLE", "bounded, excludes 0, agree -> CONTRADICTS")
  ok(all3(mkset("bounded", FALSE, "+"), q3 = -1) == "UNRESOLVED", "bounded, excludes 0, sign disagrees -> not CONTRADICTS")
  ok(all3(mkset("bounded", FALSE, "mixed")) == "UNRESOLVED", "set on both sides of 0 -> not CONTRADICTS")
  ok(all3(mkset("empty", FALSE, NA_character_)) == "INVALID", "empty primary -> INVALID, not precise")
  ok(all3(mkset(NA, NA, NA, st = "NUMERIC_FAILURE")) == "INVALID", "numeric failure -> INVALID")
  v1 <- readout_verdict(rbind(mkset("bounded", FALSE, "+"), mkset("bounded", FALSE, "+", a = "statinnaive")), rbind(mkest(1), mkest(1, "statinnaive")))
  ok(v1$label == "UNRESOLVED", "arm missing -> not CONTRADICTS")
  arset <- data.table(arm = "biochem", outcome = "FINNGEN", statistic = "AR_subset", status = "OK", set_class = "empty", includes_zero = FALSE, set_sign = NA_character_, intervals = "empty", zero_component_has_S_minimum = NA)
  v2 <- readout_verdict(rbind(mkset("unbounded", TRUE, NA_character_), arset), mkest(1))
  ok(any(grepl("never of a precise estimate", v2$lines)), "empty AR is read as heterogeneity, never precision")

  message("---- 8. full pipeline on a synthetic fixture ----")
  fx <- write_fixture(file.path(o$selftest_dir, "fx"))
  outd <- file.path(o$selftest_dir, "out")
  unlink(outd, recursive = TRUE); dir.create(outd, recursive = TRUE)
  fake_md5 <- "fixture-md5"
  gatef <- file.path(outd, "calibration_status.txt")
  opt2 <- list(out_dir = outd, cache_dir = file.path(outd, "cache"), plasmode_reps = 120L, ncores = o$ncores, nboot = 40L, grid = 201L, gamma_min = 0.05, arms = "biochem")
  mkarm <- function(ctl, nm = "biochem", req = TRUE, tr = c("LDL_C", "APOB", "TRIGLY"), dir = fx$sumstat)
    list(biochem = list(name = nm, traits = tr, files = function(t) file.path(dir, sprintf("lipids_%s.txt.gz", t)), required = req, pheno = fx$pheno, keep = NULL,
                        outcomes = list(list(id = "FINNGEN", kind = "finngen", files = fx$finngen, overlap = "none", palin = "freq"),
                                        list(id = "UKB_CAD", kind = "regenie", files = file.path(dir, "outcomes_CAD.txt.gz"), overlap = "UKB_ON_UKB_UNCORRECTED")),
                        control = ctl))
  f0 <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = fx$data$bx, bxse = fx$data$se, by = fx$data$by, byse = fx$data$sy))
  # unname: @Estimate carries lm coefficient names, and c(LDL_C = named) would rename the element.
  good <- list(outcome = "UKB_CAD", n = 151L, b = c(LDL_C = unname(f0@Estimate[1]), APOB = unname(f0@Estimate[2])), rho = fx$data$R[1, 2], condF = c(LDL_C = NA_real_)[0])
  # (a) no calibration file: must refuse before reading anything.
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(good), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("refusing to run", e), "no calibration file -> refuses")
  # (b) failed calibration, stale script, too few reps: all refuse.
  nsc <- sprintf("n_scenarios=%d", length(scenarios()))
  writeLines(c("PASSED=FALSE", sprintf("script_md5=%s", fake_md5), "n_reps=2000", nsc), gatef)
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(good), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("not passed", e), "failed calibration -> refuses")
  writeLines(c("PASSED=TRUE", "script_md5=OLDVERSION", "n_reps=2000", nsc), gatef)
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(good), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("different version", e), "stale calibration -> refuses")
  writeLines(c("PASSED=TRUE", sprintf("script_md5=%s", fake_md5), "n_reps=50", nsc), gatef)
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(good), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("fewer than 1000", e), "too few calibration reps -> refuses")
  writeLines(c("PASSED=TRUE", sprintf("script_md5=%s", fake_md5), "n_reps=2000", "n_scenarios=3"), gatef)
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(good), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("all scenarios", e), "calibration on a subset of scenarios -> refuses")
  writeLines(c("PASSED=TRUE", sprintf("script_md5=%s", fake_md5), "n_reps=2000", nsc), gatef)
  # (c) failing positive control stops the script, on a wrong count and on a wrong estimate.
  bad1 <- good; bad1$n <- 150L
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(bad1), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("POSITIVE CONTROL FAILED", e), sprintf("wrong instrument count stops the script [%s]", substr(e, 1, 200)))
  bad2 <- good; bad2$b["LDL_C"] <- bad2$b["LDL_C"] + 0.05
  e <- tryCatch({ run_real(opt2, gatef, c(.85, 1), mkarm(bad2), fake_md5); "ran" }, error = function(e) conditionMessage(e))
  ok(grepl("POSITIVE CONTROL FAILED", e), "wrong control estimate stops the script")
  # (d) a clean run, end to end.
  res <- run_real(opt2, gatef, c(.85, 1), mkarm(good), fake_md5)
  ok(file.exists(file.path(outd, "robust_contrast_sets.csv")) && file.exists(file.path(outd, "robust_contrast_verdict.txt")), "pipeline writes sets and verdict")
  fg <- res$ests[outcome == "FINNGEN"]
  ok(fg$n_instruments == 150L, sprintf("FinnGen harmonisation drops only the ambiguous palindrome (%d of 151)", fg$n_instruments))
  p_truth <- mk_prob(fx$data$bx[-7, ], fx$data$se[-7, ], fx$data$by[-7], fx$data$sy[-7], fx$data$R)
  ok(abs(fg$ivw_contrast - ivw_fit(p_truth)$th[1]) < 1e-6, "allele flips, multi-rsID cells and multi-allelic rows are harmonised exactly")
  ok(all(res$sets[outcome == "UKB_CAD", overlap_flag] == "UKB_ON_UKB_UNCORRECTED"), "UK Biobank rows carry the overlap flag")
  ok(!is.na(res$verdict$label), sprintf("verdict is labelled: %s", res$verdict$label))
  ll <- fread(file.path(outd, "robust_contrast_leave_locus.csv"))[statistic == "KOH_subset"]
  ok(nrow(ll) == 3L && all(ll$status == "OK") && identical(ll[outcome == "FINNGEN"][order(locus_removed), n_removed], c(1L, 1L, 2L)),
     sprintf("leave-locus refits run on the primary outcome and drop the planted LPA and APOE instruments (%d rows)", nrow(ll)))
  ok(any(grepl("without the LPA and APOE loci|LOCUS-SENSITIVE", res$lines)), "verdict file carries the leave-locus line")
  # (e) a secondary arm that is absent is skipped loudly, not guessed at; one with a failed control gives NA.
  arm_abs <- mkarm(good); arm_abs$statinnaive <- list(name = "statinnaive", traits = c("LDL_C", "APOB", "TRIGLY"),
    files = function(t) file.path(fx$naive, sprintf("lipids_%s.txt.gz", t)), required = FALSE, pheno = fx$pheno, keep = NULL,
    outcomes = list(), control = good)
  arm_bad <- mkarm(good)$biochem; arm_bad$name <- "nmr"; arm_bad$required <- FALSE
  arm_bad$control$n <- 62L
  arm_abs$nmr <- arm_bad
  opt3 <- opt2; opt3$out_dir <- file.path(outd, "o3"); opt3$cache_dir <- file.path(outd, "cache3")
  res3 <- run_real(opt3, gatef, c(.85, 1), arm_abs, fake_md5)
  ok(any(grepl("statinnaive skipped", res3$lines)), "absent arm is skipped with a message")
  ok(all(res3$ests[arm == "nmr", is.na(qhet_contrast)]), "arm with a failed control gives NA and no number")
  ok(any(grepl("POSITIVE CONTROL FAILED in arm nmr", res3$lines)), "failed secondary control is in the verdict file")
  ok(!anyNA(res$ests$mvmr_pkg_qhet_contrast), "the installed MVMR::qhet_mvmr runs through its real API")
  vt <- readLines(file.path(outd, "robust_contrast_verdict.txt"))
  ok(any(grepl("Primary, KOH_subset", vt)) && any(grepl("conditional F of the contrast", vt)), "verdict file carries the other sets and the instrument strength")
  # (e2) the arm configuration names the files the task defines, and the NMR route runs on them.
  fn <- write_nmr_fixture(file.path(o$selftest_dir, "fx_nmr"))
  cfg <- arm_config(list(repo = "/x", sumstat_dir = fx$sumstat, statinnaive_dir = fx$naive, nmr_dir = fn$dir, pheno = fx$pheno,
                         pheno_statinnaive = fx$pheno, keep_statinnaive = "", pheno_nmr = fn$pheno, finngen = fx$finngen, cardiogram = fn$cardiogram))
  ok(identical(basename(cfg$nmr$files("LDL_P")[22]), "nmr_chr22_LDL_P.regenie.gz") && length(cfg$nmr$files("LDL_P")) == 22L, "NMR reselection route (selftest only) is configured for 22 per-chromosome files")
  ok(identical(basename(cfg$biochem$files("APOB")), "lipids_APOB.txt.gz") && identical(cfg$biochem$control$n, 151L) && identical(cfg$statinnaive$control$n, 148L) && identical(cfg$nmr$control$n, 62L), "arm controls expect 151, 148 and 62 instruments")
  ok(abs(exp(cfg$biochem$control$b[["LDL_C"]]) - 1.81) < 0.005 && abs(exp(cfg$biochem$control$b[["APOB"]]) - 0.72) < 0.005, "biochemistry control is the published OR 1.81 and 0.72")
  ok(cfg$biochem$outcomes[[1]]$id == PRIMARY_OUT, "FinnGen is the first outcome of the primary arm")
  an <- cfg["nmr"]; names(an) <- "biochem"; an$biochem$name <- "biochem"; an$biochem$required <- TRUE; an$biochem$outcomes <- an$biochem$outcomes[1]
  an$biochem$precomputed <- NULL  # this block tests the reselection route; (e4) tests the published-file route
  nb <- fn$data; keep6 <- setdiff(seq_len(62L), 5:7)
  f6 <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = nb$bx[keep6, ], bxse = nb$se[keep6, ], by = nb$by[keep6], byse = nb$sy[keep6]))
  an$biochem$control <- list(outcome = "CARDIOGRAM", n = 59L, b = c(LDL_C_NMR = unname(f6@Estimate[1]), LDL_P = unname(f6@Estimate[2])), rho = NA_real_, condF = numeric(0))
  optn <- opt2; optn$out_dir <- file.path(outd, "on"); optn$cache_dir <- file.path(outd, "cachen")
  # This control path takes no FinnGen outcome, so the primary row is absent and the label is INVALID.
  resn <- run_real(optn, gatef, c(.85, 1), an, fake_md5)
  en <- resn$ests[outcome == "CARDIOGRAM"]
  ok(en$n_instruments == 59L, "per-chromosome exposures and CARDIoGRAM format give 59 of 62 instruments, palindromes dropped")
  ok(abs(en$ivw_contrast - (f6@Estimate[1] - f6@Estimate[2])) < 1e-6, "CARDIoGRAM allele flips are harmonised exactly")
  ok(resn$verdict$label == "INVALID", "no FinnGen row in the primary arm cannot produce a verdict")
  # (e4) the NMR arm's real route: instruments and harmonised effects fixed in the published file,
  # CARDIoGRAM only, nothing reselected and nothing dropped.
  pcf <- file.path(fn$dir, "nmr_mvmr_snps.csv")
  fwrite(data.table(SNP = sprintf("rs%d", seq_len(62L)), bxp = nb$bx[, 2], sxp = nb$se[, 2], bxc = nb$bx[, 1], sxc = nb$se[, 1],
                    by = nb$by, sy = nb$sy), pcf)
  cfgp <- arm_config(list(repo = "/x", sumstat_dir = fx$sumstat, statinnaive_dir = fx$naive, nmr_dir = fn$dir, pheno = fx$pheno,
                          pheno_statinnaive = fx$pheno, keep_statinnaive = "", pheno_nmr = fn$pheno, finngen = fx$finngen,
                          cardiogram = fn$cardiogram, nmr_instruments = pcf))
  ok(identical(cfgp$nmr$precomputed$file, pcf) && identical(cfgp$nmr$precomputed$outcome$id, "CARDIOGRAM") &&
     abs(exp(cfgp$nmr$control$b[["LDL_P"]]) - 2.47578) < 1e-5 && abs(exp(cfgp$nmr$control$b[["LDL_C_NMR"]]) - 0.686485) < 1e-6,
     "NMR arm reads the published instrument file, CARDIoGRAM only, control at the export's full precision")
  ap <- cfgp["nmr"]; names(ap) <- "biochem"; ap$biochem$name <- "biochem"; ap$biochem$required <- TRUE
  f62 <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = nb$bx, bxse = nb$se, by = nb$by, byse = nb$sy))
  ap$biochem$control <- list(outcome = "CARDIOGRAM", n = 62L, b = c(LDL_C_NMR = unname(f62@Estimate[1]), LDL_P = unname(f62@Estimate[2])),
                             rho = NA_real_, condF = numeric(0))
  optp <- opt2; optp$out_dir <- file.path(outd, "op"); optp$cache_dir <- file.path(outd, "cachep")
  resp <- run_real(optp, gatef, c(.85, 1), ap, fake_md5)
  ep <- resp$ests[outcome == "CARDIOGRAM"]
  ok(nrow(ep) == 1L && ep$status == "OK" && ep$n_instruments == 62L, "published-file route keeps all 62 instruments: nothing reselected, no palindrome dropped")
  ok(abs(ep$ivw_contrast - (f62@Estimate[1] - f62@Estimate[2])) < 1e-6, "published-file route reproduces mr_mvivw on the file exactly")
  ok(!"FINNGEN" %in% resp$ests$outcome, "published-file route has no FinnGen outcome: the file carries no alleles")
  # (e3) changed inputs must not be answered from the cache.
  cd <- file.path(outd, "cachetest"); dir.create(cd)
  fa <- file.path(cd, "g.txt"); base <- file.path(cd, "x")
  hd <- "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA"
  writeLines(c(hd, "1 10 rsA A G 0.2 9 ADD 0.1 0.01 100 20 NA"), fa)
  a1 <- extract_gws(fa, base)
  writeLines(c(hd, "1 10 rsB A G 0.2 9 ADD 0.1 0.01 100 20 NA", "1 99 rsC A G 0.2 9 ADD 0.1 0.01 100 20 NA"), fa)
  a2 <- extract_gws(fa, base)
  ok(identical(a1$ID, "rsA") && identical(a2$ID, c("rsB", "rsC")), "a changed exposure file is not served from the cache")
  b1 <- extract_ids(fa, "rsB", base); b2 <- extract_ids(fa, "rsC", base)
  ok(identical(b1$ID, "rsB") && identical(b2$ID, "rsC"), "a changed instrument list is not served from the cache")
  # (f) the selftest's own partial calibration can never open the gate of a real run.
  ok(!read_gate(file.path(dirc, "calibration_status.txt"), script_md5())$ok, "the selftest calibration does not open the gate")
  message("[selftest] PASSED")
}

# ------------------------------------------------------------
# Entry
# ------------------------------------------------------------
main <- function() {
  o <- parse_args(OptionParser(option_list = list(
    make_option("--selftest", action = "store_true", default = FALSE),
    make_option("--calibrate", action = "store_true", default = FALSE),
    make_option("--repo", type = "character", default = "/path/to/cad-genetics"),
    make_option("--out_dir", type = "character", default = ""),
    make_option("--cache_dir", type = "character", default = ""),
    make_option("--calibration_dir", type = "character", default = ""),
    make_option("--arms", type = "character", default = "biochem,statinnaive,nmr"),
    make_option("--sumstat_dir", type = "character", default = ""),
    make_option("--statinnaive_dir", type = "character", default = ""),
    make_option("--nmr_dir", type = "character", default = ""),
    make_option("--pheno", type = "character", default = ""),
    make_option("--pheno_statinnaive", type = "character", default = ""),
    make_option("--keep_statinnaive", type = "character", default = ""),
    make_option("--pheno_nmr", type = "character", default = ""),
    make_option("--nmr_instruments", type = "character", default = "", help = "harmonised NMR instruments behind Table oS3b"),
    make_option("--finngen", type = "character", default = ""),
    make_option("--cardiogram", type = "character", default = ""),
    make_option("--n_reps", type = "integer", default = 2000L),
    make_option("--plasmode_reps", type = "integer", default = 1000L),
    make_option("--nboot", type = "integer", default = 2000L),
    make_option("--grid", type = "integer", default = 2001L),
    make_option("--gamma_min", type = "double", default = 0.05, help = "minimum distortion, default as in mvmr.weakiv"),
    make_option("--ncores", type = "integer", default = 15L),
    make_option("--selftest_dir", type = "character", default = ""),
    make_option("--selftest_reps", type = "integer", default = 1000L),
    make_option("--mvmr049", type = "character", default = "")
  )))
  R <- file.path(o$repo, "results/agent1_genetics")
  G <- file.path(R, "gwas_summary_stats"); P <- file.path(R, "phenotypes")
  dflt <- function(v, x) if (nzchar(v)) v else x
  o$out_dir <- dflt(o$out_dir, file.path(o$repo, "exports/agent1/nature_round/robust_contrast"))
  o$calibration_dir <- dflt(o$calibration_dir, o$out_dir)
  o$sumstat_dir <- dflt(o$sumstat_dir, file.path(G, "merged_statinfix"))
  o$statinnaive_dir <- dflt(o$statinnaive_dir, file.path(G, "merged_statinnaive"))
  o$nmr_dir <- dflt(o$nmr_dir, file.path(G, "nmr"))
  o$pheno <- dflt(o$pheno, file.path(P, "pheno_quantitative_statinfix.txt"))
  # The statin-naive GWAS used the unrecoded phenotype and a keep list (slurm/08g).
  o$pheno_statinnaive <- dflt(o$pheno_statinnaive, file.path(P, "pheno_quantitative.txt"))
  o$keep_statinnaive <- dflt(o$keep_statinnaive, file.path(P, "regenie_keep_statinnaive.txt"))
  o$pheno_nmr <- dflt(o$pheno_nmr, file.path(P, "pheno_nmr.txt"))
  o$nmr_instruments <- dflt(o$nmr_instruments, file.path(o$repo, "exports/agent1/nmr_gencov/nmr_mvmr_snps.csv"))
  o$finngen <- dflt(o$finngen, file.path(o$repo, "data/external/finngen/finngen_R11_I9_CHD.gz"))
  o$cardiogram <- dflt(o$cardiogram, file.path(o$repo, "data/external/cardiogram/cad.add.160614.website.txt"))

  if (o$selftest) {
    if (!nzchar(o$selftest_dir)) o$selftest_dir <- tempfile("robust_selftest")
    dir.create(o$selftest_dir, recursive = TRUE, showWarnings = FALSE)
    selftest(o)
  } else if (o$calibrate) {
    r <- run_calibration(o$calibration_dir, o$n_reps, o$ncores)
    if (!r$passed) quit(status = 1)
  } else {
    run_real(o, file.path(o$calibration_dir, "calibration_status.txt"))
  }
}
if (sys.nframe() == 0L) main()

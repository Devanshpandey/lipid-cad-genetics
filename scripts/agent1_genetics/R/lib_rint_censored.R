# Rank-inverse-normal transform that keeps participants censored outside an assay's reportable range.
#
# Port of rint_censored() from Causal Atlas (atlas/phenotypes.py), so the two projects score Lp(a) identically.
# Measured values get qnorm((n_low + r - 0.5) / n), r their average rank and n = measured + censored. A censored
# block covering the bottom (top) fraction p of the ranking gets E[Z | Z in that quantile interval], -dnorm(a)/p
# with a = qnorm(p) (dnorm(b)/p with b = qnorm(p, lower.tail = FALSE)), so the scores keep mean 0.
# Returns list(z, expected_sd): the SD the scores should have, 1 less each block's within-block variance times p.
# With no censoring this is the plain rankit RINT and expected_sd is 1.
#
# Self-test: Rscript lib_rint_censored.R --selftest
rint_censored <- function(x, low = rep(FALSE, length(x)), high = rep(FALSE, length(x))) {
  stopifnot(length(x) == length(low), length(x) == length(high), !anyNA(low), !anyNA(high), !any(low & high))
  meas <- !is.na(x) & !low & !high
  n_lo <- sum(low); n_hi <- sum(high)
  n <- sum(meas) + n_lo + n_hi
  z <- rep(NA_real_, length(x))
  z[meas] <- qnorm((n_lo + rank(x[meas], ties.method = "average") - 0.5) / n)
  v <- 1
  if (n_lo > 0) {
    p <- n_lo / n; a <- qnorm(p); mu <- -dnorm(a) / p
    z[low] <- mu
    v <- v - p * (1 - a * dnorm(a) / p - mu^2)
  }
  if (n_hi > 0) {
    p <- n_hi / n; b <- qnorm(p, lower.tail = FALSE); mu <- dnorm(b) / p
    z[high] <- mu
    v <- v - p * (1 + b * dnorm(b) / p - mu^2)
  }
  list(z = z, expected_sd = sqrt(v))
}

if (sys.nframe() == 0L && "--selftest" %in% commandArgs(TRUE)) {
  set.seed(42)
  x <- round(exp(rnorm(2e5, 2.5, 1.2)), 1)
  low <- x < 3.8; high <- x > 189
  x[low | high] <- NA
  x[sample(which(!low & !high), 1e4)] <- NA       # never measured
  r <- rint_censored(x, low, high)
  stopifnot(abs(mean(r$z, na.rm = TRUE)) < 1e-3,
            abs(sd(r$z, na.rm = TRUE) - r$expected_sd) < 2e-3,
            max(r$z[!low & !high], na.rm = TRUE) < min(r$z[high]),
            min(r$z[!low & !high], na.rm = TRUE) > max(r$z[low]),
            sum(is.na(r$z)) == 1e4, length(unique(r$z[low])) == 1, length(unique(r$z[high])) == 1)
  r0 <- rint_censored(x)
  stopifnot(abs(r0$expected_sd - 1) < 1e-12, abs(mean(r0$z, na.rm = TRUE)) < 1e-3)
  cat("rint_censored selftest passed\n")
}

# Independent re-implementation of the model as coded in
# Ribosome/R/Polysome_functions_RAUC.R, written 2026-10-09 by a Claude Code
# session as an oracle for the test suite (see the handoff note alongside).
#   tau_i   = tau_c * i / (9 imax)
#   kappa_i = kappa (1 - i/imax)
#   m*_i    = mu * sum_{j>=i} m_j / tau_i   (i >= 1);  m*_0 = mu * sum_j m_j / delta
# Base R only: no limSolve.
solve_model <- function(kappa, mu, tau_c, imax, lambda = 1, delta = 1e5) {
  n <- imax + 1
  taui <- tau_c / (9 * imax) * (0:imax)
  kapi <- kappa * seq(1, 0, length.out = n)
  A <- matrix(0, n, n)
  for (i in 1:n) {
    A[i, i] <- -(kapi[i] + taui[i] + mu)
    if (i < n) { A[i + 1, i] <- kapi[i]; A[i, i + 1] <- taui[i + 1] }
  }
  m <- solve(A, c(-lambda, rep(0, imax)))
  S <- rev(cumsum(rev(m)))
  mstar <- c(mu * S[1] / delta, mu * S[-1] / taui[-1])
  list(m = m, mstar = mstar)
}

if (sys.nframe() == 0) {
  imax <- 39; tau_c <- 1
  mus <- c(0.0057762265, 0.0038508177, 0.0023104906, 0.0016503504,
           0.0008251752, 0.0004443251, 0.0003300701, 0.0002221626)
  kappas <- 10^seq(-4, log10(0.5), length.out = 60)
  cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n",
              "mu", "share_min", "share_max", "closed", "maxFull"))
  out <- list()
  for (mu in mus) {
    sh <- sapply(kappas, function(k) {
      s <- solve_model(k, mu, tau_c, imax)
      pU <- sum((0:imax) * s$m); pM <- sum((0:imax) * s$mstar)
      c(pM / (pU + pM), pU + pM)
    })
    out[[as.character(mu)]] <- sh[2, ]
    cat(sprintf("%-10.2e %-10.4f %-10.4f %-10.4f %-10.4f\n", mu,
                min(sh[1, ]), max(sh[1, ]),
                9 * mu * imax / (1 + 9 * mu * imax), max(sh[2, ])))
  }
  gmax <- max(unlist(out))
  cat("\nPeak full production normalised to the global max",
      "(compare to Fig 10 peak heights 1.000, 0.692, 0.143):\n")
  for (mu in mus)
    cat(sprintf("mu=%.2e  peak=%.3f\n", mu, max(out[[as.character(mu)]]) / gmax))
  cat("\nDecapped share with mu=5.7e-3, imax=39, for tau_c = 1, 5, 10 codons/s:\n")
  for (tc in c(1, 5, 10)) {
    r <- 9 * 5.7e-3 * imax / tc
    cat(sprintf("tau_c=%2d  share=%.3f\n", tc, r / (1 + r)))
  }
}

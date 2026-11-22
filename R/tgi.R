# Tumor size (Stein model), overall survival, and
# the phase 2 go/no-go decision. Time is in weeks.

tgi_control <- list(
  ts0_med = 70,
  ts0_sd = 0.5, # mm, log-scale SD
  d_med = 0.02,
  g_med = 0.025, # 1/week
  om_d = 0.6,
  om_g = 0.6,
  meas_sd = 0.12, # log-scale measurement error
  shape = 1.3, # Weibull shape for OS
  b_g = 0.8, # log hazard per unit log growth rate
  b_ts = 0.3 # log hazard per unit log baseline size
)

# drug-arm effects: multipliers on the median
# shrinkage and growth rates, and a hazard ratio
# acting outside the tumor
tgi_scenarios <- tibble::tribble(
  ~scenario , ~label                                   , ~d_mult , ~g_mult , ~os_direct ,
  "S0"      , "no effect"                              , 1       , 1       , 1          ,
  "S1"      , "shrinks tumors, growth unchanged"       , 5       , 1       , 1          ,
  "S2"      , "slows regrowth"                         , 1.5     , 0.6     , 1          ,
  "S3"      , "shrinks and slows regrowth"             , 5       , 0.6     , 1          ,
  "S4"      , "part of the benefit outside tumor size" , 5       , 0.85    , 0.8
)

stein_f <- function(t, d, g) exp(-d * t) + exp(g * t) - 1

# Weibull scale chosen so median control OS is 52
# weeks for a typical patient
os_lambda <- function(p = tgi_control) log(2) / 52^p$shape

sim_tgi_patients <- function(n, arm, sc, p = tgi_control) {
  drug <- arm == "drug"
  ts0 <- pmin(pmax(rlnorm(n, log(p$ts0_med), p$ts0_sd), 20), 250)
  d <- p$d_med * ifelse(drug, sc$d_mult, 1) * exp(rnorm(n, 0, p$om_d))
  g <- p$g_med * ifelse(drug, sc$g_mult, 1) * exp(rnorm(n, 0, p$om_g))
  lp <- p$b_g *
    log(g / p$g_med) +
    p$b_ts * log(ts0 / p$ts0_med) +
    ifelse(drug, log(sc$os_direct), 0)
  os <- (-log(runif(n)) / (os_lambda(p) * exp(lp)))^(1 / p$shape)
  tibble::tibble(arm = arm, ts0 = ts0, d = d, g = g, os = os)
}

# Scans every `interval` weeks from randomization
# until progression, death or data cutoff.
# Progression: 20% and 5 mm above the smallest
# size so far (RECIST-like). stop_at_pd = FALSE
# keeps scanning after progression (used only by
# the post hoc check in
# scripts/09_posthoc_scans.R).
scan_patients <- function(
  pts,
  interval,
  follow,
  p = tgi_control,
  stop_at_pd = TRUE
) {
  n <- nrow(pts)
  nmax <- floor(max(follow) / interval)
  tt <- interval * (0:nmax)
  truth <- outer(pts$ts0, rep(1, nmax + 1)) *
    stein_f(outer(rep(1, n), tt), pts$d, pts$g)
  obs <- truth * exp(matrix(rnorm(n * (nmax + 1), 0, p$meas_sd), n))
  seen <- outer(pmin(pts$os, follow), tt, ">=")
  m <- integer(n)
  pfs_event <- logical(n)
  best <- rep(0, n)
  for (i in seq_len(n)) {
    k <- max(which(seen[i, ]))
    nadir <- obs[i, 1]
    last <- 1
    if (k > 1) {
      for (j in 2:k) {
        last <- j
        ch <- obs[i, j] / obs[i, 1] - 1
        if (
          !pfs_event[i] && obs[i, j] >= 1.2 * nadir && obs[i, j] - nadir >= 5
        ) {
          pfs_event[i] <- TRUE
          if (stop_at_pd) break
        }
        if (!pfs_event[i]) {
          best[i] <- min(best[i], ch)
          nadir <- min(nadir, obs[i, j])
        }
      }
    }
    m[i] <- last - 1
  }
  list(
    y = log(obs),
    m = m,
    times = tt,
    responder = best <= -0.3,
    pfs_event = pfs_event
  )
}

sim_trial <- function(
  n_arm,
  interval,
  sc,
  p = tgi_control,
  accrual = 26,
  cutoff = 50,
  stop_at_pd = TRUE
) {
  pts <- dplyr::bind_rows(
    sim_tgi_patients(n_arm, "control", sc, p),
    sim_tgi_patients(n_arm, "drug", sc, p)
  )
  entry <- runif(2 * n_arm, 0, accrual)
  follow <- cutoff - entry
  pts <- dplyr::mutate(
    pts,
    follow = follow,
    os_time = pmin(os, follow),
    os_event = os <= follow
  )
  list(pts = pts, scans = scan_patients(pts, interval, follow, p, stop_at_pd))
}

# Grid of (log d, log g); log baseline size is
# profiled out by least squares, so a patient's
# fit needs no start values.
make_grid <- function(nd = 40, ng = 40) {
  expand.grid(
    ld = seq(log(0.002), log(0.5), length.out = nd),
    lg = seq(log(0.002), log(0.2), length.out = ng)
  )
}

# log-likelihood of each patient (rows) at each
# grid point (columns); patients with no scan
# after baseline get NA
grid_loglik <- function(y, m, times, grid, meas_sd) {
  n <- nrow(y)
  out <- matrix(NA_real_, n, nrow(grid))
  for (k in setdiff(sort(unique(m)), 0)) {
    rows <- which(m == k)
    tt <- times[seq_len(k + 1)]
    f <- log(stein_f(
      outer(rep(1, nrow(grid)), tt),
      exp(grid$ld),
      exp(grid$lg)
    ))
    yy <- y[rows, seq_len(k + 1), drop = FALSE]
    nk <- k + 1
    sse <- rowSums(yy^2) -
      2 * yy %*% t(f) +
      outer(rep(1, length(rows)), rowSums(f^2)) -
      nk * outer(rowMeans(yy), rowMeans(f), "-")^2
    out[rows, ] <- -sse / (2 * meas_sd^2)
  }
  out
}

# Marginal maximum likelihood for the arm's mean
# (log d, log g), with the spread fixed from the
# historical fit and the normal prior discretized
# on the grid. Returns the mean and each patient's
# posterior mean of log g.
grid_em <- function(ll, grid, om_d, om_g, start) {
  use <- !is.na(ll[, 1])
  ll <- ll[use, , drop = FALSE]
  shift <- apply(ll, 1, max)
  el <- exp(ll - shift)
  x <- cbind(grid$ld, grid$lg)
  om2 <- c(om_d, om_g)^2
  post <- function(mu) {
    lp <- -((x[, 1] - mu[1])^2 / om2[1] + (x[, 2] - mu[2])^2 / om2[2]) / 2
    pri <- exp(lp - max(lp))
    pri <- pri / sum(pri)
    a <- el %*% pri
    list(pri = pri, a = as.vector(a))
  }
  nll <- function(mu) -sum(log(post(mu)$a))
  grad <- function(mu) {
    p <- post(mu)
    w <- sweep(el, 2, p$pri, "*") / p$a
    m_post <- colSums(w %*% x)
    m_pri <- colSums(p$pri * x)
    -(m_post - nrow(el) * m_pri) / om2
  }
  lo <- c(min(x[, 1]), min(x[, 2]))
  hi <- c(max(x[, 1]), max(x[, 2]))
  opt <- stats::optim(
    start,
    nll,
    grad,
    method = "L-BFGS-B",
    lower = lo,
    upper = hi
  )
  p <- post(opt$par)
  w <- sweep(el, 2, p$pri, "*") / p$a
  post_lg <- rep(NA_real_, length(use))
  post_lg[use] <- as.vector(w %*% x[, 2])
  list(
    mu = opt$par,
    post_lg = post_lg,
    iter = opt$counts[[1]],
    n_used = sum(use)
  )
}

# Everything the decision rules need from one
# trial
trial_metrics <- function(tr, grid, hist) {
  pts <- tr$pts
  sc <- tr$scans
  ll <- grid_loglik(sc$y, sc$m, sc$times, grid, hist$meas_sd)
  start <- c(log(tgi_control$d_med), log(tgi_control$g_med))
  em <- lapply(c("control", "drug"), function(a) {
    i <- pts$arm == a
    grid_em(ll[i, , drop = FALSE], grid, hist$om_d, hist$om_g, start)
  })
  d_lg <- em[[2]]$mu[2] - em[[1]]$mu[2]
  cox <- survival::coxph(
    survival::Surv(os_time, os_event) ~ I(arm == "drug"),
    data = pts
  )
  tibble::tibble(
    orr_ctrl = mean(sc$responder[pts$arm == "control"]),
    orr_drug = mean(sc$responder[pts$arm == "drug"]),
    d_lg = d_lg,
    pred_hr = exp(hist$b_g * d_lg),
    os_hr = exp(unname(stats::coef(cox))),
    os_events = sum(pts$os_event)
  )
}

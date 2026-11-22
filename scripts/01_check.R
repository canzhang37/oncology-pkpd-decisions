# Rscript scripts/01_check.R
# Code checks that must pass before any result is
# read.
source("R/pk.R")
source("R/er.R")
source("R/tgi.R")
set.seed(7)
log_lines <- character()
report <- function(name, value, ok) {
  line <- sprintf("%-46s %-28s %s", name, value, if (ok) "PASS" else "FAIL")
  cat(line, "\n")
  log_lines <<- c(log_lines, line)
}

# (a) closed-form concentrations against rxode2
pk <- indiv_pk(sim_patients(5, c(100, 200)))
mod <- rxode2::rxode2({
  ka <- p_ka
  cl <- p_cl
  v <- p_v
  cp <- linCmt()
})
worst <- 0
for (i in seq_len(nrow(pk))) {
  ev <- rxode2::et(amt = pk$dose[i], ii = 24, addl = 28) |>
    rxode2::et(pk_times)
  s <- suppressMessages(rxode2::rxSolve(
    mod,
    params = c(p_ka = pk$ka[i], p_cl = pk$cl[i], p_v = pk$v[i]),
    events = ev
  ))
  ours <- conc_oral(s$time, pk$dose[i], pk$ka[i], pk$cl[i], pk$v[i])
  worst <- max(worst, max(abs(ours / s$cp - 1)))
}
report(
  "(a) PK closed form vs rxode2, max rel. diff",
  signif(worst, 2),
  worst < 1e-6
)

# (a2) steady-state records in the fitted data
# sets
pk <- indiv_pk(dplyr::bind_rows(
  sim_patients(500, 100),
  sim_patients(500, 200, 501)
))
obs <- sim_pk_obs(pk)
ev <- nm_data(pk, obs, "half")[, c("ID", "TIME", "EVID", "AMT", "SS", "II")]
names(ev) <- tolower(names(ev))
prm <- data.frame(id = pk$id, p_ka = pk$ka, p_cl = pk$cl, p_v = pk$v)
s <- as.data.frame(rxode2::rxSolve(mod, params = prm, events = ev))
full <- unlist(lapply(seq_len(nrow(pk)), function(i) {
  conc_oral(pk_times, pk$dose[i], pk$ka[i], pk$cl[i], pk$v[i])
}))
worst <- max(abs(s$cp / full - 1))
report(
  "(a2) steady-state records vs all doses, max rel.",
  signif(worst, 2),
  worst < 0.005
)

# (b) grid fit finds noise-free Stein curves
grid <- make_grid()
tt <- 8 * (0:6)
d0 <- c(0.01, 0.05, 0.2)
g0 <- c(0.01, 0.03, 0.08)
y <- t(mapply(function(d, g) log(60 * stein_f(tt, d, g)), d0, g0))
ll <- grid_loglik(y, rep(6, 3), tt, grid, 0.12)
best <- grid[apply(ll, 1, which.max), ]
step <- max(diff(unique(grid$ld))[1], diff(unique(grid$lg))[1])
err <- max(abs(best$ld - log(d0)), abs(best$lg - log(g0)))
report("(b) grid argmax vs truth, max |log error|", round(err, 3), err <= step)

# (c) grid EM recovers an arm's mean log growth
# rate
sc0 <- tgi_scenarios[1, ]
pts <- sim_tgi_patients(3000, "control", sc0)
sc <- scan_patients(pts, 6, rep(60, nrow(pts)))
ll <- grid_loglik(sc$y, sc$m, sc$times, grid, tgi_control$meas_sd)
em <- grid_em(
  ll,
  grid,
  tgi_control$om_d,
  tgi_control$om_g,
  start = c(log(0.05), log(0.05))
)
err <- abs(em$mu[2] - log(tgi_control$g_med))
report("(c) grid EM mean log g, |error|", round(err, 3), err < 0.05)

# (d) OS simulation recovers the hazard
# coefficients
pts <- sim_tgi_patients(50000, "control", sc0)
cox <- survival::coxph(survival::Surv(os) ~ log(g) + log(ts0), data = pts)
err <- max(abs(coef(cox) - c(tgi_control$b_g, tgi_control$b_ts)))
report("(d) Cox on simulated OS, max |coef error|", round(err, 3), err < 0.03)

# (e) BLQ records: M3 keeps every sample, drop
# removes BLQ
pk <- indiv_pk(sim_patients(40, 100))
obs <- sim_pk_obs(pk)
m3 <- nm_data(pk, obs, "m3")
dr <- nm_data(pk, obs, "drop")
ok <- sum(m3$CENS) == sum(obs$blq) &&
  sum(dr$EVID == 0) == sum(!obs$blq) &&
  all(m3$DV[m3$CENS == 1] == pk_truth$lloq)
report("(e) BLQ records (M3 censored, dropped)", sum(obs$blq), ok)

# (f) exposure-response gains in a very large
# trial
tr <- er_truth$small
d <- sim_er_trial(50000, tr)
g <- er_gains(d)
tg <- true_gain(tr, 2e5)
report(
  "(f) adjusted ER gain vs truth, |error|",
  round(abs(g[["adjusted"]] - tg), 4),
  abs(g[["adjusted"]] - tg) < 0.005
)
report(
  "(f) by-arm gain vs truth, |error|",
  round(abs(g[["by_arm"]] - tg), 4),
  abs(g[["by_arm"]] - tg) < 0.01
)

writeLines(log_lines, "results/check_log.txt")

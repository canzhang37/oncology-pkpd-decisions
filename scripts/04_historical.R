# Rscript scripts/04_historical.R
# B0: fit the tumor model to a historical study
# with nlmixr2, then estimate how tumor regrowth
# links to survival. The decision rule in
# 05_decisions.R uses these.
suppressMessages(library(nlmixr2))
source("R/tgi.R")
source("R/tgi_fit.R")
set.seed(4401)

earlier_drug <- tibble::tibble(
  scenario = "H",
  label = "earlier drug",
  d_mult = 3,
  g_mult = 0.7,
  os_direct = 1
)
hist_trial <- sim_trial(300, 8, earlier_drug, accrual = 26, cutoff = 130)
dat <- tgi_long(hist_trial)
cat("patients", nrow(hist_trial$pts), "scans", nrow(dat), "\n")

t0 <- Sys.time()
fit <- fit_tgi(dat, "saem")
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
th <- fit$theta
om <- sqrt(diag(fit$omega))
truth <- c(
  ld = log(tgi_control$d_med),
  lg = log(tgi_control$g_med),
  ld_drug = log(3),
  lg_drug = log(0.7),
  lsd = tgi_control$meas_sd,
  om_d = tgi_control$om_d,
  om_g = tgi_control$om_g
)
est <- c(
  th[c("ld", "lg", "ld_drug", "lg_drug", "lsd")],
  om_d = om[["eta_d"]],
  om_g = om[["eta_g"]]
)

# per-patient posterior mean log growth rate on
# the grid, with the spread and noise from the
# population fit
grid <- make_grid()
sc <- hist_trial$scans
pts <- hist_trial$pts
ll <- grid_loglik(sc$y, sc$m, sc$times, grid, est[["lsd"]])
pts$post_lg <- NA_real_
for (a in c("control", "drug")) {
  i <- pts$arm == a
  em <- grid_em(
    ll[i, , drop = FALSE],
    grid,
    est[["om_d"]],
    est[["om_g"]],
    start = c(log(0.03), log(0.03))
  )
  pts$post_lg[i] <- em$post_lg
}
pts$log_ts0 <- sc$y[, 1]
eb <- fit$eta
lg_ebe <- th[["lg"]] + th[["lg_drug"]] * (pts$arm == "drug") + eb$eta_g
ok <- !is.na(pts$post_lg)
r_ebe <- cor(lg_ebe[ok], pts$post_lg[ok])

cox <- survival::coxph(
  survival::Surv(os_time, os_event) ~ post_lg + log_ts0,
  data = pts
)
link <- coef(cox)
se <- sqrt(diag(vcov(cox)))

out <- tibble::tibble(
  parameter = c(names(est), "b_g", "b_ts"),
  truth = c(truth, tgi_control$b_g, tgi_control$b_ts),
  estimate = c(est, link[["post_lg"]], link[["log_ts0"]]),
  se = c(rep(NA, length(est)), se[["post_lg"]], se[["log_ts0"]])
)
readr::write_csv(out, "results/historical_fit.csv")
print(out)
cat(sprintf(
  "SAEM fit %.0f s; patients with a scan after baseline %d of %d; OS events %d; cor(EBE, grid posterior) %.3f\n",
  secs,
  sum(!is.na(pts$post_lg)),
  nrow(pts),
  sum(pts$os_event),
  r_ebe
))
writeLines(
  sprintf(
    "SAEM fit %.0f s, %d patients, %d scans, %d deaths, cor(EBE eta_g, grid posterior mean log g) %.3f",
    secs,
    nrow(pts),
    nrow(dat),
    sum(pts$os_event),
    r_ebe
  ),
  "results/historical_fit_log.txt"
)

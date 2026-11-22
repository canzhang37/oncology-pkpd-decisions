# Rscript scripts/05_decisions.R
# Part B: phase 2 trials in every scenario and
# design cell, with the three decision metrics for
# each trial.
source("R/tgi.R")
hf <- readr::read_csv("results/historical_fit.csv", show_col_types = FALSE)
val <- function(p) hf$estimate[hf$parameter == p]
hist <- list(
  meas_sd = val("lsd"),
  om_d = val("om_d"),
  om_g = val("om_g"),
  b_g = val("b_g")
)

reps <- 1000
cells <- tidyr::expand_grid(
  scenario = tgi_scenarios$scenario,
  n_arm = c(40, 60, 90),
  interval = c(6, 8, 12)
)
cells$seed <- 50000 + seq_len(nrow(cells))

run_cell <- function(k) {
  source("R/tgi.R")
  cl <- cells[k, ]
  set.seed(cl$seed)
  sc <- tgi_scenarios[tgi_scenarios$scenario == cl$scenario, ]
  grid <- make_grid()
  res <- lapply(seq_len(reps), function(r) {
    trial_metrics(sim_trial(cl$n_arm, cl$interval, sc), grid, hist)
  })
  dplyr::bind_rows(res) |>
    dplyr::mutate(
      scenario = cl$scenario,
      n_arm = cl$n_arm,
      interval = cl$interval,
      rep = seq_len(reps),
      .before = 1
    )
}

workers <- as.integer(Sys.getenv("PKPD_WORKERS", "8"))
cl <- parallel::makePSOCKcluster(workers)
wd <- getwd()
parallel::clusterExport(cl, c("cells", "hist", "reps", "wd"))
parallel::clusterEvalQ(cl, setwd(wd))
t0 <- Sys.time()
out <- parallel::parLapplyLB(cl, seq_len(nrow(cells)), run_cell)
parallel::stopCluster(cl)
out <- dplyr::bind_rows(out)
readr::write_csv(out, "results/decision_trials.csv.gz")
cat(
  nrow(out),
  "trials in",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
  "min\n"
)

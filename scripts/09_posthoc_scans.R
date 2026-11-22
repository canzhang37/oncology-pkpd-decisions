# Rscript scripts/09_posthoc_scans.R
# Post hoc, after P8 was not met: why did wider
# scan intervals raise the model rule's go rate in
# S2? Compares scans stopped at progression (as
# designed) with scans continued until death or
# cutoff, at 60 patients per arm.
source("R/tgi.R")
hf <- readr::read_csv("results/historical_fit.csv", show_col_types = FALSE)
val <- function(p) hf$estimate[hf$parameter == p]
hist <- list(
  meas_sd = val("lsd"),
  om_d = val("om_d"),
  om_g = val("om_g"),
  b_g = val("b_g")
)
reps <- 500
cells <- tidyr::expand_grid(
  scenario = c("S0", "S2"),
  interval = c(6, 8, 12),
  stop_at_pd = c(TRUE, FALSE)
)
cells$seed <- 90000 + seq_len(nrow(cells))

run_cell <- function(k) {
  source("R/tgi.R")
  cl <- cells[k, ]
  set.seed(cl$seed)
  sc <- tgi_scenarios[tgi_scenarios$scenario == cl$scenario, ]
  grid <- make_grid()
  res <- lapply(seq_len(reps), function(r) {
    tr <- sim_trial(60, cl$interval, sc, stop_at_pd = cl$stop_at_pd)
    ctrl <- tr$pts$arm == "control"
    dplyr::mutate(
      trial_metrics(tr, grid, hist),
      scans_per_patient = mean(tr$scans$m),
      pd_first_scan = mean(tr$scans$pfs_event[ctrl] & tr$scans$m[ctrl] == 1)
    )
  })
  dplyr::bind_rows(res) |>
    dplyr::mutate(
      scenario = cl$scenario,
      interval = cl$interval,
      stop_at_pd = cl$stop_at_pd,
      .before = 1
    )
}

workers <- as.integer(Sys.getenv("PKPD_WORKERS", "6"))
cl <- parallel::makePSOCKcluster(workers)
wd <- getwd()
parallel::clusterExport(cl, c("cells", "hist", "reps", "wd"))
parallel::clusterEvalQ(cl, setwd(wd))
out <- dplyr::bind_rows(parallel::parLapplyLB(
  cl,
  seq_len(nrow(cells)),
  run_cell
))
parallel::stopCluster(cl)

summ <- out |>
  dplyr::group_by(interval, stop_at_pd) |>
  dplyr::summarise(
    cut = quantile(pred_hr[scenario == "S0"], 0.1, type = 1),
    sd_null = sd(d_lg[scenario == "S0"]),
    mean_s2 = mean(d_lg[scenario == "S2"]),
    go_s2 = mean(pred_hr[scenario == "S2"] < cut),
    scans_per_patient = mean(scans_per_patient[scenario == "S0"]),
    pd_first_scan_ctrl = mean(pd_first_scan[scenario == "S0"]),
    .groups = "drop"
  )
readr::write_csv(summ, "results/posthoc_scans.csv")
print(summ, width = Inf)

# Rscript scripts/08_ae_ml.R
# Part C (amendment 1): a new patient's risk of a
# grade 3+ adverse event, from the
# exposure-response model on the population PK fit
# vs XGBoost and logistic regression. Training
# sets are the 20 A1 studies, rebuilt from their
# seeds; the M3 fits are refitted to get the
# empirical Bayes clearances.
source("R/pk.R")
source("R/er.R")
source("R/ae.R")
reps <- as.integer(Sys.getenv("PKPD_PK_REPS", "20"))
workers <- as.integer(Sys.getenv("PKPD_WORKERS", "4"))
test_doses <- c(100, 150, 200, 400)

run_rep <- function(r) {
  suppressMessages(library(nlmixr2))
  source("R/pk.R")
  source("R/er.R")
  source("R/ae.R")
  # same random draws as 03_poppk.R, then the new
  # outcome
  set.seed(30000 + r)
  pk <- indiv_pk(dplyr::bind_rows(
    sim_patients(60, 100),
    sim_patients(60, 200, 61)
  ))
  obs <- sim_pk_obs(pk)
  pk$resp <- rbinom(nrow(pk), 1, resp_prob(pk$cavg, pk$alb, er_truth$small))
  pk <- add_noise_covs(pk)
  pk$ae <- rbinom(nrow(pk), 1, ae_prob(pk$dose / pk$cl))
  f <- fit_pk(nm_data(pk, obs, "m3"))
  if (!f$ok) {
    return(NULL)
  }

  set.seed(80000 + r)
  new <- purrr::map(test_doses, function(dz) {
    sim_patients(5000, dz) |> add_noise_covs()
  }) |>
    dplyr::bind_rows()
  new$oracle <- oracle_risk(new)
  pr <- ae_predictions(
    pk,
    new,
    f$fixef,
    f$omega_sd[["eta_cl"]],
    f$eta_cl[as.character(pk$id)]
  )
  long <- dplyr::bind_cols(dplyr::select(new, dose, oracle), pr$p) |>
    tidyr::pivot_longer(
      c(er_model, logistic, xgboost),
      names_to = "method",
      values_to = "pred"
    ) |>
    dplyr::group_by(dose, method) |>
    dplyr::summarise(
      rmse = sqrt(mean((pred - oracle)^2)),
      mean_pred = mean(pred),
      mean_oracle = mean(oracle),
      .groups = "drop"
    )
  dplyr::mutate(
    long,
    rep = r,
    objf = f$objf,
    n_ae = sum(pk$ae),
    er_slope = pr$er_slope,
    xgb_depth = pr$xgb_depth,
    xgb_rounds = pr$xgb_rounds,
    .before = 1
  )
}

wd <- getwd()
cl <- parallel::makePSOCKcluster(workers)
parallel::clusterExport(cl, c("wd", "test_doses"))
parallel::clusterEvalQ(cl, setwd(wd))
t0 <- Sys.time()
res <- parallel::parLapplyLB(cl, seq_len(reps), run_rep)
parallel::stopCluster(cl)
out <- dplyr::bind_rows(res)
readr::write_csv(out, "results/ae_ml.csv")
cat(sprintf(
  "%d replicates, %.1f min\n",
  length(unique(out$rep)),
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
))

# Rscript scripts/03_poppk.R
# A1: population PK fits (nlmixr2 FOCEi) to
# replicate phase 2 data sets, with BLQ samples
# handled three ways. Also feeds the empirical
# Bayes clearances into the exposure-response
# check.
source("R/pk.R")
source("R/er.R")
reps <- as.integer(Sys.getenv("PKPD_PK_REPS", "20"))
workers <- as.integer(Sys.getenv("PKPD_WORKERS", "4"))

run_rep <- function(r) {
  suppressMessages(library(nlmixr2))
  source("R/pk.R")
  source("R/er.R")
  set.seed(30000 + r)
  pk <- indiv_pk(dplyr::bind_rows(
    sim_patients(60, 100),
    sim_patients(60, 200, 61)
  ))
  obs <- sim_pk_obs(pk)
  pk$resp <- rbinom(nrow(pk), 1, resp_prob(pk$cavg, pk$alb, er_truth$small))
  fits <- list()
  ests <- list()
  for (m in c("m3", "drop", "half")) {
    f <- fit_pk(nm_data(pk, obs, m))
    fits[[m]] <- f
    if (!f$ok) {
      ests[[m]] <- tibble::tibble(rep = r, method = m, ok = FALSE, msg = f$msg)
      next
    }
    ests[[m]] <- tibble::tibble(
      rep = r,
      method = m,
      ok = TRUE,
      msg = NA_character_,
      secs = f$secs,
      objf = f$objf,
      parameter = c(names(f$fixef), paste0("om_", c("ka", "cl", "v"))),
      estimate = c(unname(f$fixef), unname(f$omega_sd))
    )
  }
  # exposure-response with true and with empirical
  # Bayes clearances (pooled doses, with and without
  # albumin)
  er <- list(tibble::tibble(
    rep = r,
    exposure = "true",
    gain_pooled = er_gains(pk)[["pooled"]],
    gain_adjusted = er_gains(pk)[["adjusted"]]
  ))
  for (m in c("m3", "drop")) {
    f <- fits[[m]]
    if (!f$ok) {
      next
    }
    fe <- f$fixef
    cl_ebe <- exp(
      fe[["tcl"]] +
        0.75 * log(pk$wt / 70) +
        fe[["b_alb"]] * log(pk$alb / 40) +
        f$eta_cl[as.character(pk$id)]
    )
    g <- er_gains(pk, cl = cl_ebe)
    er[[m]] <- tibble::tibble(
      rep = r,
      exposure = paste0("ebe_", m),
      gain_pooled = g[["pooled"]],
      gain_adjusted = g[["adjusted"]]
    )
  }
  list(
    est = dplyr::bind_rows(ests),
    er = dplyr::bind_rows(er),
    blq = mean(obs$blq)
  )
}

wd <- getwd()
cl <- parallel::makePSOCKcluster(workers)
parallel::clusterExport(cl, "wd")
parallel::clusterEvalQ(cl, setwd(wd))
t0 <- Sys.time()
res <- parallel::parLapplyLB(cl, seq_len(reps), run_rep)
parallel::stopCluster(cl)

est <- dplyr::bind_rows(purrr::map(res, "est"))
readr::write_csv(est, "results/poppk_fits.csv")
er <- dplyr::bind_rows(purrr::map(res, "er"))
readr::write_csv(er, "results/poppk_er_check.csv")
cat(sprintf(
  "%d replicates, mean BLQ share %.3f, %.1f min\n",
  reps,
  mean(purrr::map_dbl(res, "blq")),
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
))

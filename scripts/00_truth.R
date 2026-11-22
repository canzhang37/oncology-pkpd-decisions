# Rscript scripts/00_truth.R
# Truths from large simulated populations: BLQ
# share, the response gain from 100 to 200 mg, and
# per scenario the response rates, medians and the
# true OS hazard ratio.
source("R/pk.R")
source("R/er.R")
source("R/tgi.R")
set.seed(20261121)

pk <- indiv_pk(dplyr::bind_rows(
  sim_patients(10000, 100),
  sim_patients(10000, 200, 10001)
))
obs <- sim_pk_obs(pk)
blq <- obs |>
  dplyr::left_join(dplyr::select(pk, id, dose), by = "id") |>
  dplyr::group_by(dose, time) |>
  dplyr::summarise(blq = mean(blq), .groups = "drop")
readr::write_csv(blq, "results/truth_blq.csv")
cat("BLQ share overall:", round(mean(obs$blq), 3), "\n")

gains <- purrr::map_dbl(er_truth, true_gain)
er <- tibble::tibble(
  truth = names(er_truth),
  orr_100 = purrr::map_dbl(er_truth, function(tr) {
    p <- indiv_pk(sim_patients(2e5, 100))
    mean(resp_prob(p$cavg, p$alb, tr))
  }),
  gain = gains
)
readr::write_csv(er, "results/truth_er.csv")
print(er)

n <- 20000
tgi <- purrr::map(seq_len(nrow(tgi_scenarios)), function(k) {
  sc <- tgi_scenarios[k, ]
  pts <- dplyr::bind_rows(
    sim_tgi_patients(n, "control", sc),
    sim_tgi_patients(n, "drug", sc)
  )
  scn <- scan_patients(pts, 8, rep(104, 2 * n))
  drug <- pts$arm == "drug"
  cox <- survival::coxph(survival::Surv(os) ~ drug, data = pts)
  tibble::tibble(
    scenario = sc$scenario,
    orr_ctrl = mean(scn$responder[!drug]),
    orr_drug = mean(scn$responder[drug]),
    os_med_ctrl = median(pts$os[!drug]),
    os_med_drug = median(pts$os[drug]),
    true_hr = exp(unname(coef(cox))),
    cond_hr = sc$g_mult^tgi_control$b_g * sc$os_direct,
    pd_by_24wk_ctrl = mean(scn$pfs_event[!drug] & scn$m[!drug] * 8 <= 24)
  )
}) |>
  dplyr::bind_rows()
readr::write_csv(tgi, "results/truth_tgi.csv")
print(tgi, width = Inf)

# Rscript scripts/02_er.R
# A2: exposure-response estimates of the gain from
# 100 to 200 mg, and the dose each estimator would
# choose.
source("R/pk.R")
source("R/er.R")
set.seed(2202)
reps <- 2000
truths <- readr::read_csv("results/truth_er.csv", show_col_types = FALSE)

out <- purrr::map(names(er_truth), function(nm) {
  tr <- er_truth[[nm]]
  purrr::map(seq_len(reps), function(r) {
    g <- er_gains(sim_er_trial(60, tr))
    tibble::tibble(truth = nm, rep = r, method = names(g), gain = unname(g))
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()
readr::write_csv(out, "results/er_trials.csv.gz")

summ <- out |>
  dplyr::left_join(
    dplyr::select(truths, truth, true_gain = gain),
    by = "truth"
  ) |>
  dplyr::group_by(truth, method) |>
  dplyr::summarise(
    true_gain = dplyr::first(true_gain),
    mean_gain = mean(gain),
    sd_gain = sd(gain),
    rel_bias = mean(gain) / dplyr::first(true_gain) - 1,
    choose_200 = mean(gain >= 0.10),
    .groups = "drop"
  )
readr::write_csv(summ, "results/er_summary.csv")
print(summ, width = Inf)

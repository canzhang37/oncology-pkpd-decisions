# Rscript scripts/06_verdicts.R
# Reads the result files and checks each rule in
# docs/design.md. Writes results/verdicts.txt and
# the summary tables.
library(dplyr, warn.conflicts = FALSE)
lines <- character()
say <- function(id, text, ok) {
  line <- sprintf("%-4s %-6s %s", id, if (ok) "MET" else "NOT MET", text)
  cat(line, "\n")
  lines <<- c(lines, line)
}

# A1: population PK under three BLQ treatments
pk_truth_log <- c(tka = log(1.0), tcl = log(10), tv = log(150), b_alb = -1.0)
fits <- readr::read_csv("results/poppk_fits.csv", show_col_types = FALSE)
a1 <- fits |>
  filter(ok, parameter %in% names(pk_truth_log)) |>
  mutate(
    truth = pk_truth_log[parameter],
    value = ifelse(parameter == "b_alb", estimate, exp(estimate)),
    tval = ifelse(parameter == "b_alb", truth, exp(truth))
  ) |>
  group_by(method, parameter) |>
  summarise(
    n_fits = n(),
    rel_bias = mean(value / tval - 1),
    mc_se = sd(value / tval) / sqrt(n()),
    rel_sd = sd(value / tval),
    .groups = "drop"
  )
fail <- fits |>
  distinct(rep, method, ok) |>
  group_by(method) |>
  summarise(fits = n(), failed = sum(!ok))
readr::write_csv(a1, "results/poppk_summary.csv")
readr::write_csv(fail, "results/poppk_failures.csv")
rb <- function(m, p) a1$rel_bias[a1$method == m & a1$parameter == p]
say(
  "P1",
  sprintf(
    "M3 relative bias CL/F %+.1f%%, V/F %+.1f%% (limit 5%%)",
    100 * rb("m3", "tcl"),
    100 * rb("m3", "tv")
  ),
  abs(rb("m3", "tcl")) < 0.05 && abs(rb("m3", "tv")) < 0.05
)
m1 <- max(abs(rb("drop", "tcl")), abs(rb("drop", "tv")))
m5 <- max(abs(c(rb("half", "tcl"), rb("half", "tv"), rb("half", "tka"))))
say(
  "P2",
  sprintf(
    "M1 largest |bias| in CL/F, V/F %.1f%%; M5 largest in CL/F, V/F, ka %.1f%% (each must exceed 5%%)",
    100 * m1,
    100 * m5
  ),
  m1 > 0.05 && m5 > 0.05
)

# A2: exposure-response and the dose choice
er <- readr::read_csv("results/er_summary.csv", show_col_types = FALSE)
g <- function(t, m, col) er[[col]][er$truth == t & er$method == m]
tg <- g("small", "pooled", "true_gain")
ok3 <- g("small", "pooled", "mean_gain") > 1.3 * tg &&
  g("small", "one_dose", "mean_gain") > g("small", "pooled", "mean_gain") &&
  abs(g("small", "adjusted", "mean_gain") - tg) < 0.015 &&
  abs(g("small", "by_arm", "mean_gain") - tg) < 0.015
say(
  "P3",
  sprintf(
    "small truth %.1f pts: one dose %.1f, pooled %.1f, adjusted %.1f, by arm %.1f",
    100 * tg,
    100 * g("small", "one_dose", "mean_gain"),
    100 * g("small", "pooled", "mean_gain"),
    100 * g("small", "adjusted", "mean_gain"),
    100 * g("small", "by_arm", "mean_gain")
  ),
  ok3
)
ok4 <- g("small", "pooled", "choose_200") >=
  1.5 * g("small", "adjusted", "choose_200") &&
  g("large", "adjusted", "choose_200") > g("large", "by_arm", "choose_200")
say(
  "P4",
  sprintf(
    "choose 200 mg: small truth pooled %.0f%% vs adjusted %.0f%%; large truth adjusted %.0f%% vs by arm %.0f%%",
    100 * g("small", "pooled", "choose_200"),
    100 * g("small", "adjusted", "choose_200"),
    100 * g("large", "adjusted", "choose_200"),
    100 * g("large", "by_arm", "choose_200")
  ),
  ok4
)

# B0: historical fit
hf <- readr::read_csv("results/historical_fit.csv", show_col_types = FALSE)
h <- function(p, col = "estimate") hf[[col]][hf$parameter == p]
ok5 <- abs(h("ld_drug") - h("ld_drug", "truth")) < 0.15 &&
  abs(h("lg_drug") - h("lg_drug", "truth")) < 0.15 &&
  abs(h("lsd") / h("lsd", "truth") - 1) < 0.2 &&
  abs(h("b_g") - 0.8) < 0.25
say(
  "P5",
  sprintf(
    "arm effect log d %.2f (truth %.2f), log g %.2f (truth %.2f), noise SD %.3f (0.12), b_g %.2f (0.8)",
    h("ld_drug"),
    h("ld_drug", "truth"),
    h("lg_drug"),
    h("lg_drug", "truth"),
    h("lsd"),
    h("b_g")
  ),
  ok5
)

# B: decisions at thresholds calibrated to 10% go
# under S0
tr <- readr::read_csv(
  "results/decision_trials.csv.gz",
  show_col_types = FALSE
) |>
  mutate(orr_diff = orr_drug - orr_ctrl)
cut <- tr |>
  filter(scenario == "S0") |>
  group_by(n_arm, interval) |>
  summarise(
    c_orr = quantile(orr_diff, 0.9, type = 1),
    c_model = quantile(pred_hr, 0.1, type = 1),
    c_os = quantile(os_hr, 0.1, type = 1),
    .groups = "drop"
  )
dec <- tr |>
  left_join(cut, by = c("n_arm", "interval")) |>
  group_by(scenario, n_arm, interval) |>
  summarise(
    go_orr = mean(orr_diff > c_orr),
    go_model = mean(pred_hr < c_model),
    go_os = mean(os_hr < c_os),
    fixed_orr = mean(orr_diff >= 0.15),
    fixed_model = mean(pred_hr <= 0.80),
    fixed_os = mean(os_hr <= 0.80),
    med_pred_hr = median(pred_hr),
    med_os_hr = median(os_hr),
    .groups = "drop"
  )
readr::write_csv(cut, "results/decision_thresholds.csv")
readr::write_csv(dec, "results/decision_summary.csv")
d <- function(s) filter(dec, scenario == s)
s1 <- d("S1")
say(
  "P6",
  sprintf(
    "S1 go: response rule %.0f%%-%.0f%%, model rule %.0f%%-%.0f%% across cells",
    100 * min(s1$go_orr),
    100 * max(s1$go_orr),
    100 * min(s1$go_model),
    100 * max(s1$go_model)
  ),
  all(s1$go_orr > 0.8) && all(s1$go_model < 0.25)
)
s2 <- d("S2")
c60 <- filter(s2, n_arm == 60, interval == 8)
say(
  "P7",
  sprintf(
    "S2 model minus OS go rate: %+.0f to %+.0f points; at 60/arm, 8 weeks %+.0f",
    100 * min(s2$go_model - s2$go_os),
    100 * max(s2$go_model - s2$go_os),
    100 * (c60$go_model - c60$go_os)
  ),
  all(s2$go_model > s2$go_os) && c60$go_model - c60$go_os >= 0.10
)
g6 <- filter(s2, n_arm == 60, interval == 6)$go_model
g12 <- filter(s2, n_arm == 60, interval == 12)$go_model
say(
  "P8",
  sprintf(
    "S2, 60/arm: model go rate %.0f%% with 6-week scans, %.0f%% with 12-week scans",
    100 * g6,
    100 * g12
  ),
  g6 - g12 >= 0.05
)
s4 <- d("S4")
med4 <- median(filter(tr, scenario == "S4")$pred_hr)
s490 <- filter(s4, n_arm == 90)
say(
  "P9",
  sprintf(
    "S4 median predicted HR %.2f; at 90/arm OS minus model go rate %+.0f to %+.0f points",
    med4,
    100 * min(s490$go_os - s490$go_model),
    100 * max(s490$go_os - s490$go_model)
  ),
  med4 > 0.80 && all(s490$go_os > s490$go_model)
)

# Part C (amendment 1): adverse-event risk for new
# patients
ae <- readr::read_csv("results/ae_ml.csv", show_col_types = FALSE)
a1_objf <- fits |>
  filter(ok, method == "m3") |>
  distinct(rep, objf)
c2 <- ae |>
  distinct(rep, objf) |>
  inner_join(a1_objf, by = "rep", suffix = c("_c", "_a1"))
say(
  "C2",
  sprintf(
    "refitted M3 objective values match A1 in %d of %d studies (max diff %.2g)",
    sum(abs(c2$objf_c - c2$objf_a1) < 1e-3),
    nrow(c2),
    max(abs(c2$objf_c - c2$objf_a1))
  ),
  all(abs(c2$objf_c - c2$objf_a1) < 1e-3) && nrow(c2) == nrow(a1_objf)
)
ae_sum <- ae |>
  group_by(dose, method) |>
  summarise(
    median_rmse = median(rmse),
    q25 = quantile(rmse, 0.25),
    q75 = quantile(rmse, 0.75),
    mean_pred = mean(mean_pred),
    mean_oracle = mean(mean_oracle),
    .groups = "drop"
  )
readr::write_csv(ae_sum, "results/ae_ml_summary.csv")
m <- function(dz, me) {
  ae_sum$median_rmse[ae_sum$dose == dz & ae_sum$method == me]
}
say(
  "P10",
  sprintf(
    "400 mg median RMSE: ER model %.1f, XGBoost %.1f, logistic %.1f points",
    100 * m(400, "er_model"),
    100 * m(400, "xgboost"),
    100 * m(400, "logistic")
  ),
  m(400, "er_model") < 0.5 * m(400, "xgboost")
)
say(
  "P11",
  sprintf(
    "median RMSE ER vs XGBoost: 100 mg %.1f vs %.1f, 200 mg %.1f vs %.1f points",
    100 * m(100, "er_model"),
    100 * m(100, "xgboost"),
    100 * m(200, "er_model"),
    100 * m(200, "xgboost")
  ),
  m(100, "er_model") <= m(100, "xgboost") &&
    m(200, "er_model") <= m(200, "xgboost")
)

writeLines(lines, "results/verdicts.txt")

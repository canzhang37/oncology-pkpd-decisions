# Rscript scripts/07_figures.R
# Figures for the readout, from the result files.
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)

blue <- "#2a78d6"
orange <- "#eb6834"
aqua <- "#1baf7a"
yellow <- "#eda100"
ink <- "#2b2b2a"
theme_set(
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      text = element_text(colour = ink),
      legend.position = "top",
      plot.title.position = "plot"
    )
)
save <- function(p, file, w, h) {
  ggsave(
    file.path("figures", file),
    p,
    width = w,
    height = h,
    dpi = 150,
    bg = "white"
  )
}

# A1: bias of the population PK estimates by BLQ
# handling
truth <- c(tka = 1.0, tcl = 10, tv = 150)
fits <- readr::read_csv("results/poppk_fits.csv", show_col_types = FALSE) |>
  filter(ok, parameter %in% names(truth)) |>
  mutate(
    rel = 100 * (exp(estimate) / truth[parameter] - 1),
    parameter = factor(
      parameter,
      c("tv", "tka", "tcl"),
      c("V/F", "ka", "CL/F")
    ),
    method = factor(
      method,
      c("m3", "half", "drop"),
      c("Censored (M3)", "Set to LLOQ/2 (M5)", "Dropped (M1)")
    )
  )
means <- fits |>
  group_by(parameter, method) |>
  summarise(rel = mean(rel), .groups = "drop")
p1 <- ggplot(fits, aes(rel, parameter, colour = method)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_point(
    position = position_dodge(width = 0.6),
    alpha = 0.35,
    size = 1.6
  ) +
  geom_point(
    data = means,
    position = position_dodge(width = 0.6),
    size = 3.6,
    shape = 21,
    fill = "white",
    stroke = 1.6
  ) +
  scale_colour_manual(values = c(blue, aqua, orange), name = NULL) +
  labs(
    x = "Relative error of the estimate (%), 20 simulated studies; ring = mean",
    y = NULL,
    title = "Dropping samples below the LLOQ biases the PK estimates"
  )
save(p1, "blq_bias.png", 7, 3.6)

# A2: estimated gain from 100 to 200 mg, small
# truth
er <- readr::read_csv("results/er_trials.csv.gz", show_col_types = FALSE) |>
  filter(truth == "small")
tg <- readr::read_csv("results/truth_er.csv", show_col_types = FALSE) |>
  filter(truth == "small") |>
  pull(gain)
lab <- c(
  adjusted = "Exposure-response, adjusted for albumin",
  by_arm = "Arm comparison",
  pooled = "Exposure-response, both doses",
  one_dose = "Exposure-response, 200 mg arm only"
)
pick <- er |>
  group_by(method) |>
  summarise(go = mean(gain >= 0.10))
er <- er |>
  left_join(pick, by = "method") |>
  mutate(
    method = factor(
      method,
      names(lab),
      sprintf(
        "%s\n(chooses 200 mg in %.0f%%)",
        lab,
        100 * pick$go[match(names(lab), pick$method)]
      )
    )
  )
p2 <- ggplot(er, aes(100 * gain, method, fill = method)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    colour = ink,
    linewidth = 0.4
  ) +
  geom_vline(xintercept = 100 * tg, colour = ink, linewidth = 0.8) +
  geom_vline(xintercept = 10, colour = ink, linetype = "dashed") +
  annotate(
    "text",
    x = 100 * tg,
    y = 4.55,
    label = sprintf("true gain %.1f", 100 * tg),
    hjust = 1.05,
    size = 3.3,
    colour = ink
  ) +
  annotate(
    "text",
    x = 10,
    y = 4.55,
    label = "choose 200 mg above 10",
    hjust = -0.05,
    size = 3.3,
    colour = ink
  ) +
  scale_fill_manual(values = c(blue, aqua, orange, yellow), guide = "none") +
  coord_cartesian(xlim = c(-15, 40)) +
  labs(
    x = "Estimated gain in response rate, 100 to 200 mg (points), 2,000 studies",
    y = NULL,
    title = "Albumin confounds exposure and response"
  )
save(p2, "er_gain.png", 7.5, 3.8)

# B: go rates at thresholds giving 10% go with no
# effect
dec <- readr::read_csv(
  "results/decision_summary.csv",
  show_col_types = FALSE
) |>
  filter(scenario %in% c("S1", "S2", "S4"), interval == 8) |>
  tidyr::pivot_longer(
    c(go_model, go_orr, go_os),
    names_to = "rule",
    values_to = "go"
  ) |>
  mutate(
    rule = factor(
      rule,
      c("go_model", "go_orr", "go_os"),
      c("Tumor regrowth model", "Response rate", "Interim OS")
    ),
    scenario = factor(
      scenario,
      c("S1", "S2", "S4"),
      c(
        "S1: shrinks tumors, no OS benefit\n(go should be rare)",
        "S2: slows regrowth, OS HR 0.71\n(go should be common)",
        "S4: OS HR 0.73, partly outside the tumor\n(go should be common)"
      )
    )
  )
p3 <- ggplot(dec, aes(n_arm, 100 * go, colour = rule)) +
  geom_hline(yintercept = 10, colour = "grey55", linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~scenario) +
  scale_colour_manual(values = c(blue, orange, aqua), name = NULL) +
  scale_x_continuous(breaks = c(40, 60, 90)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Patients per arm (scans every 8 weeks, read at week 50)",
    y = "Studies that go (%)",
    title = "Each rule goes in 10% of studies when the drug does nothing (dashed)"
  ) +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
save(p3, "go_rates.png", 9, 3.9)

# C: error of the predicted adverse-event risk by
# test dose
ae <- readr::read_csv("results/ae_ml_summary.csv", show_col_types = FALSE) |>
  mutate(
    method = factor(
      method,
      c("er_model", "xgboost", "logistic"),
      c("Exposure-response model", "XGBoost", "Logistic regression")
    ),
    dose = factor(
      dose,
      c(100, 150, 200, 400),
      c("100", "150", "200", "400\n(not studied)")
    )
  )
p4 <- ggplot(ae, aes(dose, 100 * median_rmse, colour = method)) +
  geom_linerange(
    aes(ymin = 100 * q25, ymax = 100 * q75),
    position = position_dodge(width = 0.5),
    linewidth = 0.9
  ) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  scale_colour_manual(values = c(blue, orange, aqua), name = NULL) +
  labs(
    x = "Dose of the new patients (mg); the studies gave 100 and 200 mg",
    y = "RMSE of predicted risk (points)",
    title = "Predicting a new patient's grade 3+ adverse event risk",
    subtitle = "Median and interquartile range over 20 studies of 120 patients"
  ) +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
save(p4, "ae_risk.png", 7, 3.8)
cat("figures written\n")

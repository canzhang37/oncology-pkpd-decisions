# Exposure-response for objective response and the
# 100 mg vs 200 mg dose choice.

er_truth <- list(
  small = list(b0 = qlogis(0.30), b_exp = 0.5, b_alb = 0.6),
  large = list(b0 = qlogis(0.30), b_exp = 1.2, b_alb = 0.6)
)
cavg_ref <- 100 / (pk_truth$cl * 24)

resp_prob <- function(cavg, alb, tr) {
  plogis(tr$b0 + tr$b_exp * log(cavg / cavg_ref) + tr$b_alb * (alb - 38) / 5)
}

# Population gain in response rate from 100 to 200
# mg, from a large simulated population
true_gain <- function(tr, n = 1e6) {
  pk <- indiv_pk(sim_patients(n, 100))
  cl <- pk$cl
  mean(resp_prob(200 / (cl * 24), pk$alb, tr)) -
    mean(resp_prob(100 / (cl * 24), pk$alb, tr))
}

sim_er_trial <- function(n_arm, tr) {
  pk <- indiv_pk(dplyr::bind_rows(
    sim_patients(n_arm, 100),
    sim_patients(n_arm, 200, n_arm + 1)
  ))
  pk$resp <- rbinom(nrow(pk), 1, resp_prob(pk$cavg, pk$alb, tr))
  pk
}

# Each estimator's predicted gain from 100 to 200
# mg for the trial's own patients (g-computation
# over their clearances)
er_gains <- function(d, cl = d$cl) {
  ce <- function(dose) log(dose / (cl * 24) / cavg_ref)
  gain <- function(fit, newdat) {
    mean(predict(fit, newdat(200), type = "response")) -
      mean(predict(fit, newdat(100), type = "response"))
  }
  d$lc <- log(d$dose / (cl * 24) / cavg_ref)
  nd <- function(dose) data.frame(lc = ce(dose), alb = d$alb)
  hi <- d$dose == 200
  f_one <- glm(resp ~ lc, binomial, data = d[hi, ])
  nd_one <- function(dose) data.frame(lc = ce(dose)[hi])
  f_pool <- glm(resp ~ lc, binomial, data = d)
  f_adj <- glm(resp ~ lc + alb, binomial, data = d)
  c(
    one_dose = gain(f_one, nd_one),
    pooled = gain(f_pool, nd),
    adjusted = gain(f_adj, nd),
    by_arm = mean(d$resp[hi]) - mean(d$resp[!hi])
  )
}

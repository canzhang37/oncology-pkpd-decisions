# Rscript -e
#   "testthat::test_file('tests/test_core.R')"
library(testthat)
source("../R/pk.R", chdir = TRUE)
source("../R/er.R", chdir = TRUE)
source("../R/tgi.R", chdir = TRUE)

test_that("single oral dose matches the textbook formula", {
  t <- c(0.5, 2, 10)
  k <- 10 / 150
  expected <- 100 * 1 / (150 * (1 - k)) * (exp(-k * t) - exp(-t))
  expect_equal(conc_oral(t, 100, 1, 10, 150), expected, tolerance = 1e-12)
})

test_that("concentrations approach the steady-state trough", {
  k <- 10 / 150
  ss_trough <- 100 /
    (150 * (1 - k)) *
    (exp(-k * 24) / (1 - exp(-k * 24)) - exp(-24) / (1 - exp(-24)))
  expect_equal(conc_oral(24 * 40, 100, 1, 10, 150), ss_trough, tolerance = 1e-6)
})

test_that("BLQ records are coded as designed", {
  set.seed(1)
  pk <- indiv_pk(sim_patients(10, 100))
  obs <- sim_pk_obs(pk)
  obs$blq[1:3] <- TRUE
  h <- nm_data(pk, obs, "half")
  expect_equal(sum(h$EVID == 1), 4 * 10)
  expect_true(all(h$DV[h$EVID == 0][1:3] == pk_truth$lloq / 2))
  expect_equal(sum(nm_data(pk, obs, "drop")$EVID == 0), sum(!obs$blq))
})

test_that("Stein curve starts at 1 and has the right nadir", {
  expect_equal(stein_f(0, 0.1, 0.02), 1)
  t_nadir <- log(0.1 / 0.02) / (0.1 + 0.02)
  eps <- 1e-4
  expect_lt(stein_f(t_nadir, 0.1, 0.02), stein_f(t_nadir - eps, 0.1, 0.02))
  expect_lt(stein_f(t_nadir, 0.1, 0.02), stein_f(t_nadir + eps, 0.1, 0.02))
})

test_that("grid log-likelihood skips patients without follow-up scans", {
  grid <- make_grid(5, 5)
  y <- matrix(log(c(50, 45, 44, 70, NA, NA)), 2, byrow = TRUE)
  ll <- grid_loglik(y, c(2, 0), c(0, 8, 16), grid, 0.12)
  expect_true(all(is.na(ll[2, ])))
  expect_true(all(ll[1, ] <= 0))
})

test_that("progression stops the scans and is flagged", {
  pts <- tibble::tibble(ts0 = 50, d = 0.001, g = 0.05, os = 200)
  withr::with_seed(1, {
    tgi_control2 <- tgi_control
    tgi_control2$meas_sd <- 1e-6
    sc <- scan_patients(pts, 8, 60, tgi_control2)
  })
  # the true size passes 1.2 x baseline near week 4,
  # so the first scan (week 8) is the progression
  # and the last one
  expect_true(sc$pfs_event)
  expect_equal(sc$m, 1)
})

test_that("ER estimators return the four gains", {
  set.seed(2)
  g <- er_gains(sim_er_trial(60, er_truth$small))
  expect_named(g, c("one_dose", "pooled", "adjusted", "by_arm"))
  expect_true(all(abs(g) < 1))
})

test_that("adverse-event risk with no spread in clearance is the curve itself", {
  source("../R/ae.R", chdir = TRUE)
  expect_equal(avg_over_cl(200, 10, 0, ae_prob), ae_prob(20), tolerance = 1e-12)
  expect_equal(ae_prob(20), plogis(-1.5))
})

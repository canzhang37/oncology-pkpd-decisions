# Population PK of a once-daily oral drug:
# simulation, BLQ handling and nlmixr2 fits.

pk_truth <- list(
  ka = 1.0, # 1/h
  cl = 10, # L/h at 70 kg, albumin 40 g/L
  v = 150, # L at 70 kg
  b_alb = -1.0, # exponent of albumin on CL
  om_ka = 0.4,
  om_cl = 0.3,
  om_v = 0.25,
  prop_sd = 0.15,
  add_sd = 0.01, # mg/L
  lloq = 0.2 # mg/L
)

# hours after the first dose; trough samples sit
# 30 min before the day 15 and day 29 doses
pk_times <- c(0.5, 2, 6, 335.5, 338, 671.5)

sim_patients <- function(n, dose, start_id = 1) {
  tibble::tibble(
    id = start_id + seq_len(n) - 1,
    dose = rep_len(dose, n),
    wt = pmin(pmax(rlnorm(n, log(70), 0.2), 40), 130),
    alb = pmin(pmax(rnorm(n, 38, 5), 25), 50)
  )
}

indiv_pk <- function(pts, truth = pk_truth) {
  n <- nrow(pts)
  pts |>
    dplyr::mutate(
      ka = truth$ka * exp(rnorm(n, 0, truth$om_ka)),
      cl = truth$cl *
        (wt / 70)^0.75 *
        (alb / 40)^truth$b_alb *
        exp(rnorm(n, 0, truth$om_cl)),
      v = truth$v * (wt / 70) * exp(rnorm(n, 0, truth$om_v)),
      cavg = dose / (cl * 24)
    )
}

# superposition of all doses given before time t
conc_oral <- function(t, dose, ka, cl, v, tau = 24) {
  k <- cl / v
  ndose <- ceiling(t / tau)
  out <- numeric(length(t))
  for (j in seq_along(t)) {
    td <- t[j] - tau * (seq_len(ndose[j]) - 1)
    out[j] <- sum(dose * ka / (v * (ka - k)) * (exp(-k * td) - exp(-ka * td)))
  }
  out
}

sim_pk_obs <- function(pk, truth = pk_truth, times = pk_times) {
  rows <- purrr::pmap(pk, function(id, dose, ka, cl, v, ...) {
    cc <- conc_oral(times, dose, ka, cl, v)
    eps <- sqrt(
      (truth$prop_sd * cc)^2 + truth$add_sd^2
    ) *
      rnorm(length(times))
    tibble::tibble(id = id, time = times, ipred = cc, dv = cc + eps)
  })
  dplyr::bind_rows(rows) |>
    dplyr::mutate(blq = dv < truth$lloq)
}

# NONMEM-style records: dose lines, then samples;
# method decides what happens to BLQ rows
nm_data <- function(
  pk,
  obs,
  method = c("m3", "drop", "half"),
  lloq = pk_truth$lloq
) {
  method <- match.arg(method)
  obs <- switch(
    method,
    m3 = dplyr::mutate(obs, cens = as.integer(blq), dv = ifelse(blq, lloq, dv)),
    drop = dplyr::filter(obs, !blq) |> dplyr::mutate(cens = 0L),
    half = dplyr::mutate(obs, cens = 0L, dv = ifelse(blq, lloq / 2, dv))
  )
  # by day 14 the drug is at steady state, so the
  # doses before each trough are written as
  # steady-state records (SS = 1) instead of 29
  # separate doses
  doses <- tidyr::expand_grid(
    dplyr::select(pk, id, dose),
    tibble::tibble(time = c(0, 312, 336, 648), ss = c(0L, 1L, 0L, 1L))
  ) |>
    dplyr::transmute(
      id,
      time,
      evid = 1L,
      amt = dose,
      ss,
      ii = ifelse(ss == 1L, 24, 0),
      dv = NA_real_,
      cens = 0L
    )
  samp <- obs |>
    dplyr::transmute(id, time, evid = 0L, amt = 0, ss = 0L, ii = 0, dv, cens)
  dplyr::bind_rows(doses, samp) |>
    dplyr::left_join(dplyr::select(pk, id, wt, alb), by = "id") |>
    dplyr::arrange(id, time, dplyr::desc(evid)) |>
    dplyr::rename_with(toupper, c(id, time, evid, amt, ss, ii, dv, cens))
}

pk_model <- function() {
  ini({
    tka <- log(1.3)
    tcl <- log(13)
    tv <- log(120)
    b_alb <- -0.5
    eta_ka ~ 0.25
    eta_cl ~ 0.1
    eta_v ~ 0.1
    prop_sd <- 0.2
    add_sd <- 0.02
  })
  model({
    ka <- exp(tka + eta_ka)
    cl <- exp(tcl + 0.75 * log(wt / 70) + b_alb * log(alb / 40) + eta_cl)
    v <- exp(tv + log(wt / 70) + eta_v)
    linCmt() ~ prop(prop_sd) + add(add_sd)
  })
}

fit_pk <- function(dat, est = "focei") {
  t0 <- Sys.time()
  fit <- tryCatch(
    suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
      pk_model,
      dat,
      est = est,
      control = nlmixr2est::foceiControl(print = 0, covMethod = "")
    ))),
    error = function(e) e
  )
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, msg = conditionMessage(fit), secs = secs))
  }
  fe <- fit$theta
  om <- diag(fit$omega)
  eb <- fit$eta
  list(
    ok = TRUE,
    secs = secs,
    fixef = fe,
    omega_sd = sqrt(om),
    objf = fit$objf,
    eta_cl = stats::setNames(eb$eta_cl, eb$ID)
  )
}

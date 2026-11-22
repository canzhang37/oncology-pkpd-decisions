# Population fit of the Stein model to a
# historical study, then the growth-rate to
# survival link used by the decision rule.

tgi_model <- function() {
  ini({
    lts0 <- log(60)
    ld <- log(0.03)
    lg <- log(0.03)
    ld_drug <- 0.5
    lg_drug <- -0.1
    eta_ts0 ~ 0.2
    eta_d ~ 0.3
    eta_g ~ 0.3
    lsd <- 0.2
  })
  model({
    ts0 <- exp(lts0 + eta_ts0)
    d <- exp(ld + ld_drug * drug + eta_d)
    g <- exp(lg + lg_drug * drug + eta_g)
    ts <- ts0 * (exp(-d * t) + exp(g * t) - 1)
    ts ~ lnorm(lsd)
  })
}

# long data for nlmixr2: one row per scan up to
# progression
tgi_long <- function(tr) {
  pts <- tr$pts
  sc <- tr$scans
  purrr::map(seq_len(nrow(pts)), function(i) {
    k <- sc$m[i] + 1
    tibble::tibble(
      ID = i,
      TIME = sc$times[seq_len(k)],
      DV = exp(sc$y[i, seq_len(k)]),
      drug = as.integer(pts$arm[i] == "drug")
    )
  }) |>
    dplyr::bind_rows()
}

fit_tgi <- function(dat, est = "saem") {
  ctl <- if (est == "saem") {
    nlmixr2est::saemControl(print = 0)
  } else {
    nlmixr2est::foceiControl(print = 0)
  }
  suppressMessages(suppressWarnings(
    nlmixr2est::nlmixr2(tgi_model, dat, est = est, control = ctl)
  ))
}

# Grade 3+ adverse event risk for a new patient:
# an exposure-response model built on the
# population PK fit, against XGBoost and logistic
# regression on dose and covariates.

ae_truth <- list(b0 = -1.5, b1 = 1.5, auc_ref = 20) # AUC in mg*h/L

ae_prob <- function(auc, tr = ae_truth) {
  plogis(tr$b0 + tr$b1 * log(auc / tr$auc_ref))
}

# covariates with no role in the truth
add_noise_covs <- function(d) {
  n <- nrow(d)
  dplyr::mutate(
    d,
    age = round(rnorm(n, 62, 10)),
    female = rbinom(n, 1, 0.45),
    ecog = rbinom(n, 1, 0.4),
    crcl = rlnorm(n, log(80), 0.25),
    ldh = rlnorm(n, log(250), 0.4)
  )
}

gh <- statmod::gauss.quad.prob(20, "normal")

# risk averaged over the between-patient spread of
# CL (rows: patients, columns: quadrature nodes)
avg_over_cl <- function(dose, cl_typ, om_cl, f) {
  cl <- outer(cl_typ, exp(om_cl * gh$nodes))
  as.vector(f(dose / cl) %*% gh$weights)
}

oracle_risk <- function(new) {
  cl_typ <- pk_truth$cl * (new$wt / 70)^0.75 * (new$alb / 40)^pk_truth$b_alb
  avg_over_cl(new$dose, cl_typ, pk_truth$om_cl, ae_prob)
}

feature_matrix <- function(d) {
  as.matrix(d[, c("dose", "wt", "alb", "age", "female", "ecog", "crcl", "ldh")])
}

fit_xgb <- function(d) {
  dtr <- xgboost::xgb.DMatrix(feature_matrix(d), label = d$ae)
  best <- NULL
  for (depth in 1:3) {
    prm <- list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      eta = 0.05,
      max_depth = depth,
      nthread = 1
    )
    cv <- xgboost::xgb.cv(
      params = prm,
      data = dtr,
      nrounds = 1000,
      nfold = 5,
      early_stopping_rounds = 50,
      verbose = 0
    )
    it <- cv$early_stop$best_iteration
    loss <- cv$evaluation_log$test_logloss_mean[it]
    if (is.null(best) || loss < best$loss) {
      best <- list(prm = prm, it = it, loss = loss)
    }
  }
  model <- xgboost::xgb.train(best$prm, dtr, nrounds = best$it, verbose = 0)
  list(model = model, depth = best$prm$max_depth, rounds = best$it)
}

# fe: fixed effects of the M3 fit; eta_cl: its
# empirical Bayes etas, in patient order
ae_predictions <- function(train, new, fe, om_cl, eta_cl) {
  typ <- function(d) {
    exp(fe[["tcl"]] + 0.75 * log(d$wt / 70) + fe[["b_alb"]] * log(d$alb / 40))
  }
  train$lauc <- log(train$dose / (typ(train) * exp(eta_cl)) / ae_truth$auc_ref)
  er <- glm(ae ~ lauc, binomial, data = train)
  a <- coef(er)
  p_er <- avg_over_cl(new$dose, typ(new), om_cl, function(auc) {
    plogis(a[[1]] + a[[2]] * log(auc / ae_truth$auc_ref))
  })
  lr <- glm(
    ae ~ log(dose) +
      log(wt) +
      log(alb) +
      age +
      female +
      ecog +
      log(crcl) +
      log(ldh),
    binomial,
    data = train
  )
  xg <- fit_xgb(train)
  list(
    p = tibble::tibble(
      er_model = p_er,
      logistic = unname(predict(lr, new, type = "response")),
      xgboost = predict(xg$model, xgboost::xgb.DMatrix(feature_matrix(new)))
    ),
    er_slope = a[[2]],
    xgb_depth = xg$depth,
    xgb_rounds = xg$rounds
  )
}

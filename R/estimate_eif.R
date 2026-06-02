# ============================================================
# Shared Super Learner libraries
# ============================================================

.make_Q_learners <- function(is_binary = TRUE,
                             outcome_family = NULL,
                             num_trees = 100) {
  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' is required for flexible Super Learner nuisances.")
  }

  if (is.null(outcome_family)) {
    outcome_family <- if (isTRUE(is_binary)) stats::binomial() else stats::gaussian()
  }

  if (isTRUE(is_binary)) {
    list(
      sl3::Lrnr_glm_fast$new(family = outcome_family),
      sl3::Lrnr_glmnet$new(family = "binomial"),
      sl3::Lrnr_glmnet$new(family = "binomial", alpha = 0.5),
      sl3::Lrnr_ranger$new(num.trees = num_trees),
      sl3::Lrnr_gam$new(family = outcome_family),
      sl3::Lrnr_earth$new(glm = list(family = outcome_family))
    )
  } else {
    list(
      sl3::Lrnr_glm_fast$new(family = outcome_family),
      sl3::Lrnr_glmnet$new(),
      sl3::Lrnr_ranger$new(num.trees = num_trees),
      sl3::Lrnr_gam$new(),
      sl3::Lrnr_earth$new()
    )
  }
}

.make_g_learners <- function(num_trees = 100) {
  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' is required for flexible Super Learner nuisances.")
  }

  list(
    sl3::Lrnr_glm_fast$new(family = stats::binomial()),
    sl3::Lrnr_glmnet$new(family = "binomial"),
    sl3::Lrnr_ranger$new(num.trees = num_trees)
  )
}

.make_sl <- function(learners) {
  sl3::Lrnr_sl$new(
    learners = learners,
    metalearner = sl3::Lrnr_nnls$new()
  )
}


estimate_nuisances <- function(data,
                               outcome = "y",
                               treat = "treatment",
                               covars = NULL,
                               nuisance_type = c("simple", "flexible"),
                               family_outcome = NULL,
                               clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)

  if (is.null(covars)) {
    covars <- setdiff(names(data), c(outcome, treat))
  }

  data <- as.data.frame(data)
  data <- .safe_factor_country(data)
  data[[treat]] <- as.numeric(data[[treat]])

  yvals <- data[[outcome]]
  is_binary <- all(stats::na.omit(unique(yvals)) %in% c(0, 1)) &&
    length(unique(stats::na.omit(yvals))) <= 2L

  if (is.null(family_outcome)) {
    family_outcome <- if (is_binary) stats::binomial() else stats::gaussian()
  }

  clamp_prob <- function(x) .clamp_prob(x, clamp)

  # ----------------------------
  # simple nuisance: GLM
  # ----------------------------

  if (nuisance_type == "simple") {
    g_formula <- stats::as.formula(
      paste(
        treat,
        "~",
        if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
      )
    )

    g_mod <- stats::glm(
      g_formula,
      data = data,
      family = stats::binomial()
    )

    g_hat <- stats::predict(g_mod, type = "response")

    Q_formula <- stats::as.formula(
      paste(
        outcome,
        "~",
        treat,
        if (length(covars) > 0) {
          paste("+", paste(covars, collapse = " + "))
        } else {
          ""
        }
      )
    )

    Q_mod <- stats::glm(
      Q_formula,
      data = data,
      family = family_outcome
    )

    QAW <- stats::predict(Q_mod, type = "response")

    new1 <- data
    new0 <- data
    new1[[treat]] <- 1
    new0[[treat]] <- 0

    Q1 <- stats::predict(Q_mod, newdata = new1, type = "response")
    Q0 <- stats::predict(Q_mod, newdata = new0, type = "response")

    return(list(
      Q_mod = Q_mod,
      g_mod = g_mod,
      Q1 = if (is_binary) clamp_prob(Q1) else as.numeric(Q1),
      Q0 = if (is_binary) clamp_prob(Q0) else as.numeric(Q0),
      QAW = if (is_binary) clamp_prob(QAW) else as.numeric(QAW),
      g = clamp_prob(g_hat),
      is_binary = is_binary,
      nuisance_type = nuisance_type
    ))
  }

  # ----------------------------
  # flexible nuisance: sl3 Super Learner
  # ----------------------------

  if (nuisance_type == "flexible") {
    if (!requireNamespace("sl3", quietly = TRUE)) {
      stop("Package 'sl3' is required for nuisance_type = 'flexible'.")
    }

    sl_data <- data

    Q_learners <- .make_Q_learners(
      is_binary = is_binary,
      outcome_family = family_outcome
    )
    g_learners <- .make_g_learners()

    Q_sl <- .make_sl(Q_learners)
    g_sl <- .make_sl(g_learners)

    task_g <- sl3::make_sl3_Task(
      data = sl_data,
      outcome = treat,
      covariates = covars
    )

    g_fit <- g_sl$train(task_g)
    g_hat <- as.numeric(g_fit$predict(task_g))

    Q_covars <- c(treat, covars)

    task_Q <- sl3::make_sl3_Task(
      data = sl_data,
      outcome = outcome,
      covariates = Q_covars
    )

    Q_fit <- Q_sl$train(task_Q)
    QAW <- as.numeric(Q_fit$predict(task_Q))

    new1 <- sl_data
    new0 <- sl_data
    new1[[treat]] <- 1
    new0[[treat]] <- 0

    task_Q1 <- sl3::make_sl3_Task(
      data = new1,
      outcome = outcome,
      covariates = Q_covars
    )

    task_Q0 <- sl3::make_sl3_Task(
      data = new0,
      outcome = outcome,
      covariates = Q_covars
    )

    Q1 <- as.numeric(Q_fit$predict(task_Q1))
    Q0 <- as.numeric(Q_fit$predict(task_Q0))

    return(list(
      Q_mod = Q_fit,
      g_mod = g_fit,
      Q1 = if (is_binary) clamp_prob(Q1) else as.numeric(Q1),
      Q0 = if (is_binary) clamp_prob(Q0) else as.numeric(Q0),
      QAW = if (is_binary) clamp_prob(QAW) else as.numeric(QAW),
      g = clamp_prob(g_hat),
      is_binary = is_binary,
      nuisance_type = nuisance_type
    ))
  }
}

one_step_ATE <- function(data,
                         nuisances = NULL,
                         outcome = "y",
                         treat = "treatment",
                         covars = NULL,
                         nuisance_type = c("simple", "flexible")) {
  nuisance_type <- match.arg(nuisance_type)

  if (is.null(nuisances)) {
    nuisances <- estimate_nuisances(
      data,
      outcome = outcome,
      treat = treat,
      covars = covars,
      nuisance_type = nuisance_type
    )
  }

  Y <- data[[outcome]]
  A <- data[[treat]]

  Q1 <- nuisances$Q1
  Q0 <- nuisances$Q0
  g <- .clamp_prob(nuisances$g)

  n <- nrow(data)

  psi_plug <- mean(Q1 - Q0)

  contrib <- (A / g) * (Y - Q1) -
    ((1 - A) / (1 - g)) * (Y - Q0) +
    (Q1 - Q0)

  psi_one_step <- mean(contrib)
  IC <- contrib - psi_one_step

  se <- stats::sd(IC, na.rm = TRUE) / sqrt(n)
  ci <- psi_one_step + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimator = "one-step",
    psi = as.numeric(psi_one_step),
    psi_plug = as.numeric(psi_plug),
    se = as.numeric(se),
    ci_lower = ci[1],
    ci_upper = ci[2],
    n = n,
    influence = IC,
    nuisances = nuisances
  )
}

one_step_ATE <- function(data,
                         nuisances = NULL,
                         outcome = "y",
                         treat = "treatment",
                         covars = NULL,
                         nuisance_type = c("simple", "flexible")) {
  nuisance_type <- match.arg(nuisance_type)

  if (is.null(nuisances)) {
    nuisances <- estimate_nuisances(
      data,
      outcome = outcome,
      treat = treat,
      covars = covars,
      nuisance_type = nuisance_type
    )
  }

  Y <- data[[outcome]]
  A <- data[[treat]]

  Q1 <- nuisances$Q1
  Q0 <- nuisances$Q0
  g <- .clamp_prob(nuisances$g)

  n <- nrow(data)

  psi_plug <- mean(Q1 - Q0)

  contrib <- (A / g) * (Y - Q1) -
    ((1 - A) / (1 - g)) * (Y - Q0) +
    (Q1 - Q0)

  psi_one_step <- mean(contrib)
  IC <- contrib - psi_one_step

  se <- stats::sd(IC, na.rm = TRUE) / sqrt(n)
  ci <- psi_one_step + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimator = "one-step",
    psi = as.numeric(psi_one_step),
    psi_plug = as.numeric(psi_plug),
    se = as.numeric(se),
    ci_lower = ci[1],
    ci_upper = ci[2],
    n = n,
    influence = IC,
    nuisances = nuisances
  )
}


eif_trial_onestep <- function(data,
                              outcome = "y",
                              treat = "treatment",
                              covars = c("L1", "L2"),
                              treat_s1_active = "a",
                              treat_s1_control = "sc1",
                              treat_s2_active = "b",
                              treat_s2_control = "sc2",
                              nuisance_type = c("simple", "flexible")) {
  nuisance_type <- match.arg(nuisance_type)

  dat_s1 <- subset(data, s == 1)
  dat_s2 <- subset(data, s == 2)

  if (nrow(dat_s1) == 0) stop("No data for trial s = 1")
  if (nrow(dat_s2) == 0) stop("No data for trial s = 2")

  dat_s1[[treat]] <- as.numeric(as.character(dat_s1[[treat]]) == treat_s1_active)
  dat_s2[[treat]] <- as.numeric(as.character(dat_s2[[treat]]) == treat_s2_active)

  nuis1 <- estimate_nuisances(
    dat_s1,
    outcome = outcome,
    treat = treat,
    covars = covars,
    nuisance_type = nuisance_type
  )

  res1 <- one_step_ATE(
    dat_s1,
    nuisances = nuis1,
    outcome = outcome,
    treat = treat,
    covars = covars,
    nuisance_type = nuisance_type
  )

  nuis2 <- estimate_nuisances(
    dat_s2,
    outcome = outcome,
    treat = treat,
    covars = covars,
    nuisance_type = nuisance_type
  )

  res2 <- one_step_ATE(
    dat_s2,
    nuisances = nuis2,
    outcome = outcome,
    treat = treat,
    covars = covars,
    nuisance_type = nuisance_type
  )

  ate_diff <- res1$psi - res2$psi

  se_diff <- sqrt(
    stats::var(res1$influence, na.rm = TRUE) / res1$n +
      stats::var(res2$influence, na.rm = TRUE) / res2$n
  )

  ci_lower <- ate_diff + stats::qnorm(0.025) * se_diff
  ci_upper <- ate_diff + stats::qnorm(0.975) * se_diff

  list(
    res1 = res1,
    res2 = res2,
    ate_diff = ate_diff,
    se_diff = se_diff,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n1 = res1$n,
    n2 = res2$n,
    nuisance_type = nuisance_type
  )
}



.fit_binary_glm <- function(formula, data) {
  stats::glm(formula, data = data, family = stats::binomial())
}

.predict_outcome_binary <- function(fit, newdata) {
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}

.eif_prepare_factor_covariates <- function(data, covariates, levels_map = NULL) {
  data <- as.data.frame(data)

  for (v in intersect(covariates, names(data))) {
    if (is.character(data[[v]]) || is.factor(data[[v]])) {
      if (!is.null(levels_map) && v %in% names(levels_map)) {
        data[[v]] <- factor(as.character(data[[v]]), levels = levels_map[[v]])
      } else {
        data[[v]] <- factor(as.character(data[[v]]))
      }
    }
  }

  data
}

.eif_factor_levels <- function(data, covariates) {
  out <- list()

  for (v in intersect(covariates, names(data))) {
    if (is.character(data[[v]]) || is.factor(data[[v]])) {
      out[[v]] <- levels(factor(as.character(data[[v]])))
    }
  }

  out
}

.eif_make_binary_sl <- function(learner_set = c("Q", "g")) {
  learner_set <- match.arg(learner_set)

  learners <- if (learner_set == "Q") {
    .make_Q_learners(is_binary = TRUE, outcome_family = stats::binomial())
  } else {
    .make_g_learners()
  }

  .make_sl(learners)
}

.eif_fit_binary_regression <- function(data,
                                       outcome,
                                       covariates,
                                       nuisance_type = c("simple", "flexible"),
                                       learner_set = c("Q", "g"),
                                       clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  learner_set <- match.arg(learner_set)

  data <- .safe_factor_country(as.data.frame(data))
  covariates <- intersect(covariates, names(data))

  if (!(outcome %in% names(data))) {
    stop("Outcome variable '", outcome, "' not found in training data.")
  }

  levels_map <- .eif_factor_levels(data, covariates)
  data <- .eif_prepare_factor_covariates(data, covariates, levels_map)

  if (length(unique(stats::na.omit(data[[outcome]]))) < 2) {
    stop("Binary outcome '", outcome, "' must contain both 0 and 1 in the training data.")
  }

  if (nuisance_type == "simple") {
    form <- stats::as.formula(
      paste(outcome, "~", if (length(covariates) > 0) paste(covariates, collapse = " + ") else "1")
    )

    fit <- stats::glm(
      form,
      data = data,
      family = stats::binomial()
    )

    return(list(
      predict = function(newdata) {
        newdata <- .safe_factor_country(as.data.frame(newdata))
        newdata <- .eif_prepare_factor_covariates(newdata, covariates, levels_map)

        if (!(outcome %in% names(newdata))) {
          newdata[[outcome]] <- 0
        }

        .clamp_prob(
          stats::predict(fit, newdata = newdata, type = "response"),
          clamp
        )
      }
    ))
  }

  sl_fit <- .eif_make_binary_sl(learner_set = learner_set)

  task <- sl3::make_sl3_Task(
    data = data,
    outcome = outcome,
    covariates = covariates
  )

  trained <- sl_fit$train(task)

  list(
    predict = function(newdata) {
      newdata <- .safe_factor_country(as.data.frame(newdata))
      newdata <- .eif_prepare_factor_covariates(newdata, covariates, levels_map)

      if (!(outcome %in% names(newdata))) {
        newdata[[outcome]] <- 0
      }

      task_new <- sl3::make_sl3_Task(
        data = newdata,
        outcome = outcome,
        covariates = covariates
      )

      .clamp_prob(as.numeric(trained$predict(task_new)), clamp)
    }
  )
}


eif_direct_est <- function(dataset,
                           a_val = "a",
                           b_val = "b",
                           clamp = c(1e-6, 1 - 1e-6),
                           nuisance_type = c("simple", "flexible")) {
  nuisance_type <- match.arg(nuisance_type)

  dataset <- as.data.frame(dataset)
  dataset <- .safe_factor_country(dataset)

  if (!all(c("y", "treatment", "s") %in% names(dataset))) {
    stop("dataset must contain y, treatment, s")
  }

  n <- nrow(dataset)

  data_s1 <- dataset[dataset$s == 1, , drop = FALSE]
  data_s2 <- dataset[dataset$s == 2, , drop = FALSE]

  if (nrow(data_s1) == 0) stop("no rows with s == 1")
  if (nrow(data_s2) == 0) stop("no rows with s == 2")

  Y <- dataset$y
  A <- as.character(dataset$treatment)
  S <- dataset$s

  covars <- c("L1", "L2")
  Q_covars <- c("treatment", covars)

  # Treatment mechanism within trial 1: P(A = a | S = 1, L)
  data_s1$.A_a <- as.numeric(as.character(data_s1$treatment) == a_val)
  pi_a_fit <- .eif_fit_binary_regression(
    data = data_s1,
    outcome = ".A_a",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )
  ps_a_s1 <- pi_a_fit$predict(dataset)

  # Treatment mechanism within trial 2: P(A = b | S = 2, L)
  data_s2$.A_b <- as.numeric(as.character(data_s2$treatment) == b_val)
  pi_b_fit <- .eif_fit_binary_regression(
    data = data_s2,
    outcome = ".A_b",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )
  ps_b_s2 <- pi_b_fit$predict(dataset)

  # Outcome regression in trial 1: E[Y | S = 1, A, L], evaluated at A = a
  mu_s1_fit <- .eif_fit_binary_regression(
    data = data_s1,
    outcome = "y",
    covariates = Q_covars,
    nuisance_type = nuisance_type,
    learner_set = "Q",
    clamp = clamp
  )

  pred_data_a <- dataset
  pred_data_a$treatment <- factor(a_val, levels = levels(factor(data_s1$treatment)))
  mu_a_s1 <- mu_s1_fit$predict(pred_data_a)

  # Outcome regression in trial 2: E[Y | S = 2, A, L], evaluated at A = b
  mu_s2_fit <- .eif_fit_binary_regression(
    data = data_s2,
    outcome = "y",
    covariates = Q_covars,
    nuisance_type = nuisance_type,
    learner_set = "Q",
    clamp = clamp
  )

  pred_data_b <- dataset
  pred_data_b$treatment <- factor(b_val, levels = levels(factor(data_s2$treatment)))
  mu_b_s2 <- mu_s2_fit$predict(pred_data_b)

  # Trial membership mechanism: P(S = 1 | L)
  dataset$.S1 <- as.numeric(dataset$s == 1)
  g_fit <- .eif_fit_binary_regression(
    data = dataset,
    outcome = ".S1",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )
  prob_s1 <- g_fit$predict(dataset)
  prob_s2 <- .clamp_prob(1 - prob_s1, clamp)

  I_A_a_S_1 <- as.numeric(A == a_val & S == 1)
  I_A_b_S_2 <- as.numeric(A == b_val & S == 2)
  I_S_1 <- as.numeric(S == 1)
  I_S_2 <- as.numeric(S == 2)

  term1 <- (I_A_a_S_1 / ps_a_s1) * (Y - mu_a_s1) + I_S_1 * mu_a_s1
  term2 <- (I_A_b_S_2 / ps_b_s2) * (Y - mu_b_s2) + I_S_2 * mu_b_s2
  term3 <- (I_A_a_S_1 * prob_s2) / (ps_a_s1 * prob_s1) * (Y - mu_a_s1) + I_S_2 * mu_a_s1
  term4 <- (I_A_b_S_2 * prob_s1) / (ps_b_s2 * prob_s2) * (Y - mu_b_s2) + I_S_1 * mu_b_s2

  eif_vec <- term1 - term2 + term3 - term4
  psi_hat <- mean(eif_vec, na.rm = TRUE)
  se_plug <- stats::sd(eif_vec, na.rm = TRUE) / sqrt(sum(!is.na(eif_vec)))
  ci_plug <- psi_hat + stats::qnorm(c(0.025, 0.975)) * se_plug

  list(
    eif = eif_vec,
    estimate = psi_hat,
    se_if = se_plug,
    CI_if = ci_plug,
    components = data.frame(term1 = term1, term2 = term2, term3 = term3, term4 = term4),
    nuisance_type = nuisance_type
  )
}

eif_theta_est <- function(dataset,
                          a = "a",
                          c1 = "sc1",
                          b = "b",
                          c2 = "sc2",
                          clamp = c(1e-6, 1 - 1e-6),
                          nuisance_type = c("simple", "flexible")) {
  nuisance_type <- match.arg(nuisance_type)

  dataset <- as.data.frame(dataset)
  dataset <- .safe_factor_country(dataset)

  if (!all(c("y", "treatment", "s") %in% names(dataset))) {
    stop("dataset must contain y, treatment, s")
  }

  n <- nrow(dataset)
  Y <- dataset$y
  A <- as.character(dataset$treatment)
  S <- dataset$s

  data_s1 <- dataset[S == 1, , drop = FALSE]
  data_s2 <- dataset[S == 2, , drop = FALSE]

  if (nrow(data_s1) == 0 || nrow(data_s2) == 0) {
    stop("need both trial strata")
  }

  covars <- c("L1", "L2")
  Q_covars <- c("treatment", covars)

  # Treatment mechanisms in trial 1
  data_s1$.A_a <- as.numeric(as.character(data_s1$treatment) == a)
  data_s1$.A_c1 <- as.numeric(as.character(data_s1$treatment) == c1)

  pi_a_s1_fit <- .eif_fit_binary_regression(
    data = data_s1,
    outcome = ".A_a",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )
  pi_c1_s1_fit <- .eif_fit_binary_regression(
    data = data_s1,
    outcome = ".A_c1",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )

  ps_a_s1 <- pi_a_s1_fit$predict(dataset)
  ps_c1_s1 <- pi_c1_s1_fit$predict(dataset)

  # Treatment mechanisms in trial 2
  data_s2$.A_b <- as.numeric(as.character(data_s2$treatment) == b)
  data_s2$.A_c2 <- as.numeric(as.character(data_s2$treatment) == c2)

  pi_b_s2_fit <- .eif_fit_binary_regression(
    data = data_s2,
    outcome = ".A_b",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )
  pi_c2_s2_fit <- .eif_fit_binary_regression(
    data = data_s2,
    outcome = ".A_c2",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )

  ps_b_s2 <- pi_b_s2_fit$predict(dataset)
  ps_c2_s2 <- pi_c2_s2_fit$predict(dataset)

  # Outcome regressions by trial, evaluated at each treatment value
  mu_s1_fit <- .eif_fit_binary_regression(
    data = data_s1,
    outcome = "y",
    covariates = Q_covars,
    nuisance_type = nuisance_type,
    learner_set = "Q",
    clamp = clamp
  )

  new_a <- dataset
  new_c1 <- dataset
  new_a$treatment <- factor(a, levels = levels(factor(data_s1$treatment)))
  new_c1$treatment <- factor(c1, levels = levels(factor(data_s1$treatment)))

  mu_a_s1 <- mu_s1_fit$predict(new_a)
  mu_c1_s1 <- mu_s1_fit$predict(new_c1)

  mu_s2_fit <- .eif_fit_binary_regression(
    data = data_s2,
    outcome = "y",
    covariates = Q_covars,
    nuisance_type = nuisance_type,
    learner_set = "Q",
    clamp = clamp
  )

  new_b <- dataset
  new_c2 <- dataset
  new_b$treatment <- factor(b, levels = levels(factor(data_s2$treatment)))
  new_c2$treatment <- factor(c2, levels = levels(factor(data_s2$treatment)))

  mu_b_s2 <- mu_s2_fit$predict(new_b)
  mu_c2_s2 <- mu_s2_fit$predict(new_c2)

  # Trial membership mechanism: P(S = 1 | L)
  dataset$.S1 <- as.numeric(dataset$s == 1)
  g_fit <- .eif_fit_binary_regression(
    data = dataset,
    outcome = ".S1",
    covariates = covars,
    nuisance_type = nuisance_type,
    learner_set = "g",
    clamp = clamp
  )
  prob_s1 <- g_fit$predict(dataset)
  prob_s2 <- .clamp_prob(1 - prob_s1, clamp)

  I_A_a_S_1 <- as.numeric(A == a & S == 1)
  I_A_c1_S_1 <- as.numeric(A == c1 & S == 1)
  I_A_b_S_2 <- as.numeric(A == b & S == 2)
  I_A_c2_S_2 <- as.numeric(A == c2 & S == 2)
  I_S_1 <- as.numeric(S == 1)
  I_S_2 <- as.numeric(S == 2)

  e_theta1_part1 <- (I_A_a_S_1 / ps_a_s1) * (Y - mu_a_s1) + I_S_1 * mu_a_s1
  e_theta1_part2 <- (I_A_c1_S_1 / ps_c1_s1) * (Y - mu_c1_s1) + I_S_1 * mu_c1_s1
  e_theta1_part3 <- (I_A_a_S_1 * prob_s2) / (ps_a_s1 * prob_s1) * (Y - mu_a_s1) + I_S_2 * mu_a_s1
  e_theta1_part4 <- (I_A_c1_S_1 * prob_s2) / (ps_c1_s1 * prob_s1) * (Y - mu_c1_s1) + I_S_2 * mu_c1_s1

  eif_theta1 <- e_theta1_part1 - e_theta1_part2 + e_theta1_part3 - e_theta1_part4

  e_theta2_part1 <- (I_A_b_S_2 / ps_b_s2) * (Y - mu_b_s2) + I_S_2 * mu_b_s2
  e_theta2_part2 <- (I_A_c2_S_2 / ps_c2_s2) * (Y - mu_c2_s2) + I_S_2 * mu_c2_s2
  e_theta2_part3 <- (I_A_b_S_2 * prob_s1) / (ps_b_s2 * prob_s2) * (Y - mu_b_s2) + I_S_1 * mu_b_s2
  e_theta2_part4 <- (I_A_c2_S_2 * prob_s1) / (ps_c2_s2 * prob_s2) * (Y - mu_c2_s2) + I_S_1 * mu_c2_s2

  eif_theta2 <- e_theta2_part1 - e_theta2_part2 + e_theta2_part3 - e_theta2_part4

  eif_theta <- eif_theta1 - eif_theta2
  theta_hat <- mean(eif_theta, na.rm = TRUE)
  se_theta <- stats::sd(eif_theta, na.rm = TRUE) / sqrt(sum(!is.na(eif_theta)))

  list(
    theta = theta_hat,
    eif = eif_theta,
    se = se_theta,
    nuisance_type = nuisance_type
  )
}


eif_indirect_est <- function(dataset,
                             a = "a",
                             b = "b",
                             c1 = "sc1",
                             c2 = "sc2",
                             clamp = c(1e-6, 1 - 1e-6),
                             nuisance_type = c("simple", "flexible")) {
  nuisance_type <- match.arg(nuisance_type)

  th <- eif_theta_est(
    dataset,
    a = a,
    c1 = c1,
    b = b,
    c2 = c2,
    clamp = clamp,
    nuisance_type = nuisance_type
  )

  phi_eif_res <- eif_direct_est(
    dataset,
    a_val = c1,
    b_val = c2,
    clamp = clamp,
    nuisance_type = nuisance_type
  )

  eif_indirect_vec <- th$eif + phi_eif_res$eif
  indirect_hat <- mean(eif_indirect_vec, na.rm = TRUE)
  se_indirect <- stats::sd(eif_indirect_vec, na.rm = TRUE) / sqrt(sum(!is.na(eif_indirect_vec)))
  CI_indirect <- indirect_hat + stats::qnorm(c(0.025, 0.975)) * se_indirect

  list(
    indirect = indirect_hat,
    se_if = se_indirect,
    CI_if = CI_indirect,
    eif = eif_indirect_vec,
    theta = th$theta,
    phi = phi_eif_res$estimate,
    nuisance_type = nuisance_type
  )
}

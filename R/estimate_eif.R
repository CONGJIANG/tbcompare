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
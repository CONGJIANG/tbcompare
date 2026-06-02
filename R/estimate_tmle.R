# ------------------------------------------------------------
# TMLE helpers
# ------------------------------------------------------------

.tmle_clip <- function(x, clamp = c(1e-6, 1 - 1e-6)) {
  pmin(pmax(as.numeric(x), clamp[1]), clamp[2])
}

.tmle_logit <- function(p, clamp = c(1e-6, 1 - 1e-6)) {
  stats::qlogis(.tmle_clip(p, clamp))
}

.tmle_expit <- function(x) {
  stats::plogis(x)
}

.tmle_safe_var <- function(x) {
  stats::var(x, na.rm = TRUE)
}

.tmle_safe_country <- function(data) {
  data <- as.data.frame(data)
  if ("country" %in% names(data)) {
    data$country <- factor(data$country, levels = c(1, 2, 3, 4))
  }
  data
}

.tmle_make_learners <- function(nuisance_type = c("simple", "flexible"),
                                outcome_family = stats::binomial()) {
  nuisance_type <- match.arg(nuisance_type)

  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' is required for nuisance_type = 'flexible'.", call. = FALSE)
  }

  if (nuisance_type == "simple") {
    Q_learner <- sl3::Lrnr_glm_fast$new(family = outcome_family)
    g_learner <- sl3::Lrnr_glm_fast$new(family = stats::binomial())
  } else {
    Q_learner <- .make_sl(
      .make_Q_learners(
        is_binary = identical(outcome_family$family, "binomial"),
        outcome_family = outcome_family
      )
    )

    g_learner <- .make_sl(
      .make_g_learners()
    )
  }

  list(Q = Q_learner, g = g_learner)
}

# ------------------------------------------------------------
# Core binary-treatment TMLE ATE
# ------------------------------------------------------------

tmle_ate_binary <- function(data,
                            outcome = "y",
                            treat = "A_bin",
                            covars = c("L1", "L2"),
                            nuisance_type = c("simple", "flexible"),
                            clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)

  data <- .tmle_safe_country(data)
  data[[treat]] <- as.numeric(data[[treat]])

  if (length(unique(stats::na.omit(data[[treat]]))) < 2) {
    stop("Treatment variable must contain both 0 and 1.")
  }

  Q_covars <- c(treat, covars)

  if (nuisance_type == "simple") {
    Q_formula <- stats::as.formula(
      paste(outcome, "~", paste(Q_covars, collapse = " + "))
    )

    Q_fit <- stats::glm(
      Q_formula,
      data = data,
      family = stats::binomial()
    )

    data1 <- data
    data0 <- data
    data1[[treat]] <- 1
    data0[[treat]] <- 0

    Q1 <- .tmle_clip(
      stats::predict(Q_fit, newdata = data1, type = "response"),
      clamp
    )

    Q0 <- .tmle_clip(
      stats::predict(Q_fit, newdata = data0, type = "response"),
      clamp
    )

    g_formula <- stats::as.formula(
      paste(treat, "~", paste(covars, collapse = " + "))
    )

    g_fit <- stats::glm(
      g_formula,
      data = data,
      family = stats::binomial()
    )

    g <- .tmle_clip(
      stats::predict(g_fit, type = "response"),
      clamp
    )
  } else {
    learners <- .tmle_make_learners(
      nuisance_type = nuisance_type,
      outcome_family = stats::binomial()
    )

    task_Q <- sl3::make_sl3_Task(
      data = data,
      outcome = outcome,
      covariates = Q_covars
    )

    Q_fit <- learners$Q$train(task_Q)

    data1 <- data
    data0 <- data
    data1[[treat]] <- 1
    data0[[treat]] <- 0

    task_Q1 <- sl3::make_sl3_Task(
      data = data1,
      outcome = outcome,
      covariates = Q_covars
    )

    task_Q0 <- sl3::make_sl3_Task(
      data = data0,
      outcome = outcome,
      covariates = Q_covars
    )

    Q1 <- .tmle_clip(Q_fit$predict(task_Q1), clamp)
    Q0 <- .tmle_clip(Q_fit$predict(task_Q0), clamp)

    task_g <- sl3::make_sl3_Task(
      data = data,
      outcome = treat,
      covariates = covars
    )

    g_fit <- learners$g$train(task_g)
    g <- .tmle_clip(g_fit$predict(task_g), clamp)
  }

  A <- data[[treat]]
  Y <- data[[outcome]]

  QAW <- ifelse(A == 1, Q1, Q0)
  QAW <- .tmle_clip(QAW, clamp)

  H <- A / g - (1 - A) / (1 - g)

  eps_fit <- stats::glm(
    Y ~ -1 + H,
    family = stats::binomial(),
    offset = .tmle_logit(QAW, clamp)
  )

  eps <- as.numeric(stats::coef(eps_fit)[1])

  Q1_star <- .tmle_clip(
    .tmle_expit(.tmle_logit(Q1, clamp) + eps / g),
    clamp
  )

  Q0_star <- .tmle_clip(
    .tmle_expit(.tmle_logit(Q0, clamp) - eps / (1 - g)),
    clamp
  )

  psi <- mean(Q1_star - Q0_star)

  IC <- (A / g) * (Y - Q1_star) -
    ((1 - A) / (1 - g)) * (Y - Q0_star) +
    (Q1_star - Q0_star) -
    psi

  n <- nrow(data)
  se <- stats::sd(IC, na.rm = TRUE) / sqrt(n)
  CI <- psi + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimate = as.numeric(psi),
    se = as.numeric(se),
    CI = as.numeric(CI),
    nuisance_type = nuisance_type,
    influence_curve = IC
  )
}

# ------------------------------------------------------------
# TMLE method: target trial emulation, restricted a vs b
# ------------------------------------------------------------

tmle_tte <- function(data,
                     nuisance_type = c("simple", "flexible"),
                     a_val = "a",
                     b_val = "b",
                     covars = c("L1", "L2"),
                     outcome = "y",
                     clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  data <- .tmle_safe_country(data)

  dat_ab <- data[data$treatment %in% c(a_val, b_val), , drop = FALSE]

  if (nrow(dat_ab) == 0) {
    stop("No observations in a/b dataset for TMLE TTE.")
  }

  dat_ab$A_bin <- as.numeric(as.character(dat_ab$treatment) == a_val)

  if (length(unique(stats::na.omit(dat_ab$A_bin))) < 2) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      CI = c(NA_real_, NA_real_),
      nuisance_type = nuisance_type
    ))
  }

  tmle_ate_binary(
    data = dat_ab,
    outcome = outcome,
    treat = "A_bin",
    covars = covars,
    nuisance_type = nuisance_type,
    clamp = clamp
  )
}

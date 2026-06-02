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




# ------------------------------------------------------------
# TMLE method 1: traditional contrast-of-contrasts
# ------------------------------------------------------------

tmle_traditional <- function(data,
                             nuisance_type = c("simple", "flexible"),
                             a_val = "a",
                             b_val = "b",
                             c1_val = "sc1",
                             c2_val = "sc2",
                             covars = c("L1", "L2"),
                             outcome = "y",
                             clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  data <- .tmle_safe_country(data)

  dat_s1 <- data[
    data$s == 1 & data$treatment %in% c(a_val, c1_val),
    ,
    drop = FALSE
  ]

  dat_s2 <- data[
    data$s == 2 & data$treatment %in% c(b_val, c2_val),
    ,
    drop = FALSE
  ]

  if (nrow(dat_s1) == 0) stop("No observations in trial 1 for TMLE traditional.")
  if (nrow(dat_s2) == 0) stop("No observations in trial 2 for TMLE traditional.")

  dat_s1$A_bin <- as.numeric(as.character(dat_s1$treatment) == a_val)
  dat_s2$A_bin <- as.numeric(as.character(dat_s2$treatment) == b_val)

  res1 <- tmle_ate_binary(
    data = dat_s1,
    outcome = outcome,
    treat = "A_bin",
    covars = covars,
    nuisance_type = nuisance_type,
    clamp = clamp
  )

  res2 <- tmle_ate_binary(
    data = dat_s2,
    outcome = outcome,
    treat = "A_bin",
    covars = covars,
    nuisance_type = nuisance_type,
    clamp = clamp
  )

  psi <- res1$estimate - res2$estimate
  se <- sqrt(res1$se^2 + res2$se^2)
  CI <- psi + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimate = as.numeric(psi),
    se = as.numeric(se),
    CI = as.numeric(CI),
    nuisance_type = nuisance_type
  )
}



# ------------------------------------------------------------
# Internal helpers for cross-trial TMLE direct / theta / indirect
# ------------------------------------------------------------
.align_newdata_to_glm <- function(fit, newdata) {
  newdata <- as.data.frame(newdata)

  if (!is.null(fit$xlevels)) {
    for (v in names(fit$xlevels)) {
      if (v %in% names(newdata)) {
        newdata[[v]] <- factor(newdata[[v]], levels = fit$xlevels[[v]])
      }
    }
  }

  newdata
}

.tmle_fit_binary_regression <- function(data,
                                        outcome,
                                        covariates,
                                        nuisance_type = c("simple", "flexible"),
                                        clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  data <- .tmle_safe_country(data)

  if (nuisance_type == "simple") {
    form <- stats::as.formula(
      paste(outcome, "~", paste(covariates, collapse = " + "))
    )
    data <- .tmle_safe_country(data)
    fit <- stats::glm(
      form,
      data = data,
      family = stats::binomial()
    )

    return(list(
      predict = function(newdata) {
        newdata <- .tmle_safe_country(newdata)
        newdata <- .align_newdata_to_glm(fit, newdata)

        pred <- try(
          stats::predict(fit, newdata = newdata, type = "response"),
          silent = TRUE
        )

        if (inherits(pred, "try-error")) {
          pred <- rep(mean(data[[outcome]], na.rm = TRUE), nrow(newdata))
        }

        pred <- as.numeric(pred)

        if (anyNA(pred) || any(!is.finite(pred))) {
          fallback <- mean(data[[outcome]], na.rm = TRUE)
          pred[is.na(pred) | !is.finite(pred)] <- fallback
        }

        .tmle_clip(pred, clamp)
      }
    ))
  }

  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' is required for nuisance_type = 'flexible'.")
  }

  learners <- .tmle_make_learners(
    nuisance_type = "flexible",
    outcome_family = stats::binomial()
  )

  task <- sl3::make_sl3_Task(
    data = data,
    outcome = outcome,
    covariates = covariates
  )

  fit <- learners$Q$train(task)

  list(
    predict = function(newdata) {
      newdata <- .tmle_safe_country(newdata)

      if (!(outcome %in% names(newdata))) {
        newdata[[outcome]] <- 0
      }

      task_new <- sl3::make_sl3_Task(
        data = newdata,
        outcome = outcome,
        covariates = covariates
      )

      .tmle_clip(fit$predict(task_new), clamp)
    }
  )
}

.tmle_fit_arm_mu <- function(data,
                             outcome = "y",
                             covariates = c("L1", "L2"),
                             s_value,
                             a_value,
                             nuisance_type = c("simple", "flexible"),
                             clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)

  train <- data[data$s == s_value & as.character(data$treatment) == a_value, , drop = FALSE]

  if (nrow(train) < 10) {
    stop("Too few observations in arm s=", s_value, ", treatment=", a_value)
  }

  .tmle_fit_binary_regression(
    data = train,
    outcome = outcome,
    covariates = covariates,
    nuisance_type = nuisance_type,
    clamp = clamp
  )
}

.tmle_fit_trial_g <- function(data,
                              covariates = c("L1", "L2"),
                              nuisance_type = c("simple", "flexible"),
                              clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  tmp <- data
  tmp$.S1 <- as.integer(tmp$s == 1)

  .tmle_fit_binary_regression(
    data = tmp,
    outcome = ".S1",
    covariates = covariates,
    nuisance_type = nuisance_type,
    clamp = clamp
  )
}

.tmle_fit_pi <- function(data,
                         covariates = c("L1", "L2"),
                         s_value,
                         z_value,
                         nuisance_type = c("simple", "flexible"),
                         known_pi = NULL,
                         clamp = c(1e-6, 1 - 1e-6)) {
  key <- paste0(s_value, ":", z_value)

  if (!is.null(known_pi)) {
    if (is.null(names(known_pi)) || !(key %in% names(known_pi))) {
      stop("known_pi must contain a named entry for '", key, "'.")
    }

    val <- as.numeric(known_pi[[key]])

    return(list(
      predict = function(newdata) {
        rep(.tmle_clip(val, clamp), nrow(newdata))
      }
    ))
  }

  train <- data[data$s == s_value, , drop = FALSE]
  tmp <- train
  tmp$.Az <- as.integer(as.character(tmp$treatment) == z_value)

  .tmle_fit_binary_regression(
    data = tmp,
    outcome = ".Az",
    covariates = covariates,
    nuisance_type = nuisance_type,
    clamp = clamp
  )
}

.tmle_fluctuate_binary <- function(Y_obs,
                                   mu_obs,
                                   H_obs,
                                   mu_all,
                                   H_all,
                                   clamp = c(1e-6, 1 - 1e-6)) {
  mu_obs <- .tmle_clip(mu_obs, clamp)
  mu_all <- .tmle_clip(mu_all, clamp)

  fit <- try(
    stats::glm(
      Y_obs ~ -1 + H_obs,
      family = stats::binomial(),
      offset = .tmle_logit(mu_obs, clamp),
      control = stats::glm.control(maxit = 100)
    ),
    silent = TRUE
  )

  eps <- if (inherits(fit, "try-error")) {
    0
  } else {
    as.numeric(stats::coef(fit)[1])
  }

  if (!is.finite(eps)) {
    eps <- 0
  }

  mu_star_all <- .tmle_clip(
    .tmle_expit(.tmle_logit(mu_all, clamp) + eps * H_all),
    clamp
  )
  list(
    epsilon = eps,
    mu_star_all = mu_star_all
  )
}

# ------------------------------------------------------------
# TMLE method 3: direct pooled IPD A vs B
# ------------------------------------------------------------
#' TMLE estimator of the direct effect
#'
#' Estimates the direct cross-trial contrast using targeted maximum
#' likelihood estimation (TMLE).
#'
#' The estimator combines:
#' \itemize{
#'   \item outcome regressions within each trial,
#'   \item treatment assignment mechanisms within each trial,
#'   \item a trial membership mechanism,
#'   \item a targeting step based on the efficient influence function.
#' }
#'
#' @param data A data frame containing the observed data. Must include
#'   variables `y`, `treatment`, `s`, `L1`, and `L2`.
#' @param nuisance_type Nuisance estimation strategy. Either `"simple"`
#'   or `"flexible"`.
#' @param a_val Active treatment in trial 1.
#' @param b_val Active treatment in trial 2.
#' @param covars Baseline covariates used in nuisance estimation.
#' @param outcome Outcome variable name.
#' @param known_pi Optional known treatment assignment probabilities.
#' @param trim Truncation level for estimated propensity scores.
#' @param clamp Bounds used to stabilize predicted probabilities.
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{TMLE estimate of the direct effect.}
#'   \item{se}{Influence-function-based standard error.}
#'   \item{CI}{95\% Wald confidence interval.}
#'   \item{eif}{Estimated efficient influence function values.}
#'   \item{nuisance_type}{Nuisance estimation strategy used.}
#' }
#'
#' @details
#' Outcome regressions are first estimated separately within each trial.
#' These initial estimates are then updated through a logistic fluctuation
#' step targeting the efficient influence function corresponding to the
#' direct effect parameter.
#'
#' Standard errors are obtained from the empirical variance of the estimated
#' efficient influence function.
#'
#' @export
tmle_direct <- function(data,
                        nuisance_type = c("simple", "flexible"),
                        a_val = "a",
                        b_val = "b",
                        covars = c("L1", "L2"),
                        outcome = "y",
                        known_pi = NULL,
                        trim = 0.01,
                        clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  data <- .tmle_safe_country(data)

  n <- nrow(data)

  mu1a_fit <- .tmle_fit_arm_mu(
    data = data,
    outcome = outcome,
    covariates = covars,
    s_value = 1,
    a_value = a_val,
    nuisance_type = nuisance_type,
    clamp = clamp
  )

  mu2b_fit <- .tmle_fit_arm_mu(
    data = data,
    outcome = outcome,
    covariates = covars,
    s_value = 2,
    a_value = b_val,
    nuisance_type = nuisance_type,
    clamp = clamp
  )

  g_fit <- .tmle_fit_trial_g(
    data = data,
    covariates = covars,
    nuisance_type = nuisance_type,
    clamp = clamp
  )

  pi1a_fit <- .tmle_fit_pi(
    data = data,
    covariates = covars,
    s_value = 1,
    z_value = a_val,
    nuisance_type = nuisance_type,
    known_pi = known_pi,
    clamp = clamp
  )

  pi2b_fit <- .tmle_fit_pi(
    data = data,
    covariates = covars,
    s_value = 2,
    z_value = b_val,
    nuisance_type = nuisance_type,
    known_pi = known_pi,
    clamp = clamp
  )

  mu1a_all <- mu1a_fit$predict(data)
  mu2b_all <- mu2b_fit$predict(data)

  g_all <- .tmle_clip(g_fit$predict(data), c(trim, 1 - trim))
  pi1a_all <- .tmle_clip(pi1a_fit$predict(data), c(trim, 1 - trim))
  pi2b_all <- .tmle_clip(pi2b_fit$predict(data), c(trim, 1 - trim))

  H1a_all <- 1 / (pi1a_all * g_all)
  H2b_all <- 1 / (pi2b_all * (1 - g_all))

  idx1a <- data$s == 1 & as.character(data$treatment) == a_val
  idx2b <- data$s == 2 & as.character(data$treatment) == b_val

  fl1a <- .tmle_fluctuate_binary(
    Y_obs = data[[outcome]][idx1a],
    mu_obs = mu1a_all[idx1a],
    H_obs = H1a_all[idx1a],
    mu_all = mu1a_all,
    H_all = H1a_all,
    clamp = clamp
  )

  fl2b <- .tmle_fluctuate_binary(
    Y_obs = data[[outcome]][idx2b],
    mu_obs = mu2b_all[idx2b],
    H_obs = H2b_all[idx2b],
    mu_all = mu2b_all,
    H_all = H2b_all,
    clamp = clamp
  )

  mu1a_star <- fl1a$mu_star_all
  mu2b_star <- fl2b$mu_star_all

  psi <- mean(mu1a_star - mu2b_star)

  D <- numeric(n)
  D <- D + as.numeric(idx1a) * H1a_all * (data[[outcome]] - mu1a_star)
  D <- D - as.numeric(idx2b) * H2b_all * (data[[outcome]] - mu2b_star)
  D <- D + (mu1a_star - mu2b_star) - psi

  se <- sqrt(.tmle_safe_var(D) / n)
  CI <- psi + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimate = as.numeric(psi),
    se = as.numeric(se),
    CI = as.numeric(CI),
    eif = D,
    nuisance_type = nuisance_type
  )
}


# ------------------------------------------------------------
# Internal theta TMLE for indirect decomposition
# ------------------------------------------------------------
#' TMLE estimator of the pathway-specific contrast (\eqn{\theta})
#'
#' Estimates the pathway-specific contrast
#'
#' \deqn{
#' \theta =
#' \{E[Y^a - Y^{c_1}]\}_{S=1 \rightarrow 2}
#' -
#' \{E[Y^b - Y^{c_2}]\}_{S=2 \rightarrow 1}.
#' }
#'
#' using targeted maximum likelihood estimation (TMLE).
#'
#' @param data A data frame containing the observed data.
#' @param nuisance_type Nuisance estimation strategy.
#' @param a_val Active treatment in trial 1.
#' @param b_val Active treatment in trial 2.
#' @param c1_val Standard-of-care treatment in trial 1.
#' @param c2_val Standard-of-care treatment in trial 2.
#' @param covars Baseline covariates used in nuisance estimation.
#' @param outcome Outcome variable name.
#' @param known_pi Optional known treatment probabilities.
#' @param trim Truncation level for estimated propensity scores.
#' @param clamp Bounds used to stabilize predicted probabilities.
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{TMLE estimate of \eqn{\theta}.}
#'   \item{se}{Influence-function-based standard error.}
#'   \item{CI}{95\% Wald confidence interval.}
#'   \item{eif}{Estimated efficient influence function values.}
#'   \item{nuisance_type}{Nuisance estimation strategy used.}
#' }
#'
#' @details
#' Four treatment-specific outcome regressions are estimated and targeted:
#' \eqn{\mu_{1,a}},
#' \eqn{\mu_{1,c_1}},
#' \eqn{\mu_{2,b}},
#' and
#' \eqn{\mu_{2,c_2}}.
#'
#' The resulting TMLE estimator solves the empirical efficient influence
#' function estimating equation for the parameter \eqn{\theta}.
#'
#' @export
tmle_theta <- function(data,
                       nuisance_type = c("simple", "flexible"),
                       a_val = "a",
                       b_val = "b",
                       c1_val = "sc1",
                       c2_val = "sc2",
                       covars = c("L1", "L2"),
                       outcome = "y",
                       known_pi = NULL,
                       trim = 0.01,
                       clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)
  data <- .tmle_safe_country(data)

  n <- nrow(data)

  mu1a_fit <- .tmle_fit_arm_mu(data, outcome, covars, 1, a_val, nuisance_type, clamp)
  mu1c_fit <- .tmle_fit_arm_mu(data, outcome, covars, 1, c1_val, nuisance_type, clamp)
  mu2b_fit <- .tmle_fit_arm_mu(data, outcome, covars, 2, b_val, nuisance_type, clamp)
  mu2c_fit <- .tmle_fit_arm_mu(data, outcome, covars, 2, c2_val, nuisance_type, clamp)

  g_fit <- .tmle_fit_trial_g(data, covars, nuisance_type, clamp)

  pi1a_fit <- .tmle_fit_pi(data, covars, 1, a_val, nuisance_type, known_pi, clamp)
  pi1c_fit <- .tmle_fit_pi(data, covars, 1, c1_val, nuisance_type, known_pi, clamp)
  pi2b_fit <- .tmle_fit_pi(data, covars, 2, b_val, nuisance_type, known_pi, clamp)
  pi2c_fit <- .tmle_fit_pi(data, covars, 2, c2_val, nuisance_type, known_pi, clamp)

  mu1a <- mu1a_fit$predict(data)
  mu1c <- mu1c_fit$predict(data)
  mu2b <- mu2b_fit$predict(data)
  mu2c <- mu2c_fit$predict(data)

  g_all <- .tmle_clip(g_fit$predict(data), c(trim, 1 - trim))

  pi1a <- .tmle_clip(pi1a_fit$predict(data), c(trim, 1 - trim))
  pi1c <- .tmle_clip(pi1c_fit$predict(data), c(trim, 1 - trim))
  pi2b <- .tmle_clip(pi2b_fit$predict(data), c(trim, 1 - trim))
  pi2c <- .tmle_clip(pi2c_fit$predict(data), c(trim, 1 - trim))

  H1a <- 1 / (pi1a * g_all)
  H1c <- 1 / (pi1c * g_all)
  H2b <- 1 / (pi2b * (1 - g_all))
  H2c <- 1 / (pi2c * (1 - g_all))

  idx1a <- data$s == 1 & as.character(data$treatment) == a_val
  idx1c <- data$s == 1 & as.character(data$treatment) == c1_val
  idx2b <- data$s == 2 & as.character(data$treatment) == b_val
  idx2c <- data$s == 2 & as.character(data$treatment) == c2_val

  fl1a <- .tmle_fluctuate_binary(data[[outcome]][idx1a], mu1a[idx1a], H1a[idx1a], mu1a, H1a, clamp)
  fl1c <- .tmle_fluctuate_binary(data[[outcome]][idx1c], mu1c[idx1c], H1c[idx1c], mu1c, H1c, clamp)
  fl2b <- .tmle_fluctuate_binary(data[[outcome]][idx2b], mu2b[idx2b], H2b[idx2b], mu2b, H2b, clamp)
  fl2c <- .tmle_fluctuate_binary(data[[outcome]][idx2c], mu2c[idx2c], H2c[idx2c], mu2c, H2c, clamp)

  mu1a_star <- fl1a$mu_star_all
  mu1c_star <- fl1c$mu_star_all
  mu2b_star <- fl2b$mu_star_all
  mu2c_star <- fl2c$mu_star_all

  theta <- mean((mu1a_star - mu1c_star) - (mu2b_star - mu2c_star))

  D <- numeric(n)
  D <- D + as.numeric(idx1a) * H1a * (data[[outcome]] - mu1a_star)
  D <- D - as.numeric(idx1c) * H1c * (data[[outcome]] - mu1c_star)
  D <- D - as.numeric(idx2b) * H2b * (data[[outcome]] - mu2b_star)
  D <- D + as.numeric(idx2c) * H2c * (data[[outcome]] - mu2c_star)
  D <- D + ((mu1a_star - mu1c_star) - (mu2b_star - mu2c_star)) - theta

  se <- sqrt(.tmle_safe_var(D) / n)
  CI <- theta + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimate = as.numeric(theta),
    se = as.numeric(se),
    CI = as.numeric(CI),
    eif = D,
    nuisance_type = nuisance_type
  )
}

# ------------------------------------------------------------
# TMLE method 4: indirect = theta + phi
# ------------------------------------------------------------
#' TMLE estimator of the indirect effect
#'
#' Estimates the indirect effect using the decomposition
#'
#' \deqn{
#' \delta = \theta + \phi,
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\theta} is estimated by [tmle_theta()],
#'   \item \eqn{\phi} is estimated by [tmle_direct()]
#'     evaluated at the two standard-of-care treatments.
#' }
#'
#' @param data A data frame containing the observed data.
#' @param nuisance_type Nuisance estimation strategy.
#' @param a_val Active treatment in trial 1.
#' @param b_val Active treatment in trial 2.
#' @param c1_val Standard-of-care treatment in trial 1.
#' @param c2_val Standard-of-care treatment in trial 2.
#' @param covars Baseline covariates used in nuisance estimation.
#' @param outcome Outcome variable name.
#' @param known_pi Optional known treatment probabilities.
#' @param trim Truncation level for estimated propensity scores.
#' @param clamp Bounds used to stabilize predicted probabilities.
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{TMLE estimate of the indirect effect.}
#'   \item{se}{Influence-function-based standard error.}
#'   \item{CI}{95\% Wald confidence interval.}
#'   \item{eif}{Estimated efficient influence function values for the
#'   indirect effect.}
#'   \item{components}{Estimates and inference for the \eqn{\theta}
#'   and \eqn{\phi} components.}
#'   \item{nuisance_type}{Nuisance estimation strategy used.}
#' }
#'
#' @details
#' The indirect effect estimator is constructed as
#'
#' \deqn{
#' \hat\delta =
#' \hat\theta_{\mathrm{TMLE}}
#' +
#' \hat\phi_{\mathrm{TMLE}}.
#' }
#'
#' Because both components are estimated from the same sample,
#' inference is based on the combined efficient influence function
#'
#' \deqn{
#' D_{\delta}
#' =
#' D_{\theta}
#' +
#' D_{\phi},
#' }
#'
#' which automatically accounts for the covariance between
#' \eqn{\hat\theta} and \eqn{\hat\phi}.
#'
#' @seealso [tmle_theta()], [tmle_direct()]
#'
#' @export
tmle_indirect <- function(data,
                          nuisance_type = c("simple", "flexible"),
                          a_val = "a",
                          b_val = "b",
                          c1_val = "sc1",
                          c2_val = "sc2",
                          covars = c("L1", "L2"),
                          outcome = "y",
                          known_pi = NULL,
                          trim = 0.01,
                          clamp = c(1e-6, 1 - 1e-6)) {
  nuisance_type <- match.arg(nuisance_type)

  theta_fit <- tmle_theta(
    data = data,
    nuisance_type = nuisance_type,
    a_val = a_val,
    b_val = b_val,
    c1_val = c1_val,
    c2_val = c2_val,
    covars = covars,
    outcome = outcome,
    known_pi = known_pi,
    trim = trim,
    clamp = clamp
  )

  phi_fit <- tmle_direct(
    data = data,
    nuisance_type = nuisance_type,
    a_val = c1_val,
    b_val = c2_val,
    covars = covars,
    outcome = outcome,
    known_pi = known_pi,
    trim = trim,
    clamp = clamp
  )

  # Inference for the indirect TMLE should use the EIF for the
  # combined estimand psi = theta + phi. Since theta and phi are
  # estimated on the same data, their covariance is generally nonzero.
  # Therefore, use D_indirect = D_theta + D_phi to compute the standard error.
  psi <- theta_fit$estimate + phi_fit$estimate
  eif <- theta_fit$eif + phi_fit$eif
  se <- sqrt(.tmle_safe_var(eif) / nrow(data))
  CI <- psi + stats::qnorm(c(0.025, 0.975)) * se

  list(
    estimate = as.numeric(psi),
    se = as.numeric(se),
    CI = as.numeric(CI),
    eif = eif,
    components = list(
      theta = list(
        estimate = theta_fit$estimate,
        se = theta_fit$se,
        CI = theta_fit$CI
      ),
      phi = list(
        estimate = phi_fit$estimate,
        se = phi_fit$se,
        CI = phi_fit$CI
      )
    ),
    nuisance_type = nuisance_type
  )
}

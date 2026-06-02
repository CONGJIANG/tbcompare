# ============================================================
# 2. G-computation helpers
# ============================================================

gcomp_trial <- function(data,
                        s_val,
                        treat_active,
                        treat_control,
                        covars = c("L1", "L2"),
                        outcome = "y",
                        nboot = 100,
                        quiet = TRUE) {
  dat_s <- data[data$s == s_val, , drop = FALSE]

  if (nrow(dat_s) == 0) {
    stop("No rows for that trial value.")
  }

  dat_s <- .safe_factor_country(dat_s)
  dat_s$treatment <- factor(dat_s$treatment)

  formula_str <- paste0(outcome, " ~ treatment + ", paste(covars, collapse = " + "))

  fit <- stats::glm(
    stats::as.formula(formula_str),
    data = dat_s,
    family = stats::binomial()
  )

  if (!quiet) {
    cat("Outcome model for s =", s_val, "\n")
    print(summary(fit)$coefficients)
  }

  dat_active <- dat_s
  dat_control <- dat_s

  dat_active$treatment <- factor(treat_active, levels = levels(dat_s$treatment))
  dat_control$treatment <- factor(treat_control, levels = levels(dat_s$treatment))

  p_active <- stats::predict(fit, newdata = dat_active, type = "response")
  p_control <- stats::predict(fit, newdata = dat_control, type = "response")

  eff_act <- mean(p_active)
  eff_ctr <- mean(p_control)
  ate <- mean(p_active - p_control)

  boot_est <- numeric(0)

  if (nboot > 0) {
    boot_est <- numeric(nboot)
    n <- nrow(dat_s)

    for (b in seq_len(nboot)) {
      idx <- sample.int(n, n, replace = TRUE)
      ds <- dat_s[idx, , drop = FALSE]
      ds$treatment <- factor(ds$treatment, levels = levels(dat_s$treatment))

      fitb <- try(
        stats::glm(stats::as.formula(formula_str), data = ds, family = stats::binomial()),
        silent = TRUE
      )

      if (inherits(fitb, "try-error")) {
        boot_est[b] <- NA_real_
        next
      }

      dab <- ds
      dcb <- ds
      dab$treatment <- factor(treat_active, levels = levels(dat_s$treatment))
      dcb$treatment <- factor(treat_control, levels = levels(dat_s$treatment))

      p_a <- try(stats::predict(fitb, newdata = dab, type = "response"), silent = TRUE)
      p_c <- try(stats::predict(fitb, newdata = dcb, type = "response"), silent = TRUE)

      if (inherits(p_a, "try-error") || inherits(p_c, "try-error")) {
        boot_est[b] <- NA_real_
      } else {
        boot_est[b] <- mean(p_a - p_c)
      }
    }

    boot_est <- boot_est[!is.na(boot_est)]
  }

  ci <- if (length(boot_est) > 0) {
    stats::quantile(boot_est, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    c(NA_real_, NA_real_)
  }

  se_boot <- if (length(boot_est) > 0) stats::sd(boot_est) else NA_real_

  list(
    eff_act = eff_act,
    eff_ctr = eff_ctr,
    ate = ate,
    se_boot = se_boot,
    CI = ci,
    boot_est = boot_est,
    fit = fit,
    n = nrow(dat_s)
  )
}

#' Cross-trial g-computation prediction helper
#'
#' Fits an outcome regression model within one trial and predicts the
#' counterfactual mean outcomes under two treatment values in the other trial.
#'
#' This is an internal helper used by the g-computation estimators.
#'
#' @param data A pooled trial data frame.
#' @param s_val Trial indicator value used to fit the outcome model.
#' @param treat_active Treatment label for the active regimen.
#' @param treat_control Treatment label for the comparator regimen.
#' @param covars Character vector of covariate names included in the outcome model.
#' @param outcome Name of the binary outcome variable.
#' @param quiet Logical; if `FALSE`, prints the fitted outcome model coefficients.
#'
#' @return A list with predicted mean outcomes under the active and comparator
#' treatment values in the opposite trial population.
#'
#' @export
gcomp_cross <- function(data,
                        s_val,
                        treat_active,
                        treat_control,
                        covars = c("L1", "L2"),
                        outcome = "y",
                        quiet = TRUE) {
  dat_s <- data[data$s == s_val, , drop = FALSE]

  if (nrow(dat_s) == 0) {
    stop("No rows for that trial value.")
  }

  dat_s <- .safe_factor_country(dat_s)
  dat_s$treatment <- factor(dat_s$treatment)

  formula_str <- paste0(outcome, " ~ treatment + ", paste(covars, collapse = " + "))

  fit <- stats::glm(
    stats::as.formula(formula_str),
    data = dat_s,
    family = stats::binomial()
  )

  if (!quiet) {
    cat("Outcome model for s =", s_val, "\n")
    print(summary(fit)$coefficients)
  }

  other_data <- data[data$s != s_val, , drop = FALSE]
  other_data <- .safe_factor_country(other_data)

  dat_active <- other_data
  dat_control <- other_data

  dat_active$treatment <- factor(treat_active, levels = levels(dat_s$treatment))
  dat_control$treatment <- factor(treat_control, levels = levels(dat_s$treatment))

  p_active <- stats::predict(fit, newdata = dat_active, type = "response")
  p_control <- stats::predict(fit, newdata = dat_control, type = "response")

  list(
    eff_act = mean(p_active),
    eff_ctr = mean(p_control)
  )
}

g_com_traditional <- function(data,
                              nboot = 100,
                              ncores = 1,
                              seed = NULL,
                              quiet = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  covars <- c("L1")

  res1 <- gcomp_trial(
    data = data,
    s_val = 1,
    treat_active = "a",
    treat_control = "sc1",
    covars = covars,
    outcome = "y",
    nboot = nboot,
    quiet = quiet
  )

  res2 <- gcomp_trial(
    data = data,
    s_val = 2,
    treat_active = "b",
    treat_control = "sc2",
    covars = covars,
    outcome = "y",
    nboot = nboot,
    quiet = quiet
  )

  ate_diff <- res1$ate - res2$ate

  if (length(res1$boot_est) > 0 && length(res2$boot_est) > 0) {
    min_len <- min(length(res1$boot_est), length(res2$boot_est))
    boot_diff <- res1$boot_est[seq_len(min_len)] - res2$boot_est[seq_len(min_len)]
    se_diff <- stats::sd(boot_diff)
    ci_diff <- stats::quantile(boot_diff, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    boot_diff <- numeric(0)
    se_diff <- NA_real_
    ci_diff <- c(NA_real_, NA_real_)
  }

  list(
    res1 = res1,
    res2 = res2,
    ate_diff = ate_diff,
    se_diff = se_diff,
    CI_diff = ci_diff,
    boot_est = boot_diff
  )
}

g_com_direct <- function(data,
                         eff = c("trt", "ctrl"),
                         nboot = 100,
                         ncores = 1,
                         seed = NULL,
                         quiet = TRUE) {
  eff <- match.arg(eff)
  covars <- c("L1")

  compute_effect <- function(dat) {
    res_try <- try({
      ate1 <- gcomp_trial(
        dat,
        s_val = 1,
        treat_active = "a",
        treat_control = "sc1",
        covars = covars,
        outcome = "y",
        nboot = 0,
        quiet = quiet
      )

      ate2 <- gcomp_trial(
        dat,
        s_val = 2,
        treat_active = "b",
        treat_control = "sc2",
        covars = covars,
        outcome = "y",
        nboot = 0,
        quiet = quiet
      )

      s1_cross <- gcomp_cross(
        dat,
        s_val = 2,
        treat_active = "b",
        treat_control = "sc2",
        covars = covars,
        outcome = "y",
        quiet = quiet
      )

      s2_cross <- gcomp_cross(
        dat,
        s_val = 1,
        treat_active = "a",
        treat_control = "sc1",
        covars = covars,
        outcome = "y",
        quiet = quiet
      )

      p1 <- mean(dat$s == 1)
      p2 <- 1 - p1

      if (eff == "trt") {
        (ate1$eff_act - s1_cross$eff_act) * p1 +
          (s2_cross$eff_act - ate2$eff_act) * p2
      } else {
        (ate1$eff_ctr - s1_cross$eff_ctr) * p1 +
          (s2_cross$eff_ctr - ate2$eff_ctr) * p2
      }
    }, silent = TRUE)

    if (inherits(res_try, "try-error")) NA_real_ else as.numeric(res_try)
  }

  res <- compute_effect(data)

  n <- nrow(data)
  if (!is.null(seed)) set.seed(seed)

  boot_est <- numeric(0)

  if (nboot > 0) {
    if (ncores > 1 && requireNamespace("parallel", quietly = TRUE)) {
      idx_list <- replicate(nboot, sample.int(n, n, replace = TRUE), simplify = FALSE)
      boot_est <- unlist(parallel::mclapply(
        idx_list,
        function(idx) compute_effect(data[idx, , drop = FALSE]),
        mc.cores = ncores
      ))
    } else {
      boot_est <- numeric(nboot)
      for (b in seq_len(nboot)) {
        idx <- sample.int(n, n, replace = TRUE)
        boot_est[b] <- compute_effect(data[idx, , drop = FALSE])
      }
    }
    boot_est <- boot_est[!is.na(boot_est)]
  }

  ci <- if (length(boot_est) > 0) {
    stats::quantile(boot_est, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    c(NA_real_, NA_real_)
  }

  se_boot <- if (length(boot_est) > 0) stats::sd(boot_est) else NA_real_

  list(
    eff_diff = res,
    se_boot = se_boot,
    CI = ci,
    boot_est = boot_est,
    n = n
  )
}

g_com_tte <- function(data,
                      nboot = 100,
                      seed = NULL,
                      quiet = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  datasub_trt <- data[data$treatment %in% c("a", "b"), , drop = FALSE]
  datasub_trt$s <- 3

  res <- gcomp_trial(
    data = datasub_trt,
    s_val = 3,
    treat_active = "a",
    treat_control = "b",
    covars = c("L1"),
    outcome = "y",
    nboot = nboot,
    quiet = quiet
  )

  list(
    ate = res$ate,
    se_boot = res$se_boot,
    CI = res$CI,
    boot_est = res$boot_est,
    raw = res
  )
}

g_com_indirt <- function(data,
                         nboot = 100,
                         ncores = 1,
                         seed = NULL,
                         quiet = TRUE) {
  covars <- c("L1")

  compute_theta <- function(dat) {
    ate1 <- gcomp_trial(
      dat,
      s_val = 1,
      treat_active = "a",
      treat_control = "sc1",
      covars = covars,
      outcome = "y",
      nboot = 0,
      quiet = quiet
    )

    ate2 <- gcomp_trial(
      dat,
      s_val = 2,
      treat_active = "b",
      treat_control = "sc2",
      covars = covars,
      outcome = "y",
      nboot = 0,
      quiet = quiet
    )

    s1_cross <- gcomp_cross(
      dat,
      s_val = 2,
      treat_active = "b",
      treat_control = "sc2",
      covars = covars,
      outcome = "y",
      quiet = quiet
    )

    s2_cross <- gcomp_cross(
      dat,
      s_val = 1,
      treat_active = "a",
      treat_control = "sc1",
      covars = covars,
      outcome = "y",
      quiet = quiet
    )

    p1 <- mean(dat$s == 1)
    p2 <- 1 - p1

    theta_1 <- ate1$ate * p1 + (s2_cross$eff_act - s2_cross$eff_ctr) * p2
    theta_2 <- (s1_cross$eff_act - s1_cross$eff_ctr) * p1 + ate2$ate * p2

    theta_1 - theta_2
  }

  compute_phi <- function(dat) {
    res_phi <- g_com_direct(
      data = dat,
      eff = "ctrl",
      nboot = 0,
      ncores = 1,
      seed = NULL,
      quiet = quiet
    )
    res_phi$eff_diff
  }

  compute_indirect <- function(dat) {
    res_try <- try({
      theta <- compute_theta(dat)
      phi <- compute_phi(dat)
      list(theta = theta, phi = phi, indirect = theta + phi)
    }, silent = TRUE)

    if (inherits(res_try, "try-error")) {
      list(theta = NA_real_, phi = NA_real_, indirect = NA_real_)
    } else {
      res_try
    }
  }

  est <- compute_indirect(data)
  ind_est <- est$indirect

  n <- nrow(data)
  if (!is.null(seed)) set.seed(seed)

  boot_est <- numeric(0)

  if (nboot > 0) {
    if (ncores > 1 && requireNamespace("parallel", quietly = TRUE)) {
      idx_list <- replicate(nboot, sample.int(n, n, replace = TRUE), simplify = FALSE)
      boot_list <- parallel::mclapply(
        idx_list,
        function(idx) compute_indirect(data[idx, , drop = FALSE])$indirect,
        mc.cores = ncores
      )
      boot_est <- unlist(boot_list)
    } else {
      boot_est <- numeric(nboot)
      for (b in seq_len(nboot)) {
        idx <- sample.int(n, n, replace = TRUE)
        boot_est[b] <- compute_indirect(data[idx, , drop = FALSE])$indirect
      }
    }

    boot_est <- boot_est[!is.na(boot_est)]
  }

  ci <- if (length(boot_est) > 0) {
    stats::quantile(boot_est, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    c(NA_real_, NA_real_)
  }

  se_boot <- if (length(boot_est) > 0) stats::sd(boot_est) else NA_real_

  list(
    theta = est$theta,
    phi = est$phi,
    indirect = ind_est,
    se_boot = se_boot,
    CI = ci,
    boot_est = boot_est,
    n = n
  )
}

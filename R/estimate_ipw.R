
# ============================================================
# 1. IPW unified estimator
# ============================================================
#' Unified IPW estimators for cross-trial comparisons
#'
#' Computes inverse probability weighted estimators for four cross-trial
#' comparison strategies:
#'
#' * Traditional indirect comparison
#' * Observational comparison
#' * Direct pooled IPD comparison
#' * Indirect decomposition comparison
#'
#' @param dataset A data frame containing the pooled trial data.
#' @param a_val Treatment label for regimen A.
#' @param b_val Treatment label for regimen B.
#' @param c1_val Standard-of-care comparator in trial 1.
#' @param c2_val Standard-of-care comparator in trial 2.
#' @param clamp Bounds used to truncate estimated probabilities.
#' @param compute_boot Logical; whether bootstrap inference should be computed.
#' @param B Number of bootstrap replicates.
#'
#' @return A list containing estimates, standard errors,
#' confidence intervals, and influence-function quantities for all
#' implemented comparison strategies.
#'
#' @seealso [xtrial()]
#'
#' @export

estimate_ipw_unified <- function(dataset,
                                 a_val = "a",
                                 b_val = "b",
                                 c1_val = "sc1",
                                 c2_val = "sc2",
                                 clamp = c(1e-6, 1 - 1e-6),
                                 compute_boot = FALSE,
                                 B = 200) {

  dataset <- as.data.frame(dataset)
  dataset <- .safe_factor_country(dataset)

  n <- nrow(dataset)
  Y <- dataset$y
  A <- as.character(dataset$treatment)
  S <- dataset$s

  clamp_prob <- function(x) .clamp_prob(x, clamp)

  # -----------------------------
  # Helper functions
  # -----------------------------
  ipw_term <- function(ind, y, ps, transport = 1) {
    ind * y * transport / ps
  }

  direct_pool <- function(I_left, I_right,
                          ps_left, ps_right,
                          y, delta, delta_inv, n) {

    g_left  <- I_left  * y / ps_left  * (1 + delta)
    g_right <- I_right * y / ps_right * (1 + delta_inv)

    psi <- mean(g_left - g_right)

    IF <- (g_left - g_right) - psi

    se_plugin <- stats::sd(IF, na.rm = TRUE) / sqrt(n)

    CI_plugin <- psi + stats::qnorm(c(0.025, 0.975)) * se_plugin

    list(psi = psi, IF = IF,
         se_plugin = se_plugin,
         CI_plugin = CI_plugin)
  }

  # -----------------------------
  # Propensity models
  # -----------------------------
  ds1 <- dataset[S == 1, , drop = FALSE]
  ds2 <- dataset[S == 2, , drop = FALSE]

  ds1$treat_bin <- as.numeric(ds1$treatment == a_val)
  psa_all <- clamp_prob(predict(glm(
    treat_bin ~ L1 + L2 + factor(country),
    data = ds1, family = binomial()
  ), newdata = dataset, type = "response"))

  ds1$treat_bin <- as.numeric(ds1$treatment == c1_val)
  psc1_all <- clamp_prob(predict(glm(
    treat_bin ~ L1 + L2 + factor(country),
    data = ds1, family = binomial()
  ), newdata = dataset, type = "response"))

  ds2$treat_bin <- as.numeric(ds2$treatment == b_val)
  psb_all <- clamp_prob(predict(glm(
    treat_bin ~ L1 + L2 + factor(country),
    data = ds2, family = binomial()
  ), newdata = dataset, type = "response"))

  ds2$treat_bin <- as.numeric(ds2$treatment == c2_val)
  psc2_all <- clamp_prob(predict(glm(
    treat_bin ~ L1 + L2 + factor(country),
    data = ds2, family = binomial()
  ), newdata = dataset, type = "response"))

  # Trial membership
  mod_S <- glm(I(s == 1) ~ L1 + L2 + factor(country),
               data = dataset, family = binomial())

  p_s1 <- clamp_prob(predict(mod_S, type = "response"))
  p_s2 <- 1 - p_s1

  delta <- p_s2 / p_s1
  delta_inv <- 1 / delta

  # Indicators
  I_a1  <- as.numeric(A == a_val  & S == 1)
  I_c11 <- as.numeric(A == c1_val & S == 1)
  I_b2  <- as.numeric(A == b_val  & S == 2)
  I_c22 <- as.numeric(A == c2_val & S == 2)

  # -----------------------------
  # g terms
  # -----------------------------
  g1 <- ipw_term(I_a1,  Y, psa_all)
  g2 <- ipw_term(I_c11, Y, psc1_all)
  g3 <- ipw_term(I_a1,  Y, psa_all, delta)
  g4 <- ipw_term(I_c11, Y, psc1_all, delta)

  g5 <- ipw_term(I_b2,  Y, psb_all, delta_inv)
  g6 <- ipw_term(I_c22, Y, psc2_all, delta_inv)
  g7 <- ipw_term(I_b2,  Y, psb_all)
  g8 <- ipw_term(I_c22, Y, psc2_all)

  # -----------------------------
  # 1. Traditional
  # -----------------------------
  idx_s1 <- S == 1
  idx_s2 <- S == 2

  ate1_terms <- (I_a1 * Y / psa_all) - (I_c11 * Y / psc1_all)
  ate2_terms <- (I_b2 * Y / psb_all) - (I_c22 * Y / psc2_all)

  ate1 <- mean(ate1_terms[idx_s1], na.rm = TRUE)
  ate2 <- mean(ate2_terms[idx_s2], na.rm = TRUE)

  psi_traditional <- ate1 - ate2

  IF1 <- ate1_terms[idx_s1] - ate1
  IF2 <- ate2_terms[idx_s2] - ate2

  se_traditional <- sqrt(
    stats::var(IF1, na.rm = TRUE) / sum(idx_s1) +
      stats::var(IF2, na.rm = TRUE) / sum(idx_s2)
  )

  CI_traditional <- psi_traditional +
    stats::qnorm(c(0.025, 0.975)) * se_traditional
  # -----------------------------
  # 2. TTE restricted
  # -----------------------------
  dat_ab <- dataset[A %in% c(a_val, b_val), , drop = FALSE]
  dat_ab$treat_ab <- as.numeric(dat_ab$treatment == a_val)

  if (length(unique(dat_ab$treat_ab)) < 2) {
    psi_tte_ab <- NA
    se_tte_ab <- NA
    CI_tte_ab <- c(NA, NA)
  } else {
    ps_ab <- glm(treat_ab ~ L1 + L2 + factor(country),
                 data = dat_ab, family = binomial())

    dat_ab$ps <- clamp_prob(predict(ps_ab, type = "response"))

    dat_ab$w <- ifelse(dat_ab$treat_ab == 1,
                       1 / dat_ab$ps,
                       1 / (1 - dat_ab$ps))

    mu1 <- weighted.mean(dat_ab$y[dat_ab$treat_ab == 1],
                         dat_ab$w[dat_ab$treat_ab == 1])

    mu0 <- weighted.mean(dat_ab$y[dat_ab$treat_ab == 0],
                         dat_ab$w[dat_ab$treat_ab == 0])

    psi_tte_ab <- mu1 - mu0

    IF <- with(dat_ab,
               treat_ab * w * (y - mu1) -
               (1 - treat_ab) * w * (y - mu0))

    se_tte_ab <- sd(IF) / sqrt(nrow(dat_ab))
    CI_tte_ab <- psi_tte_ab +
      qnorm(c(0.025, 0.975)) * se_tte_ab
  }

  # -----------------------------
  # 3. Direct pooled
  # -----------------------------
  dp_ab <- direct_pool(I_a1, I_b2,
                       psa_all, psb_all,
                       Y, delta, delta_inv, n)

  # -----------------------------
  # 4. Indirect
  # -----------------------------
  theta1 <- mean(g1 - g2 + g3 - g4)
  theta2 <- mean(g5 - g6 + g7 - g8)
  theta <- theta1 - theta2

  dp_c <- direct_pool(I_c11, I_c22,
                      psc1_all, psc2_all,
                      Y, delta, delta_inv, n)

  phi <- dp_c$psi
  psi_indirect <- theta + phi

  IF_theta <- ((g1 - g2 + g3 - g4) - theta1) -
              ((g5 - g6 + g7 - g8) - theta2)

  IF_indirect <- IF_theta + dp_c$IF

  se_indirect_plugin <- sd(IF_indirect) / sqrt(n)

  CI_indirect_plugin <- psi_indirect +
    qnorm(c(0.025, 0.975)) * se_indirect_plugin

  # -----------------------------
  # Bootstrap (ALL methods)
  # -----------------------------
  se_traditional_boot <- NA
  se_tte_boot         <- NA
  se_direct_boot      <- NA
  se_indirect_boot    <- NA

  if (isTRUE(compute_boot)) {

    boot_mat <- replicate(B, {

      idx <- sample(seq_len(n), replace = TRUE)

      res <- estimate_ipw_unified(
        dataset[idx, ],
        a_val, b_val, c1_val, c2_val,
        compute_boot = FALSE
      )

      c(
        res$psi_traditional,
        res$psi_tte_ab,
        res$psi_direct_pool,
        res$psi_indirect
      )
    })

    se_boot <- apply(boot_mat, 1, sd, na.rm = TRUE)

    se_traditional_boot <- se_boot[1]
    se_tte_boot         <- se_boot[2]
    se_direct_boot      <- se_boot[3]
    se_indirect_boot    <- se_boot[4]
  }

  # -----------------------------
  # Output
  # -----------------------------
  list(
    theta = theta,
    phi = phi,

    psi_traditional = psi_traditional,
    se_traditional = se_traditional,
    CI_traditional = CI_traditional,
    se_traditional_boot = se_traditional_boot,

    psi_tte_ab = psi_tte_ab,
    se_tte_ab = se_tte_ab,
    CI_tte_ab = CI_tte_ab,
    se_tte_boot = se_tte_boot,

    psi_direct_pool = dp_ab$psi,
    se_direct_pool = dp_ab$se_plugin,
    CI_direct_pool = dp_ab$CI_plugin,
    se_direct_pool_boot = se_direct_boot,

    # backward compatibility
    psi_direct_ab = dp_ab$psi,
    se_ab_direct = dp_ab$se_plugin,
    CI_ab_direct = dp_ab$CI_plugin,

    psi_indirect = psi_indirect,
    se_indirect_plugin = se_indirect_plugin,
    CI_indirect_plugin = CI_indirect_plugin,
    se_indirect_boot = se_indirect_boot,

    boot_used = compute_boot
  )
}

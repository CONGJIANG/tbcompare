#' Main interface for TB comparison estimators
#'
#' @param data A data frame containing the trial indicator, treatment, outcome,
#'   and baseline covariates.
#' @param approach Identification strategy. One of "traditional",
#'   "observational", "direct", or "indirect".
#' @param estimator Estimation method. One of "ipw", "gcomp", "onestep",
#'   or "tmle".
#' @param nboot Number of bootstrap samples for bootstrap-based estimators.
#' @param nuisance_type Nuisance estimation type. Currently mainly used for
#'   one-step/EIF and TMLE estimators; one of "simple" or "flexible".
#'
#' @return A `tbcompare` object.
#'
#' @export
compare_trials <- function(
  data,
  approach = c("traditional", "observational", "direct", "indirect"),
  estimator = c("ipw", "gcomp", "onestep", "tmle"),
  nboot = 100,
  nuisance_type = c("simple", "flexible")
) {
  approach <- match.arg(approach)
  estimator <- match.arg(estimator)
  nuisance_type <- match.arg(nuisance_type)

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (estimator == "ipw") {
    res <- estimate_ipw_unified(data)

    out <- switch(
      approach,
      traditional = list(
        estimate = res$psi_traditional,
        se = res$se_traditional,
        ci = res$CI_traditional
      ),
      observational = list(
        estimate = res$psi_tte_ab,
        se = res$se_tte_ab,
        ci = res$CI_tte_ab
      ),
      direct = list(
        estimate = res$psi_direct_pool,
        se = res$se_direct_pool,
        ci = res$CI_direct_pool
      ),
      indirect = list(
        estimate = res$psi_indirect,
        se = res$se_indirect_plugin,
        ci = res$CI_indirect_plugin
      )
    )

    return(new_tbcompare_result(
      estimate = out$estimate,
      se = out$se,
      ci = out$ci,
      approach = approach,
      estimator = estimator,
      n = nrow(data),
      raw = res
    ))
  }

  if (estimator == "gcomp") {
    fit <- switch(
      approach,
      traditional = {
        res <- g_com_traditional(data = data, nboot = nboot)
        list(
          estimate = res$ate_diff,
          se = res$se_diff,
          ci = res$CI_diff,
          raw = res
        )
      },
      observational = {
        res <- g_com_tte(data = data, nboot = nboot)
        list(
          estimate = res$ate,
          se = res$se_boot,
          ci = res$CI,
          raw = res
        )
      },
      direct = {
        res <- g_com_direct(data = data, eff = "trt", nboot = nboot)
        list(
          estimate = res$eff_diff,
          se = res$se_boot,
          ci = res$CI,
          raw = res
        )
      },
      indirect = {
        res <- g_com_indirt(data = data, nboot = nboot)
        list(
          estimate = res$indirect,
          se = res$se_boot,
          ci = res$CI,
          raw = res
        )
      }
    )

    return(new_tbcompare_result(
      estimate = fit$estimate,
      se = fit$se,
      ci = fit$ci,
      approach = approach,
      estimator = estimator,
      n = nrow(data),
      raw = fit$raw
    ))
  }

  if (estimator == "onestep") {
    if (approach == "observational") {
      dat_ab <- data[data$treatment %in% c("a", "b"), , drop = FALSE]
      dat_ab$s <- 3
      dat_ab$treatment <- as.numeric(as.character(dat_ab$treatment) == "a")

      nuis <- estimate_nuisances(
        data = dat_ab,
        outcome = "y",
        treat = "treatment",
        covars = c("L1", "L2"),
        nuisance_type = nuisance_type
      )

      res <- one_step_ATE(
        data = dat_ab,
        nuisances = nuis,
        outcome = "y",
        treat = "treatment",
        covars = c("L1", "L2"),
        nuisance_type = nuisance_type
      )

      return(new_tbcompare_result(
        estimate = res$psi,
        se = res$se,
        ci = c(res$ci_lower, res$ci_upper),
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(dat_ab),
        raw = res
      ))
    }

    if (approach == "traditional") {
      res <- eif_trial_onestep(
        data = data,
        outcome = "y",
        treat = "treatment",
        covars = c("L1", "L2", "country"),
        treat_s1_active = "a",
        treat_s1_control = "sc1",
        treat_s2_active = "b",
        treat_s2_control = "sc2",
        nuisance_type = nuisance_type
      )

      return(new_tbcompare_result(
        estimate = res$ate_diff,
        se = res$se_diff,
        ci = c(res$ci_lower, res$ci_upper),
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }

    if (approach == "direct") {
      res <- eif_direct_est(
        dataset = data,
        a_val = "a",
        b_val = "b",
        nuisance_type = nuisance_type
      )

      return(new_tbcompare_result(
        estimate = res$estimate,
        se = res$se_if,
        ci = res$CI_if,
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }

    if (approach == "indirect") {
      res <- eif_indirect_est(
        dataset = data,
        a = "a",
        b = "b",
        c1 = "sc1",
        c2 = "sc2",
        nuisance_type = nuisance_type
      )

      return(new_tbcompare_result(
        estimate = res$indirect,
        se = res$se_if,
        ci = res$CI_if,
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }
  }

  if (estimator == "tmle") {
    if (approach == "observational") {
      res <- tmle_tte(
        data = data,
        nuisance_type = nuisance_type,
        a_val = "a",
        b_val = "b",
        covars = c("L1", "L2", "country"),
        outcome = "y"
      )

      return(new_tbcompare_result(
        estimate = res$estimate,
        se = res$se,
        ci = res$CI,
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }

    if (approach == "traditional") {
      res <- tmle_traditional(
        data = data,
        nuisance_type = nuisance_type,
        a_val = "a",
        b_val = "b",
        c1_val = "sc1",
        c2_val = "sc2",
        covars = c("L1", "L2", "country"),
        outcome = "y"
      )

      return(new_tbcompare_result(
        estimate = res$estimate,
        se = res$se,
        ci = res$CI,
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }

    if (approach == "direct") {
      res <- tmle_direct(
        data = data,
        nuisance_type = nuisance_type,
        a_val = "a",
        b_val = "b",
        covars = c("L1", "L2", "country"),
        outcome = "y"
      )

      return(new_tbcompare_result(
        estimate = res$estimate,
        se = res$se,
        ci = res$CI,
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }

    if (approach == "indirect") {
      res <- tmle_indirect(
        data = data,
        nuisance_type = nuisance_type,
        a_val = "a",
        b_val = "b",
        c1_val = "sc1",
        c2_val = "sc2",
        covars = c("L1", "L2", "country"),
        outcome = "y"
      )

      return(new_tbcompare_result(
        estimate = res$estimate,
        se = res$se,
        ci = res$CI,
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ))
    }
  }

  stop("Estimator and approach combination not implemented.", call. = FALSE)
}

new_tbcompare_result <- function(estimate,
                                 se,
                                 ci,
                                 approach,
                                 estimator,
                                 n,
                                 raw,
                                 nuisance_type = NA_character_) {
  structure(
    list(
      estimate = as.numeric(estimate),
      se = as.numeric(se),
      ci = as.numeric(ci),
      approach = approach,
      estimator = estimator,
      nuisance_type = nuisance_type,
      n = n,
      raw = raw
    ),
    class = "tbcompare"
  )
}

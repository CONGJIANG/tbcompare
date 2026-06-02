#' Main interface for TB comparison estimators
#'
#' @param data Data frame.
#' @param approach One of
#'   "traditional",
#'   "observational",
#'   "direct",
#'   "indirect".
#' @param estimator One of
#'   "ipw",
#'   "gcomp",
#'   "onestep",
#'   "tmle".
#' @param nboot Number of bootstrap samples for bootstrap-based estimators.
#' @param nuisance_type Nuisance estimation type. Currently mainly used for EIF/TMLE; one of "simple" or "flexible".
#' @return A tbcompare object.
#'
#' @export
tb_compare <- function(
  data,
  approach = c("traditional", "tte", "direct", "indirect"),
  estimator = c("ipw", "gcomp", "onestep", "tmle"),
  nboot = 100,
  nuisance_type = c("simple", "flexible")
) {
  approach <- match.arg(approach)
  estimator <- match.arg(estimator)
  nuisance_type <- match.arg(nuisance_type)

  if (estimator == "ipw") {
    res <- estimate_ipw_unified(data)

    out <- switch(
      approach,
      traditional = list(
        estimate = res$psi_traditional,
        se = res$se_traditional,
        ci = res$CI_traditional
      ),
      tte = list(
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

    return(structure(
      list(
        estimate = as.numeric(out$estimate),
        se = as.numeric(out$se),
        ci = as.numeric(out$ci),
        approach = approach,
        estimator = estimator,
        n = nrow(data),
        raw = res
      ),
      class = "tbcompare"
    ))
  }
  if (estimator == "gcomp") {

  fit <- switch(
    approach,
    traditional = {
      res <- g_com_traditional(data = data)
      list(
        estimate = res$ate_diff,
        se = res$se_diff,
        ci = res$CI_diff,
        raw = res
      )
    },
    tte = {
      res <- g_com_tte(data = data, nboot = nboot)
      list(
        estimate = res$ate,
        se = res$se_boot,
        ci = res$CI,
        raw = res
      )
    },
    direct = {
      res <- g_com_direct(data = data, eff = "trt")
      list(
        estimate = res$eff_diff,
        se = res$se_boot,
        ci = res$CI,
        raw = res
      )
    },
    indirect = {
      res <- g_com_indirt(data = data)
      list(
        estimate = res$indirect,
        se = res$se_boot,
        ci = res$CI,
        raw = res
      )
    }
  )

  return(structure(
    list(
      estimate = as.numeric(fit$estimate),
      se = as.numeric(fit$se),
      ci = as.numeric(fit$ci),
      approach = approach,
      estimator = estimator,
      n = nrow(data),
      raw = fit$raw
    ),
    class = "tbcompare"
  ))
}
  if (estimator == "onestep") {

  if (approach == "tte") {

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

    return(structure(
      list(
        estimate = as.numeric(res$psi),
        se = as.numeric(res$se),
        ci = as.numeric(c(res$ci_lower, res$ci_upper)),
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(dat_ab),
        raw = res
      ),
      class = "tbcompare"
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

  return(structure(
    list(
      estimate = as.numeric(res$ate_diff),
      se = as.numeric(res$se_diff),
      ci = as.numeric(c(res$ci_lower, res$ci_upper)),
      approach = approach,
      estimator = estimator,
      nuisance_type = nuisance_type,
      n = nrow(data),
      raw = res
    ),
    class = "tbcompare"
  ))
}
    if (approach == "direct") {

  res <- eif_direct_est(
    dataset = data,
    a_val = "a",
    b_val = "b",
    nuisance_type = nuisance_type
  )

  return(structure(
    list(
      estimate = as.numeric(res$estimate),
      se = as.numeric(res$se_if),
      ci = as.numeric(res$CI_if),
      approach = approach,
      estimator = estimator,
      nuisance_type = nuisance_type,
      n = nrow(data),
      raw = res
    ),
    class = "tbcompare"
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

  return(structure(
    list(
      estimate = as.numeric(res$indirect),
      se = as.numeric(res$se_if),
      ci = as.numeric(res$CI_if),
      approach = approach,
      estimator = estimator,
      nuisance_type = nuisance_type,
      n = nrow(data),
      raw = res
    ),
    class = "tbcompare"
  ))
} 
  stop("For now, estimator = 'onestep' supports approach = 'tte', 'traditional', 'direct', or 'indirect'.")
  }
  
  if (estimator == "tmle") {

  if (approach == "tte") {

    res <- tmle_tte(
      data = data,
      nuisance_type = nuisance_type,
      a_val = "a",
      b_val = "b",
      covars = c("L1", "L2", "country"),
      outcome = "y"
    )

    return(structure(
      list(
        estimate = as.numeric(res$estimate),
        se = as.numeric(res$se),
        ci = as.numeric(res$CI),
        approach = approach,
        estimator = estimator,
        nuisance_type = nuisance_type,
        n = nrow(data),
        raw = res
      ),
      class = "tbcompare"
    ))
  }

  stop("For now, estimator = 'tmle' only supports approach = 'tte'.")
}
  stop("Estimator not yet implemented.")
}
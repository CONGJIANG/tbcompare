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
  nuisance_type = "simple"
) {
  approach <- match.arg(approach)
  estimator <- match.arg(estimator)

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
  stop("Only estimator = 'ipw' is currently implemented.")
}
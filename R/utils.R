# ============================================================
# General helpers
# ============================================================
#' @importFrom stats binomial glm predict qnorm rbinom rnorm sd weighted.mean
NULL

.clamp_prob <- function(x, clamp = c(1e-6, 1 - 1e-6)) {
  pmin(pmax(as.numeric(x), clamp[1]), clamp[2])
}

.safe_factor_country <- function(data) {
  data$country <- factor(data$country, levels = c(1, 2, 3, 4))
  data
}

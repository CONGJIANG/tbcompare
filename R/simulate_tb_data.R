#' Get TB DGP parameters
#'
#' @param dgp_choice Character. One of \code{"delta0"}, \code{"delta1"},
#'   \code{"delta2"}, or \code{"delta3"}.
#'
#' @return A list containing scenario-specific DGP parameters.
#'
#' @examples
#' get_tb_dgp_parameters("delta1")
#'
#' @export
get_tb_dgp_parameters <- function(dgp_choice = "delta1") {
  delta_c <- switch(
    as.character(dgp_choice),
    delta0 = 0.00,
    delta1 = 0.10,
    delta2 = 0.20,
    delta3 = 0.30,
    stop(sprintf("Unknown dgp_choice: %s", dgp_choice))
  )

  list(delta_c = delta_c)
}

#' Generate TB cross-trial simulation data
#'
#' This is the core data-generating mechanism for the TB cross-trial
#' comparison setting. The function generates two trials, baseline covariates,
#' country, treatment assignment, and a binary outcome.
#'
#' @param n Integer. Sample size.
#' @param beta0 Numeric. Outcome-model intercept.
#' @param betaL1 Numeric. Coefficient for baseline covariate \code{L1}.
#' @param betaL2 Numeric. Coefficient for baseline covariate \code{L2}.
#' @param betaA_a Numeric. Treatment effect for treatment \code{a}.
#' @param betaA_b Numeric. Treatment effect for treatment \code{b}.
#' @param betaA_sc1 Numeric. Treatment effect for standard care version \code{sc1}.
#' @param delta_c Numeric. Difference between \code{sc1} and \code{sc2};
#'   internally, \code{betaA_sc2 = betaA_sc1 - delta_c}.
#' @param betaC2 Numeric. Reserved country parameter for compatibility with
#'   earlier simulation scripts.
#' @param betaC3 Numeric. Reserved country parameter for compatibility with
#'   earlier simulation scripts.
#' @param betaC4 Numeric. Reserved country parameter for compatibility with
#'   earlier simulation scripts.
#'
#' @return A data frame with columns \code{s}, \code{L1}, \code{L2},
#'   \code{country}, \code{treatment}, and \code{y}.
#'
#' @examples
#' dat <- dgp_tb_trial_data(n = 100, delta_c = 0.10)
#' head(dat)
#'
#' @export
dgp_tb_trial_data <- function(
  n = 1000,
  beta0 = -0.5,
  betaL1 = 0.4,
  betaL2 = -0.3,
  betaA_a = 0.8,
  betaA_b = 0.5,
  betaA_sc1 = 0.2,
  delta_c = 0.15,
  betaC2 = 0.3,
  betaC3 = 0.5,
  betaC4 = 0.2
) {
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("n must be a positive integer.")
  }

  n <- as.integer(n)

  # derive betaA_sc2 from delta_c
  betaA_sc2 <- betaA_sc1 - delta_c

  # covariates and trial assignment
  s <- sample(1:2, n, replace = TRUE)
  L1 <- rnorm(n, mean = 1 - s, sd = 1)
  L2 <- rbinom(n, size = 1, prob = 0.35 * s)
  country <- sample(1:4, n, replace = TRUE)

  # treatment assignment by trial
  treatment <- character(n)

  treatment[s == 1] <- sample(
    c("a", "sc"),
    sum(s == 1),
    replace = TRUE
  )

  treatment[s == 2] <- sample(
    c("b", "sc"),
    sum(s == 2),
    replace = TRUE
  )

  # split standard care into sc1 / sc2 by country,
  # correlated with country but not deterministically assigned
  is_sc <- treatment == "sc"

  p_sc1 <- ifelse(country %in% c(1, 3), 0.7, 0.3)

  treatment[is_sc] <- ifelse(
    rbinom(sum(is_sc), size = 1, prob = p_sc1[is_sc]) == 1,
    "sc1",
    "sc2"
  )

  # treatment effects
  tr_map <- stats::setNames(
    c(betaA_a, betaA_b, betaA_sc1, betaA_sc2),
    c("a", "b", "sc1", "sc2")
  )

  treatment_effect <- unname(tr_map[treatment])
  treatment_effect[is.na(treatment_effect)] <- 0

  effect_mod_component <- numeric(n)

  effect_mod_component[treatment == "a"] <-
    0.4 * L1[treatment == "a"]

  effect_mod_component[treatment == "b"] <-
    0.2 * L1[treatment == "b"]

  effect_mod_component[treatment == "sc1"] <-
    0.2 * L1[treatment == "sc1"]

  effect_mod_component[treatment == "sc2"] <-
    0.1 * L1[treatment == "sc2"]

  # binary outcome model
  lp <- beta0 +
    betaL1 * L1 +
    betaL2 * L2 +
    effect_mod_component +
    treatment_effect

  prob <- stats::plogis(lp)
  y <- stats::rbinom(n, size = 1, prob = prob)

  data.frame(
    s = s,
    L1 = L1,
    L2 = L2,
    country = country,
    treatment = treatment,
    y = y
  )
}

#' Simulate TB comparison data
#'
#' User-facing wrapper for generating TB cross-trial simulation data.
#' This function uses the scenario choices \code{"delta0"}, \code{"delta1"},
#' \code{"delta2"}, and \code{"delta3"} to set the standard-care contrast.
#'
#' @param n Integer. Sample size.
#' @param dgp_choice Character. One of \code{"delta0"}, \code{"delta1"},
#'   \code{"delta2"}, or \code{"delta3"}.
#' @param seed Optional integer random seed.
#' @param ... Additional arguments passed to \code{dgp_tb_trial_data()}.
#'
#' @return A data frame with simulated TB comparison data.
#'
#' @examples
#' dat <- simulate_tb_data(n = 1000, dgp_choice = "delta1", seed = 1)
#' table(dat$s, dat$treatment)
#'
#' @export
simulate_tb_data <- function(
  n = 1000,
  dgp_choice = "delta1",
  seed = NULL,
  ...
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  pars <- get_tb_dgp_parameters(dgp_choice = dgp_choice)

  dgp_tb_trial_data(
    n = n,
    delta_c = pars$delta_c,
    ...
  )
}

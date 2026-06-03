test_that("IPW direct runs", {
  dat <- simulate_xtrial_data(
    n = 500
  )
  fit <- xtrial(
    dat,
    approach = "direct",
    estimator = "ipw"
  )
  expect_true(is.numeric(fit$estimate))
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
})
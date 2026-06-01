test_that("IPW direct runs", {
  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )
  fit <- tb_compare(
    dat,
    approach = "direct",
    estimator = "ipw"
  )
  expect_true(is.numeric(fit$estimate))
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
})
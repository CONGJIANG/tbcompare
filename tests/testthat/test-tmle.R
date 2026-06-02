test_that("TMLE TTE runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "tte",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})
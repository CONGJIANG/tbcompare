test_that("TMLE TTE runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "observational",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("TMLE traditional runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "traditional",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("TMLE direct runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "direct",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("TMLE indirect runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "indirect",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})
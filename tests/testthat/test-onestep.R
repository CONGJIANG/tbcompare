test_that("one-step tte runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "tte",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("one-step traditional runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "traditional",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("one-step direct runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "direct",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("one-step indirect runs", {

  dat <- simulate_tb_data(
    n = 500,
    seed = 1
  )

  fit <- tb_compare(
    data = dat,
    approach = "indirect",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "tbcompare")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})
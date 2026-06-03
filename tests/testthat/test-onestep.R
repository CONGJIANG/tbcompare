test_that("one-step tte runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "observational",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("one-step traditional runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "traditional",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("one-step direct runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "direct",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("one-step indirect runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "indirect",
    estimator = "onestep",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})
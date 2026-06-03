test_that("TMLE TTE runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "observational",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("TMLE traditional runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "traditional",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("TMLE direct runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "direct",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})

test_that("TMLE indirect runs", {

  dat <- simulate_xtrial_data(
    n = 500
  )

  fit <- xtrial(
    data = dat,
    approach = "indirect",
    estimator = "tmle",
    nuisance_type = "simple"
  )

  expect_s3_class(fit, "xtrial")
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$se))
  expect_length(fit$ci, 2)
})
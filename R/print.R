#' @export
print.tbcompare <- function(x, ...) {

  cat("\n")
  cat("tbcompare object\n\n")

  cat("Estimator:",
      toupper(x$estimator),
      "\n")

  cat("Approach :",
      tools::toTitleCase(x$approach),
      "\n\n")

  cat(
    sprintf(
      "Estimate : %.4f\n",
      x$estimate
    )
  )

  cat(
    sprintf(
      "SE       : %.4f\n",
      x$se
    )
  )

  cat("\n95% CI:\n")

  cat(
    sprintf(
      "(%.4f, %.4f)\n",
      x$ci[1],
      x$ci[2]
    )
  )

  invisible(x)
}


#' @export
summary.tbcompare <- function(object, ...) {

  out <- data.frame(
    estimator = object$estimator,
    approach = object$approach,
    estimate = object$estimate,
    se = object$se,
    ci_low = object$ci[1],
    ci_high = object$ci[2]
  )

  out
}

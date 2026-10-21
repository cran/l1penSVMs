#' Print Cross-Validation Results
#'
#' @param x An object of class \code{cv.l1svm}.
#' @param ... Additional arguments.
#' @return Invisibly returns the input object of class "cv.l1svm". The function prints a summary of the cross-validation results.
#' @export
print.cv.l1svm <- function(x, ...) {

  if (x$classification == "SVM") {
    cat("\nL1-penalized SVM cross-validation\n\n")
    cat("Optimal t:", x$t.opt, "\n")
    cat("Minimum CV error:", min(x$cv.err), "\n")
    cat("Number of selected variables:", sum(x$params.opt[-1] != 0), "\n")

  } else if (x$classification == "ROCSVM") {
    cat("\nL1-penalized ROC-SVM cross-validation\n\n")
    cat("Optimal t:", x$t.opt, "\n")
    cat("Minimum CV AUC-risk:", min(x$cv.AUC.risk), "\n")
    cat("Number of selected variables:", sum(x$params.opt != 0), "\n")

  } else {
    cat("\nL1-penalized FD-SVM cross-validation\n\n")
    cat("Optimal t:", x$t.opt, "\n")
    cat("Sensitivity:", x$sens.opt, "\n")
    cat("Specificity:", x$spec.opt, "\n")
    cat("Number of selected variables:", sum(x$params.opt[-1] != 0), "\n")
  }

  cat("\nCoefficients at optimal t:\n")
  print(x$params.opt)

  invisible(x)
}

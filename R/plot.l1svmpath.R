#' Plot the Entire Regularization Paths for L1-penalized SVM models
#' @description
#' Plot the whole regularization paths for L1-penalized Support Vector Machine (SVM), Receiver Operating Characteristic Support Vector Machine (ROC-SVM) and Fraud Detection Support Vector Machine (FD-SVM)
#' @param x An object of class \code{l1svmpath}.
#' @param ... Additional graphical arguments passed to plotting functions.
#' @details
#' This function is a method for the generic \code{plot()} function,
#' @return Invisibly returns the input object of class \code{l1svmpath}. The function displays the regularization paths.
#' @examples
#' library(l1penSVMs)
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p),n,p)
#' beta  <- c(2,-1,0.5,rep(0,p-3))
#' fx  <- -2 + crossprod(t(x),beta)
#' y   <- sign(fx + rnorm(n,0,1))
#'
#' # SVM
#' obj <- l1svmpath(x, y, classification="SVM")
#' plot(obj)
#' @export
plot.l1svmpath <- function(x, ...) {

  t <- na.omit(x$t)
  beta <- x$beta
  b0 <- if (!is.null(x$b0)) x$b0 else NULL

  ylim <- if (is.null(b0)) range(beta, na.rm = TRUE) else
    range(c(b0, beta), na.rm = TRUE)

  plot(t, beta[, 1], type = "n", xlim = c(0, max(t) + 0.02),
       ylim = ylim, xlab = "t", ylab = expression(beta),
       cex.lab = 1.5, ...)

  abline(v = t[-1], lwd = 0.5, col = "gray", lty = 3)
  matlines(t, beta, col = "red", lwd = 1, lty = 1)

  if (is.null(b0)) {
    legend("topleft", "coefficient", col = "red", pch = 15,
           cex = 1.2, pt.cex = 1.5, bty = "n")
  } else {
    lines(t, b0, col = "green", lwd = 1)
    legend("topleft", c("coefficient", "intercept"),
           col = c("red", "green"), pch = 15,
           cex = 1.2, pt.cex = 1.5, bty = "n")
  }

  invisible(x)
}

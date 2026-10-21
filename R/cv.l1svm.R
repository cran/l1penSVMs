#' Perform CV to choose the optimal t for three L1-penalized Support Vector Machine models
#' @description
#' This function performs k-fold cross-validation to select the optimal for the L1-penalized Support Vector Machine (SVM), Receiver Operating Characteristic Support Vector Machine (ROC-SVM) and Fraud Detection Support Vector Machine (FD-SVM) based on the regularization path algorithm.
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns); categorical variables are automatically one-hot encoded.
#' @param y Response vector with class labels in \{-1, 1\}. Labels given as \{0, 1\} or as a two-level factor/character are automatically converted to this format.
#' @param classification The classification method is to be applied. Either "SVM", "ROC-SVM" and "FD-SVM"
#' @param nfold The number of folds in the cross-validation (CV), default is \code{nfold=5} for "SVM" and "ROC-SVM" and \code{nfold=2} for "FD-SVM"
#' @param sensitivity The minimum proportion of actual fraud cases that are correctly classified as fraud for the "FD-SVM".
#' @param init.t For "SVM" and "FD-SVM", the initial regularization parameter with sufficiently small t; default is 1e-5
#' @param max.iter The maximum number of iterations (total across the entire path algorithm); default is 10000
#' @param standardize Standardize each column of the data matrix to have zero mean and unit variance
#' @param incomp.size Sample size for the incomplete U-statistic in ROC-SVM; default is the number of observations
#' @details The \code{cv.l1.svm} determines the optimal regularization parameter \eqn{t} and yields the cofficient estimates at the optimal \eqn{t}
#' @return An object \code{cv.l1svm} returns a list with following components:
#' \describe{
#' \item{t}{ A vector of t values used in K-fold CV. }
#' \item{cv.err}{A vector of the averaged CV errors for different t values.}
#' \item{cv.AUC.risk}{A vector of the averaged CV AUC-risks for different t values.}
#' \item{cv.sens}{A vector of the averaged CV sensitivies for different t values.}
#' \item{cv.spec}{A vector of the averaged CV specificies for different t values.}
#' \item{t.opt}{Optimal value of t for thee L1-penalized SVM models.}
#' \item{sens.opt}{Optimal value of sensitivity at t=t.opt for the L1-penalized FD-SVM.}
#' \item{spec.opt}{Optimal value of specificity at t=t.opt for the L1-penalized FD-SVM.}
#' \item{params.opt}{ The estimated coefficients at the optimal parameter t for three L1-penalized SVM.}
#' }
#' @importFrom graphics abline legend lines matlines par points
#' @importFrom stats na.omit rnorm sd
#' @importFrom MASS ginv
#' @references
#' Ji Zhu, Saharon Rosset, Robert Tibshirani and Trevor Hastie, "1-norm Support Vector Machines", Advances in neural information processing systems, 16.
#'
#' Hyungwoo Kim, Insuk Sohn and Seung Jun Shin, "Regularization paths of L1-penalized ROC Curve-Optimizing Support Vector Machines", Stat, 10(1), e400.
#'
#' Minhyoung Park, Hyungwoo Kim and Seung Jun Shin, "L1-penalized Fraud Detection Support Vector Machines", under
#' @author Hyungwoo Kim, MinHyung Park, Seung Jun Shin
#' @seealso \code{\link{l1.svmpath}}, \code{\link{l1.rocsvmpath}}, \code{\link{l1.fdsvmpath}}, \code{\link{plot.l1svmpath}}, \code{\link{predict.l1svmpath}}
#' @examples
#' library(l1penSVMs)
#' set.seed(2)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p),n,p)
#' beta  <- c(2,-1,0.5,rep(0,p-3))
#' fx  <- -2 + crossprod(t(x),beta)
#' y   <- sign(fx + rnorm(n,0,1))

#' # SVM
#' cv.svm <- cv.l1svm(x, y, classification="SVM")
#'
#' # ROCSVM
#' cv.rocsvm <- cv.l1svm(x, y, classification="ROCSVM")
#'
#' # FDSVM
#' cv.fdsvm <- cv.l1svm(x, y, classification="FDSVM", sensitivity=0.9)
#' @export
cv.l1svm <- function(x, y,
                     classification = c("SVM", "ROCSVM", "FDSVM"),
                     nfold = NULL, sensitivity = NULL,
                     init.t = NULL, max.iter = NULL,
                     standardize = TRUE, incomp.size = NULL) {

  # Model specification
  cl <- match.arg(classification)

  if (cl == "FDSVM" && is.null(sensitivity)) {
    stop("You must specify 'sensitivity' for FD-SVM.", call. = FALSE)
  }

  # Fit full-data path
  obj <- l1svmpath(x = x, y = y, classification = cl, init.t = init.t,
                   max.iter = max.iter, standardize = standardize,
                   incomp.size = incomp.size)

  x <- obj$x_orig
  y <- obj$y_cv
  n <- length(y)

  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  x.names <- colnames(x)

  # Candidate t-grid
  if (cl == "SVM") {
    t.grid <- c(seq(0, 5, length.out = 100),
                seq(5, 19.5, length.out = 50))
  } else {
    t.grid <- c(seq(0, 3, length.out = 100),
                seq(3, 19.5, length.out = 50))
  }

  # Number of folds
  if (is.null(nfold)) nfold <- if (cl == "FDSVM") 2L else 5L

  if (cl == "FDSVM") {
    n.pos <- sum(y == 1)
    n.neg <- sum(y == -1)

    if (nfold < 2L || nfold > min(n.pos, n.neg)) {
      stop(
        "'nfold' must be between 2 and the size of the smaller class."
      )
    }

    fold.id <- integer(n)

    for (cls in c(-1, 1)) {
      ind <- which(y == cls)

      fold.id[ind] <- sample(
        rep(seq_len(nfold), length.out = length(ind))
      )
    }

  } else {
    fold.id <- sample(
      rep(seq_len(nfold), length.out = n)
    )
  }


  # Storage
  cv.result <- matrix(NA_real_, nrow = nfold, ncol = length(t.grid))
  cv.sens <- matrix(NA_real_, nrow = nfold, ncol = length(t.grid))
  cv.spec <- matrix(NA_real_, nrow = nfold, ncol = length(t.grid))

  # Interpolate coefficient path
  interpolate_path <- function(path.t, path.beta, t.values) {

    ans <- matrix(
      NA_real_,
      nrow = ncol(path.beta),
      ncol = length(t.values)
    )

    valid <- which(
      t.values >= min(path.t) &
        t.values <= max(path.t)
    )

    for (j in valid) {
      tt <- t.values[j]

      if (tt %in% path.t) {
        ans[, j] <- path.beta[match(tt, path.t), ]
      } else {
        k <- which(path.t > tt)[1L]

        w <- (tt - path.t[k - 1L]) /
          (path.t[k] - path.t[k - 1L])

        ans[, j] <- (1 - w) * path.beta[k - 1L, ] +
          w * path.beta[k, ]
      }
    }

    ans
  }

  # Cross-validation
  for (i in seq_len(nfold)) {
    ind.tr <- which(fold.id != i)
    ind.ts <- which(fold.id == i)

    tr.x <- x[ind.tr, , drop = FALSE]
    tr.y <- y[ind.tr]
    ts.x <- x[ind.ts, , drop = FALSE]
    ts.y <- y[ind.ts]

    obj.tr <- l1svmpath(
      x = tr.x,
      y = tr.y,
      classification = cl,
      init.t = init.t,
      max.iter = max.iter,
      standardize = standardize,
      incomp.size = incomp.size,
      fraud.label.fixed = (cl == "FDSVM")
    )

    # Transform held-out data using training-fold statistics only
    if (isTRUE(obj.tr$standardize)) {
      ts.x <- sweep(ts.x, 2L, obj.tr$x.mean, FUN = "-")
      ts.x <- sweep(ts.x, 2L, obj.tr$x.sd, FUN = "/")
    }

    # SVM
    if (cl == "SVM") {
      path.t <- obj.tr$t
      path.beta <- cbind(obj.tr$b0, obj.tr$beta)

      if (length(path.t) <= 10L) next

      beta.t.val <- interpolate_path(path.t, path.beta, t.grid)
      valid <- which(colSums(is.na(beta.t.val)) == 0L)
      if (length(valid) > 0L) {
        pred.ts.y <- sign(cbind(1, ts.x) %*%beta.t.val[, valid, drop = FALSE])
        test.acc <- colMeans(pred.ts.y == ts.y)
        cv.result[i, valid] <- 1 - test.acc
      }

      # ROC-SVM
    } else if (cl == "ROCSVM") {
      path.t <- obj.tr$t
      path.beta <- obj.tr$beta

      if (length(path.t) <= 10L) next

      beta.t.val <- interpolate_path(path.t, path.beta, t.grid)

      valid <- which(colSums(is.na(beta.t.val)) == 0L)
      ind.pos <- which(ts.y == 1)
      ind.neg <- which(ts.y == -1)

      if (length(ind.pos) == 0L || length(ind.neg) == 0L) {
        warning(sprintf("Fold %d contains only one response class.", i))
        next
      }

      if (length(valid) > 0L) {
        fhat <- ts.x %*% beta.t.val[, valid, drop = FALSE]

        test.auc <- vapply(
          seq_len(ncol(fhat)),
          function(j) {

            score.pos <- fhat[ind.pos, j]
            score.neg <- fhat[ind.neg, j]

            comp <- outer(score.pos, score.neg, "-")

            mean(comp > 0) + 0.5 * mean(comp == 0)
          },
          numeric(1)
        )

        cv.result[i, valid] <- 1 - test.auc
      }


      # FD-SVM
    } else {
      path.t <- obj.tr$t
      path.beta <- cbind(obj.tr$b0, obj.tr$beta)

      if (length(path.t) <= 1L) next
      beta.t.val <- interpolate_path(path.t, path.beta, t.grid)

      valid <- which(colSums(is.na(beta.t.val)) == 0L)

      ind.pos <- which(ts.y == 1)
      ind.neg <- which(ts.y == -1)

      if (length(ind.pos) == 0L || length(ind.neg) == 0L) {
        warning(sprintf("Fold %d contains only one response class.", i))
        next
      }

      if (length(valid) > 0L) {

        pred.ts.y <- sign(
          cbind(1, ts.x) %*%
            beta.t.val[, valid, drop = FALSE]
        )

        cv.sens[i, valid] <- colMeans(
          pred.ts.y[ind.pos, , drop = FALSE] == 1
        )

        cv.spec[i, valid] <- colMeans(
          pred.ts.y[ind.neg, , drop = FALSE] == -1
        )
      }
    }
  }

  # Select optimal t
  if (cl != "FDSVM") {
    valid.ind <- which(
      colSums(is.na(cv.result)) < nfold / 2 &
        t.grid >= min(obj$t) &
        t.grid <= max(obj$t)
    )

    if (length(valid.ind) == 0L) {
      stop("No valid candidate values of 't' were obtained.",
           call. = FALSE)
    }

    t.grid <- t.grid[valid.ind]
    cv.result <- cv.result[, valid.ind, drop = FALSE]
    mean.cv <- colMeans(cv.result, na.rm = TRUE)

    opt.ind <- which.min(mean.cv)
    opt.t <- t.grid[opt.ind]

  } else {
    valid.ind <- which(
      colSums(is.na(cv.sens) | is.na(cv.spec)) == 0L &
        t.grid >= min(obj$t) &
        t.grid <= max(obj$t)
    )

    if (length(valid.ind) == 0L) {
      stop("No valid candidate values of 't' were obtained.",
           call. = FALSE)
    }

    t.grid <- t.grid[valid.ind]
    cv.sens <- cv.sens[, valid.ind, drop = FALSE]
    cv.spec <- cv.spec[, valid.ind, drop = FALSE]

    mean.sens <- colMeans(cv.sens, na.rm = TRUE)
    mean.spec <- colMeans(cv.spec, na.rm = TRUE)
    feasible.ind <- which(mean.sens >= sensitivity)

    if (length(feasible.ind) == 0L) {
      stop("There is no solution that satisfies the specified sensitivity.",
           call. = FALSE)
    }

    opt.ind <- feasible.ind[which.max(mean.spec[feasible.ind])]
    opt.t <- t.grid[opt.ind]
    opt.sens <- mean.sens[opt.ind]
    opt.spec <- mean.spec[opt.ind]
  }

  # Coefficients at optimal t
  if (cl == "ROCSVM") {
    path.beta <- obj$beta
  } else {
    path.beta <- cbind(obj$b0, obj$beta)
  }

  # Check range of t
  beta.opt <- interpolate_path(obj$t, path.beta, opt.t)
  if (anyNA(beta.opt)) {
    stop(
      "The selected t is outside the full-data regularization path.",
      call. = FALSE
    )
  }
  beta.opt <- as.numeric(beta.opt)

  if (cl == "ROCSVM") {
    names(beta.opt) <- x.names
  } else {
    names(beta.opt) <- c("(Intercept)", x.names)
  }

  # Plot and output
  if (cl == "SVM") {
    plot(t.grid, mean.cv, type = "b", xlab = "t", ylab = "CV error",
         cex.lab = 1.2, cex.axis = 1.2, col = "blue")
    abline(v = opt.t, lty = 2, col = "red")

    result <- list(classification = cl, t = t.grid,
                   cv.err = mean.cv, t.opt = opt.t,
                   params.opt = beta.opt)

  } else if (cl == "ROCSVM") {
    plot(t.grid, mean.cv, type = "b", xlab = "t", ylab = "CV AUC-risk",
         cex.lab = 1.2, cex.axis = 1.2, col = "blue")
    abline(v = opt.t, lty = 2, col = "red")

    result <- list(classification = cl, t = t.grid,
                   cv.AUC.risk = mean.cv, t.opt = opt.t,
                   params.opt = beta.opt)

  } else {
    plot(t.grid, mean.sens, type = "b", xlab = "t", ylab = "",
         ylim = c(0, 1), pch = 3, cex.lab = 1.2, col = "blue")
    points(t.grid, mean.spec, type = "b", pch = 2, col = "black")
    abline(h = sensitivity, lty = 2)
    points(opt.t, opt.sens, pch = 18, cex = 1.2, col = "red")

    legend("bottomright", legend = c("Sensitivity", "Specificity"),
           col = c("blue", "black"), pch = c(3, 2), lty = 1,
           cex = 1.2, bty = "n")

    result <- list(classification = cl, t = t.grid,
                   cv.sens = mean.sens, cv.spec = mean.spec,
                   t.opt = opt.t, sens.opt = opt.sens,
                   spec.opt = opt.spec, params.opt = beta.opt)
  }

  class(result) <- "cv.l1svm"
  return(result)

}

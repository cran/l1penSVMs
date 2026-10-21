#' Predict the class label or classification score
#'
#' @param object A fitted regularization path object of class \code{l1svmpath}.
#' @param x The test data matrix (n x p), with n rows (observations) on p variables (columns)
#' @param t The value of regularization parameter
#' @param ... Additional arguments passed to the prediction method.
#' @details
#' This function returns the predicted class labels or the fitted classification scores for L1-penalized Support Vector Machine (SVM), Receiver Operating Characteristic Support Vector Machine (ROC-SVM) and Fraud Detection Support Vector Machine (FD-SVM)
#' @return The predicted class label for the L1-penalized SVM and FDSVM,
#' or the fitted classification scores for the L1-penalized ROCSVM.
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
#' obj <- l1svmpath(x,y,classification="SVM")
#' predict(obj, x=x, t=1.5)

#' @export
predict.l1svmpath <- function(object, x, t = NULL, ...) {

  if (is.null(t)) {
    stop("The regularization parameter 't' must be specified.",
         call. = FALSE)
  }

  obj <- object
  t.val <- t

  # Coefficient path
  if (obj$classification != "ROCSVM") {
    path.beta <- cbind(obj$b0, obj$beta)
  } else {
    path.beta <- obj$beta
  }

  path.beta <- as.matrix(path.beta)

  if (any(!is.finite(t.val)) ||
      any(t.val < min(obj$t)) ||
      any(t.val > max(obj$t))) {
    stop("'t' must lie within the fitted regularization path.")
  }

  curr.b <- matrix(
    NA_real_,
    nrow = ncol(path.beta),
    ncol = length(t.val)
  )

  for (j in seq_along(t.val)) {

    tt <- t.val[j]

    if (tt %in% obj$t) {

      curr.b[, j] <- path.beta[match(tt, obj$t), ]

    } else {

      k <- which(obj$t > tt)[1L]

      w <- (tt - obj$t[k - 1L]) /
        (obj$t[k] - obj$t[k - 1L])

      curr.b[, j] <-
        (1 - w) * path.beta[k - 1L, ] +
        w * path.beta[k, ]
    }
  }

  # Check new data
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data.frame.", call. = FALSE)
  }

  if (is.data.frame(x)) {
    has.listcol <- any(vapply(x, is.list, logical(1)))
    if (has.listcol) {
      stop("'x' must not contain list-columns.", call. = FALSE)
    }
  }

  # Categorical variables from training data
  cat.vars <- object$cat.vars

  if (length(cat.vars) > 0L) {

    if (!is.data.frame(x)) {
      stop("Categorical predictors require a data.frame.")
    }

    if (!all(cat.vars %in% names(x))) {
      stop("Categorical predictors are missing in new data.")
    }

    for (nm in cat.vars) {
      vals <- as.character(x[[nm]])

      if (any(!is.na(vals) &
              !(vals %in% object$cat.levels[[nm]]))) {
        stop(sprintf("Unknown category in '%s'.", nm))
      }

      x[[nm]] <- factor(
        vals,
        levels = object$cat.levels[[nm]]
      )
    }

    x <- fastDummies::dummy_cols(
      x,
      select_columns = cat.vars,
      remove_selected_columns = TRUE,
      remove_first_dummy = TRUE
    )
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (is.null(colnames(x)) &&
      ncol(x) == length(object$x.names)) {
    colnames(x) <- object$x.names
  }

  if (!identical(colnames(x), object$x.names)) {
    if (is.null(colnames(x)) ||
        !setequal(colnames(x), object$x.names)) {
      stop("Predictor columns do not match training data.")
    }
    x <- x[, object$x.names, drop = FALSE]
  }

  if (any(!is.finite(x))) {
    stop("'x' contains non-finite values.")
  }

  # Standardize predictors
  if (isTRUE(object$standardize)) {

    x <- sweep(x, 2L, object$x.mean, FUN = "-")
    x <- sweep(x, 2L, object$x.sd, FUN = "/")

  }


  # Original response levels
  y.orig <- obj$y_orig

  if (length(unique(y.orig)) != 2L) {
    stop("The response 'y' must have exactly two classes.",
         call. = FALSE)
  }

  level.y <- levels(as.factor(y.orig))

  # Prediction
  if (obj$classification != "ROCSVM") {
    pred.y <- sign(cbind(1, x) %*% curr.b)

    prediction <- apply(pred.y, 2, function(z) {
      factor(
        ifelse(z == -1, obj$negative.class, obj$positive.class),
        levels = level.y
      )
    })


  } else {
    prediction <- round(x %*% curr.b, 4)
  }

  if (length(t.val) == 1L) {
    prediction <- drop(prediction)
  }
  return(prediction)
}

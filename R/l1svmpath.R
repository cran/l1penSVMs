#' The entire regularization path for L1-penalized Support Vector Machine models
#' @description
#' This provides the entire solution paths for L1-penalized Support Vector Machine (SVM), Receiver Operating Characteristic Support Vector Machine (ROC-SVM) and Fraud Detection Support Vector Machine (FD-SVM) with a relatively low cost compared to linear programming problem.
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns); categorical variables are automatically one-hot encoded.
#' @param y Response vector with class labels in \{-1, 1\}. Labels given as \{0, 1\} or as a two-level factor/character are automatically converted to this format.
#' @param classification The classification method is to be applied. Either "SVM", "ROCSVM" or "FDSVM"
#' @param init.t For "SVM" and "FDSVM" models, the initial regularization parameter with sufficiently small t; default is 1e-5
#' @param max.iter The maximum number of iterations (total across the entire path algorithm); default is 10000
#' @param standardize Standardize each column of the data matrix to have zero mean and unit variance
#' @param incomp.size Sample size for the incomplete U-statistic in ROC-SVM; default is the number of observations
#' @param fraud.label.fixed If TRUE, preserves the existing \{-1, +1\} class encoding for FD-SVM.
#' @details
#' The algorithm used in \code{l1.svmpath()} for "SVM" is described in "1-norm Support Vector Machines" by Ji Zhu, Saharon Rosset, Robert Tibshirani and Trevor Hastie. "ROCSVM" is described in "Regularization paths of L1-penalized ROC Curve-Optimizing Support Vector Machines" by Hyungwoo Kim, Insuk Sohn and Seung Jun Shin. "FDSVM" is described in "L1-penalized Fraud Detection Support Vector Machines" by Minhyoung Park, Hyungwoo Kim and Seung Jun Shin. For more information, see references.
#' @return A 'l1.svmpath' object is returned. It contains
#' \describe{
#' \item{t}{The sequence of regularization parameter values in the path}
#' \item{b0}{For "SVM" and "FDSVM" models, the intercept sequence of length \code{length(t)}}
#' \item{beta}{A numeric matrix with \code{length(t)} rows and \code{p} columns, containing coefficient estimates along the regularization path.}
#' \item{x}{The data matrix after dummy variable encoding}
#' \item{y}{The \code{{-1, +1}} valued response variable (For "FDSVM", Fraud +1 and Normal -1)}
#' }
#' @references
#' Ji Zhu, Saharon Rosset, Robert Tibshirani and Trevor Hastie, "1-norm Support Vector Machines", Advances in neural information processing systems, 16.
#'
#' Hyungwoo Kim, Insuk Sohn and Seung Jun Shin, "Regularization paths of L1-penalized ROC Curve-Optimizing Support Vector Machines", Stat, 10(1), e400.
#'
#' Minhyoung Park, Hyungwoo Kim and Seung Jun Shin, "L1-penalized Fraud Detection Support Vector Machines", under
#' @author Seung Jun Shin, MinHyung Park, Hyungwoo Kim
#' @seealso \code{\link{l1.svmpath}}, \code{\link{l1.rocsvmpath}}, \code{\link{l1.fdsvmpath}}, \code{\link{plot.l1svmpath}}, \code{\link{predict.l1svmpath}}
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
#'
#' # ROCSVM
#' obj <- l1svmpath(x,y,classification="ROCSVM")
#'
#' # FDSVM
#' obj <- l1svmpath(x,y,classification="FDSVM")
#' @export
#'
l1svmpath <- function(
    x,
    y,
    classification = c("SVM","ROCSVM","FDSVM"),
    init.t = NULL,
    max.iter = NULL,
    standardize = TRUE,
    incomp.size = NULL,
    fraud.label.fixed = FALSE
    )
{

  ## -----------------------------
  ## Specify SVM model
  ## -----------------------------

  cl <- match.arg(
    classification,
    choices = c("SVM", "ROCSVM", "FDSVM")
  )


  ## -----------------------------
  ## Response variable: y
  ## -----------------------------

  y0 <- y

  # Check the number of classes
  if (length(unique(y)) != 2L) {
    stop(
      "The response 'y' must have exactly two classes.",
      call. = FALSE
    )
  }

  # Check non-finite values
  if (is.numeric(y) && any(!is.finite(y))) {
    stop(
      "The response 'y' contains non-finite values.",
      call. = FALSE
    )
  }

  # Convert y to {-1, +1}
  level.y <- levels(as.factor(y))

  if (setequal(level.y, c("-1", "1"))) {
    y <- as.numeric(as.character(y))
  } else {
    y <- ifelse(as.character(y) == level.y[1L], -1, 1)
  }



  # Save encoded response for cross-validation
  y.cv <- as.numeric(y)

  if (cl == "FDSVM") {
    if (isTRUE(fraud.label.fixed)) {
      if (!setequal(unique(y.cv), c(-1, 1))) {
        stop("'y' must contain both -1 and +1.")
      }
    } else {
      fraud.y <- if (sum(y.cv) >= 0) -1 else 1
      y.cv <- ifelse(y.cv == fraud.y, 1, -1)
    }
  }


  ## -----------------------------
  ## Predictor variables: x
  ## -----------------------------

  # Check input type
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(
      "'x' must be a matrix or data.frame.",
      call. = FALSE
    )
  }

  # Reject list-columns
  if (is.data.frame(x)) {

    has.listcol <- any(
      vapply(x, is.list, logical(1))
    )

    if (has.listcol) {
      stop(
        "'x' must not contain list-columns.",
        call. = FALSE
      )
    }
  }


  ## -----------------------------
  ## Categorical variables
  ## -----------------------------
  cat.vars <- character(0)

  if (is.data.frame(x)) {

    cat.vars <- names(x)[
      vapply(
        x,
        function(z) is.character(z) || is.factor(z),
        logical(1)
      )
    ]
  }

  cat.levels <- lapply(
    x[cat.vars],
    function(z) levels(factor(z))
  )

  if (length(cat.vars) > 0L) {

    message(
      sprintf(
        "Categorical variable(s) detected and one-hot encoded: %s",
        paste(cat.vars, collapse = ", ")
      )
    )

    x <- fastDummies::dummy_cols(
      x,
      select_columns = cat.vars,
      remove_selected_columns = TRUE,
      remove_first_dummy = TRUE
    )
  }



  ## -----------------------------
  ## Convert x to numeric matrix
  ## -----------------------------

  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (any(!is.finite(x))) {
    stop(
      "'x' contains non-finite values after encoding.",
      call. = FALSE
    )
  }


  ## -----------------------------
  ## Jitter discrete predictors
  ## -----------------------------
  x_orig <- x

  jit.ind <- which(
    vapply(
      seq_len(ncol(x)),
      function(j) length(unique(x[, j])) < nrow(x) / 2,
      logical(1)
    )
  )

  if (length(jit.ind) > 0L) {

    x[,jit.ind] <- x[,jit.ind] + rnorm(nrow(x)*length(jit.ind), 0, .05)

  }


  ## -----------------------------
  ## Standardization
  ## -----------------------------
  if (isTRUE(standardize)) {

    x.mean <- colMeans(x)
    x.sd <- apply(x, 2L, sd)

    if (any(x.sd == 0)) {
      stop("'x' contains predictor(s) with zero variance.")
    }

    x <- sweep(x, 2L, x.mean, FUN = "-")
    x <- sweep(x, 2L, x.sd, FUN = "/")

  } else {

    x.mean <- rep(0, ncol(x))
    x.sd <- rep(1, ncol(x))

  }



  ## -----------------------------
  ## Default arguments
  ## -----------------------------

  if (is.null(max.iter)) {
    max.iter <- 1e4
  }

  if (is.null(init.t)) {
    init.t <- 1e-5
  }

  original.levels <- levels(as.factor(y0))

  if (cl == "FDSVM") {
    positive.class <- unique(as.character(y0[y.cv == 1]))[1L]
    negative.class <- unique(as.character(y0[y.cv == -1]))[1L]
  } else {
    positive.class <- original.levels[2L]
    negative.class <- original.levels[1L]
  }



  ## -----------------------------
  ## Fit regularization path
  ## -----------------------------

  if (cl == "SVM") {

    obj <- l1.svmpath(
      x = x,
      y = y,
      init.t = init.t,
      max.iter = max.iter
    )

    result <- list(
      t              = obj$t.seq,
      b0             = obj$b0,
      beta           = obj$beta,
      x              = obj$x,
      y              = obj$y,
      x_orig         = x_orig,
      y_orig         = y0,
      y_cv           = y.cv,
      x.mean         = x.mean,
      x.sd           = x.sd,
      cat.vars       = cat.vars,
      cat.levels     = cat.levels,
      x.names        = colnames(x_orig),
      positive.class = positive.class,
      negative.class = negative.class,
      standardize    = standardize,
      classification = cl
    )

  } else if (cl == "ROCSVM") {

    obj <- l1.rocsvmpath(
      x = x,
      y = y,
      max.iter = max.iter,
      incomp.size = incomp.size
    )

    result <- list(
      t              = na.omit(obj$t),
      beta           = na.omit(obj$beta),
      x              = obj$x,
      y              = obj$y,
      x_orig         = x_orig,
      y_orig         = y0,
      y_cv           = y.cv,
      x.mean         = x.mean,
      x.sd           = x.sd,
      cat.vars       = cat.vars,
      cat.levels     = cat.levels,
      x.names        = colnames(x_orig),
      standardize    = standardize,
      classification = cl
    )

  } else if (cl == "FDSVM") {

    obj <- l1.fdsvmpath(
      x = x,
      y = y.cv,
      init.t = init.t,
      max.iter = max.iter
    )

    result <- list(
      t              = na.omit(obj$t),
      b0             = na.omit(obj$b0),
      beta           = obj$beta,
      x              = obj$x,
      y              = obj$y,
      x_orig         = x_orig,
      y_orig         = y0,
      y_cv           = y.cv,
      x.mean         = x.mean,
      x.sd           = x.sd,
      cat.vars       = cat.vars,
      cat.levels     = cat.levels,
      x.names        = colnames(x_orig),
      positive.class = positive.class,
      negative.class = negative.class,
      standardize    = standardize,
      classification = cl
    )
  }

    structure(result, class="l1svmpath")
  }

#' Fit the L1-penalized Support Vector Machine models at a specific value of t using the the \code{lp} function
#' @description
#' It provides the estimated coefficients for L1-penalized Support Vector Machine (SVM) and Receiver Operating Characteristic Support Vector Machine (ROC-SVM) using linear programming (LP).
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns)
#' @param y The \code{{-1, 1}} valued response variable
#' @param t The value of regularization parameter
#' @param classification The classification method is to be applied. Either "SVM" (the default) or "ROCSVM"
#' @return A numeric vector containing the estimated coefficients
#' at the specified regularization parameter \code{t}.
#' For SVM, the first element is the intercept.
#' For ROC-SVM, only slope coefficients are returned.
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p),n,p)
#' beta  <- c(2,-1,0.5,rep(0,p-3))
#' fx  <- -2 + crossprod(t(x),beta)
#' y   <- sign(fx + rnorm(n,0,1))
#' obj <- lpfit.l1svm(x, y, t=1.5, classification="SVM")
#' @export

lpfit.l1svm <- function(x,
                        y,
                        t,
                        classification=c("SVM", "ROCSVM")){

  classification <- match.arg(classification)

  if (!(is.matrix(x) || is.data.frame(x)))
    stop("'x' must be a matrix or data.frame.", call. = FALSE)

  if (!is.numeric(y))
    stop("The response 'y' must be numeric or a 2-class factor/character.", call. = FALSE)
  if (any(!is.finite(y)))
    stop("y contains non-finite values after coercion.", call. = FALSE)
  if (setequal(unique(y), c(0, 1))) {
    warning("'y' contains {0, 1} labels; converting to {-1, 1}.", call. = FALSE)
    y <- ifelse(y == 0, -1, 1)
    label.y <- unique(y)
  }
  if (!setequal(unique(y), c(-1, 1)))
    stop("'y' must contain only -1 and 1.", call. = FALSE)


  # Find the initial solution
  if( classification == "SVM" ){
    n <- nrow(x)
    p <- ncol(x)
    X <- c(y) * cbind(1,x)

    r1 <- cbind(X, -X, diag(1,n))
    r2 <- diag(1, 2*(p+1)+n)
    r3 <- rbind(r1, r2)
    r4 <- rbind(r3, c(rep(c(0,rep(1,p)),2), rep(0,n)))


    f.obj <- c(rep(0,2*(p+1)),rep(1,n))
    f.con <- r4
    f.dir <- c(rep(">=",2*(n+p+1)),"<=")
    f.rhs <- c(rep(1,n),rep(0,n+2*(p+1)),t)
    par <- lp("min",f.obj,f.con,f.dir,f.rhs)$sol

    beta <- c(par[1:(p+1)] - par[(p+2):(2*(p+1))])

  } else{

    i_1 <- which(y==1)
    i_2 <- which(y==-1)
    nm <- length(i_1)*length(i_2)
    p <- ncol(x)

    # design matrix
    i_1.1 <- matrix(rep(x[i_1,],length(i_2)), ncol=length(i_1), byrow=T)
    i_1.2 <- matrix(c(i_1.1),nrow=p)
    i_2.2 <- matrix(c(rep(t(x[i_2,]),length(i_1))),nrow=p)

    M.ij <- t(i_1.2 - i_2.2)
    Xmat <- cbind(M.ij,-M.ij)

    d1 <- c(rep(0,2*p), rep(1,nm))
    dvec <- d1

    M1 <- cbind(Xmat,diag(1,nm))
    M2 <- rbind(M1, diag(1,nm+2*p))
    M3 <- rbind(M2, c(rep(1,2*p),rep(0,nm)))

    f.obj <- dvec
    f.con <- as.matrix(M3)
    f.dir <- c(rep(">=", nm), rep(">=", 2*p+nm), "<=")
    f.rhs <- c(rep(1,nm), rep(0,2*p+nm), t)

    sol  <- lp ("min", f.obj, f.con, f.dir, f.rhs)$sol
    beta <- sol[1:p]-sol[(p+1):(2*p)]

  }
  return( beta )
}






#' Finding the initial solution for L1-penalized Fraud Detection Support Vector Machine
#' @description This algorithm provides the initial solution of L1-penalized Fraud Detection Support Vector Machine (FD-SVM) for sufficiently small t.
#' @param x The data matrix after dummy variable encoding
#' @param y The \code{{-1, 1}} valued response variable (Fraud +1, normal -1)
#' @param init.t Initial regularization parameter with sufficiently small t; default is 1e-5
#' @return The intercept and coefficients at initial regularization parameter
#' @export
#' @import lpSolve
l1.fdsvm.initialization <- function(x, y, init.t) {

  if (nrow(x)!=length(y)) stop("length(y) must match nrow(x)")

  # Fraud class is already encoded as +1
  if (!setequal(unique(y), c(-1, 1))) {
    stop("'y' must contain both -1 and +1.")
  }

  id_p <- which(y == 1)
  id_n <- which(y == -1)

  xp <- x[id_p, , drop = FALSE]
  xn <- x[id_n, , drop = FALSE]

  x <- rbind(xp, xn)
  yp <- rep(1, nrow(xp))
  yn <- rep(-1, nrow(xn))
  y <- c(yp, yn)

  weight = rep_len(1,nrow(xn))
  yy_p <- matrix(rep(yp, each = ncol(xp)), nrow = nrow(xp), ncol = ncol(xp), byrow= T)
  yy_n <- matrix(rep(yn, each = ncol(xn)), nrow = nrow(xn), ncol = ncol(xn), byrow= T)


  opt = lp(direction = "min",
           objective.in = c(apply(xn, 2, sum),   # no loss in fraud class but only positive loss in majority class!
                            -apply(xn, 2, sum),
                            nrow(xn),
                            -nrow(xn)),
           const.mat = rbind(cbind(-yy_n*xn, yy_n*xn, -yn, yn),
                             cbind(yy_p*xp, -yy_p*xp, yp, -yp),
                             diag(2*(ncol(xp)+1)),
                             c(rep(-1, 2* ncol(x)), rep(0, 2))),
           const.dir = c(rep(">=", (nrow(x) +2*(ncol(x)+1) + 1))), # delta s is equal for strong duality
           const.rhs = c(rep(-1,nrow(xn)),rep(1,nrow(xp)), rep(0,(2*(ncol(xn)+1))), -init.t)
  )
  sol <- opt$solution

  u = opt$solution[seq(1L,length.out=ncol(x))]
  v = opt$solution[seq(ncol(x)+1,length.out=ncol(x))]
  b = opt$solution[2*ncol(x)+1] - opt$solution[2*ncol(x)+2]
  w = c(b,u-v)
  list(w)
}

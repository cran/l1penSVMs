#' Finding the initial solution for L1-penalized SVM
#' @description This algorithm provides the initial solution of L1-penalized SVM for sufficiently small t.
#' @param x The data matrix after dummy variable encoding
#' @param y The \code{{-1, 1}} valued response variable
#' @param init.t Initial regularization parameter with sufficiently small t; default is 1e-5
#' @return The intercept and coefficients  at initial regularization parameter
#' @export
#' @import lpSolve
l1.svm.initialization <- function(x,y,init.t) {

  if (nrow(x)!=length(y)) stop("length(y) must match nrow(x)")

  plus <- if (sum(y) >= 0) 1 else -1
  id_p <- c(y==plus)
  id_n <- c(y!=plus)
  xp <- x[id_p,,drop=FALSE]
  xn <- x[id_n,,drop=FALSE]
  x <- rbind(xp, xn)
  yp <- rep(1, nrow(xp))
  yn <- rep(-1, nrow(xn))
  y <- c(yp,yn)

  weight = rep_len(1,nrow(xp))
  yy_p <- matrix(rep(yp, each = ncol(xp)), nrow = nrow(xp), ncol = ncol(xp), byrow= T)
  yy_n <- matrix(rep(yn, each = ncol(xn)), nrow = nrow(xn), ncol = ncol(xn), byrow= T)


  opt = lp(direction = "min",
           objective.in = c(apply(xn, 2, sum),
                            -apply(xn, 2, sum),
                            nrow(xn),
                            -nrow(xn),
                            weight),
           const.mat = rbind(cbind(yy_p*xp, -yy_p*xp, yp, -yp, diag(nrow(xp))),
                             cbind(-yy_n*xn, yy_n*xn, -yn, yn, matrix(rep(0, len=nrow(xn)*nrow(xp)), nrow = nrow(xn))),
                             diag(2*(ncol(xp)+1)+nrow(xp)),
                             c(rep(-1, 2* ncol(xp)), rep(0, 2+nrow(xp)))),
           const.dir = c(rep(">=", (nrow(x) +2*(ncol(x)+1) + nrow(xp)) ), ">="), # delta s is equal for strong duality
           const.rhs = c(rep(1,nrow(xp)),rep(-1,nrow(xn)), rep(0,(2*(ncol(xp)+1)+nrow(xp))), -init.t)
  )
  sol <- opt$solution

  u = opt$solution[seq(1L,length.out=ncol(x))]
  v = opt$solution[seq(ncol(x)+1,length.out=ncol(x))]
  b = opt$solution[2*ncol(x)+1] - opt$solution[2*ncol(x)+2]

  w = c(b,u-v)

  if(plus == 1) {
    list(w)
  } else {
    list(-w)
  }

}

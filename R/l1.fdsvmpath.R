#' Fit the entire solution paths for L1-penalized Fraud Detection Support Vector Machine
#' @description
#' This algorithm computes the entire solution paths for L1-penalized FDSVM with a relatively low cost compared to linear programming problem.
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns)
#' @param y The \code{{-1, 1}} valued response variable (Fraud +1, normal -1)
#' @param init.t Initial regularization parameter with sufficiently small t; default is 1e-5
#' @param max.iter The maximum number of iterations (total across the entire path algorithm); default is 3000
#' @return A 'l1.fdsvmpath' object is returned. It contains
#' \describe{
#' \item{t}{The sequence of regularization parameter values in the path}
#' \item{b0}{The intercept sequence of length \code{length(t)}}
#' \item{beta}{A \code{p*length(t)} matrix of coefficients.}
#' }
#' @references
#' Minhyoung Park, Hyungwoo Kim and Seung Jun Shin, "L1-penalized Fraud Detection Support Vector Machines", under
#' @author Seung Jun Shin, MinHyung Park, Hyungwoo Kim
#' @seealso \code{\link{l1.svmpath}}, \code{\link{l1.fdsvm.initialization}}, \code{\link{plot.l1svmpath}}
#' @examples
#' library(l1penSVMs)
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p),n,p)
#' beta  <- c(2,-1,0.5,rep(0,p-3))
#' fx  <- -2 + crossprod(t(x),beta)
#' y   <- sign(fx + rnorm(n,0,1))
#' obj <- l1svmpath(x,y,classification="FDSVM")
#' @export
l1.fdsvmpath <- function(x, y, init.t, max.iter){


  if (!setequal(unique(y), c(-1, 1))) {
    stop("'y' must contain both -1 and +1.")
  }

  id_p <- which(y == 1)
  id_n <- which(y == -1)

  xp <- x[id_p,,drop=FALSE]
  xn <- x[id_n,,drop=FALSE]
  x <- rbind(xp, xn)
  yp <- rep(1, nrow(xp))
  yn <- rep(-1, nrow(xn))
  y <- c(yp,yn)

  n <- nrow(x)
  p <- ncol(x)

  beta.init <- l1.fdsvm.initialization(x,y,init.t)[[1]]
  beta.init

  # initial w
  eps <- 1e-8

  w <- beta.init[-1] / init.t
  w[(abs(w) < eps)] <- 0
  w.0 <- (beta.init[1] + 1) / init.t
  beta <- c(-1, rep(0, p))
  S <- which(w != 0)     # index of on-zero beta
  t <- 0

  # init.set
  E <- seq(y)[abs(1- y * cbind(1, x) %*% beta.init) < eps]
  R <- seq(y)[(y * cbind(1, x) %*% beta.init - 1) > eps]
  L <- seq(y)[(y * cbind(1, x) %*% beta.init - 1) < -eps]

  temp.beta2 <- matrix(NA, max.iter, p+1)
  temp.t2 <- rep(NA, max.iter)


  for(z in 1:max.iter){

    ### Step1

    if(length(S) == 1) {

      # hit Elbow candidate for finding new joint
      fx.candidate <- cbind(1, x) %*% beta   ### All the majority points at E, thus no need to have delta to enter E
      t.cand.E <- (1-y*fx.candidate) / (y*(cbind(1, x) %*% c(w.0, w)))
      t.cand.n.E <- t.cand.E[-E]

      if (length(t.cand.n.E[t.cand.n.E > eps]) == 0){
        break
      }

      t.E <- min(t.cand.n.E[t.cand.n.E > eps]) # delta.s must non-negative and larger than eps

      delt.t  <- t.E

      if (delt.t == Inf){
        break
      }

    } else {

      # hit Elbow candidate for finding new joint
      fx.candidate <- cbind(1, x) %*% beta
      t.cand.E <- (1-y*fx.candidate) / (y*(cbind(1, x) %*% c(w.0, w)))
      t.cand.n.E <- t.cand.E[-E]


      if (length(t.cand.n.E[t.cand.n.E > eps]) == 0){
        break
      }

      t.E <- min(t.cand.n.E[t.cand.n.E > eps]) # delta.s must non-negative and larger than eps


      # basis leaves
      t.cand.S <- -(beta[-1] / w)
      t.cand.S.ind <- which(t.cand.S > eps)

      if(length(t.cand.S.ind) != 0){
        t.S <- min(t.cand.S[t.cand.S.ind])
      } else{
        t.S <- Inf
      }


      # finding minimum delta s
      delt.t <- min(t.E, t.S)

      if (delt.t == Inf){
        break
      }

    }

    if(t>20) break
    if( max(abs(beta)) > 1e3 ) break


    # save delta.s and beta
    temp.beta2[z,] <- beta
    temp.t2[z] <- t


    # Update delta.s and beta
    t <- t + delt.t
    beta[1] <- beta[1] + delt.t*w.0  # beta.0 hat old
    beta[-1] <- beta[-1] + delt.t*w  # beta.j hat old
    beta[abs(beta) < eps] <- 0  # replace zero beta if less than eps

    if (delt.t == 0){
      break
    }

    # Update sets
    if(delt.t == t.E){

      E.ask <- unique(c(E, which(t.cand.E == t.E)[1]))
      L.ask <- setdiff(L, E.ask)
      S.ask <- S


    } else{
      E.ask <- E
      L.ask <- L
      S.ask <- setdiff(S, which(t.cand.S == t.S))
    }

    ### violate KKT condition
    if( length(L.ask) == 0) break
    if( length(S.ask)+1 != length(E.ask)) break

    ### Step2 (A2 in thesis)
    beta.j <- beta[-1]
    not.S <- setdiff(1:p,S.ask) ### j candidate number of variable
    A2 <- matrix(numeric(),length(not.S), p+3)
    #colnames(step2) <- c("jstar", "loss.delta", "w.0", rep("w.j", p))

    if(length(not.S) > 0){
      for (j in 1:length(not.S)){

        temp.A <- cbind(1, x[ c(E.ask) , c(S.ask, not.S[j])])
        ### w.jstar(+)
        A.p <- rbind(temp.A, c(0, sign(beta.j[S.ask]), 1))
        b <- c(rep(0,length(E.ask)),1)

        w.A2.p <- rep(0, p+1)
        w.A2.p[c(1, c(S.ask, not.S[j])+1)] <- solve(A.p + eps*diag(rep(1, nrow(A.p))), b)

        loss.A2.p <- -sum(y[L.ask]*(cbind(1,x[L.ask, c(S.ask,  not.S[j]), drop=F])) %*% w.A2.p[c(1,c(S.ask,  not.S[j])+1)])
        loss.A2.p

        if(w.A2.p[not.S[j]+1] < 0){
          loss.A2.p <- 1
        }


        ### w.jstar(-)
        A.n <- rbind(temp.A, c(0, sign(beta.j[S.ask]), -1))

        w.A2.n <- rep(0, p+1)
        w.A2.n[c(1, c(S.ask, not.S[j])+1)] <- solve(A.n + eps*diag(rep(1, nrow(A.n))),b)

        loss.A2.n <- -sum(y[L.ask]*(cbind(1,x[L.ask, c(S.ask,  not.S[j]), drop=F])) %*% w.A2.n[c(1,c(S.ask,  not.S[j])+1)])
        loss.A2.n

        if(w.A2.n[not.S[j]+1] > 0){
          loss.A2.n <- 1
        }


        ### u.jstar
        if(loss.A2.p < loss.A2.n) {
          w.A2 <- w.A2.p
          loss.A2 <- loss.A2.p
        } else {
          w.A2 <- w.A2.n
          loss.A2 <- loss.A2.n
        }

        A2[j,] <- c(not.S[j], loss.A2, w.A2)
      }
    } else {
      A2 <- matrix(rep(0,p+3),1, p+3)
    }


    A2

    loss.A2 <- A2[,2]
    minstar.A2 <- which.min(loss.A2)
    A2.candidate <- A2[minstar.A2, ]
    A2.candidate


    ### Step3 (A1)
    A1 <- matrix(numeric(),length(E.ask), p+3)
    #colnames(step3) <- c("iprime", "loss.delta", "u.0", rep("u.j", p))

    for (i in 1:length(E.ask)){

      temp.A <- cbind(1, x[ E.ask[-i] , S.ask])
      A <- rbind(temp.A, c(0, sign(beta.j[S.ask])))
      b <- c(rep(0,length(E.ask)-1),1)

      w.A1 <- rep(0, p+1)
      w.A1[c(1, c(S.ask)+1)] <- solve(A + eps*diag(rep(1, nrow(A))), b)

      # E on the fraud class must not to move to L
      if(y[E.ask[i]] == 1
         && (y[E.ask[i]]*(cbind(1,x[E.ask[i], S.ask, drop = FALSE]) %*% w.A1[c(1, S.ask+1)])) < 0) {
        loss.A1 <- 1
      } else {
        loss.A1 <- ifelse( (y[E.ask[i]]*(cbind(1,x[E.ask[i], S.ask, drop = FALSE]) %*% w.A1[c(1, S.ask+1)])) < 0,  # Elbow to Left condition. If each E have same loss value, leave the first existing E
                           -sum(y[c(E.ask[i], L.ask)]*(cbind(1, x[c(E.ask[i], L.ask), S.ask, drop=F])) %*% w.A1[c(1, S.ask+1)]),
                           -sum(y[c(L.ask)]*(cbind(1, x[L.ask, S.ask, drop=F])) %*% w.A1[c(1, S.ask+1)]))
      }

      A1[i,] <- c(E.ask[i], loss.A1, w.A1)
    }
    A1

    loss.A1 <- A1[,2]
    minstar.A1 <- which.min(loss.A1)

    A1.candidate <- A1[minstar.A1, ]
    A1.candidate

    ### Step4 Compare loss (A1) and (A2)

    loss.tot <- min(min(A1.candidate[2]), min(A2.candidate[2]))

    if (loss.tot >= 0)
      break
    if (loss.tot == A2.candidate[2]) {

      ### In (A2)
      # Find minimum
      A2.ind <- A2.candidate[1]

      # update derivative
      w.0 <- A2.candidate[3]
      w <- A2.candidate[-c(1:3)]

      # update set
      E <- E.ask
      L <- L.ask
      S <- c(S.ask, A2.ind)

    } else {

      ### In (A1)
      # Find minimum
      A1.ind <- A1.candidate[1]

      # update derivative
      w.0 <- A1.candidate[3]
      w <- A1.candidate[-c(1:3)]

      # update set
      E <- E.ask[E.ask != A1.ind]
      # Determine L from E
      if ( (y[A1.candidate[1]] * cbind(1, x[A1.candidate[1],,drop = FALSE]) %*% c(w.0,w)) < 0){  # Elbow to Left condition
        L <- c(L.ask, A1.ind)
      } else {
        L <- L.ask
      }
      S <- S.ask

    }
  }


  ### Go back to step1
  temp.t <- na.omit(temp.t2)
  temp.beta <- na.omit(temp.beta2)




  result <- structure(list(b0 = temp.beta[,1],
                           beta = temp.beta[,-1, drop = FALSE],
                           t = temp.t,
                           x=x, y=y), class="l1.fdsvmpath")

  return(result)
}

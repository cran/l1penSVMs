#' Fit the entire solution paths for L1-penalized Support Vector Machine
#' @description
#' This algorithm computes the entire solution paths for L1-penalized  Support Vector Machine (SVM)
#' @param x The data matrix after dummy variable encoding
#' @param y The \code{{-1, 1}} valued response variable
#' @param init.t The initial regularization parameter with sufficiently small t; default is 1e-5
#' @param max.iter The maximum number of iterations (total across the entire path algorithm); default is 3000.
#' @return An object of class \code{l1.svmpath}. It contains
#' \describe{
#' \item{b0}{A vector of intercept estimates for different values of t }
#' \item{beta}{ A \code{p*length(t)} matrix of coefficients}
#' \item{t}{The sequence of regularization parameter values in the path}
#' \item{x}{The data matrix after dummy variable encoding }
#' \item{y}{The \code{{-1, 1}} valued response variable }
#' }
#' @references
#' Hyungwoo Kim, Insuk Sohn and Seung Jun Shin, "Regularization paths of L1-penalized ROC Curve-Optimizing Support Vector Machines", Stat, 10(1), e400.
#' @author Seung Jun Shin, MinHyung Park, Hyungwoo Kim
#' @seealso \code{\link{l1svmpath}}, \code{\link{plot.l1svmpath}}
#' @examples
#' library(l1penSVMs)
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p),n,p)
#' beta  <- c(2,-1,0.5,rep(0,p-3))
#' fx  <- -2 + crossprod(t(x),beta)
#' y   <- sign(fx + rnorm(n,0,1))
#' obj <- l1svmpath(x,y,classification="SVM")
#' @export

l1.svmpath <- function(x, y, init.t, max.iter){

  x <- as.matrix(x)
  y <- y
  n <- nrow(x)
  p <- ncol(x)

  if( is.null(init.t) ) init.t <- 1e-5
  beta.init = lpfit.l1svm(x, y, init.t, classification="SVM")

  if( sum(abs(abs(beta.init[-1]))) < init.t/2 ) beta.init <- lpfit.l1svm(x, y, t=init.t, classification="SVM")

  S <- which(beta.init[-1] != 0)
  S.true <- ifelse(length(S)==1, TRUE, prod(sapply(apply(x[,S], 2, unique), length) == n) == 0)

  if( S.true ){
     x <- as.matrix(x)
     beta.init <- lpfit.l1svm(x, y, init.t, classification="SVM")
  }

  # initial w
  w <- beta.init[-1] / init.t
  w.0 <- (beta.init[1] - 1) / init.t
  beta <- c(1 , rep(0, p))
  S <- which(w != 0)     # index of on-zero beta
  t <- 0
  eps <- 1e-10

  E <- seq(y)[abs(1- y * cbind(1, x) %*% beta.init) < eps]
  R <- seq(y)[(1- y * cbind(1, x) %*% beta.init) < -eps]
  L <- seq(y)[(1- y * cbind(1, x) %*% beta.init) > eps]

  temp.beta <- matrix(NA, max.iter, p+1)
  temp.t <- rep(NA, max.iter)


  for(z in 1:max.iter){

    if(length(S) == 1) {

      # hit Elbow candidate for finding new joint
      fx.candidate <- cbind(1, x) %*% beta   ### All the majority points at E, thus no need to have delta to enter E
      t.cand.E <- (1-y*fx.candidate) / (y*(cbind(1, x) %*% c(w.0, w)))
      t.cand.n.E <- t.cand.E[-E]

      # if( sum(t.cand.n.E > eps)==0)
      t.E <- min(t.cand.n.E[t.cand.n.E > eps]) # delta.s must non-negative and larger than eps

      delt.t <- t.E

    } else {

      # hit Elbow candidate for finding new joint
      fx.candidate <- cbind(1, x) %*% beta   ### All the majority points at E, thus no need to have delta to enter E
      t.cand.E <- (1-y*fx.candidate) / (y*(cbind(1, x) %*% c(w.0, w)))
      t.cand.n.E <- t.cand.E[-E]

      t.E <- min(t.cand.n.E[t.cand.n.E > eps]) # delta.s must non-negative and larger than eps

      delt.t <- t.E

      # basis leaves
      t.cand.S <- -(beta[-1] / w)
      t.cand.S.ind <- which(t.cand.S > eps)


      if(length(t.cand.S.ind) != 0){
        t.S <- min(t.cand.S[t.cand.S.ind])
      } else{
        t.S <- Inf
      }

      # stop iteration case
      # if(z>10){
      #   if( abs(temp.t[z-1]-temp.t[z-10])<1e-19 ){
      #     t.cand.n.E <- t.cand.E[-E]
      #     t.E <- min(t.cand.n.E[t.cand.n.E>1e-15])
      #   }
      # }

      # finding minimum delta s
      delt.t <- min(t.E, t.S)

      if (delt.t == Inf){
        break
      }
    }

    # save delta_s and beta
    if(z == 1) {
      temp.beta[z,] <- beta.init
      temp.t[z] <- init.t
    } else {
      temp.beta[z,] <- beta
      temp.t[z] <- t
    }


    # Update delta_s and beta
    t <- t + delt.t
    beta[1] <- beta[1] + delt.t*w.0  # beta_0 hat old
    beta[-1] <- beta[-1] + delt.t*w  # beta_j hat old
    beta[abs(beta) < eps] <- 0  # replace zero beta if less than eps

    if (delt.t == 0){
      break
    }
    if (t>20){ break}

    # Update sets
    if(delt.t == t.E){

      add.E <- which(t.cand.E == t.E)
      E.ask <- unique(c(E, add.E))
      L.ask <- setdiff(L, E.ask)
      S.ask <- S

    } else{
      E.ask <- E
      L.ask <- L
      S.ask <- setdiff(S, which(t.cand.S == t.S))
    }

    ### violate KKT condition
    if( length(L.ask) == 0) break
    if( length(S.ask)+1 != length(E.ask) ) break

    ### Step2 (A2 in thesis)
    beta.j <- beta[-1]
    not.S <- setdiff(1:p,S.ask) ### j candidate number of variable
    A2 <- matrix(numeric(),length(not.S), p+3)

    if(length(not.S) >0 ){
      for (j in 1:length(not.S)){

        temp.A <- cbind(1, x[ c(E.ask) , c(S.ask, not.S[j])])
        ### w_jstar(+)
        A.p <- rbind(temp.A, c(0, sign(beta.j[S.ask]), 1))
        b <- c(rep(0,length(E.ask)),1)

        w.A2.p <- rep(0, p+1)
        w.A2.p[c(1, c(S.ask, not.S[j])+1)] <- solve(A.p + eps*diag(rep(1, nrow(A.p))), b)

        loss.A2.p <- -sum(y[L.ask]*(cbind(1,x[L.ask, c(S.ask,  not.S[j]), drop=F])) %*% w.A2.p[c(1,c(S.ask,  not.S[j])+1)])
        loss.A2.p

        if(w.A2.p[not.S[j]+1] < 0){
          loss.A2.p <- 1
        }

        ### w_jstar(-)

        A.n <- rbind(temp.A, c(0, sign(beta.j[S.ask]), -1))

        w.A2.n <- rep(0, p+1)
        w.A2.n[c(1, c(S.ask, not.S[j])+1)] <- solve(A.n + eps*diag(rep(1, nrow(A.n))),b)

        loss.A2.n <- -sum(y[L.ask]*(cbind(1,x[L.ask, c(S.ask,  not.S[j]), drop=F])) %*% w.A2.n[c(1,c(S.ask,  not.S[j])+1)])
        loss.A2.n

        if(w.A2.n[not.S[j]+1] > 0){
          loss.A2.n <- 1
        }

        ### u_jstar

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
    loss.A2 <- A2[,2]
    minstar.A2 <- which.min(loss.A2)
    A2.candidate <- A2[minstar.A2, ]
    A2.candidate

    ### Step3 (A1)
    A1 <- matrix(numeric(),length(E.ask), p+3)

    for (i in 1:length(E.ask)){

      temp.A <- cbind(1, x[ E.ask[-i] , S.ask])
      A <- rbind(temp.A, c(0, sign(beta.j[S.ask])))
      b <- c(rep(0,length(E.ask)-1),1)

      w.A1 <- rep(0, p+1)
      w.A1[c(1, c(S.ask)+1)] <- solve(A + eps*diag(rep(1, nrow(A))), b, tol=1e-16)


      loss.A1 <- ifelse((y[E.ask[i]]*(cbind(1,x[E.ask[i], S.ask, drop = FALSE]) %*% w.A1[c(1, S.ask+1)])) < 0,  # Elbow to Left condition. If each E have same loss value, leave the first existing E
                        -sum(y[c(E.ask[i], L.ask)]*(cbind(1, x[c(E.ask[i], L.ask), S.ask, drop=F])) %*% w.A1[c(1, S.ask+1)]),
                        -sum(y[c(L.ask)]*(cbind(1, x[L.ask, S.ask, drop=F])) %*% w.A1[c(1, S.ask+1)]))

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
  temp.t <- na.omit(temp.t)
  temp.beta <- na.omit(temp.beta)

  ### Remove the unusual value at the last t
  t.length <- length(temp.t)
  if( max(abs(temp.beta[t.length,])) > 1e4 ){
    temp.t <- temp.t[-t.length]
    temp.beta <- temp.beta[-t.length, , drop = FALSE]
  }

  temp.beta[abs(temp.beta) < 1e-6] <- 0

  structure(list(b0 = temp.beta[,1], beta = temp.beta[,2:(p+1), drop = FALSE],
                 t.seq = temp.t, x=x, y=y), class="l1.svmpath")

}

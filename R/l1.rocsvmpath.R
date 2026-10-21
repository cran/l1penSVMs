#' Fit the entire solution paths for L1-penalized Receiver Operating Characteristic Support Vector Machine
#' @description
#' This algorithm computes the entire solution paths for L1-penalized Receiver Operating Characteristic Support Vector Machine (ROC-SVM)
#' @param x The data matrix after dummy variable encoding
#' @param y The \code{{-1, 1}} valued response variable
#' @param max.iter The maximum number of iterations (total across the entire path algorithm); default is 3000.
#' @param incomp.size Sample Size for the incomplete U-statistic in ROC-SVM; default is 2n.
#' @return A 'l1.rocsvmpath' object is returned. It contains
#' \describe{
#' \item{beta}{A \code{p*length(t)} matrix of coefficients.}
#' \item{t}{The sequence of regularization parameter values in the path}
#' \item{x}{ The data matrix after dummy variable encoding }
#' \item{y}{ The \code{{-1, 1}} valued response variable }
#' }
#' @references
#' Hyungwoo Kim, Insuk Sohn and Seung Jun Shin, "Regularization paths of L1-penalized ROC Curve-Optimizing Support Vector Machines", Stat, 10(1), e400.
#' @author Seung Jun Shin, MinHyung Park, Hyungwoo Kim
#' @seealso \code{\link{l1.svmpath}}, \code{\link{plot.l1svmpath}}
#' @examples
#' library(l1penSVMs)
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n*p),n,p)
#' beta  <- c(2,-1,0.5,rep(0,p-3))
#' fx  <- -2 + crossprod(t(x),beta)
#' y   <- sign(fx + rnorm(n,0,1))
#' obj <- l1svmpath(x,y,classification="ROCSVM")
#' @export

l1.rocsvmpath <- function(x, y, max.iter, incomp.size){

  x <- as.matrix(x)
  y <- y
  n <- nrow(x)
  p <- ncol(x)

  i_1 <- which(y==1)
  i_2 <- which(y==-1)
  nm  <- c(length(i_1)*length(i_2))

  if( is.null(max.iter) ) max.iter = 3000

  # Sampling the training data for Incomplete U-statistic
  B <- length(i_1) * length(i_2)

  if (is.null(incomp.size)) {
    incomp.size <- min(2L * n, B)
  }

  if (length(incomp.size) != 1L ||
      !is.finite(incomp.size) ||
      incomp.size < 1 ||
      incomp.size > B ||
      incomp.size != floor(incomp.size)) {
    stop("'incomp.size' must be an integer between 1 and B.")
  }



  #### Initialization
  samp_ind <- sample.int(B, size = incomp.size)

  #### Initialization
  # Design matrix
  i_1.1 <- matrix(rep(x[i_1,],length(i_2)), ncol=length(i_1), byrow=T)
  i_1.2 <- matrix(c(i_1.1),nrow=p)
  i_2.2 <- matrix(c(rep(t(x[i_2,]),length(i_1))),nrow=p)
  x.ij0 <- t(i_1.2 - i_2.2)
  x.ij <- x.ij0[samp_ind, ] # Use n pairs for approximation


  # Find init variable and set
  init.opt <- apply(x.ij,2,sum)
  init.ind <- which.max(abs(init.opt))
  init.ind.sign <- sign(init.opt)[init.ind]


  S <- init.ind
  if( n != length(unique(x[,S])) ){

    S = c(which(sapply(apply(x, 2, unique), length)< n/3),S)
    x[,S] <- x[,S] + rnorm(nrow(x),0,0.01)
    # x <- as.matrix(scale(x))
    x <- as.matrix(x)

    # Design matrix
    i_1.1 <- matrix(rep(x[i_1,],length(i_2)), ncol=length(i_1), byrow=T)
    i_1.2 <- matrix(c(i_1.1),nrow=p)
    i_2.2 <- matrix(c(rep(t(x[i_2,]),length(i_1))),nrow=p)
    x.ij <- t(i_1.2 - i_2.2)[samp_ind, ]
  }

  # Find derivative
  w <- rep(0,p)
  w[init.ind] <- init.ind.sign

  t <- 0
  E <- NULL
  L <- 1:nrow(x.ij)
  S <- init.ind
  beta <- rep(0,p)

  temp.beta <- matrix(NA, max.iter, p)
  temp.t   <- rep(NA, max.iter)



  # Iteration
  for(z in 1:max.iter){

    if(length(S) == 1){

      t.cand.E <- (1-x.ij%*%beta) / (x.ij%*%w)
      t.E <- min(t.cand.E[t.cand.E>0])
      delt.t  <- t.E

    } else{

      # Event 1. Non-elbow hits elbow
      t.cand.E   <- (1-x.ij%*%beta)  / (x.ij%*%w)
      t.cand.n.E <- t.cand.E[-E]

      t.E <- if(length(t.cand.n.E[t.cand.n.E>1e-8]) == 0) break
      t.E <- min(na.omit(t.cand.n.E[t.cand.n.E>1e-8]))

      # Event 2, S becomes not.S
      t.cand.S   <- -beta/w
      t.cand.S.ind <- which(t.cand.S>0)

      if(length(t.cand.S.ind) != 0){
        t.S <- min(t.cand.S[t.cand.S.ind])
      } else{
        t.S <- 1e4
      }

      # Find t
      delt.t <- min(t.S,t.E)
    }


    # save t and beta
    temp.beta[z,] <- beta
    temp.t[z]   <- t

    # Update t and beta
    t <- t + delt.t
    beta <- beta + delt.t*w
    #if( (delt.t == t.E) & (sum(t.cand.E == t.E) != 1)  ){


    # Update sets
    if(delt.t == t.E){

      E.ask <- unique(c(E, which(t.cand.E == t.E) ))
      L.ask <- setdiff(L, E.ask)
      S.ask <- S

    } else{
      E.ask <- E
      L.ask <- L
      S.ask <- setdiff(S, which(t.cand.S == t.S))
    }

    if( length(L.ask) == 0 ){
      L.ask <- which((1-x.ij%*%beta)>0)
    }

    if( length(L.ask) == 0) break
    #if( length(S.ask) != length(E.ask)) break

    ### Take action  (A1) - (A2)
    not.S <- setdiff(1:p,S.ask)


    # Action (A1)
    loss.A1 <- rep(NA, length(S.ask))
    w.A1 <- matrix(NA, length(E.ask), length(S.ask))


    if(length(E.ask) == 1){

      w.A1 <- sign(beta[S.ask])
      loss.A1 <- ifelse( x.ij[E.ask,S.ask]%*%w.A1 < 0, -sum(x.ij[c(E.ask,L.ask), S.ask] * w.A1),
                         -sum(x.ij[L.ask,S.ask] * w.A1) )

    } else{

      for(j in 1:length(E.ask)){
        ind1 <- E.ask[-j]
        ind2 <- E.ask[j]

        e1 <- x.ij[ind1, S.ask]
        e2 <- rbind(e1,sign(beta[S.ask]))

        # Find solution
        sol.e <- tryCatch(solve(e2, c(rep(0,length(E.ask)-1),1)), error = function(t) NA)
        if(is.na(sol.e[1])) break

        w.A1[j,] <- sol.e
        loss.A1[j] <- ifelse( (x.ij[ind2,S.ask] %*% sol.e) < 0,  -sum(x.ij[c(ind2,L.ask), S.ask] %*% sol.e),
                              -sum(x.ij[L.ask, S.ask]%*%sol.e))
      }
    }
    if(sum(is.na(loss.A1)) != 0) break
    # action (A2)
    loss.A2 <- rep(NA, length(not.S))
    w.A2 <- matrix(NA, length(not.S), length(E.ask)+1)
    iter <- 1

    for(k in not.S){

      f1 <- x.ij[E.ask, c(S.ask,k)]
      f2 <- rbind(f1, c(sign(beta)[S.ask], 1))

      g1 <- x.ij[E.ask, c(S.ask,k)]
      g2 <- rbind(g1, c(sign(beta)[S.ask],-1))

      sol.f <- tryCatch(solve(f2, c(rep(0,length(E.ask)),1)), error = function(t) ginv(f2) %*% c(rep(0, length(E.ask)), 1))
      sol.g <- sol.g <- tryCatch(
        solve(g2, c(rep(0, length(E.ask)), 1)),
        error = function(e) MASS::ginv(g2) %*%
          c(rep(0, length(E.ask)), 1)
      )

      loss.f <- -sum( x.ij[L.ask, c(S.ask,k)] %*% sol.f )
      loss.g <- -sum( x.ij[L.ask, c(S.ask,k)] %*% sol.g )

      if(sol.f[length(sol.f)]<0) loss.f = 0
      if(sol.g[length(sol.g)]>0) loss.g = 0

      w.A2[iter,] <- ifelse(rep(loss.f<loss.g, length(sol.f)), sol.f, sol.g)
      loss.A2[iter] <- min(loss.f,loss.g)

      iter <- iter + 1
    }

    # Compare loss (A1) and (A2)
    loss.tot <- c(loss.A1, loss.A2)

    # If all positive, terminate.
    if( (sum(loss.tot>=1e-8) == length(loss.tot)) | (t>20) )
      break

    ### Select Action and Update sets
    if( min(loss.tot) %in% loss.A2 ){

      ### In (A2)
      # Find minimum
      A2.ind <- not.S[which(min(loss.tot) == loss.A2)]

      # update set
      E <- E.ask
      L <- L.ask
      S <- c(S.ask, A2.ind)

      # update derivative
      w.S <- w.A2[which(min(loss.tot) == loss.A2),]
      w <- rep(0,p)
      w[S] <- w.S

    } else{

      ### In (A1)
      # Find minimum
      A1.ind <- which( min(loss.tot) == loss.A1 )[1]

      # update set
      E <- E.ask[-A1.ind]

      # Determine L from E
      if( length(E.ask) == 1){
        Move.E <- x.ij[E.ask[A1.ind],S.ask] * w.A1
      } else{
        Move.E <- x.ij[E.ask[A1.ind],S.ask] %*% w.A1[A1.ind,]
      }
      if( Move.E<0 ){
        L <- c(L.ask, E.ask[A1.ind])
      } else{
        L <- L.ask
      }
      S <- S.ask


      # update derivative
      w.S <- if(length(E) == 0){
        w.A1
      } else{
        w.A1[A1.ind,]
      }
      w <- rep(0,p)
      w[S] <- w.S
    }
  }

  result <- structure(list(beta = temp.beta,
                           t = temp.t,
                           x=x, y=y), class="l1.rocsvmpath")
  return(result)
}

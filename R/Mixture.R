#' EM algorithm for Gaussian Mixture Models
#'
#' This function allows you to do Model Selection for Gaussian Mixture Models.
#' @param X The input data matrix.
#' @param M Maximum number of mixture components.
#' @param epsilon This parameter is used in penalty function, default 1e-6.
#' @param n.iter This parameter is the number of iteration for EM algorithm, default 200.
#' @param tol This parameter is the tolerance for likelihood function, default 1e-5.
#' @keywords Gaussian mixture
#' @export
#' @examples library(MASS)
#' @examples M = 10
#' @examples n <- 900
#' @examples S1 <- matrix(c(0.65,0.7794,0.7794,1.55),nrow=2,byrow=TRUE)
#' @examples mu1 <- c(-1,1)
#' @examples S2 <- matrix(c(0.65,-0.7794,-0.7794,1.55),nrow=2,byrow=TRUE)
#' @examples mu2 <- c(1,1)
#' @examples S3 <- matrix(c(2,0,0,0.2),nrow=2,byrow=TRUE)
#' @examples mu3 <- c(0,-sqrt(2))
#' @examples mu <- rbind(mu1,mu2,mu3)
#' @examples S <- list(S1,S2,S3)
#' @examples prop <- c(1/3,1/3,1/3)
#' @examples val1 <- mvrnorm(n*prop[1],mu=mu1,Sigma=S1)
#' @examples val2 <- mvrnorm(n*prop[2],mu=mu2,Sigma=S2)
#' @examples val3 <- mvrnorm(n*prop[3],mu=mu3,Sigma=S3)
#' @examples allval <- rbind(val1,val2,val3)
#' @examples X <- allval[sample(n,n),]
#' @examples lambda = sqrt(log(n)/n)/12
#' @examples res = gmm_EM(X, M, epsilon = 1e-6, n.iter=1000, tol=1e-5)
#' @examples k.hat = res$KA
#' @examples prop.hat = res$PI
#' @examples mu.hat = res$MEAN
#' @examples Sigma.hat = res$VAR


gmm_EM = function(X, M, epsilon, n.iter, tol){
# EM algorithm for Gaussian Mixture Model
# X -- data
# M -- maximum no. of mixtures
library(mvtnorm)
n = dim(X)[1]
p = dim(X)[2]
Df = p^2 + p*3/2 + 1

########### K-means method
res.kmeans = gmm_kmeans(X, M)
pi.hat = res.kmeans$PI
MU = res.kmeans$MEAN
Sigma = res.kmeans$VAR

######### EM algorithm
R = matrix(0, nrow = n, ncol = M)
MU.TEMP = list()
Sigma.TEMP = list()
for (i in 1:M){
  MU.TEMP[[i]] = MU
  Sigma.TEMP[[i]] = Sigma
}
loglik = rep(0, n.iter)
loglik.diff = tol + 1
i = 1
for (iter in 1:n.iter){
  # E-Step
  for (j in 1:M) {R[,j] = dmvnorm(X, MU[j, ], Sigma[[j]])}
  R = R * matrix(pi.hat, nrow = n, ncol = M, byrow = TRUE)
  temp = R * log(R)
  temp[!is.finite(temp)] = 0
  loglik[iter] = sum(temp / rowSums(R)) - penalty(X, pi.hat,epsilon)
  loglik.diff = loglik[iter] - loglik[iter-1]
  if((iter > 1) && (abs(loglik.diff) < tol)) break
  R = R / rowSums(R)
  # M-Step
  pi.hat = (colMeans(R) - lambda*Df) / (1 - M*lambda*Df)
  m.res = (pi.hat > 0)
  M.res = sum(m.res)
  if((M.res > 0) && (M.res <= M)){
    pi.hat = pi.hat[m.res]
    R = R[ , m.res]
    MU = MU[m.res, ]
    Sigma = Sigma[m.res]
   for (k in 1: M.res){
    X.temp = X - matrix(MU[k, ], nrow = n, ncol = p, byrow = TRUE)
    X.temp = X.temp * matrix(sqrt(R[,k]), nrow = n, ncol = p, byrow = FALSE)
    Sigma[[k]] = t(X.temp) %*% X.temp / sum(R[,k])
    MU[k, ] = t(R[,k]) %*% X / sum(R[,k])
    }
  }
  i = i + 1
  MU.TEMP[[i]] = MU
  Sigma.TEMP[[i]] = Sigma
  M = M.res
}
ind <- (pi.hat>0)
k = sum(ind)
prop = pi.hat[ind];
prop = prop/sum(prop)
return(list(KA = k, PI = prop, MEAN = MU, VAR = Sigma, loglik = loglik[1:(iter-1)]))
}            #MEAN.TEMP = MU.TEMP[1:i], VAR.TEMP = Sigma.TEMP[1:i]))


# BIC = function(X,lambda){
#   n = dim(X)[1]
#   p = dim(X)[2]
#   Df = p^2 + p*3/2 + 1
#   R = matrix(0, nrow = n, ncol = M)
#   for (j in 1:M) {
#     R[,j] = dmvnorm(X, MU[j, ], Sigma[[j]])
#     }
#   R = R * matrix(pi.hat, nrow = n, ncol = M, byrow = TRUE)
#
# }

#' K-means function for Gaussian Mixture Models
#'
#' This function allows you to use k-means clustering method to give the initial values for EM algorithm.
#' @param X The input data matrix.
#' @param K The initial number of mixing components.
#' @keywords K-means clustering

gmm_kmeans  = function (X, K){
  res.kmeans = kmeans(X, M)
  pi.hat = res.kmeans$size /n
  MU = res.kmeans$centers
  Sigma = list()
  for (i in 1:K){
    Sigma[[i]] = var(X[res.kmeans$cluster == i, ])
     if(det(Sigma[[i]]) <= 0){
       Sigma[[i]] = Sigma[[i]] + tol*diag(rep(1, p))
     }
  }
  return(list(PI = pi.hat, MEAN = MU, VAR = Sigma))
}



#' Penalty function for Gaussian Mixture Models
#'
#' This function allows you to calculate the penalty function for Gaussian mixture model.
#' @param X The input data matrix.
#' @param prop The mixing probability.
#' @param epsilon This parameter is used in penalty function, default 1e-6.
#' @keywords Penalty function

penalty = function(X, prop, epsilon){
  n = dim(X)[1]
  p = dim(X)[2]
  Df = p^2 + p*3/2 + 1
  epsilon = 1e-6
  lambda = sqrt(log(n)/n)/12
  penalty = n * lambda * Df * sum(log((epsilon + prop) / epsilon))
}















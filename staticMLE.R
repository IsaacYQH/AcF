staticMLE <- function(start = c(min(Q)-1,0.14,5), dim = 2, c = seq(0.01,0.1,0.002), alpha = seq(0.5,10,0.1)){
  # start: starting point of MLE iterations
  # dim: dimension of MLE parameters. 
  #  dim can be 2 by letting partial derivative with respect to sigma be 0
  source("staticMLE_funs.R")
  if(dim == 2){
    ######################1. Choose Appropriate Starting Point######################
    # choose best alpha by setting first partial derivative of alpha to zero
    mu = min(Q)-c[1]*log(length(Q))^(-1/mean(alpha))
    max = c(mu[1],alpha[1])
    #1.1 find the grid having the largest likelihood
    for (i in alpha){
      for (j in mu){
        if(logL_stat_2d(c(j,i)) > logL_stat_2d(c(max[1],max[2]))){
          max[1] = j
          max[2] = i
        }
      }
    }
    alpha = max[2]
    mu = max[1]
    #1.2 change initial value when hessian matrix is not negative-definite
    if(sum(eigen(Hess_stat_2d(c(mu, alpha)))$value>0) > 0){
      for (i in alpha) {
        mu = max(Q)-log(length(Q))^(-1/i)
        sigma = (length(Q)*sum((Q-mu)^-i)^-1)^(1/i)
        if(sum(eigen(Hess_stat_2d(c(mu, i)))$value<0) == 3){
          alpha = i
          break
        }
      }
    }
    #################################2. Static MLE##################################
    library(maxLik)
    A=matrix(c(-1,0,0,1),nrow=2);B=c(min(Q),0)
    r = maxLik(logL_stat_2d, start = c(mu=mu, alpha=alpha), 
               grad = grad_stat_2d, hess = Hess_stat_2d,
               constraints = list(ineqA=A, ineqB=B))
    sum.r = summary(r);sum.r
    
    #############################3. store results###################################
    mu_0 <<- sum.r[["estimate"]][1,1]
    alpha_1 <<- sum.r[["estimate"]][2,1]
    sigma_1 <<- (length(Q)*sum((Q-mu_0)^-alpha_1)^-1)^(1/alpha_1)
    H = Hess_stat_3d(c(mu_0, sigma_1, alpha_1))
    sd = c(diag(solve(-H))); sd = c(mu=sd[1],sigma=sd[2],alpha=sd[3])
    
    ##############################4. validation####################################
    message("Gradient of log-likelihhod:")
    print(grad_stat_3d(c(mu_0, sigma_1, alpha_1)))
    message("Number of positive eigenvalues of log-likelihhod:")
    print(sum(eigen(Hess_stat_3d(c(mu_0, sigma_1, alpha_1)))$value>0))
    
    return(data.frame("Estimate" = c(mu_0, sigma_1, alpha_1), "Std.error" = sd))
  }else if(dim == 3){
    A=matrix(c(-1,0,0,0,1,0,0,0,1), nrow=3, byrow = T);B=c(min(Q),0,0)
    r = maxLik(logL_stat_3d, start = start, grad = grad_stat_3d, 
               hess = Hess_stat_3d, constraints = list(ineqA=A, ineqB=B),
               iterlim=1000,reltol=1e-10)
    sum.r = summary(r)
    mu_0 <<- sum.r[["estimate"]][1,1]
    sigma_1 <<- sum.r[["estimate"]][2,1]
    alpha_1 <<- sum.r[["estimate"]][3,1]
    return(summary(r))
  }
}
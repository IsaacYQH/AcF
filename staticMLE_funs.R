logL_stat_3d = function(params){
  mu = params[1]
  sigma = params[2]
  alpha = params[3]
  n = length(Q)
  return(log(alpha)+alpha*log(sigma)-(alpha+1)*sum(log(Q-mu))/n
         -sum(((Q-mu)/sigma)^-alpha)/n)
}

grad_stat_3d <- function(params){
  mu = params[1]
  sigma = params[2]
  alpha = params[3]
  n = length(Q)
  return(c((alpha+1)*sum(1/(Q-mu))-alpha*sigma^alpha*sum((Q-mu)^(-alpha-1)),
           n*alpha/sigma-alpha*sigma^(alpha-1)*sum((Q-mu)^-alpha),
           n/alpha+n*log(sigma)-sum(log(Q-mu))+sigma^alpha*sum(log((Q-mu)/sigma)*(Q-mu)^-alpha))/n)
}

Hess_stat_3d <- function(params){
  mu = params[1]
  sigma = params[2]
  alpha = params[3]
  n = length(Q)
  result = matrix(0,3,3)
  result[1,1] = ((alpha+1)*sum((Q - mu)^-2)
                 - alpha*(alpha+1)*sigma^alpha*sum((Q-mu)^(-alpha-2)))/n
  result[2,2] = -alpha/sigma^2 - alpha*(alpha-1)/sigma^2*sum((Q-mu)^-alpha)*sigma^alpha/n
  result[3,3] = -1/alpha^2 - sum((Q-mu)^-alpha*(log((Q-mu)/sigma))^2)*sigma^alpha/n
  result[1,2] = result[2,1] = (-alpha^2/sigma^2)*sum((Q-mu)^(-alpha-1))*sigma^(alpha+1)/n
  result[1,3] = result[3,1] = (sum(1/(Q-mu)) - sigma^alpha*sum((Q-mu)^(-alpha-1))
                               + alpha*sigma^alpha*sum((Q-mu)^(-alpha-1)*log((Q-mu)/sigma)))/n
  result[2,3] = result[3,2] = (n/sigma - sigma^(alpha-1)*sum((Q-mu)^-alpha)
                               + alpha*sigma^(alpha-1)*sum((Q-mu)^-alpha*log((Q-mu)/sigma)))/n
  return(result)
}


logL_stat_2d <- function(params2){
  mu = params2[1]
  alpha = params2[2]
  n = length(Q)
  sigma = (n/sum((Q-mu)^-alpha))^(1/alpha)
  return(log(alpha)+alpha*log(sigma)-(alpha+1)*sum(log(Q-mu))/n-sum(((Q-mu)/sigma)^-alpha)/n)
}

grad_stat_2d <- function(params2){
  mu = params2[1]
  alpha = params2[2]
  n = length(Q)
  sigma = (n/sum((Q-mu)^-alpha))^(1/alpha)
  return(c(((alpha+1)*sum(1/(Q-mu))-alpha*sigma^alpha*sum((Q-mu)^(-alpha-1)))/n,
           (n/alpha+n*log(sigma)-sum(log(Q-mu))+sigma^alpha*sum(log((Q-mu)/sigma)*(Q-mu)^-alpha)))/n)
}

Hess_stat_2d <- function(params2){
  mu = params2[1]
  alpha = params2[2]
  sigma = (length(Q)*sum((Q-mu)^-alpha)^-1)^(1/alpha)
  return(Hess_stat_3d(c(mu, sigma, alpha))[-2,-2])
}
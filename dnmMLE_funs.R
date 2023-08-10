sigma_gen <- function(params, Q, sigma_1){
  # params = c("b0","b1","b2","b3")
  sigma_t <- c(sigma_1)
  T = length(Q)
  for (t in 2:T) {
    sigma_t[t] <- exp(params[1] + params[2] * log(sigma_t[t-1]) - 
                        params[3] * exp(-params[4] * Q[t-1]))
  }
  return(sigma_t)
}
alpha_gen <- function(params, Q, alpha_1){
  # params = c("g0","g1","g2","g3")
  alpha_t <- c(alpha_1)
  T = length(Q)
  for (t in 2:T) {
    alpha_t[t] <- exp(params[1] + params[2] * log(alpha_t[t-1]) + 
                        params[3] * exp(-params[4] * Q[t-1]))
  }
  return(alpha_t)
}

sigma_deriv_beta <- function(params, Q, index, sigma_t, sigma_1, t){
  # params=c("b0","b1","b2","b3")
  # index=0,1,2,3
  b0 = params["b0"]
  b1 = params["b1"]
  b2 = params["b2"]
  b3 = params["b3"]
  if(index==0){
    return(sigma_t*sum(sapply(b1, function(x) x = x^(0:(t-2)))))
  }else if(index==1){# special case that equals to 0 when t==2
    if(t==2){
      return(0)
    }else{
      return(
        sigma_t*(b0*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))))
                 - b2*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-b3*Q[(t-2):1])))
                 + (t-1)*b1^(t-2)*log(sigma_1)))
    }
  }else if(index==2){
    ####看到这！！！
    return(
      -sigma_t*(sum(sapply(b1, function(x) x = x^(0:(t-2))*exp(-b3*Q[(t-1):1])))))
  }else if(index==3){
    return(
      sigma_t*b2*sum(sapply(b1, function(x) 
        x = (x^(0:(t-2)))*Q[(t-1):1]*exp(-b3*Q[(t-1):1]))))
  }else{
    stop("Wrong Index!")
  }
}
alpha_deriv_gamma <- function(params, Q, index, alpha_t, alpha_1, t){
  # params=c("g0","g1","g2","g3")
  # index=0,1,2,3
  g0 = params["g0"]
  g1 = params["g1"]
  g2 = params["g2"]
  g3 = params["g3"]
  if(index==0){
    return(alpha_t*sum(sapply(g1, function(x) x = x^(0:(t-2)))))
  }else if(index==1){# special case that equals to 0 when t==2
    if(t==2){
      return(0)
    }else{
      return(
        alpha_t*(g0*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3))))
                 + g2*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-g3*Q[(t-2):1])))
                 + (t-1)*g1^(t-2)*log(alpha_1)))
    }
  }else if(index==2){
    return(
      alpha_t*(sum(sapply(g1, function(x) x = x^(0:(t-2))*exp(-g3*Q[(t-1):1])))))
  }else if(index==3){
    return(
      -alpha_t*g2*sum(sapply(g1, function(x) 
        x = (x^(0:(t-2)))*Q[(t-1):1]*exp(-g3*Q[(t-1):1]))))
  }else{
    stop("Wrong Index!")
  }
}
sigma_deriv2_beta <- function(params, Q, index, sigma_t, sigma_1, t){
  # params=c("b0","b1","b2","b3")
  # index=c(i,j)
  b0 = params["b0"]
  b1 = params["b1"]
  b2 = params["b2"]
  b3 = params["b3"]
  index = sort(index)
  if(index[1]==0){
    if(index[2]==1){
      return(
        sum(sapply(b1, function(x) x = x^(0:(t-2))))*
          sigma_deriv_beta(params, Q, 1, sigma_t, sigma_1, t)
        + sigma_t*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3)))))
    }else{
      return(
        sum(sapply(b1, function(x) x = x^(0:(t-2))))*
          sigma_deriv_beta(params, Q, index[2], sigma_t, sigma_1, t))
    }
  }else if(index[1]==1){
    if(index[2]==1){# special case that equals to 0 when t==3
      if(t==3){
        return(0)
      }else{
        return(
          sigma_deriv_beta(params, Q, 1, sigma_t, sigma_1, t)^2/sigma_t
          # (b0*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))))
          #   - b2*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-b3*Q[(t-2):1])))
          #   + (t-1)*b1^(t-2)*log(sigma_1))*sigma_deriv_beta(params, Q, 1, sigma_t, sigma_1, t)
          + sigma_t*(
            b0*sum(sapply(b1, function(x) x = (2:(t-2))*(1:(t-3))*x^(0:(t-4))))
            - b2*sum(sapply(b1, function(x) x = (2:(t-2))*(1:(t-3))*x^(0:(t-4))*exp(-b3*Q[(t-3):1])))
            + (t-1)*(t-2)*b1^(t-3)*log(sigma_1)))
      }
      
    }else if(index[2]==2){
      return(
        sigma_deriv_beta(params, Q, 1, sigma_t, sigma_1, t)/sigma_t
        *sigma_deriv_beta(params, Q, 2, sigma_t, sigma_1, t)
        # (b0*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3)))) - 
        #       b2*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-b3*Q[(t-2):1])))+
        #       (t-1)*b1^(t-2)*log(sigma_1))*sigma_deriv_beta(params, Q, 2, sigma_t, sigma_1, t)
           -sigma_t*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-b3*Q[(t-2):1]))))
      
    }else if(index[2]==3){
      return(
        sigma_deriv_beta(params, Q, 1, sigma_t, sigma_1, t)/sigma_t
        *sigma_deriv_beta(params, Q, 3, sigma_t, sigma_1, t)
        # (b0*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3)))) - 
        #       b2*sum(sapply(b1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-b3*Q[(t-2):1])))+
        #       (t-1)*b1^(t-2)*log(sigma_1))*sigma_deriv_beta(params, Q, 3, sigma_t, sigma_1, t)
           +sigma_t*b2*sum(sapply(b1, function(x)
             x = (1:(t-2))*Q[(t-2):1]*x^(0:(t-3))*exp(-b3*Q[(t-2):1]))))
    }else{
      stop("Wrong Index!")
    }
  }else if(index[1]==2){
    if(index[2]==2){
      return(
        sigma_deriv_beta(params, Q, 2, sigma_t, sigma_1, t)^2/sigma_t)
    }else if(index[2]==3){
      return(
        sigma_deriv_beta(params, Q, 3, sigma_t, sigma_1, t)*
          sigma_deriv_beta(params, Q, 2, sigma_t, sigma_1, t)/sigma_t
        # -sigma_deriv_beta(params, Q, 3, sigma_t, sigma_1, t)*
        #   sum(sapply(b1, function(x) x = x^(0:(t-2))*exp(-b3*Q[(t-1):1])))
        +sigma_t*sum(sapply(b1, function(x) x = Q[(t-1):1]*x^(0:(t-2))*exp(-b3*Q[(t-1):1]))))
    }else{
      stop("Wrong Index!")
    }
  }else if(index[1]==3){
    return(
      sigma_deriv_beta(params, Q, 3, sigma_t, sigma_1, t)^2/sigma_t
      # sigma_deriv_beta(params, Q, 3, sigma_t, sigma_1, t)*
      #   b2*sum(sapply(b1, function(x) x = Q[(t-1):1]*x^(0:(t-2))*exp(-b3*Q[(t-1):1])))
      -sigma_t*b2*sum(sapply(b1, function(x) 
        x = (Q[(t-1):1])^2*x^(0:(t-2))*exp(-b3*Q[(t-1):1]))))
  }else{
    stop("Wrong Index!")
  }
}
###改到这里！！！
alpha_deriv2_gamma <- function(params, Q, index, alpha_t, alpha_1, t){
  # params=c("g0","g1","g2","g3")
  # index=c(i,j)
  g0 = params["g0"]
  g1 = params["g1"]
  g2 = params["g2"]
  g3 = params["g3"]
  index = sort(index)
  if(index[1]==0){
    if(index[2]==1){
      return(
        sum(sapply(g1, function(x) x = x^(0:(t-2))))*
          alpha_deriv_gamma(params, Q, 1, alpha_t, alpha_1, t)+
          alpha_t*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3)))))
    }else{
      return(
        sum(sapply(g1, function(x) x = x^(0:(t-2))))*
          alpha_deriv_gamma(params, Q, index[2], alpha_t, alpha_1, t))
    }
    
  }else if(index[1]==1){
    if(index[2]==1){# special case that equals to 0 when t==3
      if(t==3){
        return(0)
      }else{
        return(
          alpha_deriv_gamma(params, Q, 1, alpha_t, alpha_1, t)^2/alpha_t
          # (g0*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3)))) + 
          #    g2*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-g3*Q[(t-2):1])))+
          #    (t-1)*g1^(t-2)*log(alpha_1))*alpha_deriv_gamma(params, Q, 1, alpha_t, alpha_1, t)
          +alpha_t*(
            g0*sum(sapply(g1, function(x) x = (2:(t-2))*(1:(t-3))*x^(0:(t-4))))
            +g2*sum(sapply(g1, function(x) x = (2:(t-2))*(1:(t-3))*x^(0:(t-4))*exp(-g3*Q[(t-3):1])))
            +(t-1)*(t-2)*g1^(t-3)*log(alpha_1)))
      }
      
    }else if(index[2]==2){
      return(
        alpha_deriv_gamma(params, Q, 1, alpha_t, alpha_1, t)/alpha_t
        *alpha_deriv_gamma(params, Q, 2, alpha_t, alpha_1, t)
        # (g0*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3)))) + 
        #    g2*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-g3*Q[(t-2):1])))+
        #    (t-1)*g1^(t-2)*log(alpha_1))*alpha_deriv_gamma(params, Q, 2, alpha_t, alpha_1, t)
        +alpha_t*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-g3*Q[(t-2):1]))))
      
    }else if(index[2]==3){
      return(
        alpha_deriv_gamma(params, Q, 1, alpha_t, alpha_1, t)/alpha_t
        *alpha_deriv_gamma(params, Q, 3, alpha_t, alpha_1, t)
        # (g0*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3)))) + 
        #    g2*sum(sapply(g1, function(x) x = (1:(t-2))*x^(0:(t-3))*exp(-g3*Q[(t-2):1])))+
        #    (t-1)*g1^(t-2)*log(alpha_1))*alpha_deriv_gamma(params, Q, 3, alpha_t, alpha_1, t)
        -alpha_t*g2*sum(sapply(g1, function(x)
          x = (1:(t-2))*Q[(t-2):1]*x^(0:(t-3))*exp(-g3*Q[(t-2):1]))))
    }else{
      stop("Wrong Index!")
    }
  }else if(index[1]==2){
    if(index[2]==2){
      return(
        alpha_deriv_gamma(params, Q, 2, alpha_t, alpha_1, t)^2/alpha_t)
    }else if(index[2]==3){
      return(
        alpha_deriv_gamma(params, Q, 3, alpha_t, alpha_1, t)*
          alpha_deriv_gamma(params, Q, 2, alpha_t, alpha_1, t)/alpha_t
        # alpha_deriv_gamma(params, Q, 3, alpha_t, alpha_1, t)*
        #   sum(sapply(g1, function(x) x = x^(0:(t-2))*exp(-g3*Q[(t-1):1])))
        -alpha_t*sum(sapply(g1, function(x) x = Q[(t-1):1]*x^(0:(t-2))*exp(-g3*Q[(t-1):1]))))
    }else{
      stop("Wrong Index!")
    }
  }else if(index[1]==3){
    return(
      (-alpha_deriv_gamma(params, Q, 3, alpha_t, alpha_1, t))^2/alpha_t
      # -alpha_deriv_gamma(params, Q, 3, alpha_t, alpha_1, t)*
      #   g2*sum(sapply(g1, function(x) x = Q[(t-1):1]*x^(0:(t-2))*exp(-g3*Q[(t-1):1])))
      +alpha_t*g2*sum(sapply(g1, function(x) 
        x = (Q[(t-1):1])^2*x^(0:(t-2))*exp(-g3*Q[(t-1):1]))))
  }else{
    stop("Wrong Index!")
  }
}

firstPD_beta <- function(params, Q, sigma_t, alpha_t){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  L = c(rep(0,n))
  Ln = c(rep(0,n))
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 2:T) {
    L_sigma = alpha_t[t]/sigma_t[t] - 
      alpha_t[t]/sigma_t[t]*((Q[t]-mu)/sigma_t[t])^-alpha_t[t]
    for (i in 0:(n-1)) {
      Ln[i+1] = L_sigma * sigma_deriv_beta(params[1:n], Q, i, sigma_t[t], sigma_t[1], t)
    }
    L = Ln + L
  }
  return(L)
}
firstPD_gamma <- function(params, Q, sigma_t, alpha_t){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  L = c(rep(0,n))
  Ln = c(rep(0,n))
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 2:T) {
    L_alpha = 1/alpha_t[t] - 
      log((Q[t]-mu)/sigma_t[t]) +
      ((Q[t]-mu)/sigma_t[t])^-alpha_t[t]*log((Q[t]-mu)/sigma_t[t])
    for (i in 0:(n-1)) {
      Ln[i+1] = L_alpha * alpha_deriv_gamma(params[(n+1):(2*n)], Q, i, alpha_t[t],
                                            alpha_t[1], t)
    }
    L = Ln + L
  }
  return(L)
}
firstPD_m_dnm <- function(mu, Q, sigma_t, alpha_t){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  L = 0
  T = length(Q)
  for (t in 1:T) {
    Ln = (alpha_t[t]+1)/(Q[t] - mu) - 
      (alpha_t[t]/sigma_t[t])*((Q[t] - mu)/sigma_t[t])^(-alpha_t[t]-1)
    L = Ln + L
  }
  return(L)
}
secondPD_mm_dnm <- function(mu, Q, sigma_t, alpha_t){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  L = 0
  T = length(Q)
  for (t in 1:T) {
    L = (alpha_t[t]+1)/(Q[t] - mu)^2 - 
      alpha_t[t]*(alpha_t[t]+1)*sigma_t[t]^alpha_t[t]*(Q[t]-mu)^(-alpha_t[t]-2) + L
  }
  return(L)
}
secondPD_mb_dnm <- function(params, Q, sigma_t, alpha_t, index){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  L = 0
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 2:T) {
    L_mu = -alpha_t[t]^2/sigma_t[t]^2*((Q[t]-mu)/sigma_t[t])^(-alpha_t[t]-1)
    Ln = L_mu * sigma_deriv_beta(params[1:n], Q, index, sigma_t[t], sigma_t[1], t)
    L = Ln + L
  }
  return(L)
}
secondPD_mg_dnm <- function(params, Q, sigma_t, alpha_t, index){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  L = 0
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 2:T) {
    L_mu = 1/(Q[t]-mu) - sigma_t[t]^alpha_t[t]*(Q[t]-mu)^(-alpha_t[t]-1) +
      alpha_t[t]*sigma_t[t]^alpha_t[t]*(Q[t]-mu)^(-alpha_t[t]-1)*log((Q[t]-mu)/sigma_t[t])
    Ln = L_mu * alpha_deriv_gamma(params[(n+1):(2*n)], Q, index, alpha_t[t], alpha_t[1], t)
    L = Ln + L
  }
  return(L)
}
secondPD_bb_dnm <- function(params, Q, sigma_t, alpha_t, index){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  L = 0
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 3:T) {
    L = (alpha_t[t]/sigma_t[t] - alpha_t[t]/sigma_t[t]*((Q[t]-mu)/sigma_t[t])^-alpha_t[t])*
      sigma_deriv2_beta(params[1:n], Q, index, sigma_t[t], sigma_t[1], t) +
      (-alpha_t[t]/sigma_t[t]^2-alpha_t[t]*(alpha_t[t]-1)/sigma_t[t]^2*
      ((Q[t]-mu)/sigma_t[t])^-alpha_t[t])*
      sigma_deriv_beta(params[1:n], Q, index[1], sigma_t[t], sigma_t[1], t)*
      sigma_deriv_beta(params[1:n], Q, index[2], sigma_t[t], sigma_t[1], t) + L
  }
  return(L)
}
secondPD_gg_dnm <- function(params, Q, sigma_t, alpha_t, index){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  L = 0
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 3:T) {
    L = (1/alpha_t[t]-log((Q[t]-mu)/sigma_t[t])+((Q[t]-mu)/sigma_t[t])^-alpha_t[t]*
      log((Q[t]-mu)/sigma_t[t]))*alpha_deriv2_gamma(params[(n+1):(2*n)], Q, index, alpha_t[t], alpha_t[1], t)
      + (-1/alpha_t[t]^2-((Q[t]-mu)/sigma_t[t])^-alpha_t[t]*(log((Q[t]-mu)/sigma_t[t]))^2)*
      # + (-1/alpha_t[t]^2+((Q[t]-mu)/sigma_t[t])^-alpha_t[t]*(log((Q[t]-mu)/sigma_t[t]))^2)*
      alpha_deriv_gamma(params[(n+1):(2*n)], Q, index[1], alpha_t[t], alpha_t[1], t)*
      alpha_deriv_gamma(params[(n+1):(2*n)], Q, index[2], alpha_t[t], alpha_t[1], t) + L
  }
  return(L)
}
secondPD_bg_dnm <- function(params, Q, sigma_t, alpha_t, index){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  # index = c(i,j), i for beta, j for gamma
  mu = params["mu"]
  n = (length(params)-1)/2
  L = 0
  T = length(Q)
  # Note that L1=0, start from t=2
  for (t in 2:T) {
    L = (1/sigma_t[t] - 1/sigma_t[t]*((Q[t]-mu)/sigma_t[t])^-alpha_t[t] + 
           alpha_t[t]/sigma_t[t]*((Q[t]-mu)/sigma_t[t])^-alpha_t[t]*log((Q[t]-mu)/sigma_t[t]))*
      sigma_deriv_beta(params[1:n], Q, index[1], sigma_t[t], sigma_t[1], t)*
      alpha_deriv_gamma(params[(n+1):(2*n)], Q, index[2], alpha_t[t], alpha_t[1], t) + L
  }
  return(L)
}

logL <- function(params){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  sigma_t = sigma_gen(params[1:n], Q, sigma_1)
  alpha_t = alpha_gen(params[(n+1):(2*n)], Q, alpha_1)
  L = sum(log(alpha_t) + alpha_t*log(sigma_t)
      - (alpha_t + 1)*log(Q - mu)
      - sigma_t^alpha_t * (Q - mu)^-alpha_t)
  return(L/length(Q))
}

grad_logL <- function(params){
  # params = c("b0","b1","b2","b3","g0","g1","g2","g3","mu")
  n = (length(params)-1)/2
  sigma_t = sigma_gen(params[1:n], Q, sigma_1)
  alpha_t = alpha_gen(params[(n+1):(2*n)], Q, alpha_1)
  return(c(firstPD_beta(params, Q, sigma_t, alpha_t),
           firstPD_gamma(params, Q, sigma_t, alpha_t),
           firstPD_m_dnm(params["mu"], Q, sigma_t, alpha_t)))
}

Hess_logL <- function(params){
  n = length(params)
  nb = ng = (length(params)-1)/2
  Hess = matrix(0,n,n)
  sigma_t = sigma_gen(params[1:nb], Q, sigma_1)
  alpha_t = alpha_gen(params[(nb+1):(2*nb)], Q, alpha_1)
  # Partial Derivative betabeta
  for (i in 1:nb) {
    for (j in 1:i) {
      # note that the indexes are from 0 to nb-1,
      # but the indexes of the matrix are from 1 to nb
      Hess[i,j] = Hess[j,i] = secondPD_bb_dnm(params, Q, sigma_t, alpha_t, c(i-1,j-1))
    }
  }
  # Partial Derivative gammagamma
  for (i in 1:ng) {
    for (j in 1:i) {
      Hess[i+nb,j+nb] = Hess[j+nb,i+nb] = secondPD_gg_dnm(params, Q, sigma_t, alpha_t, c(i-1,j-1))
    }
  }
  # Partial Derivative mumu
  Hess[n,n] = secondPD_mm_dnm(params["mu"], Q, sigma_t, alpha_t)
  # Partial Derivative betagamma
  for (i in 1:nb) {
    for (j in 1:ng) {
      Hess[i,j+nb] = Hess[j+nb,i] = secondPD_bg_dnm(params, Q, sigma_t, alpha_t, c(i-1,j-1))
    }
  }
  # Partial Derivative mubeta
  for (i in 1:nb) {
    Hess[i,n] = Hess[n,i] = secondPD_mb_dnm(params, Q, sigma_t, alpha_t, i-1)
  }
  # Partial Derivative mugamma
  for (i in 1:ng) {
    Hess[i+nb,n] = Hess[n,i+nb] = secondPD_mg_dnm(params, Q, sigma_t, alpha_t, i-1)
  }
  return(Hess)
}

Hess_bb <- function(params){
  n = length(params)
  nb = ng = (length(params)-1)/2
  Hess = matrix(0,4,4)
  sigma_t = sigma_gen(params[1:nb], Q, sigma_1)
  alpha_t = alpha_gen(params[(nb+1):(2*nb)], Q, alpha_1)
  # Partial Derivative betabeta
  for (i in 1:nb) {
    for (j in 1:i) {
      # note that the indexes are from 0 to nb-1,
      # but the indexes of the matrix are from 1 to nb
      Hess[i,j] = Hess[j,i] = secondPD_bb_dnm(params, Q, sigma_t, alpha_t, c(i-1,j-1))
    }
  }
  # return(Hess)
  return(Hess+diag(-1000,nrow=9))
}
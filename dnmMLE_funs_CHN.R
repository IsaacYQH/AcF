sigma_gen_CHN <- function(params, Q, limit10, limit20, sigma_1){
  # params = c("b0","b1","b2","b3","b4","b5")
  sigma_t <- c(sigma_1)
  T = length(Q)
  for (t in 2:T) {
    sigma_t[t] <- exp(params[1] + params[2] * log(sigma_t[t-1])
                      - params[3] * exp(-params[4] * Q[t-1])
                      + params[5] * limit10[t-1] + params[6] * limit20[t-1])
  }
  return(sigma_t)
}
alpha_gen_CHN <- function(params, Q, limit10, limit20, alpha_1){
  # params = c("g0","g1","g2","g3","g4","g5")
  alpha_t <- c(alpha_1)
  T = length(Q)
  for (t in 2:T) {
    alpha_t[t] <- exp(params[1] + params[2] * log(alpha_t[t-1])
                      + params[3] * exp(-params[4] * Q[t-1])
                      - params[5] * limit10[t-1] - params[6] * limit20[t-1])
  }
  return(alpha_t)
}

logL_CHN <- function(params){
  # params=c("b0","b1","b2","b3","b4","b5","g0","g1","g2","g3","g4","g5","mu")
  mu = params["mu"]
  n = (length(params)-1)/2
  sigma_t = sigma_gen_CHN(params[1:n], Q, limit10, limit20, sigma_1)
  alpha_t = alpha_gen_CHN(params[(n+1):(2*n)], Q, limit10, limit20, alpha_1)
  L = sum(log(alpha_t) + alpha_t*log(sigma_t)
          - (alpha_t + 1)*log(Q - mu)
          - sigma_t^alpha_t * (Q - mu)^-alpha_t)
  return(L/length(Q))
}
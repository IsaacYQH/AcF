source("dnmMLE_funs.R")
# public tool grids
grid = seq(1,10,4)

# constraints: A*theta + B > 0
A = matrix(c(rep(0,1),1,rep(0,7),
             rep(0,2),1,rep(0,6),
             rep(0,3),1,rep(0,5),
             rep(0,5),1,rep(0,3),
             rep(0,6),1,rep(0,2),
             rep(0,7),1,rep(0,1),
             0,-1,rep(0,7),
             rep(0,5),-1,rep(0,3),
             rep(0,8),-1), nrow = 9, byrow = T)
B = c(rep(0,6),rep(1,2),min(Q))

library(foreach)
library(doParallel)

# Register a backend for parallel computing, where workers indicates how many cores to use for computing.
# Usually set to the number of available CPU cores minus 1, but adjusted on a case-by-case basis.
registerDoParallel(workers = detectCores() - 3)

max = -Inf
r_max = NULL
start_max = NULL
b0b2 = apply(X = expand.grid(grid/40-0.125, grid^(5/2)/100), FUN = list, MARGIN = 1)
start.time = Sys.time()
# start initial searching for best starting point
results <- foreach(first2para=b0b2, .combine='rbind', .multicombine=TRUE, .packages = "maxLik") %dopar% {
  local_results <- list()
  for (b3 in grid) {
    for (g0 in grid/40-0.125) {
      for (g2 in grid^(5/2)/100) {
        for (g3 in grid) {
          
          start=c(unname(unlist(first2para))[1],0.9,unname(unlist(first2para))[2],b3,g0,0.9,g2,g3,mu=mu_0)
          r0 = try(maxLik(logLik = logL,
                          start=start,
                          iterlim=100,
                          reltol=1e-5,
                          constraints=list(ineqA=A,ineqB=B),
                          fixed=c("mu")))
          if(is(r0,"try-error")){
            next
          }else{
            # save starting points
            r0[["start"]]=start
            local_results <- append(local_results, list(r0))
          }
          
        }}}
  }
  local_results
}
end.time = Sys.time()
print(end.time-start.time)
# Find the maximum value and update r_max and max
for (result in results) {
  if (!is.null(result$maximum)){
    if (result$maximum > max) {
      r_max = result
      max = result$maximum
    }
  }
};summary(r_max)

r0 = r_max
start = r_max$start
# iterate until every standard error is not Nan
for (i in 1:100) {
  r0=maxLik(logLik = logL,
            start=start,
            iterlim=10000,
            reltol=1e-8,
            constraints=list(ineqA=A,ineqB=B),
            # do not fix mu
  )
  sum.r0 = summary(r0)
  start = r0$estimate
  if (!is.na(sum(sum.r0$estimate[,"Std. error"]))) {break}
}
# ignore the warnings produces by failing attempts
if (i == 100) {message("Standard errors are unavailable")}

message("Check positive-definiteness (show eigenvalues):")
eigen(r0$hessian)$values

# plots the evolution of alpha and sigma
plot(sigma_gen(r0$estimate[1:4],Q,sigma_1),
     type = "l", main = "sigma for Original AcF")
abline(sigma_1)
plot(alpha_gen(r0$estimate[5:8],Q,alpha_1),
     type = "l", main = "alpha for Original AcF")
abline(alpha_1)
source("dnmMLE_funs_CHN.R")
# start=c(b0=-0.050,b1=0.96,b2=0.051,b3=7,g0=0.1,g1=0.89,g2=0.33,g3=8,mu=mu_0)
grid = seq(1,10,4)

# start=c(b0=-0.050,b1=0.96,b2=0.051,b3=7,g0=-0.068,g1=0.89,g2=0.33,g3=5.33,mu=-0.069)
# constraints: A*theta + B > 0
A = matrix(c(rep(0,1),1,rep(0,11),
             rep(0,2),1,rep(0,10),
             rep(0,3),1,rep(0,9),
             rep(0,4),1,rep(0,8),
             rep(0,5),1,rep(0,7),
             rep(0,7),1,rep(0,5),
             rep(0,8),1,rep(0,4),
             rep(0,9),1,rep(0,3),
             rep(0,10),1,rep(0,2),
             rep(0,11),1,rep(0,1),
             0,-1,rep(0,11),
             rep(0,7),-1,rep(0,5),
             rep(0,12),-1), nrow = 13, byrow = T)
B = c(rep(0,10),rep(1,2),min(Q))



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
results <- foreach(first2para=b0b2, .combine='rbind', .multicombine=TRUE, .packages = "maxLik") %dopar% {
  local_results <- list()
  for (b3 in grid) {
    for (g0 in grid/40-0.125) {
      for (g2 in grid^(5/2)/100) {
        for (g3 in grid) {
          
          start=c(unname(unlist(first2para))[1],0.9,unname(unlist(first2para))[2],b3,0.1,0.1,g0,0.9,g2,g3,0.1,0.1,mu=mu_0)
          r0 = try(maxLik(logLik = logL_CHN,
                          start=start,
                          iterlim=100,
                          reltol=1e-5,
                          constraints=list(ineqA=A,ineqB=B),
                          fixed=c("mu")))
          if(is(r0,"try-error")){
            next
          }else{
            r0[["start"]]=start
            local_results <- append(local_results, list(r0))
          }
          
        }}}
  }
  local_results
}
end.time = Sys.time()
print(end.time-start.time)
# 查找最大值并更新r_max和max
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
for (i in 1:100) {
  r0=maxLik(logLik = logL_CHN,
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

plot(sigma_gen_CHN(r0$estimate[1:6],Q,limit10,limit20,sigma_1),type = "l")
plot(alpha_gen_CHN(r0$estimate[7:12],Q,limit10,limit20,alpha_1),type = "l")
a = alpha_gen_CHN(r0$estimate[7:12],Q,limit10,limit20,alpha_1)

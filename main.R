set.seed(1)
# mode: "static", "dynamic_normal", "dynamic_China"
mode = "dynamic_China"
stat_dim = 2
# Data Preparation
source("Data_Pre.R")
data <- Data_Pre(simu = F, mode = mode)
Q <- data$Q
source("staticMLE.R")
# mu_0, sigma_1, alpha_1 can be obtained be staticMLE()

if(mode == "static"){
  rm(data)
  initPara <- staticMLE(dim=stat_dim)
}else if (mode == "dynamic_normal"){
  rm(data)
  source("staticMLE.R")
  initPara <- staticMLE(dim=stat_dim)
  source("dnmMLE.R")
}else if (mode == "dynamic_China"){
  initPara <- staticMLE(dim=stat_dim, alpha = seq(1,50,1), c = seq(1,150,1))
  limit10 <- data$limit10
  limit20 <- data$limit20
  rm(data)
  source("dnmMLE_CHN.R")
}

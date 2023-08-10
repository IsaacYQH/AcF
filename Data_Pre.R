# Data Preparation
# similation
Data_Pre <- function(simu = T, mode = "static"){
  if(simu){
    source("gen_simu.R")
    simu <- gen_simu(MAX_T=200,seed=100)
    Q <<- simu$Q
    return(list(Q=Q))
  }else{ # Real Data
    # Read Data
    library(readxl)
    if(mode == "dynamic_China"){
      ret <- as.matrix(read_excel("Dataset_all_PRC_firms.xlsx", sheet = "Return"))
      mt <- as.matrix(read_excel("Dataset_all_PRC_firms.xlsx", sheet = "Market Type"))
      ts <- as.matrix(read_excel("Dataset_all_PRC_firms.xlsx", sheet = "Trading Status"))
      date <- read_excel("Dataset_all_PRC_firms.xlsx", sheet = "Dates", col_names = F)
      
      # check data
      if (nrow(ret)!=nrow(mt) || nrow(ret)!=nrow(ts) || nrow(ts)!=nrow(mt))
        stop("Uncomfortable Dimension!")
      # Only Considerate Shenzhen Stock Exchange
      SZret <- ret[,1:2870]
      SZmt <- mt[,1:2870]
      SZts <- ts[,1:2870]
      rm(ret,mt,ts)
      # Prepare covariates
      T <- nrow(SZret)
      limit10 <- c()
      limit20 <- c()
      Q <- c()
      for (t in 1:T) {
        # make sure:
        #  1. Return value is significantly far from 0.1 or 0.2
        #  2. Grouping by market type
        #  3. Trading Status are normal (avoiding adding more limits levels)
        temp10 = ((SZret[t,]>0 & (0.1-SZret[t,])>1e-2) | (SZret[t,]<=0 
            & (0.1+SZret[t,])>1e-2)) & (SZmt[t,]==4 | SZmt[t,]==8) & SZts[t,]==1
        temp20 = ((SZret[t,]>0 & (0.2-SZret[t,])>1e-2) | (SZret[t,]<=0 
            & (0.2+SZret[t,])>1e-2)) & SZmt[t,]==16 & SZts[t,]==1
        limit10[t] = sum(subset(temp10, !is.na(temp10)))
        limit20[t] = sum(subset(temp20, !is.na(temp20)))
        Q[t] <- max(
          max(subset(SZret[t,temp10],!is.na(SZret[t,temp10]))),
          max(subset(SZret[t,temp20],!is.na(SZret[t,temp20])))
        )
      }
      # to make the inferrence on parameters fair
      limit10 = min(Q) + (limit10 - min(limit10))*(max(Q)-min(Q))/(max(limit10)-min(limit10))
      limit20 = min(Q) + (limit20 - min(limit20))*(max(Q)-min(Q))/(max(limit20)-min(limit20))
      return(list(Q=Q, limit10=limit10, limit20=limit20))
    }else if(mode == "static" || mode == "dynamic_normal"){
      ret <- as.matrix(read_excel("Dataset_all_PRC_firms.xlsx", sheet = "Return"))
      # Only Considerate Shenzhen Stock Exchange
      SZret <- ret[,1:2870]
      rm(ret)
      # Prepare covariates
      T <- nrow(SZret)
      Q <- c()
      for (t in 1:T) {Q[t] <- max(subset(SZret[t,],!is.na(SZret[t,])))}
      return(list(Q=Q))
    }else{stop("Wrong mode!")}
  }
}

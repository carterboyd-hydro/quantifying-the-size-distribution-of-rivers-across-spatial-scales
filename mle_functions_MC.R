#Author: Carter Boyd
#Date: 04/26/2024
#Purpose: Functions to run maximum likelihood analysis to fit log-normal and Pareto
#distributions to river width data, including Monte Carlo uncertainty analysis

###Load libraries
library(emdbook)


###Functions

#Pareto functions
#Function to find mode
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Fit Pareto distribution to data with MLE:
pareto.mle <- function(x, xm=get_mode(x)){
  a = length(x)/(sum(log(x))-length(x)*log(xm))
  a_se = sqrt((a^2)/length(x))
  return( list(xm = xm, a = a, a_se = a_se))
}  

#Compute the KS statistic, and use parametric bootstrap to estimate the p-value.
pareto.test <- function(x, xm=get_mode(x), B=2.5e2){
  require(stats)
  require(EnvStats)
  
  a = pareto.mle(x, xm)
  
  # KS statistic
  D = ks.test(x, function(q) ppareto(q, a$xm, a$a))$statistic
  
  # Estimating p value with parametric bootstrap
  n = length(x)
  emp.D = numeric(B)
  for(b in 1:B){
    xx = rpareto(n, a$xm, a$a);
    aa = pareto.mle(xx)
    emp.D[b] = ks.test(xx, function(q) ppareto(q, aa$xm, aa$a))$statistic
  }
  
  return(list(xm = a$xm, a = a$a, D = D, p = sum(emp.D > D)/B))
}

#Create a Pareto probability density function
pareto_pdf <- function(x, lambda = 1, k = 1){
  density <- (k*(lambda^k)) / (x^(k + 1))
  return(density)
}

#Create a Pareto probability density function for plotting in log-space
pareto_pdf_log <- function(x, lambda = 1, k = 1, bw, n){
  density <- (k*(lambda^k)) / ((x)^(k + 1)) * bw * n
  return(density)
}

#Create a Pareto probability density function for plotting in log-space, scaled to river length instead of count
pareto_pdf_log_scaled <- function(x, lambda = 1, k = 1, bw, n, scale_factor){
  density <- (k*(lambda^k)) / ((x)^(k + 1)) * bw * n * scale_factor
  return(density)
}

#Function to test fit of lognormal and Pareto distributions
test_fits = function(data, pareto = F){
  require(stats)
  require(EnvStats)
  require(MASS)
  require(transport)
  
  # Create a vector of width data
  w <- data$Width
  
  # Find mode of data for Pareto fit
  minW <- get_mode(w)
  # print(paste0("Mode width is: ", minW, " m"))
  
  
  ## Custom graphing parameters: 
  # Histogram binning interval (m):
  # int = 2*IQR(w) / length(w)^(1/3)
  # print(paste0('Calculated binwidth: ', int))
  
  # Set up table to store all the fitted distribution and GOF statistics info:
  if(pareto == T){
    names = c("Stream_Order", "N",
            "meanln", "sdlog", "lnWass", #"lnKS_D", "lnKS_p",
            "meanln_sd", "sdlog_sd",
            "pXmin", "pAlpha", "pWass", #"pKS_D", "pKS_p",
            "pAlpha_se"
    )
  }else{
    names = c("Stream_Order", "N",
              "meanln", "sdlog", "lnWass", #"lnKS_D", "lnKS_p",
              "meanln_sd", "sdlog_sd")
  }
  statTab = data.frame(array(NA, c(1, length(names))))
  names(statTab) = names
  
  ## Distribution fitting with maximum likelihood estimation:
  # print('Fitting lognormal distribution...')
  lntest = fitdistr(w, "log-normal")
  lnfit = lntest$estimate
  lnsd = lntest$sd
  
  if(pareto == T){
    # print('Fitting pareto distribution...')
    parfit = pareto.mle(w[w>minW], minW)
  }
  
  #### Quantify goodness of fit:
  # print('Starting goodness of fit tests...')
  
  ### Wasserstein 1D distance metric for goodness of fit
  ## Lognormal
  # print('Performing lognormal Wasserstein GOF test...')
  # Create a synthetic distribution based on the best fit lognormal function
  ln_syntheticDist <- rlnorm(n = length(w), meanlog = lnfit[1], sdlog = lnfit[2])
  # Calculate Wasserstein distance metric
  ln_wass <- wasserstein1d(w, ln_syntheticDist, p = 1)
  
  if(pareto == T){
    ## Pareto
    # print('Performing Pareto Wasserstein GOF test...')
    # Create a synthetic distribution based on the best fit Pareto function
    par_syntheticDist <- EnvStats::rpareto(n = length(w[w>minW]), location = parfit[[1]], shape = parfit[[2]])
    # Calculate Wasserstein distance metric
    par_wass <- wasserstein1d(w[w>minW], par_syntheticDist, p = 1)
  }
  
  # #### Two sided One sample KS GOF test:
  # jw = jitter(w) # to remove ties
  # 
  # ### Lognormal
  # print('Performing lognormal KS GOF test...')
  # lnks = ks.test(jw, "plnorm", lnfit[1], lnfit[2], alternative="two.sided", simulate.p.value = T, B = 250)
  # 
  # if(pareto == T){
  #   ### Pareto
  #   print('Performing Pareto KS GOF test...')
  #   parks = pareto.test(jw[jw>minW], minW) #modeW
  # }

  
  ## Create a histogram to write out number of bins using specified binwidth
  # breaks = seq(0, (max(w)+int), int)
  # h = hist(w, breaks, plot=F)
  
  ## Add distribution fit and GOF statistics to a table:
  if(pareto == T){
    statTab[1, ] = c(NA, length(w),
                    lnfit[[1]], lnfit[[2]], ln_wass,# lnks$statistic[[1]], lnks$p.value[[1]][[1]],
                    lnsd[[1]], lnsd[[2]],
                    parfit[[1]], parfit[[2]], par_wass, #parks$D[[1]], parks$p[[1]],
                    parfit[[3]]
    )
  }else{
    statTab[1, ] = c(NA, length(w),
                     lnfit[[1]], lnfit[[2]], ln_wass, #lnks$statistic[[1]], lnks$p.value[[1]][[1]],
                     lnsd[[1]], lnsd[[2]])
  }
  
  ## Create an output table
  oTab = as.data.frame(t(statTab))
  return(oTab)
}

#Author: Carter Boyd
#Date: 04/26/2024
#Purpose: Fit statistical distributions to river width observations by stream order
#and generate synthetic widths covering the Mississippi River Basin, including
#Monte Carlo uncertainty analysis

#Load necessary libraries
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(EnvStats)
library(emdbook)
library(MASS)
library(grid)
library(scales)

###Load Data

#Load river width data
all_widths <- read.csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_04_26_MC_Analysis/Data/all_widths.csv")
all_widths_sampled <- read.csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_04_26_MC_Analysis/Data/all_widths_sampled.csv")


###Functions

#Load functions for maximum likelihood estimation
source("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_04_26_MC_Analysis/mle_functions_MC.R")

#Create a function for testing the best fits by stream order
mle_by_order <- function(width_data, csvOut = NA, orders, parOption = F){
  #Create an empty dataframe for export
  if(parOption == T){
    ensembleDF <- data.frame(array(NA,c(11 #15
                                        , 0)))
  }else{
    ensembleDF <- data.frame(array(NA, c(7 #9
                                         , 0)))
  }
  
  #For loop to iterate over stream orders
  for (i in orders){
    # print(i)
    filtered_widths <- filter(width_data, Order == i)
    tempDF <- test_fits(filtered_widths, parOption)
    tempDF[1,1] <- orders[i]
    ensembleDF <- cbind(ensembleDF, tempDF)
  }
  colnames(ensembleDF) <- c(orders)
  
  return(ensembleDF)
}

#Create a function for testing the best fits by basin
mle_by_basin <- function(width_data, csvOut = NA, basins, pareto = F){
  #Create an empty dataframe for export
  if(pareto == T){  
    ensembleDF <- data.frame(array(NA,c(11 #15
                                        , 0)))
  }else{
    ensembleDF <- data.frame(array(NA, c(7 #9
                                         , 0)))
  }
  
  #For loop to iterate over stream orders
  for (i in basins){
    print(i)
    filtered_widths <- filter(width_data, Basin == i)
    tempDF <- test_fits(filtered_widths)
    ensembleDF <- cbind(ensembleDF, tempDF)
  }
  
  colnames(ensembleDF) <- basins
  
  return(ensembleDF)
}

#Function to find average widths by stream order
avg_width_by_order <- function(width_data, csvOut = NA, orders){
  #Create an empty dataframe for export
  ensembleDF <- data.frame(array(NA,c(0, 2)))
  
  #For loop to iterate over stream orders
  for (i in orders){
    print(i)
    filtered_widths <- filter(width_data, Order == i)
    avgWidth <- mean(filtered_widths$Width)
    medianWidth <- median(filtered_widths$Width)
    ensembleDF <- rbind(ensembleDF, c(orders[i], avgWidth, medianWidth))
  }
  colnames(ensembleDF) <- c("Order", "Avg_Width_m", "Median_Width_m")
  return(ensembleDF)
}

#Define functions to find number and average length by order
find_number <- function(a, b, order){
  number <- a*(b^order)
  return(number)
}

find_length <- function(c, d, order){
  length <- c*(d^order)
  return(length)
}

#Function to find number and average length of streams of each order
find_length_by_order <- function(orders, a, b, c, d){
  #Create an empty dataframe for export
  orderDF <- data.frame(matrix(ncol = 3, nrow = 0))
  
  #For loop to iterate over stream orders
  for (i in orders){
    print(i)
    number <- find_number(a, b, i)
    length <- find_length(c, d, i)
    tempDF <- c(i, number, length)
    orderDF <- rbind(orderDF, tempDF)
  }
  colnames(orderDF) <- c("Order", "Number", "Avg_Length_m")
  orderDF <- mutate(orderDF, Total_Length_m = Number*Avg_Length_m)
  return(orderDF)
}

#Function to estimate RSSA using complex method
complex_rssa_estimate <- function(statDF, lengthDF, orders){
  areaDF <- data.frame(matrix(ncol = 2 #4
                              , nrow = 0))
  for(i in orders){
    # print(i)
    distribution_params <- statDF[,i]
    length <- lengthDF[i,2] * lengthDF[i,3]
    count <- round(length/30)
    synthetic_widths <- rlnorm(n = count,
                               meanlog = distribution_params[3],
                               sdlog = distribution_params[4]) 
    # synthetic_widths_high <- rlnorm(n = count,
    #                                 meanlog = distribution_params[5] +distribution_params[15]
    #                                 ,
    #                                 sdlog = distribution_params[6] +distribution_params[16]
    # )
    # synthetic_widths_low <- rlnorm(n = count,
    #                                meanlog = distribution_params[5] -distribution_params[15]
    #                                ,
    #                                sdlog = distribution_params[6] -distribution_params[16]
    # )
    
    synthetic_widths <- as.data.frame(synthetic_widths)
    # synthetic_widths_high <- as.data.frame(synthetic_widths_high)
    # synthetic_widths_low <- as.data.frame(synthetic_widths_low)
    colnames(synthetic_widths) <- c("Width_m")
    # colnames(synthetic_widths_high) <- c("Width_m_high")
    # colnames(synthetic_widths_low) <- c("Width_m_low")
    synthetic_widths <- mutate(synthetic_widths, Surf_Area_km2 = Width_m * 30 * 0.000001)
    # synthetic_widths_high <- mutate(synthetic_widths_high, Surf_Area_km2_high = Width_m_high * 30 * 0.000001)
    # synthetic_widths_low <- mutate(synthetic_widths_low, Surf_Area_km2_low = Width_m_low * 30 * 0.000001)
    order_area <- sum(synthetic_widths$Surf_Area_km2)
    # order_area_high <- sum(synthetic_widths_high$Surf_Area_km2_high)
    # order_area_low <- sum(synthetic_widths_low$Surf_Area_km2_low)
    tempDF <- c(i, order_area #, order_area_high, order_area_low
                )
    areaDF <- rbind(areaDF, tempDF)
  }
  colnames(areaDF) <- c("Order", "Surf_Area_km2" #, "Surf_Area_km2_high", "Surf_Area_km2_low"
                        )
  return(areaDF)
}

#Function to estimate RSSA using complex method and nested MC
complex_rssa_estimate_nested <- function(statDF, lengthDF, orders){
  areaDF <- data.frame(matrix(ncol = 2 #4
                              , nrow = 0))
  for(i in orders){
    # print(i)
    length <- lengthDF[i,2] * lengthDF[i,3]
    count <- round(length/30)
    meanlog_index <- 2*i
    sdlog_index <- meanlog_index + 1
    synthetic_widths <- rlnorm(n = count,
                               meanlog = statDF[1,meanlog_index],
                               sdlog = statDF[1,sdlog_index]) 
    # synthetic_widths_high <- rlnorm(n = count,
    #                                 meanlog = distribution_params[5] +distribution_params[15]
    #                                 ,
    #                                 sdlog = distribution_params[6] +distribution_params[16]
    # )
    # synthetic_widths_low <- rlnorm(n = count,
    #                                meanlog = distribution_params[5] -distribution_params[15]
    #                                ,
    #                                sdlog = distribution_params[6] -distribution_params[16]
    # )
    
    synthetic_widths <- as.data.frame(synthetic_widths)
    # synthetic_widths_high <- as.data.frame(synthetic_widths_high)
    # synthetic_widths_low <- as.data.frame(synthetic_widths_low)
    colnames(synthetic_widths) <- c("Width_m")
    # colnames(synthetic_widths_high) <- c("Width_m_high")
    # colnames(synthetic_widths_low) <- c("Width_m_low")
    synthetic_widths <- mutate(synthetic_widths, Surf_Area_km2 = Width_m * 30 * 0.000001)
    # synthetic_widths_high <- mutate(synthetic_widths_high, Surf_Area_km2_high = Width_m_high * 30 * 0.000001)
    # synthetic_widths_low <- mutate(synthetic_widths_low, Surf_Area_km2_low = Width_m_low * 30 * 0.000001)
    order_area <- sum(synthetic_widths$Surf_Area_km2)
    # order_area_high <- sum(synthetic_widths_high$Surf_Area_km2_high)
    # order_area_low <- sum(synthetic_widths_low$Surf_Area_km2_low)
    tempDF <- c(i, order_area #, order_area_high, order_area_low
    )
    areaDF <- rbind(areaDF, tempDF)
  }
  colnames(areaDF) <- c("Order", "Surf_Area_km2" #, "Surf_Area_km2_high", "Surf_Area_km2_low"
  )
  return(areaDF)
}

#Function to return synthetic datasets for plotting a histogram
generate_synth_datasets <- function(statDF, lengthDF, orders){
  synthDF <- data.frame(matrix(ncol = 2, nrow = 0))
  for(i in orders){
    print(i)
    distribution_params <- statDF[,i]
    length <- lengthDF[i,2] * lengthDF[i,3]
    count <- round(length/30)
    print('Best fit: lognormal')
    synthetic_widths <- rlnorm(n = count,
                               meanlog = distribution_params[3],
                               sdlog = distribution_params[4]) 
    synthetic_widths <- as.data.frame(synthetic_widths)
    synthetic_widths <- mutate(synthetic_widths, Order = i)
    colnames(synthetic_widths) <- c("Width_m", "Order")
    synthDF <- rbind(synthDF, synthetic_widths)
  }
  colnames(synthDF) <- c("Width", "Order")
  return(synthDF)
}

#Function to calculate percent of land surface area from RSSA
percent_river <- function(RSSA){
  #Land Area (km2) of Mississippi Basin
  mississippi_land_area <- 3287987.2594242636
  RSSA_percent <- 100*RSSA/mississippi_land_area
  return(RSSA_percent)
}

#Function to estimate surface area, including Monte Carlo analysis
RSSA_estimation_MC <- function(w, order_list, lengthDF, nRun = 1, pareto = F, numSD = 1){
  
  # get start time:
  old <- Sys.time()
  
  #Create large table to store outputs of each Monte Carlo run
  if(pareto == T){
    ensembleNames <- c("nRun",
                       "numSD",
                       "nRemoved",
                       "grwlAdditive",
                       "s2Additive",
                       "naipAdditive",
                       "fieldAdditive",
                       "o1_lnMean",
                       "o1_lnSD",
                       "o1_lnW",
                       # "o1_lnKSd",
                       # "o1_lnKSp",
                       "o1_pXmin",
                       "o1_pAlpha",
                       "o1_pW",
                       # "o1_pKSd",
                       # "o1_pKSp",
                       "o1_RSSA_km2",
                       "o2_lnMean",
                       "o2_lnSD",
                       "o2_lnW",
                       # "o2_lnKSd",
                       # "o2_lnKSp",
                       "o2_pXmin",
                       "o2_pAlpha",
                       "o2_pW",
                       # "o2_pKSd",
                       # "o2_pKSp",
                       "o2_RSSA_km2",
                       "o3_lnMean",
                       "o3_lnSD",
                       "o3_lnW",
                       # "o3_lnKSd",
                       # "o3_lnKSp",
                       "o3_pXmin",
                       "o3_pAlpha",
                       "o3_pW",
                       # "o3_pKSd",
                       # "o3_pKSp",
                       
                       "o3_RSSA_km2",
                       "o4_lnMean",
                       "o4_lnSD",
                       "o4_lnW",
                       # "o4_lnKSd",
                       # "o4_lnKSp",
                       "o4_pXmin",
                       "o4_pAlpha",
                       "o4_pW",
                       # "o4_pKSd",
                       # "o4_pKSp",
                       "o4_RSSA_km2",
                       "o5_lnMean",
                       "o5_lnSD",
                       "o5_lnW",
                       # "o5_lnKSd",
                       # "o5_lnKSp",
                       "o5_pXmin",
                       "o5_pAlpha",
                       "o5_pW",
                       # "o5_pKSd",
                       # "o5_pKSp",
                       "o5_RSSA_km2",
                       "o6_lnMean",
                       "o6_lnSD",
                       "o6_lnW",
                       # "o6_lnKSd",
                       # "o6_lnKSp",
                       "o6_pXmin",
                       "o6_pAlpha",
                       "o6_pW",
                       # "o6_pKSd",
                       # "o6_pKSp",
                       "o6_RSSA_km2",
                       "o7_lnMean",
                       "o7_lnSD",
                       "o7_lnW",
                       # "o7_lnKSd",
                       # "o7_lnKSp",
                       "o7_pXmin",
                       "o7_pAlpha",
                       "o7_pW",
                       # "o7_pKSd",
                       # "o7_pKSp",
                       "o7_RSSA_km2",
                       "o8_lnMean",
                       "o8_lnSD",
                       "o8_lnW",
                       # "o8_lnKSd",
                       # "o8_lnKSp",
                       "o8_pXmin",
                       "o8_pAlpha",
                       "o8_pW",
                       # "o8_pKSd",
                       # "o8_pKSp",
                       "o8_RSSA_km2",
                       "o9_lnMean",
                       "o9_lnSD",
                       "o9_lnW",
                       # "o9_lnKSd",
                       # "o9_lnKSp",
                       "o9_pXmin",
                       "o9_pAlpha",
                       "o9_pW",
                       # "o9_pKSd",
                       # "o9_pKSp",
                       "o9_RSSA_km2",
                       "o10_lnMean",
                       "o10_lnSD",
                       "o10_lnW",
                       # "o10_lnKSd",
                       # "o10_lnKSp",
                       "o10_pXmin",
                       "o10_pAlpha",
                       "o10_pW",
                       # "o10_pKSd",
                       # "o10_pKSp",
                       "o10_RSSA_km2",
                       "o11_lnMean",
                       "o11_lnSD",
                       "o11_lnW",
                       # "o11_lnKSd",
                       # "o11_lnKSp",
                       "o11_pXmin",
                       "o11_pAlpha",
                       "o11_pW",
                       # "o11_pKSd",
                       # "o11_pKSp",
                       "o11_RSSA_km2",
                       "o12_lnMean",
                       "o12_lnSD",
                       "o12_lnW",
                       # "o12_lnKSd",
                       # "o12_lnKSp",
                       "o12_pXmin",
                       "o12_pAlpha",
                       "o12_pW",
                       # "o12_pKSd",
                       # "o12_pKSp",
                       "o12_RSSA_km2",
                       "o13_lnMean",
                       "o13_lnSD",
                       "o13_lnW",
                       # "o13_lnKSd",
                       # "o13_lnKSp",
                       "o13_pXmin",
                       "o13_pAlpha",
                       "o13_pW",
                       # "o13_pKSd",
                       # "o13_pKSp",
                       "o13_RSSA_km2",
                       "RSSA_km2",
                       "RSSA_perc"
    )
  }
  else{
    ensembleNames <- c("nRun",
                       "numSD",
                       "nRemoved",
                       "grwlAdditive",
                       "s2Additive",
                       "naipAdditive",
                       "fieldAdditive",
                       "o1_lnMean",
                       "o1_lnSD",
                       "o1_lnW",
                       # "o1_lnKSd",
                       # "o1_lnKSp",
                       "o1_RSSA_km2",
                       "o2_lnMean",
                       "o2_lnSD",
                       "o2_lnW",
                       # "o2_lnKSd",
                       # "o2_lnKSp",
                       "o2_RSSA_km2",
                       "o3_lnMean",
                       "o3_lnSD",
                       "o3_lnW",
                       # "o3_lnKSd",
                       # "o3_lnKSp",
                       "o3_RSSA_km2",
                       "o4_lnMean",
                       "o4_lnSD",
                       "o4_lnW",
                       # "o4_lnKSd",
                       # "o4_lnKSp",
                       "o4_RSSA_km2",
                       "o5_lnMean",
                       "o5_lnSD",
                       "o5_lnW",
                       # "o5_lnKSd",
                       # "o5_lnKSp",
                       "o5_RSSA_km2",
                       "o6_lnMean",
                       "o6_lnSD",
                       "o6_lnW",
                       # "o6_lnKSd",
                       # "o6_lnKSp",
                       "o6_RSSA_km2",
                       "o7_lnMean",
                       "o7_lnSD",
                       "o7_lnW",
                       # "o7_lnKSd",
                       # "o7_lnKSp",
                       "o7_RSSA_km2",
                       "o8_lnMean",
                       "o8_lnSD",
                       "o8_lnW",
                       # "o8_lnKSd",
                       # "o8_lnKSp",
                       "o8_RSSA_km2",
                       "o9_lnMean",
                       "o9_lnSD",
                       "o9_lnW",
                       # "o9_lnKSd",
                       # "o9_lnKSp",
                       "o9_RSSA_km2",
                       "o10_lnMean",
                       "o10_lnSD",
                       "o10_lnW",
                       # "o10_lnKSd",
                       # "o10_lnKSp",
                       "o10_RSSA_km2",
                       "o11_lnMean",
                       "o11_lnSD",
                       "o11_lnW",
                       # "o11_lnKSd",
                       # "o11_lnKSp",
                       "o11_RSSA_km2",
                       "o12_lnMean",
                       "o12_lnSD",
                       "o12_lnW",
                       # "o12_lnKSd",
                       # "o12_lnKSp",
                       "o12_RSSA_km2",
                       "o13_lnMean",
                       "o13_lnSD",
                       "o13_lnW",
                       # "o13_lnKSd",
                       # "o13_lnKSp",
                       "o13_RSSA_km2",
                       "RSSA_km2",
                       "RSSA_perc"
    )
  }
  ensembleTable <- data.frame(array(NA, c(nRun, length(ensembleNames))))
  names(ensembleTable) <- ensembleNames
  
  
  #Create MC variables to propagate uncertainty of width measurements (from width residuals)
  if(nRun > 1){
    grwlUC <- 43.14193*numSD
    s2UC <- 37.26073*numSD
    naipUC <- 1.392238*numSD
    fieldUC <- 0.01*numSD
    
    grwlAdditive_list <- sort(rnorm(n = nRun,
                                      mean = 1,
                                      sd = grwlUC))
    s2Additive_list <- sort(rnorm(n = nRun,
                                      mean = 1,
                                      sd = s2UC))
    naipAdditive_list <- sort(rnorm(n = nRun,
                                    mean = 1,
                                    sd = naipUC))
    fieldAdditive_list <- sort(rnorm(n = nRun,
                                     mean = 0,
                                     sd = fieldUC))
  }
  else{
    grwlAdditive_list <- 0
    s2Additive_list <- 0
    naipAdditive_list <- 0
    fieldAdditive_list <- 0
  }
  
  widthMCvars <- data.frame(nRun = 1:nRun,
                            grwlAdditive = grwlAdditive_list,
                            s2Additive = s2Additive_list,
                            naipAdditive = naipAdditive_list,
                            fieldAdditive = fieldAdditive_list)
  for(i in 1:nRun){
    print(paste("Run:", i, "of", nRun))
    w_temp_1 <- dplyr::mutate(w, Width = case_when(Basin == "SSVC" ~ (Width + widthMCvars[i,5]),
                                                   Basin == "SVC" ~ (Width + widthMCvars[i,4]),
                                                   Basin == "Platte" ~ (Width + widthMCvars[i,3]),
                                                   Basin == "Mississippi" ~ (Width + widthMCvars[i,2])))
    w_temp <- dplyr::filter(w_temp_1, Width > 0)
    nRemoved <- nrow(w_temp_1) - nrow(w_temp)
    
    #Run the mle_by_order function on width data
    orderDF <- mle_by_order(width_data = w_temp,
                            orders = order_list,
                            parOption = pareto)
    
    
    #Run function to calculate complex estimate of RSSA by order
    RSSA <- complex_rssa_estimate(orderDF, lengthDF, order_list)
    
    #Estimate total RSSA using complex estimates by order
    estimate_RSSA <- sum(RSSA[,2])
    estimate_RSSA_perc <- percent_river(estimate_RSSA)
    print(paste("RSSA:", estimate_RSSA_perc, "% of land surface area"))
    # estimate_RSSA_high <- sum(RSSA[,3])
    # estimate_RSSA_low <- sum(RSSA[,4])
    
    #If not a MC run, remove all rows but the first one
    if(nRun == 1){ensembleTable <- data.frame(
      Order1 = c(orderDF[3,1], orderDF[4,1], orderDF[5,1], RSSA[1,2]),
      Order2 = c(orderDF[3,2], orderDF[4,2], orderDF[5,2], RSSA[2,2]),
      Order3 = c(orderDF[3,3], orderDF[4,3], orderDF[5,3], RSSA[3,2]),
      Order4 = c(orderDF[3,4], orderDF[4,4], orderDF[5,4], RSSA[4,2]),
      Order5 = c(orderDF[3,5], orderDF[4,5], orderDF[5,5], RSSA[5,2]),
      Order6 = c(orderDF[3,6], orderDF[4,6], orderDF[5,6], RSSA[6,2]),
      Order7 = c(orderDF[3,7], orderDF[4,7], orderDF[5,7], RSSA[7,2]),
      Order8 = c(orderDF[3,8], orderDF[4,8], orderDF[5,8], RSSA[8,2]),
      Order9 = c(orderDF[3,9], orderDF[4,9], orderDF[5,9], RSSA[9,2]),
      Order10 = c(orderDF[3,10], orderDF[4,10], orderDF[5,10], RSSA[10,2]),
      Order11 = c(orderDF[3,11], orderDF[4,11], orderDF[5,11], RSSA[11,2]),
      Order12 = c(orderDF[3,12], orderDF[4,12], orderDF[5,12], RSSA[12,2]),
      Order13 = c(orderDF[3,13], orderDF[4,13], orderDF[5,13], RSSA[13,2])
      )
    return(ensembleTable)
    }
    else{
      
      #Write data to ensemble table
      if(pareto == T){
        ensembleTable[i,] <- as.vector(c(
          i,
          numSD,
          nRemoved,
          widthMCvars[i,2],
          widthMCvars[i,3],
          widthMCvars[i,4],
          widthMCvars[i,5],
          orderDF[3,1],
          orderDF[4,1],
          orderDF[5,1],
          # orderDF[6,1],
          # orderDF[7,1],
          orderDF[10,1],
          orderDF[11,1],
          orderDF[12,1],
          # orderDF[13,1],
          # orderDF[14,1],
          RSSA[1,2],
          orderDF[3,2],
          orderDF[4,2],
          orderDF[5,2],
          # orderDF[6,2],
          # orderDF[7,2],
          orderDF[10,2],
          orderDF[11,2],
          orderDF[12,2],
          # orderDF[13,2],
          # orderDF[14,2],
          RSSA[2,2],
          orderDF[3,3],
          orderDF[4,3],
          orderDF[5,3],
          # orderDF[6,3],
          # orderDF[7,3],
          orderDF[10,3],
          orderDF[11,3],
          orderDF[12,3],
          # orderDF[13,3],
          # orderDF[14,3],
          RSSA[3,2],
          orderDF[3,4],
          orderDF[4,4],
          orderDF[5,4],
          # orderDF[6,4],
          # orderDF[7,4],
          orderDF[10,4],
          orderDF[11,4],
          orderDF[12,4],
          # orderDF[13,4],
          # orderDF[14,4],
          RSSA[4,2],
          orderDF[3,5],
          orderDF[4,5],
          orderDF[5,5],
          # orderDF[6,5],
          # orderDF[7,5],
          orderDF[10,5],
          orderDF[11,5],
          orderDF[12,5],
          # orderDF[13,5],
          # orderDF[14,5],
          RSSA[5,2],
          orderDF[3,6],
          orderDF[4,6],
          orderDF[5,6],
          # orderDF[6,6],
          # orderDF[7,6],
          orderDF[10,6],
          orderDF[11,6],
          orderDF[12,6],
          # orderDF[13,6],
          # orderDF[14,6],
          RSSA[6,2],
          orderDF[3,7],
          orderDF[4,7],
          orderDF[5,7],
          # orderDF[6,7],
          # orderDF[7,7],
          orderDF[10,7],
          orderDF[11,7],
          orderDF[12,7],
          # orderDF[13,7],
          # orderDF[14,7],
          RSSA[7,2],
          orderDF[3,8],
          orderDF[4,8],
          orderDF[5,8],
          # orderDF[6,8],
          # orderDF[7,8],
          orderDF[10,8],
          orderDF[11,8],
          orderDF[12,8],
          # orderDF[13,8],
          # orderDF[14,8],
          RSSA[8,2],
          orderDF[3,9],
          orderDF[4,9],
          orderDF[5,9],
          # orderDF[6,9],
          # orderDF[7,9],
          orderDF[10,9],
          orderDF[11,9],
          orderDF[12,9],
          # orderDF[13,9],
          # orderDF[14,9],
          RSSA[9,2],
          orderDF[3,10],
          orderDF[4,10],
          orderDF[5,10],
          # orderDF[6,10],
          # orderDF[7,10],
          orderDF[10,10],
          orderDF[11,10],
          orderDF[12,10],
          # orderDF[13,10],
          # orderDF[14,10],
          RSSA[10,2],
          orderDF[3,11],
          orderDF[4,11],
          orderDF[5,11],
          # orderDF[6,11],
          # orderDF[7,11],
          orderDF[10,11],
          orderDF[11,11],
          orderDF[12,11],
          # orderDF[13,11],
          # orderDF[14,11],
          RSSA[11,2],
          orderDF[3,12],
          orderDF[4,12],
          orderDF[5,12],
          # orderDF[6,12],
          # orderDF[7,12],
          orderDF[10,12],
          orderDF[11,12],
          orderDF[12,12],
          # orderDF[13,12],
          # orderDF[14,12],
          RSSA[12,2],
          orderDF[3,13],
          orderDF[4,13],
          orderDF[5,13],
          # orderDF[6,13],
          # orderDF[7,13],
          orderDF[10,13],
          orderDF[11,13],
          orderDF[12,13],
          # orderDF[13,13],
          # orderDF[14,13],
          RSSA[13,2],
          estimate_RSSA,
          estimate_RSSA_perc
        ))
      }
      else{
        ensembleTable[i,] <- as.vector(c(
          i,
          numSD,
          nRemoved,
          widthMCvars[i,2],
          widthMCvars[i,3],
          widthMCvars[i,4],
          widthMCvars[i,5],
          orderDF[3,1],
          orderDF[4,1],
          orderDF[5,1],
          # orderDF[6,1],
          # orderDF[7,1],
          RSSA[1,2],
          orderDF[3,2],
          orderDF[4,2],
          orderDF[5,2],
          # orderDF[6,2],
          # orderDF[7,2],
          RSSA[2,2],
          orderDF[3,3],
          orderDF[4,3],
          orderDF[5,3],
          # orderDF[6,3],
          # orderDF[7,3],
          RSSA[3,2],
          orderDF[3,4],
          orderDF[4,4],
          orderDF[5,4],
          # orderDF[6,4],
          # orderDF[7,4],
          RSSA[4,2],
          orderDF[3,5],
          orderDF[4,5],
          orderDF[5,5],
          # orderDF[6,5],
          # orderDF[7,5],
          RSSA[5,2],
          orderDF[3,6],
          orderDF[4,6],
          orderDF[5,6],
          # orderDF[6,6],
          # orderDF[7,6],
          RSSA[6,2],
          orderDF[3,7],
          orderDF[4,7],
          orderDF[5,7],
          # orderDF[6,7],
          # orderDF[7,7],
          RSSA[7,2],
          orderDF[3,8],
          orderDF[4,8],
          orderDF[5,8],
          # orderDF[6,8],
          # orderDF[7,8],
          RSSA[8,2],
          orderDF[3,9],
          orderDF[4,9],
          orderDF[5,9],
          # orderDF[6,9],
          # orderDF[7,9],
          RSSA[9,2],
          orderDF[3,10],
          orderDF[4,10],
          orderDF[5,10],
          # orderDF[6,10],
          # orderDF[7,10],
          RSSA[10,2],
          orderDF[3,11],
          orderDF[4,11],
          orderDF[5,11],
          # orderDF[6,11],
          # orderDF[7,11],
          RSSA[11,2],
          orderDF[3,12],
          orderDF[4,12],
          orderDF[5,12],
          # orderDF[6,12],
          # orderDF[7,12],
          RSSA[12,2],
          orderDF[3,13],
          orderDF[4,13],
          orderDF[5,13],
          # orderDF[6,13],
          # orderDF[7,13],
          RSSA[13,2],
          estimate_RSSA,
          estimate_RSSA_perc
        ))
      }
      new <- Sys.time() - old
      print(new)
    }
  }
  return(ensembleTable)
}

#Function to estimate surface area, including nested Monte Carlo analysis (stochastic distribution of fitted parameters)
# RSSA_estimation_MC_nested <- function(w, order_list, lengthDF, nRun = 1, pareto = F, numSD = 1){
# 
#   # get start time:
#   old <- Sys.time()
# 
#   #Create large table to store outputs of each Monte Carlo run
#   if(pareto == T){
#     ensembleNames <- c("nRun",
#                        "numSD",
#                        "nRemoved",
#                        "rsMultiplier",
#                        "fieldAdditive",
#                        "o1_lnMean",
#                        "o1_lnSD",
#                        "o1_lnW",
#                        # "o1_lnKSd",
#                        # "o1_lnKSp",
#                        "o1_pXmin",
#                        "o1_pAlpha",
#                        "o1_pW",
#                        # "o1_pKSd",
#                        # "o1_pKSp",
#                        "o1_RSSA_km2",
#                        "o2_lnMean",
#                        "o2_lnSD",
#                        "o2_lnW",
#                        # "o2_lnKSd",
#                        # "o2_lnKSp",
#                        "o2_pXmin",
#                        "o2_pAlpha",
#                        "o2_pW",
#                        # "o2_pKSd",
#                        # "o2_pKSp",
#                        "o2_RSSA_km2",
#                        "o3_lnMean",
#                        "o3_lnSD",
#                        "o3_lnW",
#                        # "o3_lnKSd",
#                        # "o3_lnKSp",
#                        "o3_pXmin",
#                        "o3_pAlpha",
#                        "o3_pW",
#                        # "o3_pKSd",
#                        # "o3_pKSp",
# 
#                        "o3_RSSA_km2",
#                        "o4_lnMean",
#                        "o4_lnSD",
#                        "o4_lnW",
#                        # "o4_lnKSd",
#                        # "o4_lnKSp",
#                        "o4_pXmin",
#                        "o4_pAlpha",
#                        "o4_pW",
#                        # "o4_pKSd",
#                        # "o4_pKSp",
#                        "o4_RSSA_km2",
#                        "o5_lnMean",
#                        "o5_lnSD",
#                        "o5_lnW",
#                        # "o5_lnKSd",
#                        # "o5_lnKSp",
#                        "o5_pXmin",
#                        "o5_pAlpha",
#                        "o5_pW",
#                        # "o5_pKSd",
#                        # "o5_pKSp",
#                        "o5_RSSA_km2",
#                        "o6_lnMean",
#                        "o6_lnSD",
#                        "o6_lnW",
#                        # "o6_lnKSd",
#                        # "o6_lnKSp",
#                        "o6_pXmin",
#                        "o6_pAlpha",
#                        "o6_pW",
#                        # "o6_pKSd",
#                        # "o6_pKSp",
#                        "o6_RSSA_km2",
#                        "o7_lnMean",
#                        "o7_lnSD",
#                        "o7_lnW",
#                        # "o7_lnKSd",
#                        # "o7_lnKSp",
#                        "o7_pXmin",
#                        "o7_pAlpha",
#                        "o7_pW",
#                        # "o7_pKSd",
#                        # "o7_pKSp",
#                        "o7_RSSA_km2",
#                        "o8_lnMean",
#                        "o8_lnSD",
#                        "o8_lnW",
#                        # "o8_lnKSd",
#                        # "o8_lnKSp",
#                        "o8_pXmin",
#                        "o8_pAlpha",
#                        "o8_pW",
#                        # "o8_pKSd",
#                        # "o8_pKSp",
#                        "o8_RSSA_km2",
#                        "o9_lnMean",
#                        "o9_lnSD",
#                        "o9_lnW",
#                        # "o9_lnKSd",
#                        # "o9_lnKSp",
#                        "o9_pXmin",
#                        "o9_pAlpha",
#                        "o9_pW",
#                        # "o9_pKSd",
#                        # "o9_pKSp",
#                        "o9_RSSA_km2",
#                        "o10_lnMean",
#                        "o10_lnSD",
#                        "o10_lnW",
#                        # "o10_lnKSd",
#                        # "o10_lnKSp",
#                        "o10_pXmin",
#                        "o10_pAlpha",
#                        "o10_pW",
#                        # "o10_pKSd",
#                        # "o10_pKSp",
#                        "o10_RSSA_km2",
#                        "o11_lnMean",
#                        "o11_lnSD",
#                        "o11_lnW",
#                        # "o11_lnKSd",
#                        # "o11_lnKSp",
#                        "o11_pXmin",
#                        "o11_pAlpha",
#                        "o11_pW",
#                        # "o11_pKSd",
#                        # "o11_pKSp",
#                        "o11_RSSA_km2",
#                        "o12_lnMean",
#                        "o12_lnSD",
#                        "o12_lnW",
#                        # "o12_lnKSd",
#                        # "o12_lnKSp",
#                        "o12_pXmin",
#                        "o12_pAlpha",
#                        "o12_pW",
#                        # "o12_pKSd",
#                        # "o12_pKSp",
#                        "o12_RSSA_km2",
#                        "o13_lnMean",
#                        "o13_lnSD",
#                        "o13_lnW",
#                        # "o13_lnKSd",
#                        # "o13_lnKSp",
#                        "o13_pXmin",
#                        "o13_pAlpha",
#                        "o13_pW",
#                        # "o13_pKSd",
#                        # "o13_pKSp",
#                        "o13_RSSA_km2",
#                        "RSSA_km2",
#                        "RSSA_perc"
#     )
#   }
#   else{
#     ensembleNames <- c("nRun",
#                        "numSD",
#                        "nRemoved",
#                        "rsMultiplier",
#                        "fieldAdditive",
#                        "o1_lnMean",
#                        "o1_lnSD",
#                        "o1_lnW",
#                        # "o1_lnKSd",
#                        # "o1_lnKSp",
#                        "o1_RSSA_km2",
#                        "o2_lnMean",
#                        "o2_lnSD",
#                        "o2_lnW",
#                        # "o2_lnKSd",
#                        # "o2_lnKSp",
#                        "o2_RSSA_km2",
#                        "o3_lnMean",
#                        "o3_lnSD",
#                        "o3_lnW",
#                        # "o3_lnKSd",
#                        # "o3_lnKSp",
#                        "o3_RSSA_km2",
#                        "o4_lnMean",
#                        "o4_lnSD",
#                        "o4_lnW",
#                        # "o4_lnKSd",
#                        # "o4_lnKSp",
#                        "o4_RSSA_km2",
#                        "o5_lnMean",
#                        "o5_lnSD",
#                        "o5_lnW",
#                        # "o5_lnKSd",
#                        # "o5_lnKSp",
#                        "o5_RSSA_km2",
#                        "o6_lnMean",
#                        "o6_lnSD",
#                        "o6_lnW",
#                        # "o6_lnKSd",
#                        # "o6_lnKSp",
#                        "o6_RSSA_km2",
#                        "o7_lnMean",
#                        "o7_lnSD",
#                        "o7_lnW",
#                        # "o7_lnKSd",
#                        # "o7_lnKSp",
#                        "o7_RSSA_km2",
#                        "o8_lnMean",
#                        "o8_lnSD",
#                        "o8_lnW",
#                        # "o8_lnKSd",
#                        # "o8_lnKSp",
#                        "o8_RSSA_km2",
#                        "o9_lnMean",
#                        "o9_lnSD",
#                        "o9_lnW",
#                        # "o9_lnKSd",
#                        # "o9_lnKSp",
#                        "o9_RSSA_km2",
#                        "o10_lnMean",
#                        "o10_lnSD",
#                        "o10_lnW",
#                        # "o10_lnKSd",
#                        # "o10_lnKSp",
#                        "o10_RSSA_km2",
#                        "o11_lnMean",
#                        "o11_lnSD",
#                        "o11_lnW",
#                        # "o11_lnKSd",
#                        # "o11_lnKSp",
#                        "o11_RSSA_km2",
#                        "o12_lnMean",
#                        "o12_lnSD",
#                        "o12_lnW",
#                        # "o12_lnKSd",
#                        # "o12_lnKSp",
#                        "o12_RSSA_km2",
#                        "o13_lnMean",
#                        "o13_lnSD",
#                        "o13_lnW",
#                        # "o13_lnKSd",
#                        # "o13_lnKSp",
#                        "o13_RSSA_km2",
#                        "RSSA_km2",
#                        "RSSA_perc"
#     )
#   }
#   ensembleTable <- data.frame(array(NA, c(nRun*nRun*nRun, length(ensembleNames))))
#   names(ensembleTable) <- ensembleNames
# 
#   #Create a count variable to keep track of iteration
#   count <- 1
#   #Create MC variables to propagate uncertainty of width measurements
#   if(nRun > 1){
#     rsUC <- 0.04*numSD
#     fieldUC <- 0.01*numSD
# 
#     rsMultiplier_list <- sort(rnorm(n = nRun,
#                                     mean = 1,
#                                     sd = rsUC))
#     fieldAdditive_list <- sort(rnorm(n = nRun,
#                                      mean = 0,
#                                      sd = fieldUC))
#   }
#   else{
#     rsMultiplier_list <- 1
#     fieldAdditive_list <- 0
#   }
# 
#   widthMCvars <- data.frame(nRun = 1:nRun,
#                             rsMultiplier = rsMultiplier_list,
#                             fieldAdditive = fieldAdditive_list)
#   for(i in 1:nRun){
#     w_temp_1 <- dplyr::mutate(w, Width = case_when(Basin == "SVC" | Basin == "Platte" | Basin == "Mississippi" ~ Width*widthMCvars[i,2],
#                                                  Basin == "SSVC" ~ (Width + widthMCvars[i,3])))
#     w_temp <- dplyr::filter(w_temp_1, Width > 0)
# 
#     nRemoved <- nrow(w_temp1) - nrow(w_temp)
# 
#     #Run the mle_by_order function on width data
#     orderDF <- mle_by_order(width_data = w_temp,
#                             orders = order_list,
#                             parOption = pareto)
#     if(nRun > 1){
#       o1_lnMeanUC <- orderDF[6,1]*numSD
#       o1_lnSDUC <- orderDF[7,1]*numSD
#       o2_lnMeanUC <- orderDF[6,2]*numSD
#       o2_lnSDUC <- orderDF[7,2]*numSD
#       o3_lnMeanUC <- orderDF[6,3]*numSD
#       o3_lnSDUC <- orderDF[7,3]*numSD
#       o4_lnMeanUC <- orderDF[6,4]*numSD
#       o4_lnSDUC <- orderDF[7,4]*numSD
#       o5_lnMeanUC <- orderDF[6,5]*numSD
#       o5_lnSDUC <- orderDF[7,5]*numSD
#       o6_lnMeanUC <- orderDF[6,6]*numSD
#       o6_lnSDUC <- orderDF[7,6]*numSD
#       o7_lnMeanUC <- orderDF[6,7]*numSD
#       o7_lnSDUC <- orderDF[7,7]*numSD
#       o8_lnMeanUC <- orderDF[6,8]*numSD
#       o8_lnSDUC <- orderDF[7,8]*numSD
#       o9_lnMeanUC <- orderDF[6,9]*numSD
#       o9_lnSDUC <- orderDF[7,9]*numSD
#       o10_lnMeanUC <- orderDF[6,10]*numSD
#       o10_lnSDUC <- orderDF[7,10]*numSD
#       o11_lnMeanUC <- orderDF[6,11]*numSD
#       o11_lnSDUC <- orderDF[7,11]*numSD
#       o12_lnMeanUC <- orderDF[6,12]*numSD
#       o12_lnSDUC <- orderDF[7,12]*numSD
#       o13_lnMeanUC <- orderDF[6,13]*numSD
#       o13_lnSDUC <- orderDF[7,13]*numSD
# 
#       o1_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,1],
#                                    sd = o1_lnMeanUC))
#       o1_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,1],
#                                  sd = o1_lnSDUC))
#       o2_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,2],
#                                    sd = o2_lnMeanUC))
#       o2_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,2],
#                                  sd = o2_lnSDUC))
#       o3_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,3],
#                                    sd = o3_lnMeanUC))
#       o3_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,3],
#                                  sd = o3_lnSDUC))
#       o4_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,4],
#                                    sd = o4_lnMeanUC))
#       o4_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,4],
#                                  sd = o4_lnSDUC))
#       o5_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,5],
#                                    sd = o5_lnMeanUC))
#       o5_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,5],
#                                  sd = o5_lnSDUC))
#       o6_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,6],
#                                    sd = o6_lnMeanUC))
#       o6_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,6],
#                                  sd = o6_lnSDUC))
#       o7_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,7],
#                                    sd = o7_lnMeanUC))
#       o7_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,7],
#                                  sd = o7_lnSDUC))
#       o8_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,8],
#                                    sd = o8_lnMeanUC))
#       o8_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,8],
#                                  sd = o8_lnSDUC))
#       o9_lnMean_list <- sort(rnorm(n = nRun,
#                                    mean = orderDF[3,9],
#                                    sd = o9_lnMeanUC))
#       o9_lnSD_list <- sort(rnorm(n = nRun,
#                                  mean = orderDF[4,9],
#                                  sd = o9_lnSDUC))
#       o10_lnMean_list <- sort(rnorm(n = nRun,
#                                     mean = orderDF[3,10],
#                                     sd = o10_lnMeanUC))
#       o10_lnSD_list <- sort(rnorm(n = nRun,
#                                   mean = orderDF[4,10],
#                                   sd = o10_lnSDUC))
#       o11_lnMean_list <- sort(rnorm(n = nRun,
#                                     mean = orderDF[3,11],
#                                     sd = o11_lnMeanUC))
#       o11_lnSD_list <- sort(rnorm(n = nRun,
#                                   mean = orderDF[4,11],
#                                   sd = o11_lnSDUC))
#       o12_lnMean_list <- sort(rnorm(n = nRun,
#                                     mean = orderDF[3,12],
#                                     sd = o12_lnMeanUC))
#       o12_lnSD_list <- sort(rnorm(n = nRun,
#                                   mean = orderDF[4,12],
#                                   sd = o12_lnSDUC))
#       o13_lnMean_list <- sort(rnorm(n = nRun,
#                                     mean = orderDF[3,13],
#                                     sd = o13_lnMeanUC))
#       o13_lnSD_list <- sort(rnorm(n = nRun,
#                                   mean = orderDF[4,13],
#                                   sd = o13_lnSDUC))
#     }
#     else{
#       o1_lnMean_list <- orderDF[3,1]
#       o1_lnSD_list <- orderDF[4,1]
#       o2_lnMean_list <- orderDF[3,2]
#       o2_lnSD_list <- orderDF[4,2]
#       o3_lnMean_list <- orderDF[3,3]
#       o3_lnSD_list <- orderDF[4,3]
#       o4_lnMean_list <- orderDF[3,4]
#       o4_lnSD_list <- orderDF[4,4]
#       o5_lnMean_list <- orderDF[3,5]
#       o5_lnSD_list <- orderDF[4,5]
#       o6_lnMean_list <- orderDF[3,6]
#       o6_lnSD_list <- orderDF[4,6]
#       o7_lnMean_list <- orderDF[3,7]
#       o7_lnSD_list <- orderDF[4,7]
#       o8_lnMean_list <- orderDF[3,8]
#       o8_lnSD_list <- orderDF[4,8]
#       o9_lnMean_list <- orderDF[3,9]
#       o9_lnSD_list <- orderDF[4,9]
#       o10_lnMean_list <- orderDF[3,10]
#       o10_lnSD_list <- orderDF[4,10]
#       o11_lnMean_list <- orderDF[3,11]
#       o11_lnSD_list <- orderDF[4,11]
#       o12_lnMean_list <- orderDF[3,12]
#       o12_lnSD_list <- orderDF[4,12]
#       o13_lnMean_list <- orderDF[3,13]
#       o13_lnSD_list <- orderDF[4,13]
#     }
#     distMCvars <- data.frame(nRun = 1:nRun,
#                              o1_lnMean = o1_lnMean_list,
#                              o1_lnSD = o1_lnSD_list,
#                              o2_lnMean = o2_lnMean_list,
#                              o2_lnSD = o2_lnSD_list,
#                              o3_lnMean = o3_lnMean_list,
#                              o3_lnSD = o3_lnSD_list,
#                              o4_lnMean = o4_lnMean_list,
#                              o4_lnSD = o4_lnSD_list,
#                              o5_lnMean = o5_lnMean_list,
#                              o5_lnSD = o5_lnSD_list,
#                              o6_lnMean = o6_lnMean_list,
#                              o6_lnSD = o6_lnSD_list,
#                              o7_lnMean = o7_lnMean_list,
#                              o7_lnSD = o7_lnSD_list,
#                              o8_lnMean = o8_lnMean_list,
#                              o8_lnSD = o8_lnSD_list,
#                              o9_lnMean = o9_lnMean_list,
#                              o9_lnSD = o9_lnSD_list,
#                              o10_lnMean = o10_lnMean_list,
#                              o10_lnSD = o10_lnSD_list,
#                              o11_lnMean = o11_lnMean_list,
#                              o11_lnSD = o11_lnSD_list,
#                              o12_lnMean = o12_lnMean_list,
#                              o12_lnSD = o12_lnSD_list,
#                              o13_lnMean = o13_lnMean_list,
#                              o13_lnSD = o13_lnSD_list
#     )
#     for(j in 1:nRun){
#       distribution_params <- distMCvars[j,]
# 
#       for(k in 1:nRun){
#         print(paste("Run:", count, "of", nRun*nRun*nRun))
# 
#         #Run function to calculate complex estimate of RSSA by order
#         RSSA <- complex_rssa_estimate_nested(distribution_params, lengthDF, order_list)
# 
#         #Estimate total RSSA using complex estimates by order
#         estimate_RSSA <- sum(RSSA[,2])
#         estimate_RSSA_perc <- percent_river(estimate_RSSA)
#         print(paste("RSSA:", estimate_RSSA_perc, "% of land surface area"))
#         # estimate_RSSA_high <- sum(RSSA[,3])
#         # estimate_RSSA_low <- sum(RSSA[,4])
# 
#         #If not a MC run, remove all rows but the first one
#         if(nRun == 1){
#           ensembleTable <- ensembleTable[1,]
#           iteration = 1
#         }else{
#             iteration = i+j+k
#           }
# 
#         #Write data to ensemble table
#         if(pareto == T){
#           ensembleTable[count,] <- as.vector(c(
#             count,
#             numSD,
#             nRemoved,
#             widthMCvars[i,2],
#             widthMCvars[i,3],
#             orderDF[3,1],
#             orderDF[4,1],
#             orderDF[5,1],
#             # orderDF[6,1],
#             # orderDF[7,1],
#             orderDF[10,1],
#             orderDF[11,1],
#             orderDF[12,1],
#             # orderDF[13,1],
#             # orderDF[14,1],
#             RSSA[1,2],
#             orderDF[3,2],
#             orderDF[4,2],
#             orderDF[5,2],
#             # orderDF[6,2],
#             # orderDF[7,2],
#             orderDF[10,2],
#             orderDF[11,2],
#             orderDF[12,2],
#             # orderDF[13,2],
#             # orderDF[14,2],
#             RSSA[2,2],
#             orderDF[3,3],
#             orderDF[4,3],
#             orderDF[5,3],
#             # orderDF[6,3],
#             # orderDF[7,3],
#             orderDF[10,3],
#             orderDF[11,3],
#             orderDF[12,3],
#             # orderDF[13,3],
#             # orderDF[14,3],
#             RSSA[3,2],
#             orderDF[3,4],
#             orderDF[4,4],
#             orderDF[5,4],
#             # orderDF[6,4],
#             # orderDF[7,4],
#             orderDF[10,4],
#             orderDF[11,4],
#             orderDF[12,4],
#             # orderDF[13,4],
#             # orderDF[14,4],
#             RSSA[4,2],
#             orderDF[3,5],
#             orderDF[4,5],
#             orderDF[5,5],
#             # orderDF[6,5],
#             # orderDF[7,5],
#             orderDF[10,5],
#             orderDF[11,5],
#             orderDF[12,5],
#             # orderDF[13,5],
#             # orderDF[14,5],
#             RSSA[5,2],
#             orderDF[3,6],
#             orderDF[4,6],
#             orderDF[5,6],
#             # orderDF[6,6],
#             # orderDF[7,6],
#             orderDF[10,6],
#             orderDF[11,6],
#             orderDF[12,6],
#             # orderDF[13,6],
#             # orderDF[14,6],
#             RSSA[6,2],
#             orderDF[3,7],
#             orderDF[4,7],
#             orderDF[5,7],
#             # orderDF[6,7],
#             # orderDF[7,7],
#             orderDF[10,7],
#             orderDF[11,7],
#             orderDF[12,7],
#             # orderDF[13,7],
#             # orderDF[14,7],
#             RSSA[7,2],
#             orderDF[3,8],
#             orderDF[4,8],
#             orderDF[5,8],
#             # orderDF[6,8],
#             # orderDF[7,8],
#             orderDF[10,8],
#             orderDF[11,8],
#             orderDF[12,8],
#             # orderDF[13,8],
#             # orderDF[14,8],
#             RSSA[8,2],
#             orderDF[3,9],
#             orderDF[4,9],
#             orderDF[5,9],
#             # orderDF[6,9],
#             # orderDF[7,9],
#             orderDF[10,9],
#             orderDF[11,9],
#             orderDF[12,9],
#             # orderDF[13,9],
#             # orderDF[14,9],
#             RSSA[9,2],
#             orderDF[3,10],
#             orderDF[4,10],
#             orderDF[5,10],
#             # orderDF[6,10],
#             # orderDF[7,10],
#             orderDF[10,10],
#             orderDF[11,10],
#             orderDF[12,10],
#             # orderDF[13,10],
#             # orderDF[14,10],
#             RSSA[10,2],
#             orderDF[3,11],
#             orderDF[4,11],
#             orderDF[5,11],
#             # orderDF[6,11],
#             # orderDF[7,11],
#             orderDF[10,11],
#             orderDF[11,11],
#             orderDF[12,11],
#             # orderDF[13,11],
#             # orderDF[14,11],
#             RSSA[11,2],
#             orderDF[3,12],
#             orderDF[4,12],
#             orderDF[5,12],
#             # orderDF[6,12],
#             # orderDF[7,12],
#             orderDF[10,12],
#             orderDF[11,12],
#             orderDF[12,12],
#             # orderDF[13,12],
#             # orderDF[14,12],
#             RSSA[12,2],
#             orderDF[3,13],
#             orderDF[4,13],
#             orderDF[5,13],
#             # orderDF[6,13],
#             # orderDF[7,13],
#             orderDF[10,13],
#             orderDF[11,13],
#             orderDF[12,13],
#             # orderDF[13,13],
#             # orderDF[14,13],
#             RSSA[13,2],
#             estimate_RSSA,
#             estimate_RSSA_perc
#           ))
#         }
#         else{
#           ensembleTable[count,] <- as.vector(c(
#             count,
#             numSD,
#             nRemoved,
#             widthMCvars[i,2],
#             widthMCvars[i,3],
#             orderDF[3,1],
#             orderDF[4,1],
#             orderDF[5,1],
#             # orderDF[6,1],
#             # orderDF[7,1],
#             RSSA[1,2],
#             orderDF[3,2],
#             orderDF[4,2],
#             orderDF[5,2],
#             # orderDF[6,2],
#             # orderDF[7,2],
#             RSSA[2,2],
#             orderDF[3,3],
#             orderDF[4,3],
#             orderDF[5,3],
#             # orderDF[6,3],
#             # orderDF[7,3],
#             RSSA[3,2],
#             orderDF[3,4],
#             orderDF[4,4],
#             orderDF[5,4],
#             # orderDF[6,4],
#             # orderDF[7,4],
#             RSSA[4,2],
#             orderDF[3,5],
#             orderDF[4,5],
#             orderDF[5,5],
#             # orderDF[6,5],
#             # orderDF[7,5],
#             RSSA[5,2],
#             orderDF[3,6],
#             orderDF[4,6],
#             orderDF[5,6],
#             # orderDF[6,6],
#             # orderDF[7,6],
#             RSSA[6,2],
#             orderDF[3,7],
#             orderDF[4,7],
#             orderDF[5,7],
#             # orderDF[6,7],
#             # orderDF[7,7],
#             RSSA[7,2],
#             orderDF[3,8],
#             orderDF[4,8],
#             orderDF[5,8],
#             # orderDF[6,8],
#             # orderDF[7,8],
#             RSSA[8,2],
#             orderDF[3,9],
#             orderDF[4,9],
#             orderDF[5,9],
#             # orderDF[6,9],
#             # orderDF[7,9],
#             RSSA[9,2],
#             orderDF[3,10],
#             orderDF[4,10],
#             orderDF[5,10],
#             # orderDF[6,10],
#             # orderDF[7,10],
#             RSSA[10,2],
#             orderDF[3,11],
#             orderDF[4,11],
#             orderDF[5,11],
#             # orderDF[6,11],
#             # orderDF[7,11],
#             RSSA[11,2],
#             orderDF[3,12],
#             orderDF[4,12],
#             orderDF[5,12],
#             # orderDF[6,12],
#             # orderDF[7,12],
#             RSSA[12,2],
#             orderDF[3,13],
#             orderDF[4,13],
#             orderDF[5,13],
#             # orderDF[6,13],
#             # orderDF[7,13],
#             RSSA[13,2],
#             estimate_RSSA,
#             estimate_RSSA_perc
#           ))
#         }
#         new <- Sys.time() - old
#         print(new)
#         count <- count + 1
#       }
#     }
#   }
#   return(ensembleTable)
# }

#Function to estimate surface area, including nested Monte Carlo analysis (run distribution fitting N times)
# RSSA_estimation_MC_nested_V2 <- function(w, order_list, lengthDF, nRun = 1, pareto = F){
#   
#   # get start time:
#   old <- Sys.time()
#   
#   #Create large table to store outputs of each Monte Carlo run
#   if(pareto == T){
#     ensembleNames <- c("nRun",
#                        "rsMultiplier",
#                        "fieldAdditive",
#                        "o1_lnMean",
#                        "o1_lnSD",
#                        "o1_lnW",
#                        # "o1_lnKSd",
#                        # "o1_lnKSp",
#                        "o1_pXmin",
#                        "o1_pAlpha",
#                        "o1_pW",
#                        # "o1_pKSd",
#                        # "o1_pKSp",
#                        "o1_RSSA_km2",
#                        "o2_lnMean",
#                        "o2_lnSD",
#                        "o2_lnW",
#                        # "o2_lnKSd",
#                        # "o2_lnKSp",
#                        "o2_pXmin",
#                        "o2_pAlpha",
#                        "o2_pW",
#                        # "o2_pKSd",
#                        # "o2_pKSp",
#                        "o2_RSSA_km2",
#                        "o3_lnMean",
#                        "o3_lnSD",
#                        "o3_lnW",
#                        # "o3_lnKSd",
#                        # "o3_lnKSp",
#                        "o3_pXmin",
#                        "o3_pAlpha",
#                        "o3_pW",
#                        # "o3_pKSd",
#                        # "o3_pKSp",
#                        
#                        "o3_RSSA_km2",
#                        "o4_lnMean",
#                        "o4_lnSD",
#                        "o4_lnW",
#                        # "o4_lnKSd",
#                        # "o4_lnKSp",
#                        "o4_pXmin",
#                        "o4_pAlpha",
#                        "o4_pW",
#                        # "o4_pKSd",
#                        # "o4_pKSp",
#                        "o4_RSSA_km2",
#                        "o5_lnMean",
#                        "o5_lnSD",
#                        "o5_lnW",
#                        # "o5_lnKSd",
#                        # "o5_lnKSp",
#                        "o5_pXmin",
#                        "o5_pAlpha",
#                        "o5_pW",
#                        # "o5_pKSd",
#                        # "o5_pKSp",
#                        "o5_RSSA_km2",
#                        "o6_lnMean",
#                        "o6_lnSD",
#                        "o6_lnW",
#                        # "o6_lnKSd",
#                        # "o6_lnKSp",
#                        "o6_pXmin",
#                        "o6_pAlpha",
#                        "o6_pW",
#                        # "o6_pKSd",
#                        # "o6_pKSp",
#                        "o6_RSSA_km2",
#                        "o7_lnMean",
#                        "o7_lnSD",
#                        "o7_lnW",
#                        # "o7_lnKSd",
#                        # "o7_lnKSp",
#                        "o7_pXmin",
#                        "o7_pAlpha",
#                        "o7_pW",
#                        # "o7_pKSd",
#                        # "o7_pKSp",
#                        "o7_RSSA_km2",
#                        "o8_lnMean",
#                        "o8_lnSD",
#                        "o8_lnW",
#                        # "o8_lnKSd",
#                        # "o8_lnKSp",
#                        "o8_pXmin",
#                        "o8_pAlpha",
#                        "o8_pW",
#                        # "o8_pKSd",
#                        # "o8_pKSp",
#                        "o8_RSSA_km2",
#                        "o9_lnMean",
#                        "o9_lnSD",
#                        "o9_lnW",
#                        # "o9_lnKSd",
#                        # "o9_lnKSp",
#                        "o9_pXmin",
#                        "o9_pAlpha",
#                        "o9_pW",
#                        # "o9_pKSd",
#                        # "o9_pKSp",
#                        "o9_RSSA_km2",
#                        "o10_lnMean",
#                        "o10_lnSD",
#                        "o10_lnW",
#                        # "o10_lnKSd",
#                        # "o10_lnKSp",
#                        "o10_pXmin",
#                        "o10_pAlpha",
#                        "o10_pW",
#                        # "o10_pKSd",
#                        # "o10_pKSp",
#                        "o10_RSSA_km2",
#                        "o11_lnMean",
#                        "o11_lnSD",
#                        "o11_lnW",
#                        # "o11_lnKSd",
#                        # "o11_lnKSp",
#                        "o11_pXmin",
#                        "o11_pAlpha",
#                        "o11_pW",
#                        # "o11_pKSd",
#                        # "o11_pKSp",
#                        "o11_RSSA_km2",
#                        "o12_lnMean",
#                        "o12_lnSD",
#                        "o12_lnW",
#                        # "o12_lnKSd",
#                        # "o12_lnKSp",
#                        "o12_pXmin",
#                        "o12_pAlpha",
#                        "o12_pW",
#                        # "o12_pKSd",
#                        # "o12_pKSp",
#                        "o12_RSSA_km2",
#                        "o13_lnMean",
#                        "o13_lnSD",
#                        "o13_lnW",
#                        # "o13_lnKSd",
#                        # "o13_lnKSp",
#                        "o13_pXmin",
#                        "o13_pAlpha",
#                        "o13_pW",
#                        # "o13_pKSd",
#                        # "o13_pKSp",
#                        "o13_RSSA_km2",
#                        "RSSA_km2",
#                        "RSSA_perc"
#     )
#   }
#   else{
#     ensembleNames <- c("nRun",
#                        "rsMultiplier",
#                        "fieldAdditive",
#                        "o1_lnMean",
#                        "o1_lnSD",
#                        "o1_lnW",
#                        # "o1_lnKSd",
#                        # "o1_lnKSp",
#                        "o1_RSSA_km2",
#                        "o2_lnMean",
#                        "o2_lnSD",
#                        "o2_lnW",
#                        # "o2_lnKSd",
#                        # "o2_lnKSp",
#                        "o2_RSSA_km2",
#                        "o3_lnMean",
#                        "o3_lnSD",
#                        "o3_lnW",
#                        # "o3_lnKSd",
#                        # "o3_lnKSp",
#                        "o3_RSSA_km2",
#                        "o4_lnMean",
#                        "o4_lnSD",
#                        "o4_lnW",
#                        # "o4_lnKSd",
#                        # "o4_lnKSp",
#                        "o4_RSSA_km2",
#                        "o5_lnMean",
#                        "o5_lnSD",
#                        "o5_lnW",
#                        # "o5_lnKSd",
#                        # "o5_lnKSp",
#                        "o5_RSSA_km2",
#                        "o6_lnMean",
#                        "o6_lnSD",
#                        "o6_lnW",
#                        # "o6_lnKSd",
#                        # "o6_lnKSp",
#                        "o6_RSSA_km2",
#                        "o7_lnMean",
#                        "o7_lnSD",
#                        "o7_lnW",
#                        # "o7_lnKSd",
#                        # "o7_lnKSp",
#                        "o7_RSSA_km2",
#                        "o8_lnMean",
#                        "o8_lnSD",
#                        "o8_lnW",
#                        # "o8_lnKSd",
#                        # "o8_lnKSp",
#                        "o8_RSSA_km2",
#                        "o9_lnMean",
#                        "o9_lnSD",
#                        "o9_lnW",
#                        # "o9_lnKSd",
#                        # "o9_lnKSp",
#                        "o9_RSSA_km2",
#                        "o10_lnMean",
#                        "o10_lnSD",
#                        "o10_lnW",
#                        # "o10_lnKSd",
#                        # "o10_lnKSp",
#                        "o10_RSSA_km2",
#                        "o11_lnMean",
#                        "o11_lnSD",
#                        "o11_lnW",
#                        # "o11_lnKSd",
#                        # "o11_lnKSp",
#                        "o11_RSSA_km2",
#                        "o12_lnMean",
#                        "o12_lnSD",
#                        "o12_lnW",
#                        # "o12_lnKSd",
#                        # "o12_lnKSp",
#                        "o12_RSSA_km2",
#                        "o13_lnMean",
#                        "o13_lnSD",
#                        "o13_lnW",
#                        # "o13_lnKSd",
#                        # "o13_lnKSp",
#                        "o13_RSSA_km2",
#                        "RSSA_km2",
#                        "RSSA_perc"
#     )
#   }
#   ensembleTable <- data.frame(array(NA, c(nRun*nRun*nRun, length(ensembleNames))))
#   names(ensembleTable) <- ensembleNames
#   
#   #Create count variable for keeping track of iteration
#   count <- 1
#   #Create MC variables to propagate uncertainty of width measurements
#   if(nRun > 1){
#     rsUC <- 0.04
#     fieldUC <- 0.01
#     
#     rsMultiplier_list <- sort(rnorm(n = nRun,
#                                     mean = 1,
#                                     sd = rsUC))
#     fieldAdditive_list <- sort(rnorm(n = nRun,
#                                      mean = 0,
#                                      sd = fieldUC))
#   }else{
#     rsMultiplier_list <- 1
#     fieldAdditive_list <- 0
#   }
#   
#   widthMCvars <- data.frame(nRun = 1:nRun,
#                             rsMultiplier = rsMultiplier_list,
#                             fieldAdditive = fieldAdditive_list)
#   for(i in 1:nRun){
#     w_temp <- dplyr::mutate(w, Width = case_when(Basin == "SVC" | Basin == "Platte" | Basin == "Mississippi" ~ Width*widthMCvars[i,2],
#                                                  Basin == "SSVC" ~ (Width + widthMCvars[i,3])))
#     w_temp <- dplyr::filter(w_temp, Width > 0)
#     
#     #Run the mle_by_order function on width data
#     # orderDF <- mle_by_order(width_data = w_temp,
#     #                         orders = order_list,
#     #                         parOption = pareto)
#     # if(nRun > 1){
#     #   o1_lnMeanUC <- orderDF[6,1]
#     #   o1_lnSDUC <- orderDF[7,1]
#     #   o2_lnMeanUC <- orderDF[6,2]
#     #   o2_lnSDUC <- orderDF[7,2]
#     #   o3_lnMeanUC <- orderDF[6,3]
#     #   o3_lnSDUC <- orderDF[7,3]
#     #   o4_lnMeanUC <- orderDF[6,4]
#     #   o4_lnSDUC <- orderDF[7,4]
#     #   o5_lnMeanUC <- orderDF[6,5]
#     #   o5_lnSDUC <- orderDF[7,5]
#     #   o6_lnMeanUC <- orderDF[6,6]
#     #   o6_lnSDUC <- orderDF[7,6]
#     #   o7_lnMeanUC <- orderDF[6,7]
#     #   o7_lnSDUC <- orderDF[7,7]
#     #   o8_lnMeanUC <- orderDF[6,8]
#     #   o8_lnSDUC <- orderDF[7,8]
#     #   o9_lnMeanUC <- orderDF[6,9]
#     #   o9_lnSDUC <- orderDF[7,9]
#     #   o10_lnMeanUC <- orderDF[6,10]
#     #   o10_lnSDUC <- orderDF[7,10]
#     #   o11_lnMeanUC <- orderDF[6,11]
#     #   o11_lnSDUC <- orderDF[7,11]
#     #   o12_lnMeanUC <- orderDF[6,12]
#     #   o12_lnSDUC <- orderDF[7,12]
#     #   o13_lnMeanUC <- orderDF[6,13]
#     #   o13_lnSDUC <- orderDF[7,13]
#     #   
#     #   o1_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,1],
#     #                                sd = o1_lnMeanUC))
#     #   o1_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,1],
#     #                              sd = o1_lnSDUC))
#     #   o2_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,2],
#     #                                sd = o2_lnMeanUC))
#     #   o2_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,2],
#     #                              sd = o2_lnSDUC))        
#     #   o3_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,3],
#     #                                sd = o3_lnMeanUC))
#     #   o3_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,3],
#     #                              sd = o3_lnSDUC))        
#     #   o4_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,4],
#     #                                sd = o4_lnMeanUC))
#     #   o4_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,4],
#     #                              sd = o4_lnSDUC))        
#     #   o5_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,5],
#     #                                sd = o5_lnMeanUC))
#     #   o5_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,5],
#     #                              sd = o5_lnSDUC))        
#     #   o6_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,6],
#     #                                sd = o6_lnMeanUC))
#     #   o6_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,6],
#     #                              sd = o6_lnSDUC))        
#     #   o7_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,7],
#     #                                sd = o7_lnMeanUC))
#     #   o7_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,7],
#     #                              sd = o7_lnSDUC))        
#     #   o8_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,8],
#     #                                sd = o8_lnMeanUC))
#     #   o8_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,8],
#     #                              sd = o8_lnSDUC))        
#     #   o9_lnMean_list <- sort(rnorm(n = nRun,
#     #                                mean = orderDF[3,9],
#     #                                sd = o9_lnMeanUC))
#     #   o9_lnSD_list <- sort(rnorm(n = nRun,
#     #                              mean = orderDF[4,9],
#     #                              sd = o9_lnSDUC))        
#     #   o10_lnMean_list <- sort(rnorm(n = nRun,
#     #                                 mean = orderDF[3,10],
#     #                                 sd = o10_lnMeanUC))
#     #   o10_lnSD_list <- sort(rnorm(n = nRun,
#     #                               mean = orderDF[4,10],
#     #                               sd = o10_lnSDUC))        
#     #   o11_lnMean_list <- sort(rnorm(n = nRun,
#     #                                 mean = orderDF[3,11],
#     #                                 sd = o11_lnMeanUC))
#     #   o11_lnSD_list <- sort(rnorm(n = nRun,
#     #                               mean = orderDF[4,11],
#     #                               sd = o11_lnSDUC))        
#     #   o12_lnMean_list <- sort(rnorm(n = nRun,
#     #                                 mean = orderDF[3,12],
#     #                                 sd = o12_lnMeanUC))
#     #   o12_lnSD_list <- sort(rnorm(n = nRun,
#     #                               mean = orderDF[4,12],
#     #                               sd = o12_lnSDUC))        
#     #   o13_lnMean_list <- sort(rnorm(n = nRun,
#     #                                 mean = orderDF[3,13],
#     #                                 sd = o13_lnMeanUC))
#     #   o13_lnSD_list <- sort(rnorm(n = nRun,
#     #                               mean = orderDF[4,13],
#     #                               sd = o13_lnSDUC))
#     # }
#     # else{
#     #   o1_lnMean_list <- orderDF[3,1]
#     #   o1_lnSD_list <- orderDF[4,1]
#     #   o2_lnMean_list <- orderDF[3,2]
#     #   o2_lnSD_list <- orderDF[4,2]
#     #   o3_lnMean_list <- orderDF[3,3]
#     #   o3_lnSD_list <- orderDF[4,3]
#     #   o4_lnMean_list <- orderDF[3,4]
#     #   o4_lnSD_list <- orderDF[4,4]
#     #   o5_lnMean_list <- orderDF[3,5]
#     #   o5_lnSD_list <- orderDF[4,5]
#     #   o6_lnMean_list <- orderDF[3,6]
#     #   o6_lnSD_list <- orderDF[4,6]
#     #   o7_lnMean_list <- orderDF[3,7]
#     #   o7_lnSD_list <- orderDF[4,7]
#     #   o8_lnMean_list <- orderDF[3,8]
#     #   o8_lnSD_list <- orderDF[4,8]
#     #   o9_lnMean_list <- orderDF[3,9]
#     #   o9_lnSD_list <- orderDF[4,9]
#     #   o10_lnMean_list <- orderDF[3,10]
#     #   o10_lnSD_list <- orderDF[4,10]
#     #   o11_lnMean_list <- orderDF[3,11]
#     #   o11_lnSD_list <- orderDF[4,11]
#     #   o12_lnMean_list <- orderDF[3,12]
#     #   o12_lnSD_list <- orderDF[4,12]
#     #   o13_lnMean_list <- orderDF[3,13]
#     #   o13_lnSD_list <- orderDF[4,13]
#     # }
#     # distMCvars <- data.frame(nRun = 1:nRun,
#     #                          o1_lnMean = o1_lnMean_list,
#     #                          o1_lnSD = o1_lnSD_list,
#     #                          o2_lnMean = o2_lnMean_list,
#     #                          o2_lnSD = o2_lnSD_list,
#     #                          o3_lnMean = o3_lnMean_list,
#     #                          o3_lnSD = o3_lnSD_list,
#     #                          o4_lnMean = o4_lnMean_list,
#     #                          o4_lnSD = o4_lnSD_list,
#     #                          o5_lnMean = o5_lnMean_list,
#     #                          o5_lnSD = o5_lnSD_list,
#     #                          o6_lnMean = o6_lnMean_list,
#     #                          o6_lnSD = o6_lnSD_list,
#     #                          o7_lnMean = o7_lnMean_list,
#     #                          o7_lnSD = o7_lnSD_list,
#     #                          o8_lnMean = o8_lnMean_list,
#     #                          o8_lnSD = o8_lnSD_list,
#     #                          o9_lnMean = o9_lnMean_list,
#     #                          o9_lnSD = o9_lnSD_list,
#     #                          o10_lnMean = o10_lnMean_list,
#     #                          o10_lnSD = o10_lnSD_list,
#     #                          o11_lnMean = o11_lnMean_list,
#     #                          o11_lnSD = o11_lnSD_list,
#     #                          o12_lnMean = o12_lnMean_list,
#     #                          o12_lnSD = o12_lnSD_list,
#     #                          o13_lnMean = o13_lnMean_list,
#     #                          o13_lnSD = o13_lnSD_list
#     # )
#     for(j in 1:nRun){
#       # distribution_params <- distMCvars[j,]
#       orderDF <- mle_by_order(width_data = w_temp,
#                               orders = order_list,
#                               parOption = pareto)
#       for(k in 1:nRun){
#         print(paste("Run:", count, "of", nRun*nRun*nRun))
#         
#         #Run function to calculate complex estimate of RSSA by order
#         RSSA <- complex_rssa_estimate(orderDF, lengthDF, order_list)
#         
#         #Estimate total RSSA using complex estimates by order
#         estimate_RSSA <- sum(RSSA[,2])
#         estimate_RSSA_perc <- percent_river(estimate_RSSA)
#         print(paste("RSSA:", estimate_RSSA_perc, "% of land surface area"))
#         # estimate_RSSA_high <- sum(RSSA[,3])
#         # estimate_RSSA_low <- sum(RSSA[,4])
#         
#         #If not a MC run, remove all rows but the first one
#         if(nRun == 1){
#           ensembleTable <- ensembleTable[1,]
#           iteration = 1
#         }else{
#           iteration = i+j+k
#         }
#         
#         #Write data to ensemble table
#         if(pareto == T){
#           ensembleTable[count,] <- as.vector(c(
#             count,
#             widthMCvars[i,2],
#             widthMCvars[i,3],
#             orderDF[3,1],
#             orderDF[4,1],
#             orderDF[5,1],
#             # orderDF[6,1],
#             # orderDF[7,1],
#             orderDF[10,1],
#             orderDF[11,1],
#             orderDF[12,1],
#             # orderDF[13,1],
#             # orderDF[14,1],
#             RSSA[1,2],
#             orderDF[3,2],
#             orderDF[4,2],
#             orderDF[5,2],
#             # orderDF[6,2],
#             # orderDF[7,2],
#             orderDF[10,2],
#             orderDF[11,2],
#             orderDF[12,2],
#             # orderDF[13,2],
#             # orderDF[14,2],
#             RSSA[2,2],
#             orderDF[3,3],
#             orderDF[4,3],
#             orderDF[5,3],
#             # orderDF[6,3],
#             # orderDF[7,3],
#             orderDF[10,3],
#             orderDF[11,3],
#             orderDF[12,3],
#             # orderDF[13,3],
#             # orderDF[14,3],
#             RSSA[3,2],
#             orderDF[3,4],
#             orderDF[4,4],
#             orderDF[5,4],
#             # orderDF[6,4],
#             # orderDF[7,4],
#             orderDF[10,4],
#             orderDF[11,4],
#             orderDF[12,4],
#             # orderDF[13,4],
#             # orderDF[14,4],
#             RSSA[4,2],
#             orderDF[3,5],
#             orderDF[4,5],
#             orderDF[5,5],
#             # orderDF[6,5],
#             # orderDF[7,5],
#             orderDF[10,5],
#             orderDF[11,5],
#             orderDF[12,5],
#             # orderDF[13,5],
#             # orderDF[14,5],
#             RSSA[5,2],
#             orderDF[3,6],
#             orderDF[4,6],
#             orderDF[5,6],
#             # orderDF[6,6],
#             # orderDF[7,6],
#             orderDF[10,6],
#             orderDF[11,6],
#             orderDF[12,6],
#             # orderDF[13,6],
#             # orderDF[14,6],
#             RSSA[6,2],
#             orderDF[3,7],
#             orderDF[4,7],
#             orderDF[5,7],
#             # orderDF[6,7],
#             # orderDF[7,7],
#             orderDF[10,7],
#             orderDF[11,7],
#             orderDF[12,7],
#             # orderDF[13,7],
#             # orderDF[14,7],
#             RSSA[7,2],
#             orderDF[3,8],
#             orderDF[4,8],
#             orderDF[5,8],
#             # orderDF[6,8],
#             # orderDF[7,8],
#             orderDF[10,8],
#             orderDF[11,8],
#             orderDF[12,8],
#             # orderDF[13,8],
#             # orderDF[14,8],
#             RSSA[8,2],
#             orderDF[3,9],
#             orderDF[4,9],
#             orderDF[5,9],
#             # orderDF[6,9],
#             # orderDF[7,9],
#             orderDF[10,9],
#             orderDF[11,9],
#             orderDF[12,9],
#             # orderDF[13,9],
#             # orderDF[14,9],
#             RSSA[9,2],
#             orderDF[3,10],
#             orderDF[4,10],
#             orderDF[5,10],
#             # orderDF[6,10],
#             # orderDF[7,10],
#             orderDF[10,10],
#             orderDF[11,10],
#             orderDF[12,10],
#             # orderDF[13,10],
#             # orderDF[14,10],
#             RSSA[10,2],
#             orderDF[3,11],
#             orderDF[4,11],
#             orderDF[5,11],
#             # orderDF[6,11],
#             # orderDF[7,11],
#             orderDF[10,11],
#             orderDF[11,11],
#             orderDF[12,11],
#             # orderDF[13,11],
#             # orderDF[14,11],
#             RSSA[11,2],
#             orderDF[3,12],
#             orderDF[4,12],
#             orderDF[5,12],
#             # orderDF[6,12],
#             # orderDF[7,12],
#             orderDF[10,12],
#             orderDF[11,12],
#             orderDF[12,12],
#             # orderDF[13,12],
#             # orderDF[14,12],
#             RSSA[12,2],
#             orderDF[3,13],
#             orderDF[4,13],
#             orderDF[5,13],
#             # orderDF[6,13],
#             # orderDF[7,13],
#             orderDF[10,13],
#             orderDF[11,13],
#             orderDF[12,13],
#             # orderDF[13,13],
#             # orderDF[14,13],
#             RSSA[13,2],
#             estimate_RSSA,
#             estimate_RSSA_perc
#           ))
#         }else{
#           ensembleTable[count,] <- as.vector(c(
#             count,
#             widthMCvars[i,2],
#             widthMCvars[i,3],
#             orderDF[3,1],
#             orderDF[4,1],
#             orderDF[5,1],
#             # orderDF[6,1],
#             # orderDF[7,1],
#             RSSA[1,2],
#             orderDF[3,2],
#             orderDF[4,2],
#             orderDF[5,2],
#             # orderDF[6,2],
#             # orderDF[7,2],
#             RSSA[2,2],
#             orderDF[3,3],
#             orderDF[4,3],
#             orderDF[5,3],
#             # orderDF[6,3],
#             # orderDF[7,3],
#             RSSA[3,2],
#             orderDF[3,4],
#             orderDF[4,4],
#             orderDF[5,4],
#             # orderDF[6,4],
#             # orderDF[7,4],
#             RSSA[4,2],
#             orderDF[3,5],
#             orderDF[4,5],
#             orderDF[5,5],
#             # orderDF[6,5],
#             # orderDF[7,5],
#             RSSA[5,2],
#             orderDF[3,6],
#             orderDF[4,6],
#             orderDF[5,6],
#             # orderDF[6,6],
#             # orderDF[7,6],
#             RSSA[6,2],
#             orderDF[3,7],
#             orderDF[4,7],
#             orderDF[5,7],
#             # orderDF[6,7],
#             # orderDF[7,7],
#             RSSA[7,2],
#             orderDF[3,8],
#             orderDF[4,8],
#             orderDF[5,8],
#             # orderDF[6,8],
#             # orderDF[7,8],
#             RSSA[8,2],
#             orderDF[3,9],
#             orderDF[4,9],
#             orderDF[5,9],
#             # orderDF[6,9],
#             # orderDF[7,9],
#             RSSA[9,2],
#             orderDF[3,10],
#             orderDF[4,10],
#             orderDF[5,10],
#             # orderDF[6,10],
#             # orderDF[7,10],
#             RSSA[10,2],
#             orderDF[3,11],
#             orderDF[4,11],
#             orderDF[5,11],
#             # orderDF[6,11],
#             # orderDF[7,11],
#             RSSA[11,2],
#             orderDF[3,12],
#             orderDF[4,12],
#             orderDF[5,12],
#             # orderDF[6,12],
#             # orderDF[7,12],
#             RSSA[12,2],
#             orderDF[3,13],
#             orderDF[4,13],
#             orderDF[5,13],
#             # orderDF[6,13],
#             # orderDF[7,13],
#             RSSA[13,2],
#             estimate_RSSA,
#             estimate_RSSA_perc
#           ))
#         }
#         new <- Sys.time() - old
#         print(new)
#         count <- count + 1
#       }
#     }
#   }
#   return(ensembleTable)
# }

###Variables

#Color palettes
colors <- brewer.pal(12, "Paired")
colors[13] <- "#5C4033"
binwidth <- 0.05


###Analysis

#Create a sorted list of all stream orders and basins in the dataset to iterate over
order_list <- sort(unique(all_widths_sampled[,2]))
basin_list <- c("SSVC", "SVC", "Platte", "Mississippi")

#Find mean width by stream order
avg_widths <- avg_width_by_order(width_data = all_widths_sampled,
                                 orders = order_list)

#Run the mle_by_basin function on width data
basinDF <- mle_by_basin(width_data = all_widths,
                                       basins = basin_list)
write.csv(basinDF, "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_04_26_MC_Analysis/Output/statDF_basin.csv")

#Find length of streams by order using Horton scaling principles
#Length (m) of 13th order rivers from NHDPlus HR
mississippi_13_length <- 1878054

#List of 13th order river length
rivers_13 <- c(mississippi_13_length)
#Average length of 13th order rivers
avg_13_length <- mean(rivers_13)

#Number of 13th order rivers
num_13 <- length(rivers_13)

#Length (m) of 12th order rivers from NHDPlus HR
missouri_12_length <- 2531746
arkansas_12_length <- 599298
ohio_12_length <- 210251
mississippi_12_length <- 35166

#List of 12th order river length
rivers_12 <- c(missouri_12_length, arkansas_12_length, ohio_12_length, mississippi_12_length)
#Average length of 12th order rivers
avg_12_length <- mean(rivers_12)

#Number of 12th order rivers
num_12 <- length(rivers_12)

#Find the number ratio and length ratio based on 12th and 13th order streams in the Mississippi Basin
number_ratio <- num_12/num_13
length_ratio <- avg_13_length/avg_12_length

#Variables "b" and "d" can be approximated from the number and length ratios
b <- 1/number_ratio
d <- length_ratio
#Variables "a" and "c" can be calculated using the actual length and number of streams, here using 13th order streams
a <- num_13/(number_ratio^(-13))
c <- avg_13_length/(length_ratio^(13))

#Run function to find average lengths and counts for each stream order
lengthDF <- find_length_by_order(order_list, a, b, c, d)




#non MC run
nonMCtable <- RSSA_estimation_MC(w = all_widths_sampled,
                                 order_list = order_list,
                                 lengthDF = lengthDF,
                                 nRun = 1,
                                 pareto = F)

#MC run
MCtable <- RSSA_estimation_MC(w = all_widths_sampled,
                              order_list = order_list,
                              lengthDF = lengthDF,
                              nRun = 1000,
                              pareto = F)
MCtable_RSSA_mean <- mean(MCtable$RSSA_km2)
MCtable_RSSA_sd <- sd(MCtable$RSSA_km2)
# 
# #non MC nested run
# nonMCnested_table <- RSSA_estimation_MC_nested(w = all_widths_sampled,
#                                                order_list = order_list,
#                                                lengthDF = lengthDF,
#                                                nRun = 1,
#                                                pareto = F)
# 
# #MC nested run
# MCnested_table <- RSSA_estimation_MC_nested(w = all_widths_sampled,
#                                                order_list = order_list,
#                                                lengthDF = lengthDF,
#                                                nRun = 10,
#                                                pareto = F)
# MCnested_table_RSSA_mean <- mean(MCnested_table$RSSA_km2)
# MCnested_table_RSSA_sd <- sd(MCnested_table$RSSA_km2)

# #non MC run (2 sd)
# nonMCtable_2sd <- RSSA_estimation_MC(w = all_widths_sampled,
#                                  order_list = order_list,
#                                  lengthDF = lengthDF,
#                                  nRun = 1,
#                                  pareto = F,
#                                  numSD = 2)
# 
# #MC run (2 sd)
# MCtable_2sd <- RSSA_estimation_MC(w = all_widths_sampled,
#                               order_list = order_list,
#                               lengthDF = lengthDF,
#                               nRun = 1000,
#                               pareto = F,
#                               numSD = 2)
# MCtable_RSSA_mean_2sd <- mean(MCtable_2sd$RSSA_km2)
# MCtable_RSSA_sd_2sd <- sd(MCtable_2sd$RSSA_km2)

# #non MC nested run (2 sd)
# nonMCnested_table_2sd <- RSSA_estimation_MC_nested(w = all_widths_sampled,
#                                                order_list = order_list,
#                                                lengthDF = lengthDF,
#                                                nRun = 1,
#                                                pareto = F,
#                                                numSD = 2)
# 
# #MC nested run (2 sd)
# MCnested_tabl_2sd <- RSSA_estimation_MC_nested(w = all_widths_sampled,
#                                             order_list = order_list,
#                                             lengthDF = lengthDF,
#                                             nRun = 10,
#                                             pareto = F,
#                                             numSD = 2)
# MCnested_table_RSSA_mean_2sd <- mean(MCnested_table_2sd$RSSA_km2)
# MCnested_table_RSSA_sd_2sd <- sd(MCnested_table_2sd$RSSA_km2)

# #non MC nested V2 run
# nonMCnested_V2_table <- RSSA_estimation_MC_nested_V2(w = all_widths_sampled,
#                                                order_list = order_list,
#                                                lengthDF = lengthDF,
#                                                nRun = 1,
#                                                pareto = F)
# 
# #MC nested V2 run
# MCnested_V2_table <- RSSA_estimation_MC_nested_V2(w = all_widths_sampled,
#                                             order_list = order_list,
#                                             lengthDF = lengthDF,
#                                             nRun = 10,
#                                             pareto = F)
# MCnested_V2_table_RSSA_mean <- mean(MCnested_V2_table$RSSA_km2)
# MCnested_V2_table_RSSA_sd <- sd(MCnested_V2_table$RSSA_km2)

#Find synthetic widths
# synth_widths <- generate_synth_datasets(orderDF, lengthDF, order_list)
# write.csv(synth_widths, "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_04_26_MC_Analysis/Output/synthetic_widths.csv")


#Allen and Pavelsky, 2018 Estimate of RSA
#Extrapolate the Pareto fit found in the Mississippi Basin by Allen and Pavelsky, 2018 to the total length of streams
grwl_mississippi_widths <- read.csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_04_26_MC_Analysis/Data/grwl_centerlines.csv")
grwl_mississippi_widths <- dplyr::select(grwl_mississippi_widths, width_m)
colnames(grwl_mississippi_widths) <- c("Width")
grwl_mississippi_widths <- filter(grwl_mississippi_widths, Width >= 90)
grwl_mississippi_widths <- dplyr::mutate(grwl_mississippi_widths, Surf_Area_km2 = Width*30*0.000001)
grwlObs_area <- sum(grwl_mississippi_widths$Surf_Area_km2)

firstOrder_width_mean <- 0.321
firstOrder_width_sd <- 0.077
firstOrder_width_range <- c(firstOrder_width_mean-firstOrder_width_sd, firstOrder_width_mean, firstOrder_width_mean+firstOrder_width_sd)
alpha <- 0.96
minObs_width <- 90
maxObs_width <- max(grwl_mississippi_widths$Width)
extrapIntegral <- -(alpha*minObs_width*(firstOrder_width_mean/minObs_width)^alpha)/(alpha-1) + (alpha*firstOrder_width_range*(firstOrder_width_mean/firstOrder_width_range)^alpha)/(alpha-1)
obsIntegral <- -(alpha*maxObs_width*(firstOrder_width_mean/maxObs_width)^alpha)/(alpha-1) + (alpha*firstOrder_width_range*(firstOrder_width_mean/firstOrder_width_range)^alpha)/(alpha-1)
extrap2ObsRatio <- extrapIntegral/obsIntegral
grwlExtrap_area <- grwlObs_area*extrap2ObsRatio+grwlObs_area
grwlExtrap_area_perc <- percent_river(grwlExtrap_area)


###Plotting
#Figure showing bar chart of RSSA by stream order
barchartDF <- data.frame(Order = c(1:13),
                         RSSA = c(nonMCtable[1,11],
                                  nonMCtable[1,15],
                                  nonMCtable[1,19],
                                  nonMCtable[1,23],
                                  nonMCtable[1,27],
                                  nonMCtable[1,31],
                                  nonMCtable[1,35],
                                  nonMCtable[1,39],
                                  nonMCtable[1,43],
                                  nonMCtable[1,47],
                                  nonMCtable[1,51],
                                  nonMCtable[1,55],
                                  nonMCtable[1,59]),
                         RSSA_sd = c(sd(MCtable[,11]),
                                      sd(MCtable[,15]),
                                      sd(MCtable[,19]),
                                      sd(MCtable[,23]),
                                      sd(MCtable[,27]),
                                      sd(MCtable[,31]),
                                      sd(MCtable[,35]),
                                      sd(MCtable[,39]),
                                      sd(MCtable[,43]),
                                      sd(MCtable[,47]),
                                      sd(MCtable[,51]),
                                      sd(MCtable[,55]),
                                      sd(MCtable[,59])),
                         RSSA_2sd = c(2*sd(MCtable[,11]),
                                     2*sd(MCtable[,15]),
                                     2*sd(MCtable[,19]),
                                     2*sd(MCtable[,23]),
                                     2*sd(MCtable[,27]),
                                     2*sd(MCtable[,31]),
                                     2*sd(MCtable[,35]),
                                     2*sd(MCtable[,39]),
                                     2*sd(MCtable[,43]),
                                     2*sd(MCtable[,47]),
                                     2*sd(MCtable[,51]),
                                     2*sd(MCtable[,55]),
                                     2*sd(MCtable[,59])
                           
                         ))
ggplot(data = barchartDF, aes(x = as.factor(Order), y = RSSA))+
  geom_bar(aes(fill = as.factor(Order)),alpha = 1, stat = "identity")+
  geom_errorbar(data = barchartDF, aes(ymin = RSSA-RSSA_2sd, ymax = RSSA+RSSA_2sd),
                color = "black", width  = 0.4, linewidth = 1.3)+
  scale_fill_manual(values = colors)+
  guides(fill = guide_legend(title = "Stream Order", ncol = 2))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),
                   breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000),
                   labels = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000))+
  labs(x = "Stream Order",
       y = "River Surface Area (km2)")+
  theme_classic()
  
#(Fig 3 relationship) Plot Average Width vs Standard Deviation for each stream order
avgWidth_sd <- data.frame(Order = c(1:13), Avg_Width_m = NA, SD = NA)
for(i in order_list){
  tempWidths <- filter(all_widths_sampled, Order == i)
  avgWidth <- mean(tempWidths$Width)
  sd <- sd(tempWidths$Width)
  avgWidth_sd[i,2] <- avgWidth
  avgWidth_sd[i,3] <- sd
}
avgWidth_sd_lm <- lm(SD~Avg_Width_m, avgWidth_sd)
linear_fit_function <- function(x){
  coef(avgWidth_sd_lm)[["Avg_Width_m"]]*x+coef(avgWidth_sd_lm)[["(Intercept)"]]
}
linear_fitDF <- data.frame(x = seq(0.001, 1000))
linear_fitDF$y <- linear_fit_function(linear_fitDF$x)
avgWidth_vs_StDv <- ggplot(avgWidth_sd, aes(x = Avg_Width_m, y = SD, color = as.factor(Order)))+
  geom_point(size = 3,
             alpha = 1,
             show.legend = T)+
  geom_abline(slope = coef(avgWidth_sd_lm)[["Avg_Width_m"]],
              intercept = coef(avgWidth_sd_lm)[["(Intercept)"]],
              linewidth = 1)+
  # geom_line(data = linear_fitDF, aes(x=x,y=y), color = "black", linewidth = 1)+
  scale_color_manual(values = colors)+
  annotate("text", label = bquote(y == .(round(coef(avgWidth_sd_lm)[["Avg_Width_m"]], digits = 2))*x + .(round(coef(avgWidth_sd_lm)[["(Intercept)"]], digits = 2))~", "~R^2 == .(round(summary(avgWidth_sd_lm)$r.squared, digits = 2))),
           x = 550, y = 65,
           family = "sans",
           size = 8)+
  # annotate("text", label = paste0("y = ", round(coef(avgWidth_sd_lm)[["Avg_Width_m"]], digits = 4), " x + ", round(coef(avgWidth_sd_lm)[["(Intercept)"]], digits = 2), ", ", bquote(R^2~" = "), round(summary(avgWidth_sd_lm)$r.squared, digits = 2)),
  #          x = 550, y = 50,
  #          family = "sans",
  #          size = 8)+
  guides(color = guide_legend(title = "Stream Order", ncol = 2))+
  xlab("Mean Width (m)")+
  ylab("Standard Deviation")+
  # scale_x_log10()+
  # scale_y_log10()+
  # scale_x_continuous(expand = c(0,0))+
  # scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 28))

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/avgWidth_vs_StDv.svg",
       plot = avgWidth_vs_StDv,
       width = 9.10,
       height = 6.02,
       units = "in")

#(Figure 3) Function to plot density plots and fitted distributions for measured streams of each order
plot_density_and_dist_order <- function(widths, statDF, order, binwidth = 1, ypos1 = 1, ypos2 = 1, ypos3 = 1, xlim = 2000){
  measured_streams <- filter(widths, Order == order)
  # # Clip top 5% of measurements
  # measured_streams <- subset(measured_streams, Width <= quantile(Width, prob = 1-(2/100)))
  # print('Best fit: lognormal')
  lnmean <- round(statDF[1,order],3)
  lnsd <- round(statDF[2,order],3)
  lnwass <- round(statDF[3, order],3)
  temp1 <- paste0("\u03bc", " = ", lnmean)
  temp2 <- paste0("\u03c3", " = ", lnsd)
  temp3 <- paste0("W = ", lnwass)
  ggplot(measured_streams, aes(x = Width))+
    geom_histogram(aes(y=after_stat(density)),
                   fill = colors[order], color = colors[order],
                   binwidth = binwidth,
    )+
    # geom_density(fill = colors[order], color = colors[order],
    #   position = "identity",adjust = 2
    #   , alpha = 0.6)+
    stat_function(fun = dlnorm, args = list(meanlog = statDF[1,order], sdlog = statDF[2,order]),
                  color = "black", linewidth = 1.25)+
    xlab("Width (m)")+
    ylab("Density")+
    annotate("text",
             x = 0.7*xlim, y = ypos1,
             family= "sans",
             size = 10,
             label = temp1
             , parse = F)+
    annotate("text",
             x = 0.7*xlim, y = ypos2,
             family= "sans",
             size = 10,
             label = temp2
             , parse = F)+
    annotate("text",
             x = 0.7*xlim, y = ypos3,
             family= "sans",
             size = 10,
             label = temp3
             , parse = F)+
    scale_x_continuous(expand = c(0,0), limits = c(0, xlim))+
    scale_y_continuous(expand = c(0,0))+
    # ggtitle(paste0("Density Plot and Fitted Log-normal Distribution for Stream Order ", order))+
    theme_classic()+
    theme(text = element_text(family = "sans", size = 28),
          plot.margin=unit(c(.2,.7,.2,.2),"cm"),
          axis.text=element_text(size=18),
          axis.title=element_text(size=28))
}

#Run the mle_by_order function on width data
# orderDF <- mle_by_order(width_data = all_widths_sampled,
#                         orders = order_list,
#                         parOption = F)

order1_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 1, 0.1524,
                                           ypos1 = 1.05,
                                           ypos2 = 0.9,
                                           ypos3 = 0.75,
                              xlim = 2)
order2_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 2, 0.175,
                                           ypos1 = 0.635,
                                           ypos2 = 0.55,
                                           ypos3 = 0.465,
                              xlim = 4)
order3_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 3, 0.27,
                                           ypos1 = 0.34,
                                           ypos2 = 0.29,
                                           ypos3 = 0.24,
                              xlim = 7)
order4_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 4, 1.27,
                                           ypos1 = 0.14,
                                           ypos2 = 0.118,
                                           ypos3 = 0.096,
                              xlim = 19)
order5_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 5, 1.499,
                                           ypos1 = 0.123,
                                           ypos2 = 0.104,
                                           ypos3 = 0.085,
                              xlim = 25)
order6_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 6, 2.5,
                                           ypos1 = 0.068,
                                           ypos2 = 0.059,
                                           ypos3 = 0.05,
                              xlim = 30)
order7_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 7, 3.0,
                                           ypos1 = 0.058,
                                           ypos2 = 0.050,
                                           ypos3 = 0.042,
                              xlim = 55)
order8_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 8, 14.6,
                                           ypos1 = 0.0126,
                                           ypos2 = 0.0108,
                                           ypos3 = 0.009,
                              xlim = 250)
order9_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 9, 26.3,
                                           ypos1 = 0.009,
                                           ypos2 = 0.0078,
                                           ypos3 = 0.0066,
                              xlim = 300)
order10_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 10, 12,
                                            ypos1 = 0.0099,
                                            ypos2 = 0.0085,
                                            ypos3 = 0.0071,
                              xlim = 300)
order11_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 11, 40,
                                            ypos1 = 0.00123,
                                            ypos2 = 0.00106,
                                            ypos3 = 0.00089,
                              xlim = 2050)
order12_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 12, 94,
                                            ypos1 = 0.00135,
                                            ypos2 = 0.00116,
                                            ypos3 = 0.00097,
                              xlim = 2050)
order13_hist <-
  plot_density_and_dist_order(all_widths_sampled, nonMCtable, 13, 197,
                                            ypos1 = 0.00087,
                                            ypos2 = 0.00075,
                                            ypos3 = 0.00063,
                              xlim = 3050)

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order1_hist.svg",
       plot = order1_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order2_hist.svg",
       plot = order2_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order3_hist.svg",
       plot = order3_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order4_hist.svg",
       plot = order4_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order5_hist.svg",
       plot = order5_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order6_hist.svg",
       plot = order6_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order7_hist.svg",
       plot = order7_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order8_hist.svg",
       plot = order8_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order9_hist.svg",
       plot = order9_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order10_hist.svg",
       plot = order10_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order11_hist.svg",
       plot = order11_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order12_hist.svg",
       plot = order12_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/Fig 3 Fits by Order/order13_hist.svg",
       plot = order13_hist,
       width = 9.10,
       height = 6.02,
       units = "in")

#(Figure 2, e-h) Function to plot density plots and fitted distributions for measured streams in each basins
plot_density_and_dist_basin <- function(widths, basin, xmax = 4003, binwidth = 1){
  measured_streams <- filter(widths, Basin == basin)
  ggplot(measured_streams, aes(x = Width))+
    geom_histogram(aes(y=after_stat(density)),
                   fill = "grey50", color = "grey30",
                   binwidth = binwidth)+
    # stat_function(fun = pareto_pdf, args = list(lambda = statDF[10,index], k = statDF[11,index]),
    #               color = "black", linewidth = 1.25)+
    xlab("Width (m)")+
    ylab("Density")+
    # xlim(c(0,500))+
    # ylim(c(0, 0.025))+
    # ggtitle(paste0("Density Plot and Fitted Pareto Distribution for ", basin, " Basin"))+
    scale_x_continuous(expand = c(0,0), limits = c(0,xmax))+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    theme(text = element_text(family = "sans", size = 22),
          plot.margin=unit(c(.2,.7,.2,.2),"cm"))
}
mississippi_hist <- plot_density_and_dist_basin(all_widths, "Mississippi", xmax = 1000, binwidth = 15.81)
platte_hist <- plot_density_and_dist_basin(all_widths, "Platte", xmax = 500, binwidth = 6.73)
svc_hist <- plot_density_and_dist_basin(all_widths, "SVC", xmax = 100, binwidth = 2.5)
ssvc_hist <- plot_density_and_dist_basin(all_widths, "SSVC", xmax = 15, binwidth = 0.15)

ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 1 Width Map and Raw Data/mississippi_hist.svg",
       plot = mississippi_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 1 Width Map and Raw Data/platte_hist.svg",
       plot = platte_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 1 Width Map and Raw Data/svc_hist.svg",
       plot = svc_hist,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 1 Width Map and Raw Data/ssvc_hist.svg",
       plot = ssvc_hist,
       width = 9.10,
       height = 6.02,
       units = "in")

#(Figure 4) Plot a stacked histogram showing all synthetic widths by order
#large scale
big_stack <- ggplot(synth_widths, aes(x = Width, fill = as.factor(Order)))+
  geom_histogram(position = "stack"
                 ,binwidth = binwidth
                 # ,bins = 100
                 ,show.legend = F
  )+
  scale_fill_manual(values = colors)+
  #GRWL Pareto Fit with no extrapolation
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = 90, k = 0.96, bw = binwidth, n = nrow(grwl_mississippi_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "black", linewidth = 1.25, xlim = c(log10(90),log10(20000)), show.legend = F)+
  #GRWL Pareto Fit with extrapolation
  stat_function(fun = function(x){pareto_pdf_log(x,lambda = 0.321, k = 0.96, bw = binwidth, n = nrow(synth_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
                color = "black", linewidth = 1.25, xlim = c(log10(0.321),log10(20000)), show.legend = F)+
  #Pareto Fit to synthetic data
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = synthDF[10,2], k = synthDF[11,2], bw = binwidth, n = npareto) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "#840032", linewidth = 1.25, xlim = c(log10(0.321),log10(20000)), show.legend = F)+
  # stat_function(fun = function(x){dlnorm(x,meanlog = synthDF[5,2], sdlog = synthDF[6,2])*binwidth*nrow(synth_widths) * x * log(10)},
  #               color = "#002642", linewidth = 1.25, show.legend = F)+
  scale_x_continuous(breaks = c(0,0.01,0.1,1,10,100,1000,10000)
                     , trans = "log10"
                     , expand = c(0,0)
                     , labels = label_log())+
  scale_y_continuous(breaks = c(0,2000000,4000000,6000000,8000000,10000000)
                     # breaks = c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
                     # , trans = "log10"
                     ,expand = c(0,0), limits = c(0,10000000)
                     # ,labels = label_log()
  )+
  # geom_vline(xintercept = 20, color = "red", linewidth = 1.25)+
  
  # scale_x_log10(expand = c(0,0))+
  # scale_y_log10(expand = c(0,0))+
  # coord_trans(x = "log10")+
  xlab("Width (m)")+
  ylab("Count")+
  # guides(fill = guide_legend(title = "Stream Order", ncol = 2))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 20))

# synth_widths_bigOrder <- filter(synth_widths, Order >= 8)
# colors_bigOrder <- colors[8:13]

#medium scale
medium_stack <- ggplot(synth_widths, aes(x = Width, fill = as.factor(Order)))+
  geom_histogram(position = "stack"
                 ,binwidth = binwidth
                 # ,bins = 100
                 ,show.legend = F
  )+
  scale_fill_manual(values = colors)+
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = 90, k = 0.96, bw = binwidth, n = nrow(grwl_mississippi_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "black", linewidth = 1.25, xlim = c(log10(90),log10(20000)), show.legend = F)+
  #GRWL Pareto Fit with extrapolation
  stat_function(fun = function(x){pareto_pdf_log(x,lambda = 0.321, k = 0.96, bw = binwidth, n = nrow(synth_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
                color = "black", linewidth = 1.25, xlim = c(log10(0.321),log10(20000)), show.legend = F)+
  #Pareto Fit to synthetic data
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = synthDF[10,2], k = synthDF[11,2], bw = binwidth, n = npareto) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "#840032", linewidth = 1.25, xlim = c(log10(0.321),log10(20000)), show.legend = F)+
  # stat_function(fun = function(x){dlnorm(x,meanlog = synthDF[5,2], sdlog = synthDF[6,2])*binwidth*nrow(synth_widths) * x * log(10)},
  #               color = "#002642", linewidth = 1.25, show.legend = F)+
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = getmode(synth_widths$Width), k = synthDF[11,2], bw = binwidth, n = nrow(synth_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "black", linewidth = 1.25, xlim = c(log10(0.32),log10(6000)), show.legend = F)+
  scale_x_continuous(breaks = c(0,0.01,0.1,1,10,100,1000,10000)
                     , trans = "log10"
                     , expand = c(0,0)
                     #, limits = c(10, 6000)
                     , labels = label_log())+
  scale_y_continuous(breaks = c(0,400000,800000,1200000,1600000,2000000)
                     # breaks = c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
                     # , trans = "log10"
                     ,expand = c(0,0), #limits = c(0,150000)
  )+
  # geom_vline(xintercept = 20, color = "red", linewidth = 1.25)+
  
  # scale_x_log10(expand = c(0,0))+
  # scale_y_log10(expand = c(0,0))+
  # coord_trans(x = "log10")+
  coord_cartesian(xlim = c(1,2500), ylim = c(0,2000000))+
  xlab("Width (m)")+
  ylab("Count")+
  # guides(fill = guide_legend(title = "Stream Order", ncol = 2))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 20))

#small scale
small_stack <- ggplot(synth_widths, aes(x = Width, fill = as.factor(Order)))+
  geom_histogram(position = "stack"
                 ,binwidth = binwidth
                 # ,bins = 100
                 ,show.legend = F
  )+
  scale_fill_manual(values = colors)+
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = 90, k = 0.96, bw = binwidth, n = nrow(grwl_mississippi_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "black", linewidth = 1.25, xlim = c(log10(90),log10(20000)), show.legend = F)+
  #GRWL Pareto Fit with extrapolation
  stat_function(fun = function(x){pareto_pdf_log(x,lambda = 0.321, k = 0.96, bw = binwidth, n = nrow(synth_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
                color = "black", linewidth = 1.25, xlim = c(log10(0.321),log10(20000)), show.legend = F)+
  #Pareto Fit to synthetic data
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = synthDF[10,2], k = synthDF[11,2], bw = binwidth, n = npareto) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "#840032", linewidth = 1.25, xlim = c(log10(0.321),log10(20000)), show.legend = F)+
  # stat_function(fun = function(x){dlnorm(x,meanlog = synthDF[5,2], sdlog = synthDF[6,2])*binwidth*nrow(synth_widths) * x * log(10)},
  #               color = "#002642", linewidth = 1.25, show.legend = F)+
  # stat_function(fun = function(x){pareto_pdf_log(x,lambda = getmode(synth_widths$Width), k = synthDF[11,2], bw = binwidth, n = nrow(synth_widths)) * x * log(10)}, #args = list(lambda = 90, k = 0.96, bw = binwidth, n = nrow(synth_widths_sampled_lessSO_full)),
  #               color = "black", linewidth = 1.25, xlim = c(log10(0.32),log10(6000)), show.legend = F)+
  scale_x_continuous(breaks = c(0,0.01,0.1,1,10,100,1000,10000)
                     , trans = "log10"
                     , expand = c(0,0)
                     #, limits = c(10, 6000)
                     , labels = label_log())+
  scale_y_continuous(breaks = c(0,30000,60000,90000,120000,150000)
                     # breaks = c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
                     # , trans = "log10"
                     ,expand = c(0,0), #limits = c(0,150000)
  )+
  # scale_x_log10(expand = c(0,0))+
  # scale_y_log10(expand = c(0,0))+
  # coord_trans(x = "log10")+
  coord_cartesian(xlim = c(5,20000), ylim = c(0,100000))+
  xlab("Width (m)")+
  ylab("Count")+
  guides(fill = guide_legend(title = "Stream Order", ncol = 2))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 20),
        legend.position = c(0.7,0.7))

ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 4 Synth Hist/big_stack.svg",
       plot = big_stack,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 4 Synth Hist/medium_stack.svg",
       plot = medium_stack,
       width = 9.10,
       height = 6.02,
       units = "in")
ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig 4 Synth Hist/small_stack.pdf",
       plot = small_stack,
       width = 9.10,
       height = 6.02,
       units = "in")


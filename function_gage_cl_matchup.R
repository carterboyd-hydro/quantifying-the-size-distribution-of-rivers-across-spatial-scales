#Author: Carter Boyd
#Date: 11/08/2023
#Purpose: Validate RivWidth measurements at GRWL, S2, and NAIP scales by comparing to USGS reported widths at gages

library(tidyverse)
library(fields)
library(NISTunits)
library(ggpmisc)
library(mblm)
library(geomtextpath)
library(scales)
library(rlist)
library(Metrics)
library(MASS)
library(moments)


# Define a function that calculates the geodesic distance between two points specified by radian latitude/longitude using the Haversine formula (hf)
gcd.hf <- function(lon1, lat1, lon2, lat2) {
  lon1 <- NISTdegTOradian(lon1)
  lat1 <- NISTdegTOradian(lat1)
  lon2 <- NISTdegTOradian(lon2)
  lat2 <- NISTdegTOradian(lat2)
  R <- 6371 # Earth mean radius [km]
  delta.lon <- (lon2 - lon1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.lon/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

#Define a function that returns the gage dataframe with new 'measured_width_m' and 'width_residual_m' attributes for each gage
match_widths = function(gages, centerlines, number_of_points){
  gages$measured_width_m <- NA
  gages$width_residual_m <- NA
  for (i in 1:nrow(gages)){
    print(i)
    gage_lon <- gages[i,]$lon
    gage_lat <- gages[i,]$lat
    close_centerlines <- filter(centerlines, (lon <= (gage_lon + 0.01)) & (lon >= (gage_lon - 0.01)) & (lat <= (gage_lat + 0.01)) & (lat >= (gage_lat - 0.01)))
    close_centerlines$dist_to_gage <- NA
    for (j in 1:nrow(close_centerlines)){
      cl_lon <- close_centerlines[j,]$lon
      cl_lat <- close_centerlines[j,]$lat
      #if ((cl_lon <= (gage_lon + 0.01)&
      #     (cl_lon >= gage_lon - 0.01))&
      #    (cl_lat <= (gage_lat + 0.01)&
      #    cl_lat >= (gage_lat - 0.01))){
      #  dist <- gcd.hf(gage_lon, gage_lat, cl_lon, cl_lat)
      #  centerlines[j,]$dist_to_gage <- dist
      #  print(dist)
      #}
      dist <- gcd.hf(gage_lon, gage_lat, cl_lon, cl_lat)
      close_centerlines[j,]$dist_to_gage <- dist
    }
    close_centerlines <- mutate(close_centerlines, rank = rank(dist_to_gage)) %>% 
      arrange(rank)
    #print(close_centerlines)

    mean_measured_width <- mean(close_centerlines[1:number_of_points,]$width)
    print(mean_measured_width)
    gages[i,]$measured_width_m <- mean_measured_width
    gages[i,]$width_residual_m <- gages[i,]$measured_width_m - gages[i,]$mean_width
  }
  return(gages)
}

#Define a function that returns the gage dataframe with new 'measured_width_m' and 'width_residual_m' attributes for each gage for the NAIP corrected widths
match_widths_corrected = function(gages, centerlines, number_of_points){
  gages$measured_width_m <- NA
  gages$width_residual_m <- NA
  for (i in 1:nrow(gages)){
    print(i)
    gage_lon <- gages[i,]$lon
    gage_lat <- gages[i,]$lat
    close_centerlines <- filter(centerlines, (lon <= (gage_lon + 0.01)) & (lon >= (gage_lon - 0.01)) & (lat <= (gage_lat + 0.01)) & (lat >= (gage_lat - 0.01)))
    close_centerlines$dist_to_gage <- NA
    for (j in 1:nrow(close_centerlines)){
      cl_lon <- close_centerlines[j,]$lon
      cl_lat <- close_centerlines[j,]$lat
      #if ((cl_lon <= (gage_lon + 0.01)&
      #     (cl_lon >= gage_lon - 0.01))&
      #    (cl_lat <= (gage_lat + 0.01)&
      #    cl_lat >= (gage_lat - 0.01))){
      #  dist <- gcd.hf(gage_lon, gage_lat, cl_lon, cl_lat)
      #  centerlines[j,]$dist_to_gage <- dist
      #  print(dist)
      #}
      dist <- gcd.hf(gage_lon, gage_lat, cl_lon, cl_lat)
      close_centerlines[j,]$dist_to_gage <- dist
    }
    close_centerlines <- mutate(close_centerlines, rank = rank(dist_to_gage)) %>% 
      arrange(rank)
    #print(close_centerlines)
    
    mean_measured_width <- mean(close_centerlines[1:number_of_points,]$width_corrected)
    print(mean_measured_width)
    gages[i,]$measured_width_m <- mean_measured_width
    gages[i,]$width_residual_m <- gages[i,]$measured_width_m - gages[i,]$mean_width
  }
  return(gages)
}

#Create a function to use Theil-Sen estimation for regression plots
# sen <- function(..., weights = NULL) {
#   mblm::mblm(...)
# }

#Identify gages to remove for each dataset
grwl_bad_gages <- c(3049676, 3049819, 3062500, 3071000, 3082500, 3190400, 3219590, 3260015, 3332555, 3411500, 3413000, 3422500, 3426310, 344878100, 3548500, 3573500, 3589500, 5247500, 5341550, 5344490, 5344500, 5360000, 5400000, 5487500, 5488110, 5533500, 5533600, 6037100, 6038500, 6039200, 6065500, 6224000, 6225000, 6226000, 6453000, 6646600, 6681998, 6684499, 6687499, 6687500, 6690500, 6714360, 6720990, 6766498, 6766499, 6767998, 6767999, 6768025, 6768035, 6790000, 6792500, 6828490, 6850200, 6851500, 6892495, 6920500, 6922450, 7047810, 7059000, 7099973, 7110000, 7130500, 7134100, 7138075, 7148140, 7157940, 7188007, 7191500, 7263450, 7263451, 7290650, 7331600, 7358000, 7359002, 7372200, 394839104570300, 414757087490401, 6818494, 6461150, 4570300, 6465000, 7312700, 6465500, 7177600, 6172000, 6463000, 6796500, 6778000, 6824500, 5411500, 7156900)
s2_bad_gages <- c(6678000, 6680700, 6681998, 6681999, 6682300, 6687499, 6687500, 6690500, 6767998, 6767999, 6768025, 6770500, 6792500, 6794000, 6799445, 6799450, 6799460, 6799500, 6803495, 6803500, 6803513, 6805000, 394839104570300, 6778000)

#Read in gage files
grwl_gages <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/grwl_gages.csv") %>% 
  filter(!(site_no %in% grwl_bad_gages)) %>% 
  rename(lat = dec_lat_va, lon = dec_long_va)
s2_gages <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/s2_gages.csv") %>% 
  filter(!(site_no %in% s2_bad_gages)) %>% 
  rename(lat = dec_lat_va, lon = dec_long_va)
naip_gages <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/naip_gages.csv") %>% 
  rename(lat = dec_lat_va, lon = dec_long_va)

#Read in centerline point files and remove lakes, reservoirs, and canals
grwl_cl_pts <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/grwl_cl_pts.csv") %>% 
  filter(lakeFlag == 0) %>% 
  rename(width = width_m)
s2_cl_pts <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/s2_cl_pts.csv") %>% 
  filter(reservoir == 0 & canal == 0)
naip_cl_pts <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/naip_cl_pts.csv") %>%
  filter(reservoir == 0 & canal ==0)
naip_cl_pts_corrected <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_12_19_Correcting_Dry_NAIP/naip_cl_pts_corrected_exponential.csv") %>% 
  filter(reservoir == 0 & canal == 0)

#Use function to match in situ and measured widths
naip_matched_widths <- match_widths(naip_gages, naip_cl_pts, 10) %>% 
  mutate(scale = 'NAIP') %>% 
  mutate(percent_error = (abs(width_residual_m) / mean_width) * 100) %>% 
  rename(gage_width_m = mean_width) #%>% 
#filter(percent_error <= 100)
naip_matched_widths_corrected <- match_widths_corrected(naip_gages, naip_cl_pts_corrected, 10) %>% 
  mutate(scale = 'NAIP') %>% 
  mutate(percent_error = (abs(width_residual_m) / mean_width) * 100) %>% 
  rename(gage_width_m = mean_width) #%>% 
  #filter(percent_error <= 100)
s2_matched_widths <- match_widths(s2_gages, s2_cl_pts, 10) %>% 
  mutate(scale = 'Sentinel2') %>% 
  mutate(percent_error = (abs(width_residual_m) / mean_width) * 100) %>% 
  rename(gage_width_m = mean_width) #%>% 
  #filter(percent_error <= 100)
grwl_matched_widths <- match_widths(grwl_gages, grwl_cl_pts, 10) %>% 
  mutate(scale = 'GRWL') %>% 
  mutate(percent_error = (abs(width_residual_m) / mean_width) * 100) %>% 
  rename(gage_width_m = mean_width) #%>% 
  #filter(percent_error <= 100)

#Join all gage records for plotting
all_matched_widths <- rbind(grwl_matched_widths, s2_matched_widths, naip_matched_widths_corrected)

# Plot all matched width/discharge data points and three linear regression model fits
# Find linear regressions at NAIP, S2, and GRWL scales
naipCorrelationCoeff <- cor.test(naip_matched_widths$gage_width_m, naip_matched_widths$measured_width_m)
naipLM <- lm(measured_width_m ~ gage_width_m, data = naip_matched_widths)
naipLMcoefficients <- naipLM$coefficients
naipInSitu <- c(gage_width_m = naip_matched_widths$gage_width_m)
naipPredicted <- c(width = predict(naipLM, data.frame(gage_width_m = naip_matched_widths$gage_width_m)))
naipR2 <- round(summary(naipLM)$r.squared, 2)
naipRMSE <- round(sqrt(mean(naipLM$residuals^2)),2)
naipMAE <- round(mae(naip_matched_widths$measured_width_m, predict(naipLM)),2)
naipBias <- round(bias(naipInSitu, naipPredicted),2)
naipSlope <- round(naipLMcoefficients[2], 2)
naipInt <- round(naipLMcoefficients[1], 2)

naipCorrelationCoeff_corr <- cor.test(naip_matched_widths_corrected$gage_width_m, naip_matched_widths_corrected$measured_width_m)
naipLM_corr <- lm(measured_width_m ~ gage_width_m, data = naip_matched_widths_corrected)
naipLMcoefficients_corr <- naipLM_corr$coefficients
naipInSitu_corr <- c(gage_width_m = naip_matched_widths_corrected$gage_width_m)
naipPredicted_corr <- c(width = predict(naipLM_corr, data.frame(gage_width_m = naip_matched_widths_corrected$gage_width_m)))
naipR2_corr <- round(summary(naipLM_corr)$r.squared, 2)
naipRMSE_corr <- round(sqrt(mean(naipLM_corr$residuals^2)),2)
naipMAE_corr <- round(mae(naip_matched_widths_corrected$measured_width_m, predict(naipLM_corr)),2)
naipBias_corr <- round(bias(naipInSitu_corr, naipPredicted_corr),2)
naipSlope_corr <- round(naipLMcoefficients_corr[2], 2)
naipInt_corr <- round(naipLMcoefficients_corr[1], 2)
  
s2CorrelationCoeff <- cor.test(s2_matched_widths$gage_width_m, s2_matched_widths$measured_width_m)
s2LM <- lm(measured_width_m ~ gage_width_m, data = s2_matched_widths)
s2LMcoefficients <- s2LM$coefficients
s2InSitu <- c(gage_width_m = s2_matched_widths$gage_width_m)
s2Predicted <- c(width = predict(s2LM, data.frame(gage_width_m = s2_matched_widths$gage_width_m)))
s2R2 <- round(summary(s2LM)$r.squared, 2)
s2RMSE <- round(sqrt(mean(s2LM$residuals^2)),1)
s2MAE <- round(mae(s2_matched_widths$measured_width_m, predict(s2LM)),1)
s2Bias <- round(bias(s2InSitu, s2Predicted),2)
s2Slope <- round(s2LMcoefficients[2], 2)
s2Int <- round(s2LMcoefficients[1], 2)

grwlCorrelationCoeff <- cor.test(grwl_matched_widths$gage_width_m, grwl_matched_widths$measured_width_m)
grwlLM <- lm(measured_width_m ~ gage_width_m, data = grwl_matched_widths)
grwlLMcoefficients <- grwlLM$coefficients
grwlInSitu <- c(gage_width_m = grwl_matched_widths$gage_width_m)
grwlPredicted <- c(width = predict(grwlLM, data.frame(gage_width_m = grwl_matched_widths$gage_width_m)))
grwlR2 <- round(summary(grwlLM)$r.squared, 2)
grwlRMSE <- round(sqrt(mean(grwlLM$residuals^2)),1)
grwlMAE <- round(mae(grwl_matched_widths$measured_width_m, predict(grwlLM)),1)
grwlBias <- round(bias(grwlInSitu, grwlPredicted),2)
grwlSlope <- round(grwlLMcoefficients[2], 2)
grwlInt <- round(grwlLMcoefficients[1], 2)

# Color palette
hokies <- c("#FFA85B", "#e36414", "#9a031e", "#5f0f40")
primary <- c("#1c3144", "#d00000", "#ffba08", "#52a447")

#Plot validation with uncorrected NAIP data
pdf(file = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Figure S1 Validation Uncorrected/validation_plot_uncorrected.pdf",
    width = 5.6,
    height = 5.6
    # ,units = 'in'
    ,bg = 'white'
    # ,res = 600
)

par(ps = 12, mar = c(4,4,1,1))

# Plot
# windowsFonts(Times = windowsFont("Times New Roman"))
plot(x = grwl_matched_widths$gage_width_m, y = grwl_matched_widths$measured_width_m,
     axes=F, main='', 
     xlab="In Situ Width (m)", ylab="Remote Sensing Measured Width (m)",
     ylim=c(0, 1200), xlim=c(0,1200), 
     col= scales::alpha(primary[1], 0.5), cex.lab=1,
     pch = 16, bty = 'l',
     family = "sans")
points(x = s2_matched_widths$gage_width_m, y = s2_matched_widths$measured_width_m,
       col= scales::alpha(primary[2], 0.6), cex.lab=1,
       pch = 16
)

tx <- seq(0, 1200, length = 7)
axis(1, tx, c(0,200,400,600,800,1000,1200), lwd = 0.7, cex.axis = 1, family = "sans")
ty <- seq(0, 1200, length = 7)
axis(2, ty, c(0,200,400,600,800,1000,1200), lwd = 0.7, cex.axis = 1, family = "sans")

# 1:1 line
lines(x = c(0, 1200), y = c(0,1200), col = "#C1C1C0", type = 'l', lwd = 3.5, lty = 1)

# NAIP Linear Regression
naipLM_yValues <- naipLMcoefficients[1] + naipLMcoefficients[2]*tx
lines(x = tx, y = naipLM_yValues, pch = 0, col = primary[3], lty = 2, lwd = 3.5)

# s2 Linear Regression
s2LM_yValues <- s2LMcoefficients[1] + s2LMcoefficients[2]*tx
lines(x = tx, y = s2LM_yValues, pch = 1, col = primary[2], lty = 1, lwd = 3.5)

# GRWL Linear Regression
grwltx <- seq(0, 1200, length = 7)
grwlLM_yValues <- grwlLMcoefficients[1] + grwlLMcoefficients[2]*grwltx
lines(x = grwltx, y = grwlLM_yValues, pch = 4, col = primary[1], lty = 3, lwd = 3.5)

points(x = naip_matched_widths$gage_width_m, y = naip_matched_widths$measured_width_m,
       col= scales::alpha(primary[3], 1), cex.lab=1,
       pch = 16
)

#Rectangles for legend
# rect(
#   0,
#   850,
#   650,
#   1190,
#   border = 'black', lwd = 2
# )
# rect(
#   550,
#   0,
#   1190,
#   240,
#   border = 'black', lwd = 2
# )

# Create legend text and symbols
lines(x = c(0, 100), y = c(1150,1150), type = 'l', col = primary[1], lty = 3, lwd = 3.5)
lines(x = c(0, 100), y = c(950, 950), type = 'l', col = primary[2], lty = 1, lwd = 3.5)
lines(x = c(500, 600), y = c(190, 190), type = 'l', col = primary[3], lty = 2, lwd = 3.5)
lines(x = c(500, 600), y = c(40, 40), type = 'l', col = "#C1C1C0", lty = 1, lwd = 3.5)

text(x = 100, y = 1200, bquote(paste(bold('Mississippi River:')
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 1145, bquote(paste('y = ', .(grwlSlope), 'x + ', .(grwlInt),'; R'^'2'*' = 0.90'#, .(grwlR2)
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 1090, bquote(paste('Bias = ', .(grwlBias), ' m; ', 'RMSE = ', '43.0', ' m'
)), pos = 4, family = "sans", cex=0.8)

text(x = 100, y = 1000, bquote(paste(bold('Platte River:')
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 945, bquote(paste('y = ', .(s2Slope), 'x - ',.(abs(s2Int)),'; R'^'2'*' = 0.80'#, .(s2R2)
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 890, bquote(paste('Bias = ', .(s2Bias), ' m; ','RMSE = ', .(s2RMSE), ' m'
)), pos = 4, family = "sans", cex=0.8)

text(x = 600, y = 240, bquote(paste(bold('St. Vrain Creek:')
)), pos = 4, family = "sans", cex=0.8)
text(x = 600, y = 185, bquote(paste('y = ', .(naipSlope), 'x - ',.(abs(naipInt)),'; R'^'2'*' = ', .(naipR2)
)), pos = 4, family = "sans", cex=0.8)
text(x = 600, y = 130, bquote(paste('Bias = ', .(naipBias), ' m; ', 'RMSE = ', .(naipRMSE), ' m'
)), pos = 4, family = "sans", cex=0.8)

text(x = 600, y = 35, bquote(paste(bold("Reference Line:"), " y = x"
)), pos = 4, family = "sans", cex=0.8)

dev.off()

#Plot validation with corrected NAIP data
pdf(file = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Figure 2 Validation Corrected/validation_plot.pdf",
    width = 5.6,
    height = 5.6,
    # units = 'in',
    bg = 'white'
    # res = 600
)
par(ps = 12, mar = c(4,4,1,1))

# Plot
# windowsFonts(Times = windowsFont("Times New Roman"))
plot(x = grwl_matched_widths$gage_width_m, y = grwl_matched_widths$measured_width_m,
     axes=F, main='', 
     xlab="In Situ Width (m)", ylab="Remote Sensing Measured Width (m)",
     ylim=c(0, 1200), xlim=c(0,1200), 
     col= scales::alpha(primary[1], 0.5), cex.lab=1,
     pch = 16, bty = 'l',
     family = "sans")
points(x = s2_matched_widths$gage_width_m, y = s2_matched_widths$measured_width_m,
       col= scales::alpha(primary[2], 0.6), cex.lab=1,
       pch = 16
)

tx <- seq(0, 1200, length = 7)
axis(1, tx, c(0,200,400,600,800,1000,1200), lwd = 0.7, cex.axis = 1, family = "sans")
ty <- seq(0, 1200, length = 7)
axis(2, ty, c(0,200,400,600,800,1000,1200), lwd = 0.7, cex.axis = 1, family = "sans")

# 1:1 line
lines(x = c(0, 1200), y = c(0,1200), col = "#C1C1C0", type = 'l', lwd = 3.5, lty = 1)

# NAIP Linear Regression
naipLM_yValues_corr <- naipLMcoefficients_corr[1] + naipLMcoefficients_corr[2]*tx
lines(x = tx, y = naipLM_yValues_corr, pch = 0, col = primary[3], lty = 2, lwd = 3.5)

# s2 Linear Regression
s2LM_yValues <- s2LMcoefficients[1] + s2LMcoefficients[2]*tx
lines(x = tx, y = s2LM_yValues, pch = 1, col = primary[2], lty = 1, lwd = 3.5)

# GRWL Linear Regression
grwltx <- seq(0, 1200, length = 7)
grwlLM_yValues <- grwlLMcoefficients[1] + grwlLMcoefficients[2]*grwltx
lines(x = grwltx, y = grwlLM_yValues, pch = 4, col = primary[1], lty = 3, lwd = 3.5)

points(x = naip_matched_widths_corrected$gage_width_m, y = naip_matched_widths_corrected$measured_width_m,
       col= scales::alpha(primary[3], 1), cex.lab=1,
       pch = 16
)

#Rectangles for legend
# rect(
#   0,
#   850,
#   650,
#   1190,
#   border = 'black', lwd = 2
# )
# rect(
#   550,
#   0,
#   1190,
#   240,
#   border = 'black', lwd = 2
# )

# Create legend text and symbols
lines(x = c(0, 100), y = c(1150,1150), type = 'l', col = primary[1], lty = 3, lwd = 3.5)
lines(x = c(0, 100), y = c(950, 950), type = 'l', col = primary[2], lty = 1, lwd = 3.5)
lines(x = c(500, 600), y = c(190, 190), type = 'l', col = primary[3], lty = 2, lwd = 3.5)
lines(x = c(500, 600), y = c(40, 40), type = 'l', col = "#C1C1C0", lty = 1, lwd = 3.5)

text(x = 100, y = 1200, bquote(paste(bold('Mississippi River:')
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 1145, bquote(paste('y = ', .(grwlSlope), 'x + ', .(grwlInt),'; R'^'2'*' = 0.90'#, .(grwlR2)
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 1090, bquote(paste('Bias = ', .(grwlBias), ' m; ', 'RMSE = ', '43.0', ' m'
)), pos = 4, family = "sans", cex=0.8)

text(x = 100, y = 1000, bquote(paste(bold('Platte River:')
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 945, bquote(paste('y = ', .(s2Slope), 'x - ',.(abs(s2Int)),'; R'^'2'*' = 0.80'#, .(s2R2)
)), pos = 4, family = "sans", cex=0.8)
text(x = 100, y = 890, bquote(paste('Bias = ', .(s2Bias), ' m; ','RMSE = ', .(s2RMSE), ' m'
)), pos = 4, family = "sans", cex=0.8)

text(x = 600, y = 240, bquote(paste(bold('Corrected St. Vrain Creek:')
)), pos = 4, family = "sans", cex=0.8)
text(x = 600, y = 185, bquote(paste('y = ', .(naipSlope_corr), 'x - ',.(abs(naipInt_corr)),'; R'^'2'*' = ', .(naipR2_corr)
)), pos = 4, family = "sans", cex=0.8)
text(x = 600, y = 130, bquote(paste('Bias = ', .(naipBias_corr), ' m; ', 'RMSE = ', .(naipRMSE_corr), ' m'
)), pos = 4, family = "sans", cex=0.8)

text(x = 600, y = 35, bquote(paste(bold("Reference Line:"), " y = x"
)), pos = 4, family = "sans", cex=0.8)

dev.off()
######

ggplot(all_matched_widths, mapping = aes(x=gage_width_m, y=measured_width_m, color = scale, shape = scale))+
  geom_point()+
  stat_smooth(method = sen, se = FALSE, linewidth = 1.25)+
  geom_textabline(slope = 1, intercept = 0, color = 'black', linewidth = 1.25, label='slope = 1')+
  stat_poly_eq(parse=T, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula=y~x)+
  ggtitle('Validating Remote Sensing Measurements')+
  xlab('In Situ Width [m]')+
  ylab('Remote Sensing Measured Width [m]')+
  theme_classic()

#####

ggplot(all_matched_widths, aes(x=width_residual_m))+
  geom_histogram(aes(y=..density..), binwidth=5)+
  stat_function(color='red', fun = dnorm, args = list(mean = mean(all_matched_widths$width_residual_m), sd = sd(all_matched_widths$width_residual_m)))+
  geom_vline(xintercept = mean(all_matched_widths$width_residual_m), color = 'red')+
  geom_vline(xintercept = (mean(all_matched_widths$width_residual_m) - stdev(all_matched_widths$width_residual_m)), linetype="dashed", color='red')+
  geom_vline(xintercept = (mean(all_matched_widths$width_residual_m) + stdev(all_matched_widths$width_residual_m)), linetype="dashed", color='red')+
  xlab('Measured Width - In Situ Width (m)')+
  ylab('density')+
  ggtitle('Histogram of width residuals')+
  geom_text(x=200, y=0.016, label= paste0('\u03bc', " = 6.16 m"), color = 'red')+
  geom_text(x=200, y=0.015, label= paste0('\u03c3', " = 45.62 m"), color = 'red')+
  theme_classic()

ggplot(naip_matched_widths, mapping = aes(x=gage_width_m, y=measured_width_m))+
  geom_point()+
  stat_smooth(method = 'lm', se = FALSE, linewidth = 1.25, formula = y ~ x)+
  geom_textabline(slope = 1, intercept = 0, color = 'black', linewidth = 1.25, label = 'slope = 1', hjust = 0.6, vjust = -0.2)+
  stat_poly_eq(parse=T, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula=y~x)+
  ggtitle('Validating Corrected NAIP Measurements')+
  xlab('In Situ Width [m]')+
  ylab('Remote Sensing Measured Width [m]')+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_classic()
#Export all joined gage records to check for outliers
write_csv(all_matched_widths, "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/dirty_joined_gages.csv")

residual_dist <- fitdistr(all_matched_widths$width_residual_m, "normal")
grwl_residual_dist <- fitdistr(grwl_matched_widths$width_residual_m, "normal")
s2_residual_dist <- fitdistr(s2_matched_widths$width_residual_m, "normal")
naip_residual_dist <- fitdistr(naip_matched_widths$width_residual_m, "normal")
naip_residual_dist_corr <- fitdistr(naip_matched_widths_corrected$width_residual_m, "normal")

#Load field repeat measurements
field_validation <- read.csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/Field Validation/ssvc_fieldwork_validation.csv")
field_validation_means <- c(mean(field_validation[,1]),mean(field_validation[,2]),
                            mean(field_validation[,3]),mean(field_validation[,4]),
                            mean(field_validation[,5]),mean(field_validation[,6]),
                            mean(field_validation[,7]),mean(field_validation[,8]),
                            mean(field_validation[,9]),mean(field_validation[,10]),
                            mean(field_validation[,11]),mean(field_validation[,12]),
                            mean(field_validation[,13]),mean(field_validation[,14]),
                            mean(field_validation[,15]),mean(field_validation[,16]),
                            mean(field_validation[,17]),mean(field_validation[,18]),
                            mean(field_validation[,19]),mean(field_validation[,20]))
field_residuals <- c(field_validation[,1]-field_validation_means[1],
                     field_validation[,2]-field_validation_means[2],
                     field_validation[,3]-field_validation_means[3],
                     field_validation[,4]-field_validation_means[4],
                     field_validation[,5]-field_validation_means[5],
                     field_validation[,6]-field_validation_means[6],
                     field_validation[,7]-field_validation_means[7],
                     field_validation[,8]-field_validation_means[8],
                     field_validation[,9]-field_validation_means[9],
                     field_validation[,10]-field_validation_means[10],
                     field_validation[,11]-field_validation_means[11],
                     field_validation[,12]-field_validation_means[12],
                     field_validation[,13]-field_validation_means[13],
                     field_validation[,14]-field_validation_means[14],
                     field_validation[,15]-field_validation_means[15],
                     field_validation[,16]-field_validation_means[16],
                     field_validation[,17]-field_validation_means[17],
                     field_validation[,18]-field_validation_means[18],
                     field_validation[,19]-field_validation_means[19],
                     field_validation[,20]-field_validation_means[20])

field_residual_dist <- fitdistr(field_residuals, "normal")

field_residuals_df <- data.frame(field_residuals, scale = "Field")

#Plot distribution of GRWL width residuals
width_residuals_grwl <- ggplot()+
  geom_histogram(data = grwl_matched_widths, aes(x = width_residual_m#, fill = as.factor(scale)
                                                 ),
                 position = "identity",
                 binwidth = 10,
                 alpha = 1,
                 fill = primary[1],
                 color = primary[1])+
  #scale_fill_manual(values = c(hokies[4], hokies[3], hokies[2]), breaks = c("GRWL", "Sentinel2", "NAIP"))+
  # geom_histogram(data = grwl_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[4],
  #                color = "black")+
  # geom_histogram(data = s2_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[3],
  #                color = "black")+
  # geom_histogram(data = naip_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[2],
  #                color = "black")+
  # stat_function(fun = function(x){dnorm(x,mean = residual_dist$estimate[1], sd = residual_dist$estimate[2])*10*nrow(all_matched_widths)},
  #               color = "black", linewidth = 1.25, show.legend = F)+
  stat_function(data = grwl_matched_widths, fun = function(x){dnorm(x,mean = grwl_residual_dist$estimate[1], sd = grwl_residual_dist$estimate[2])*10*nrow(grwl_matched_widths)},
                color = "black", linewidth = 1, show.legend = F)+
  # stat_function(data = s2_matched_widths, fun = function(x){dnorm(x,mean = s2_residual_dist$estimate[1], sd = s2_residual_dist$estimate[2])*10*nrow(s2_matched_widths)},
  #               color = hokies[3], linewidth = 1.25, show.legend = F)+
  # stat_function(data = naip_matched_widths, fun = function(x){dnorm(x,mean = naip_residual_dist$estimate[1], sd = naip_residual_dist$estimate[2])*10*nrow(naip_matched_widths)},
  #               color = hokies[2], linewidth = 1.25, show.legend = F)+
  # stat_function(data = field_residuals_df, fun = function(x){dnorm(x, mean = field_residual_dist$estimate[1], sd = field_residual_dist$estimate[2])*10*nrow(field_residuals_df)},
  #               color = hokies[1], linewidth = 1.25, show.legend = F)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Count")+
  xlab("Width Residual (m)")+
  annotate("text",
           x = 200,
           y = 80,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03bc", " = ", round(grwl_residual_dist$estimate[1], 1), " m"),
           color = "black"
           , parse = F)+
  annotate("text",
           x = 200,
           y = 70,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03c3", " = ", round(grwl_residual_dist$estimate[2], 1), " m"),
           color = "black"
           , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 80,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03bc", " = ", round(s2_residual_dist$estimate[1], 3), " m"),
  #          color = hokies[3]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 70,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03c3", " = ", round(s2_residual_dist$estimate[2], 3), " m"),
  #          color = hokies[3]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 50,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03bc", " = ", round(naip_residual_dist$estimate[1], 3), " m"),
  #          color = hokies[2]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 40,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03c3", " = ", round(naip_residual_dist$estimate[2], 3), " m"),
  #          color = hokies[2]
  #          , parse = F)+
  guides(fill = guide_legend(title = "Data Source"))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 12))

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 6 Figures/width_residuals_grwl.svg",
       plot = width_residuals_grwl,
       width = 3,
       height = 3,
       units = "in")

#Plot distribution of S2 width residuals
width_residuals_s2 <- ggplot()+
  geom_histogram(data = s2_matched_widths, aes(x = width_residual_m#, fill = as.factor(scale)
  ),
  position = "identity",
  binwidth = 7,
  alpha = 1,
  fill = primary[2],
  color = primary[2])+
  #scale_fill_manual(values = c(hokies[4], hokies[3], hokies[2]), breaks = c("GRWL", "Sentinel2", "NAIP"))+
  # geom_histogram(data = grwl_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[4],
  #                color = "black")+
  # geom_histogram(data = s2_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[3],
  #                color = "black")+
  # geom_histogram(data = naip_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[2],
  #                color = "black")+
  # stat_function(fun = function(x){dnorm(x,mean = residual_dist$estimate[1], sd = residual_dist$estimate[2])*10*nrow(all_matched_widths)},
  #               color = "black", linewidth = 1.25, show.legend = F)+
  # stat_function(data = grwl_matched_widths, fun = function(x){dnorm(x,mean = grwl_residual_dist$estimate[1], sd = grwl_residual_dist$estimate[2])*10*nrow(grwl_matched_widths)},
  #               color = primary[1], linewidth = 1.25, show.legend = F)+
  stat_function(data = s2_matched_widths, fun = function(x){dnorm(x,mean = s2_residual_dist$estimate[1], sd = s2_residual_dist$estimate[2])*7*nrow(s2_matched_widths)},
                color = "black", linewidth = 1, show.legend = F)+
  # stat_function(data = naip_matched_widths, fun = function(x){dnorm(x,mean = naip_residual_dist$estimate[1], sd = naip_residual_dist$estimate[2])*10*nrow(naip_matched_widths)},
  #               color = hokies[2], linewidth = 1.25, show.legend = F)+
  # stat_function(data = field_residuals_df, fun = function(x){dnorm(x, mean = field_residual_dist$estimate[1], sd = field_residual_dist$estimate[2])*10*nrow(field_residuals_df)},
  #               color = hokies[1], linewidth = 1.25, show.legend = F)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Count")+
  xlab("Width Residual (m)")+
  # annotate("text",
  #          x = 200,
  #          y = 80,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03bc", " = ", round(grwl_residual_dist$estimate[1], 3), " m"),
  #          color = primary[1]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 70,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03c3", " = ", round(grwl_residual_dist$estimate[2], 3), " m"),
  #          color = primary[1]
  #          , parse = F)+
  annotate("text",
           x = 120,
           y = 7.25,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03bc", " = ", round(s2_residual_dist$estimate[1], 1), " m"),
           color = "black"
           , parse = F)+
  annotate("text",
           x = 120,
           y = 6.25,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03c3", " = ", round(s2_residual_dist$estimate[2], 1), " m"),
           color = "black"
           , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 50,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03bc", " = ", round(naip_residual_dist$estimate[1], 3), " m"),
  #          color = hokies[2]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 40,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03c3", " = ", round(naip_residual_dist$estimate[2], 3), " m"),
  #          color = hokies[2]
  #          , parse = F)+
  guides(fill = guide_legend(title = "Data Source"))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 12))

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 6 Figures/width_residuals_s2.svg",
       plot = width_residuals_s2,
       width = 3,
       height = 3,
       units = "in")

#Plot distribution of NAIP width residuals
width_residuals_naip <- ggplot()+
  geom_histogram(data = naip_matched_widths_corrected, aes(x = width_residual_m#, fill = as.factor(scale)
  ),
  position = "identity",
  binwidth = 1.5,
  alpha = 1,
  fill = primary[3],
  color = primary[3])+
  #scale_fill_manual(values = c(hokies[4], hokies[3], hokies[2]), breaks = c("GRWL", "Sentinel2", "NAIP"))+
  # geom_histogram(data = grwl_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[4],
  #                color = "black")+
  # geom_histogram(data = s2_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[3],
  #                color = "black")+
  # geom_histogram(data = naip_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[2],
  #                color = "black")+
  # stat_function(fun = function(x){dnorm(x,mean = residual_dist$estimate[1], sd = residual_dist$estimate[2])*10*nrow(all_matched_widths)},
  #               color = "black", linewidth = 1.25, show.legend = F)+
  # stat_function(data = grwl_matched_widths, fun = function(x){dnorm(x,mean = grwl_residual_dist$estimate[1], sd = grwl_residual_dist$estimate[2])*10*nrow(grwl_matched_widths)},
  #               color = primary[1], linewidth = 1.25, show.legend = F)+
  # stat_function(data = s2_matched_widths, fun = function(x){dnorm(x,mean = s2_residual_dist$estimate[1], sd = s2_residual_dist$estimate[2])*10*nrow(s2_matched_widths)},
  #               color = primary[2], linewidth = 1.25, show.legend = F)+
  stat_function(data = naip_matched_widths_corrected, fun = function(x){dnorm(x,mean = naip_residual_dist_corr$estimate[1], sd = naip_residual_dist_corr$estimate[2])*1.5*nrow(naip_matched_widths_corrected)},
                color = "black", linewidth = 1, show.legend = F)+
  # stat_function(data = field_residuals_df, fun = function(x){dnorm(x, mean = field_residual_dist$estimate[1], sd = field_residual_dist$estimate[2])*10*nrow(field_residuals_df)},
  #               color = hokies[1], linewidth = 1.25, show.legend = F)+
  scale_x_continuous(expand = c(0,0),
                     limits = c(-6,6))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Count")+
  xlab("Width Residual (m)")+
  # annotate("text",
  #          x = 200,
  #          y = 80,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03bc", " = ", round(grwl_residual_dist$estimate[1], 3), " m"),
  #          color = primary[1]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 70,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03c3", " = ", round(grwl_residual_dist$estimate[2], 3), " m"),
  #          color = primary[1]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 80,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03bc", " = ", round(s2_residual_dist$estimate[1], 3), " m"),
  #          color = primary[2]
  #          , parse = F)+
  # annotate("text",
  #          x = 200,
  #          y = 70,
  #          family= "sans",
  #          size = 8,
  #          label = paste0("\u03c3", " = ", round(s2_residual_dist$estimate[2], 3), " m"),
  #          color = primary[2]
  #          , parse = F)+
  annotate("text",
           x = 3,
           y = 1.3,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03bc", " = ", round(naip_residual_dist_corr$estimate[1], 1), " m"),
           color = "black"#primary[3]
           , parse = F)+
  annotate("text",
           x = 3,
           y = 1.15,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03c3", " = ", round(naip_residual_dist_corr$estimate[2], 1), " m"),
           color = "black"#primary[3]
           , parse = F)+
  guides(fill = guide_legend(title = "Data Source"))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 12))

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 6 Figures/width_residuals_naip.svg",
       plot = width_residuals_naip,
       width = 3,
       height = 3,
       units = "in")

#Plot distribution of field width residuals
width_residuals_field <- ggplot()+
  geom_histogram(data = field_residuals_df, aes(x = field_residuals#, fill = as.factor(scale)
                                                ),
                 position = "identity",
                 binwidth = 0.01,
                 alpha = 1,
                 fill = primary[4],
                 color = primary[4])+
  # scale_fill_manual(values = c(primary[4]), breaks = c("Field"))+
  # geom_histogram(data = grwl_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[4],
  #                color = "black")+
  # geom_histogram(data = s2_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[3],
  #                color = "black")+
  # geom_histogram(data = naip_matched_widths, aes(x = width_residual_m),
  #                binwidth = 10,
  #                alpha = 0.5,
  #                fill = hokies[2],
  #                color = "black")+
  # stat_function(fun = function(x){dnorm(x,mean = residual_dist$estimate[1], sd = residual_dist$estimate[2])*10*nrow(all_matched_widths)},
  #               color = "black", linewidth = 1.25, show.legend = F)+
  # stat_function(data = grwl_matched_widths, fun = function(x){dnorm(x,mean = grwl_residual_dist$estimate[1], sd = grwl_residual_dist$estimate[2])*10*nrow(grwl_matched_widths)},
  #               color = hokies[4], linewidth = 1.25, show.legend = F)+
  # stat_function(data = s2_matched_widths, fun = function(x){dnorm(x,mean = s2_residual_dist$estimate[1], sd = s2_residual_dist$estimate[2])*10*nrow(s2_matched_widths)},
  #               color = hokies[3], linewidth = 1.25, show.legend = F)+
  # stat_function(data = naip_matched_widths, fun = function(x){dnorm(x,mean = naip_residual_dist$estimate[1], sd = naip_residual_dist$estimate[2])*10*nrow(naip_matched_widths)},
  #               color = hokies[2], linewidth = 1.25, show.legend = F)+
  stat_function(data = field_residuals_df, fun = function(x){dnorm(x, mean = field_residual_dist$estimate[1], sd = field_residual_dist$estimate[2])*0.01*nrow(field_residuals_df)},
                color = "black", linewidth = 1, show.legend = F)+
  scale_x_continuous(expand = c(0,0),
                     limits = c(-0.1,0.2))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Count")+
  xlab("Width Residual (m)")+
  annotate("text",
           x = 0.12,
           y = 24.8,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03bc", " = ", "0.00"#round(field_residual_dist$estimate[1], 2)
                          , " m"),
           color = "black"
           , parse = F)+
  annotate("text",
           x = 0.12,
           y = 21.8,
           family= "sans",
           size = 10/.pt,
           label = paste0("\u03c3", " = ", round(field_residual_dist$estimate[2], 2), " m"),
           color = "black"
           , parse = F)+
  # guides(fill = guide_legend(title = "Data Source"))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 12))

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 6 Figures/width_residuals_field.svg",
       plot = width_residuals_field,
       width = 3,
       height = 3,
       units = "in")

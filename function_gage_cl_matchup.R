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
# naip_cl_pts <- read_csv("C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/naip_cl_pts.csv") %>% 
#   filter(reservoir == 0 & canal ==0)
naip_cl_pts_corrected <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_12_19_Correcting_Dry_NAIP/naip_cl_pts_corrected_exponential.csv") %>% 
  filter(reservoir == 0 & canal == 0)

#Use function to match in situ and measured widths
naip_matched_widths <- match_widths_corrected(naip_gages, naip_cl_pts_corrected, 10) %>% 
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
all_matched_widths <- rbind(grwl_matched_widths, s2_matched_widths, naip_matched_widths)

# Plot all matched width/discharge data points and three linear regression model fits
# Find linear regressions at NAIP, S2, and GRWL scales
naipCorrelationCoeff <- cor.test(naip_matched_widths$gage_width_m, naip_matched_widths$measured_width_m)
naipLM <- lm(measured_width_m ~ gage_width_m, data = naip_matched_widths)
naipLMcoefficients <- naipLM$coefficients
naipR2 <- round(summary(naipLM)$r.squared, 3)
naipRMSE <- round(sqrt(mean(naipLM$residuals^2)),3)
naipSlope <- round(naipLMcoefficients[2], 2)
naipInt <- round(naipLMcoefficients[1], 2)
  
s2CorrelationCoeff <- cor.test(s2_matched_widths$gage_width_m, s2_matched_widths$measured_width_m)
s2LM <- lm(measured_width_m ~ gage_width_m, data = s2_matched_widths)
s2LMcoefficients <- s2LM$coefficients
s2R2 <- round(summary(s2LM)$r.squared, 3)
s2RMSE <- round(sqrt(mean(s2LM$residuals^2)),3)
s2Slope <- round(s2LMcoefficients[2], 2)
s2Int <- round(s2LMcoefficients[1], 2)

grwlCorrelationCoeff <- cor.test(grwl_matched_widths$gage_width_m, grwl_matched_widths$measured_width_m)
grwlLM <- lm(measured_width_m ~ gage_width_m, data = grwl_matched_widths)
grwlLMcoefficients <- grwlLM$coefficients
grwlR2 <- round(summary(grwlLM)$r.squared, 3)
grwlRMSE <- round(sqrt(mean(grwlLM$residuals^2)),3)
grwlSlope <- round(grwlLMcoefficients[2], 2)
grwlInt <- round(grwlLMcoefficients[1], 2)

# Color palette
hokies <- c("#FFA85B", "#e36414", "#9a031e", "#5f0f40")

pdf("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 5 Figures/validation_plot.pdf",
    width = 11.25,
    height = 8.84,
    #units = 'in',
    bg = 'white'
    #res = 600
    )

# Plot
# windowsFonts(Times = windowsFont("Times New Roman"))
plot(x = grwl_matched_widths$gage_width_m, y = grwl_matched_widths$measured_width_m,
     axes=F, main='', 
     xlab="In Situ Width (m)", ylab="Remote Sensing Measured Width (m)",
     ylim=c(0, 1100), xlim=c(0,1100), 
     col= scales::alpha("#5f0f40", 0.3), cex.lab=1.5,
     pch = 16, bty = 'l',
     family = "sans")
points(x = s2_matched_widths$gage_width_m, y = s2_matched_widths$measured_width_m,
     # axes=F, main='', 
     # xlab="In Situ Width (m)", ylab="Remote Sensing Measured Width (m)",
     # ylim=c(0, 1100), xlim=c(0,1100), 
     col= scales::alpha("#9a031e", 0.3), cex.lab=1.25,
     pch = 16#, bty = 'l',
     # family = "sans",
)
points(x = naip_matched_widths$gage_width_m, y = naip_matched_widths$measured_width_m,
       # axes=F, main='', 
       # xlab="In Situ Width (m)", ylab="Remote Sensing Measured Width (m)",
       # ylim=c(0, 1100), xlim=c(0,1100), 
       col= scales::alpha("#e36414", 0.3), cex.lab=1.25,
       pch = 16#, bty = 'l',
       # family = "sans",
)
tx <- seq(0, 1100, length = 12)
axis(1, tx, c(0,100,200,300,400,500,600,700,800,900,1000,1100), lwd = 0.7, cex.axis = 1.5, family = "sans")
ty <- seq(0, 1100, length = 12)
axis(2, ty, c(0,100,200,300,400,500,600,700,800,900,1000,1100), lwd = 0.7, cex.axis = 1.5, family = "sans")

# 1:1 line
lines(x = c(0, 1100), y = c(0,1100), col = "#C1C1C0", type = 'l', lwd = 3.5, lty = 1)

# NAIP Linear Regression
naipLM_yValues <- naipLMcoefficients[1] + naipLMcoefficients[2]*tx
lines(x = tx, y = naipLM_yValues, pch = 0, col = hokies[2], lty = 2, lwd = 3.5)

# s2 Linear Regression
s2LM_yValues <- s2LMcoefficients[1] + s2LMcoefficients[2]*tx
lines(x = tx, y = s2LM_yValues, pch = 1, col = hokies[3], lty = 1, lwd = 3.5)

# GRWL Linear Regression
grwltx <- seq(0, 1100, length = 11.5)
grwlLM_yValues <- grwlLMcoefficients[1] + grwlLMcoefficients[2]*grwltx
lines(x = grwltx, y = grwlLM_yValues, pch = 4, col = hokies[4], lty = 3, lwd = 3.5)

#Rectangle for legend
rect(
  510,
  0,
  1100,
  400,
  border = 'black', lwd = 2
)

# Create legend text and symbols
lines(x = c(530, 585), y = c(320,320), type = 'l', col = hokies[4], lty = 3, lwd = 3.5)
lines(x = c(530, 585), y = c(240, 240), type = 'l', col = hokies[3], lty = 1, lwd = 3.5)
lines(x = c(530, 585), y = c(160, 160), type = 'l', col = hokies[2], lty = 2, lwd = 3.5)
lines(x = c(530, 585), y = c(80, 80), type = 'l', col = "#C1C1C0", lty = 1, lwd = 3.5)

# legend(650, 400, legend = c(bquote(paste('Mississippi River: y = ', .(grwlSlope), 'x + ', .(grwlInt), '; R'^'2'*' = 0.90')),
#                             bquote(paste('Platte River: y = ', .(s2Slope), 'x - 7.51','; R'^'2'*' = 0.80')),
#                             bquote(paste('St. Vrain Creek: y = ', .(naipSlope), 'x - 0.89','; R'^'2'*' = 0.99'))), 
#        bty = "n")

text(x = 600, y = 325, bquote(paste('Mississippi River: y = ', .(grwlSlope), 'x + ', .(grwlInt),
                                    '; RMSE = ', .(grwlRMSE), ' m'#'; R'^'2'*' = 0.90' #, .(grwlR2)
                                    )), pos = 4, family = "sans")
text(x = 600, y = 245, bquote(paste('Platte River: y = ', .(s2Slope), 'x - 7.51', #.(s2Int),
                                    '; RMSE = ', .(s2RMSE), ' m'#'; R'^'2'*' = 0.80' #, .(s2R2)
                                    )), pos = 4, family = "sans")
text(x = 600, y = 165, bquote(paste('St. Vrain Creek: y = ', .(naipSlope), 'x - 1.11', #.(naipInt),
                                    '; RMSE = ', .(naipRMSE), ' m'#'; R'^'2'*' = 0.99' #, .(naipR2)
                                   )), pos = 4, family = "sans")
text(x = 600, y = 85, "Reference Line: y = x", pos = 4, family = "sans")

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

#Plot distribution of width residuals
width_residuals <- ggplot()+
  geom_histogram(data = all_matched_widths, aes(x = width_residual_m, fill = as.factor(scale)),
                 position = "identity",
                 binwidth = 10,
                 alpha = 0.5,
                 color = "black")+
  scale_fill_manual(values = c(hokies[4], hokies[3], hokies[2]), breaks = c("GRWL", "Sentinel2", "NAIP"))+
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
                color = hokies[4], linewidth = 1.25, show.legend = F)+
  stat_function(data = s2_matched_widths, fun = function(x){dnorm(x,mean = s2_residual_dist$estimate[1], sd = s2_residual_dist$estimate[2])*10*nrow(s2_matched_widths)},
                color = hokies[3], linewidth = 1.25, show.legend = F)+
  stat_function(data = naip_matched_widths, fun = function(x){dnorm(x,mean = naip_residual_dist$estimate[1], sd = naip_residual_dist$estimate[2])*10*nrow(naip_matched_widths)},
                color = hokies[2], linewidth = 1.25, show.legend = F)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Count")+
  xlab("Width Residual (m)")+
  annotate("text",
           x = 200,
           y = 110,
           family= "sans",
           size = 7,
           label = paste0("\u03bc", " = ", round(grwl_residual_dist$estimate[1], 3), " m"),
           color = hokies[4]
           , parse = F)+
  annotate("text",
           x = 200,
           y = 100,
           family= "sans",
           size = 7,
           label = paste0("\u03c3", " = ", round(grwl_residual_dist$estimate[2], 3), " m"),
           color = hokies[4]
           , parse = F)+
  annotate("text",
           x = 200,
           y = 80,
           family= "sans",
           size = 7,
           label = paste0("\u03bc", " = ", round(s2_residual_dist$estimate[1], 3), " m"),
           color = hokies[3]
           , parse = F)+
  annotate("text",
           x = 200,
           y = 70,
           family= "sans",
           size = 7,
           label = paste0("\u03c3", " = ", round(s2_residual_dist$estimate[2], 3), " m"),
           color = hokies[3]
           , parse = F)+
  annotate("text",
           x = 200,
           y = 50,
           family= "sans",
           size = 7,
           label = paste0("\u03bc", " = ", round(naip_residual_dist$estimate[1], 3), " m"),
           color = hokies[2]
           , parse = F)+
  annotate("text",
           x = 200,
           y = 40,
           family= "sans",
           size = 7,
           label = paste0("\u03c3", " = ", round(naip_residual_dist$estimate[2], 3), " m"),
           color = hokies[2]
           , parse = F)+
  guides(fill = guide_legend(title = "Data Source"))+
  theme_classic()+
  theme(text = element_text(family = "sans", size = 20))

ggsave(filename = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 4 Figures/Fig A2 Width Residuals/width_residuals.svg",
       plot = width_residuals,
       width = 9.10,
       height = 6.02,
       units = "in")

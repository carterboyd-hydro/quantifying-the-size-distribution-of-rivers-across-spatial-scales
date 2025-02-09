#Author: Carter Boyd
#Date: 12/22/2023
#Purpose: create a regional rating curve to correct dry NAIP widths

#Attach libraries
library(dataRetrieval)
library(tidyverse)
library(patchwork)
library(igraph)
library(poweRlaw)

#Identify gages in the NAIP study area
gages <- c('06730400','06726900','06727500','06730200','06724970','06730500','06725450','06730525','06721500')

#Create empty dataframe to append flow data
flowDF <- data.frame()

for (i in gages){
  dv <- readNWISdv(
    siteNumbers = i,
    parameterCd = '00060'
  ) %>% renameNWISColumns() %>% 
  mutate(percentile=percent_rank(Flow)*100)
  flowDF <- rbind(flowDF, dv)
}

flowDF$site_no <- as.numeric(flowDF$site_no)
flowDF_old <- filter(flowDF, Date <= "2023-12-22")

#Group flowDF by site number
flowDF_grouped <- group_by(flowDF, site_no) %>% 
  filter(site_no == '6730525')
  

ggplot(flowDF_grouped, aes(Date, percentile, color = site_no))+
  geom_line()+
  theme_classic()

#Plot and compare percentile rank for the period of NAIP imagery, S2 imagery, and fieldwork
flowDF_NAIP <-  filter(flowDF, Date >= as.Date('2021-07-26') & Date <= as.Date('2021-08-25'))
print(paste0('NAIP Period mean percentile flow: ', mean(flowDF_NAIP$percentile)))
meanLowPercFlow <- mean(flowDF_NAIP$percentile)
percFlowDuringNAIP <-  ggplot(flowDF_NAIP, aes(Date, percentile
             # , color = site_no
             ))+
  geom_point()+
  geom_hline(yintercept = mean(flowDF_NAIP$percentile))+
  ggtitle('NAIP Period')+
  ylab('Percentile Flow')+
  theme_classic()

flowDF_S2 <-  filter(flowDF, Date >= as.Date('2019-07-01') & Date <= as.Date('2019-09-01'))
print(paste0('S2 Period mean percentile flow: ', mean(flowDF_S2$percentile)))
percFlowDuringS2 <-  ggplot(flowDF_S2, aes(Date, percentile
             # , color = site_no
             ))+
  geom_point()+
  geom_hline(yintercept = mean(flowDF_S2$percentile))+
  ggtitle('S2 Period')+
  ylab('Percentile Flow')+
  #geom_text(paste0("Mean Percentile Q: ", mean(flowDF_S2$percentile)))+
  theme_classic()

flowDF_Field <-  filter(flowDF, Date >= as.Date('2023-07-07') & Date <= as.Date('2023-08-06'))
print(paste0('Field Period mean percentile flow: ', mean(flowDF_Field$percentile)))
percFlowDuringField <-  ggplot(flowDF_Field, aes(Date, percentile
             # , color = site_no
             ))+
  geom_point()+
  geom_hline(yintercept = mean(flowDF_Field$percentile))+
  ggtitle('Fieldwork Period')+
  ylab('Percentile Flow')+
  #geom_text(paste0("Mean Percentile Q: ", mean(flowDF_Field$percentile)))+
  theme_classic()

percFlowDuringField + percFlowDuringNAIP + percFlowDuringS2

#Read in in_situ_width_and_discharge file
#Filter in situ data to only include NAIP gages
in_situ <- read.csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/in_situ_width_and_discharge.csv") %>% 
  filter(site_no %in% c(6730400,6726900,6727500,6730200,6724970,6730500,6725450,6730525,6721500))

in_situ$Date <- ymd(in_situ$Date)

#Join the calculated percent rank flow to in situ data by date
in_situ_percentile <- left_join(in_situ, flowDF, by=c('Date', 'site_no')) %>%
  na.omit() #%>% 
  #select(q_cms, w_m, V4, Date, site_no, percentile) #%>% 
  #(Optional) Filter to a specific range of percentile flow
  #filter(percentile < 90)


gages <- c('6730400','6726900','6727500','6730200','6724970','6730500','6725450','6730525','6721500')

#Power law linear correction factor (using igraph package)
# columns <- c("site_no", "alpha", "xmin", "logLik", "KS.stat", "KS.p")
# plFits_byGage <- data.frame(matrix(nrow = 0, ncol = length(columns)))
# 
# for (i in gages){
#   tempDF <- filter(in_situ_percentile, site_no == i) %>% 
#     dplyr::select(w_m)
#   pl <- fit_power_law(x = unlist(tempDF),
#                       force.continuous = T,
#                       implementation = "plfit")
#   plFits_byGage <- rbind(plFits_byGage, c(i, pl$alpha, pl$xmin, pl$logLik, pl$KS.stat, pl$KS.p))
# }
# colnames(plFits_byGage) <- columns

#Exponential linear correction factor (using nls() method)
singleGageDF <- filter(in_situ_percentile, site_no == gages[3]) %>% dplyr::select(c("w_m", "percentile"))
x <- singleGageDF$percentile
y <- singleGageDF$w_m
df <- data.frame(x,y)
s <- seq(from = 0, to = 100, length = 33)
m <- nls(y ~ I(coeff*exp(exp*x)), data = df, start = list(coeff = 2, exp = 0.01), trace = T)
# m <- nls(y ~ I(coeff*x^exp), data = df, start = list(coeff = 1, exp = 1)
#          , trace = T)

plot(x, y, main = "Single Gage Exponential Rating Curve",
     xlab = "Percentile Flow", ylab = "Width (m)")
lines(s, stats::predict(m, list(x = s)), col = "green")

#Put the above in a for loop that writes the exponent, coefficient, and R2 (at least) into a df
columns <- c("site_no", "exponent", "coefficient", "r_squared")
expFits_byGage <- data.frame(matrix(nrow=0, ncol=length(columns)))

for (i in gages){
  tempDF <- filter(in_situ_percentile, site_no == i) %>% dplyr::select(c("w_m", "percentile"))
  x <- tempDF$percentile
  y <- tempDF$w_m
  df <- data.frame(x,y)
  m <- nls(y ~ I(coeff*exp(exp*x)), data = df, start = list(coeff = 2, exp = 0.01), trace = T)
  RSS <- sum(residuals(m)^2)
  TSS <- sum((y-mean(y))^2)
  r_squared <- 1 - (RSS/TSS)
  expFits_byGage <- rbind(expFits_byGage, c(i, coef(m)[2], coef(m)[1], r_squared))
}
colnames(expFits_byGage) <- columns

mean_r2 <- mean(as.numeric(expFits_byGage$r_squared))

#Code to view all paired width/discharge measurement in NAIP study area
plot(in_situ_percentile$percentile, in_situ_percentile$w_m,
     col = factor(in_situ_percentile$site_no),
     main = "Average Exponential Gage Fit",
     xlab = "Percentile Flow (%)",
     ylab = "Width (m)",
     pch = 16)

legend("topleft",
       title = "Site Number",
       legend = levels(factor(in_situ_percentile$site_no)),
       pch = 16,
       col = factor(levels(factor(in_situ_percentile$site_no))),
       ncol = 3)

#Code to view average power law fit
avg_exponent <- mean(as.numeric(expFits_byGage$exponent))
avg_coefficient <- mean(as.numeric(expFits_byGage$coefficient))

curve(avg_coefficient*exp(x*avg_exponent), 1, 100, add=T, col = "black", lwd = 5,
      lty = 1)

text(75, 50, expression(w==5.56**(0.0096*Q[perc])))

meanHighPercFlow <- mean(c(mean(flowDF_Field$percentile), mean(flowDF_S2$percentile)))
width_atHighFlow <- avg_coefficient*exp(meanHighPercFlow*avg_exponent)
width_atLowFlow <- avg_coefficient*exp(meanLowPercFlow*avg_exponent)
multiplicativeCorrFactor <- width_atHighFlow/width_atLowFlow



#Power law linear correction factor (using log-log linear fit)
columns <- c("site_no", "exponent", "coefficient")
plFits_byGage <- data.frame(matrix(nrow = 0, ncol = length(columns)))

for (i in gages){
  tempDF <- filter(in_situ_percentile, site_no == i)
  pl <- lm(log(tempDF$w_m) ~ log(tempDF$percentile))
  plFits_byGage <- rbind(plFits_byGage, c(i, pl$coefficients[2], exp(pl$coefficients[1])))
}
colnames(plFits_byGage) <- columns

#Code to view the power law fits at each gage
#Simply change the index from 1 to 9 to view each gage
index <- 5
singleGage <- filter(in_situ_percentile, site_no == gages[index])
plot(singleGage$percentile, singleGage$w_m,
     main = "Width vs Percentile Discharge at NAIP gages",
     xlab = "Percentile Flow (%)",
     ylab = "Width (m)")
curve(as.numeric(plFits_byGage$coefficient[index])*x^as.numeric(plFits_byGage$exponent[index]), 1, 100, add=T, col = "blue")

#Code to view all paired width/discharge measurement in NAIP study area
plot(in_situ_percentile$percentile, in_situ_percentile$w_m,
     main = "Average Power Law Gage Fit",
     xlab = "Percentile Flow (%)",
     ylab = "Width (m)")

#Code to view average power law fit
avg_exponent <- mean(as.numeric(plFits_byGage$exponent))
avg_coefficient <- mean(as.numeric(plFits_byGage$coefficient))

curve(avg_coefficient*x^avg_exponent, 1, 100, add=T, col = "blue")
#Code to view all fitted power laws
for (i in 1:9){
  curve(as.numeric(plFits_byGage$coefficient[i])*x^as.numeric(plFits_byGage$exponent[i]), 1, 100, add = T, col = "blue")
}

#Use average power law fit to find percentage change in width between dry and wet periods
width_atHighFlow <- avg_coefficient*meanHighPercFlow^avg_exponent
width_atLowFlow <- avg_coefficient*meanLowPercFlow^avg_exponent
multiplicativeCorrFactor <- width_atHighFlow/width_atLowFlow

filter(in_situ_percentile, site_no == gages[9]) %>% 
ggplot(aes(percentile, w_m, color = as.character(site_no)))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method='lm', formula= y~x)+
  ggtitle(paste0("Testing Power Law fit"))+
  theme_classic()
#Additive correction factor
# lmFits_byGage <- data.frame()
# 
# for (i in gages){
#   tempDF <- filter(in_situ_percentile, site_no == i)
#   lm <- lm(w_m ~ percentile, data = tempDF)
#   lmCoefficients <- lm$coefficients
#   r2 <- as.numeric(round(summary(lm)$r.squared, 3))
#   slope <- as.numeric(round(lmCoefficients[2],2))
#   int <- as.numeric(round(lmCoefficients[1],2))
#   lmFits_byGage <- rbind(lmFits_byGage, c(i, slope, int, r2))
# }
# 
# #Create linear model fits by gage with no percentile flow threshold
# 
# # lmFits_byGage <- rename(lmFits_byGage, 
# #        site_no = "X.6730400.",
# #        slope = "X.0.02.",
# #        int = "X.2.77.",
# #        r2 = "X.0.221.")
# 
# #Create linear model fits by gage with a 90 percentile flow threshold
# lmFits_byGage <- rename(lmFits_byGage, 
#                         site_no = "X.6730400.",
#                         slope = "X.0.01.",
#                         int = "X.3.29.",
#                         r2 = "X.0.125.")
# 
# lmFits_byGage <- mutate(lmFits_byGage,
#                         site_no = as.numeric(lmFits_byGage$site_no),
#                         slope = as.numeric(lmFits_byGage$slope),
#                         int = as.numeric(lmFits_byGage$int),
#                         r2 = as.numeric(lmFits_byGage$r2))
# 
# #Plot percent rank flow vs width for all width measurements at all gages
# #filter(in_situ_percentile, site_no == gages[9]) %>% 
# ggplot(in_situ_percentile, aes(percentile, w_m, color = as.character(site_no)))+
#   geom_point()+
#   #geom_smooth(method = 'lm')+
#   geom_abline(aes(intercept = lmFits_byGage[1,3], slope = lmFits_byGage[1,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[2,3], slope = lmFits_byGage[2,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[3,3], slope = lmFits_byGage[3,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[4,3], slope = lmFits_byGage[4,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[5,3], slope = lmFits_byGage[5,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[6,3], slope = lmFits_byGage[6,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[7,3], slope = lmFits_byGage[7,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[8,3], slope = lmFits_byGage[8,2]))+
#   geom_abline(aes(intercept = lmFits_byGage[9,3], slope = lmFits_byGage[9,2]))+
#   theme_classic()
# 
# naipCorrelationCoeff <- cor.test(in_situ_percentile$percentile, in_situ_percentile$w_m)
# naipLM <- lm(w_m ~ percentile, data = in_situ_percentile)
# naipLMcoefficients <- naipLM$coefficients
# naipR2 <- round(summary(naipLM)$r.squared, 3)
# naipSlope <- round(naipLMcoefficients[2], 2)
# naipInt <- round(naipLMcoefficients[1], 2)
# 
# #Find mean of regression slopes for all gages
# meanSlope <- mean(lmFits_byGage$slope)
# 
# #Find mean of mean percentile flow for S2 and Field periods of observation
# meanHighPercFlow <- mean(c(mean(flowDF_Field$percentile), mean(flowDF_S2$percentile)))
# 
# #Find the difference in mean percentile flow between NAIP and S2/Field periods of observation
# diffPercFlow <- meanHighPercFlow - mean(flowDF_NAIP$percentile)
# 
# #Create a correction factor by multiplying the difference in mean percentile flow
# #by the mean of regression slopes for all gages
# corrFactor <- diffPercFlow*meanSlope

#Read in NAIP centerlines (points, not lines)
naip_cl_pts <- read_csv("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_10_24_Validation/import data/naip_cl_pts.csv") %>% 
  filter(reservoir == 0 & canal ==0)

#Create a new column for additive corrected widths in the NAIP centerline file
# naip_cl_pts_corrected <- mutate(naip_cl_pts, width_corrected = width+corrFactor)

#Create a new column for multiplicative corrected widths in the NAIP centerline file
naip_cl_pts_corrected <- mutate(naip_cl_pts, width_corrected = width*multiplicativeCorrFactor)

#Read out corrected NAIP centerline file
write.csv(naip_cl_pts_corrected, file = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_12_19_Correcting_Dry_NAIP/naip_cl_pts_corrected_exponential.csv")


#Calculate parameters of rating curve using MLE or similar to optimize (minimize)
#the sum square error between observations and rating curve
#(look for a consistent relationship between width and discharge)
#(might need to calculate percent rank of width too???)


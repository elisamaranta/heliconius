# Environmental analysis of hybrid zone 
# Elisa Mogollon, September 2024

library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(raster)
library(sp)
library(data.table)
library(mapview)
library(gridExtra)

setwd("C:/Users/USER/OneDrive - CACTUS/Documents/OneDrive/EES Year 2/IRT3/Climate data")

## For CHELSA data
# I downloaded raster files as .tif files directly from CHELSA website and then read the raster into R. 

#Temp
chelsa_1 <- raster("CHELSA_bio1_1981-2010_V.2.1.tif")
#Rainfall
chelsa_12 <- raster("CHELSA_bio12_1981-2010_V.2.1.tif")
#Light intensity (rsds)
rsds_mean <- raster("CHELSA_rsds_1981-2010_mean_V.2.1.tif")
#Net primary productivity
npp <- raster("CHELSA_npp_1981-2010_V.2.1.tif")

#r <- raster("CHELSA_swb_2018_V.2.1.tif")
#r <- r[[c(1,12)]]

plot(chelsa_1)
plot(chelsa_12)
plot(rsds_mean)
summary(r)

## For woldclim data
# Following https://gis.stackexchange.com/questions/227585/using-r-to-extract-data-from-worldclim

r <- getData("worldclim",var="bio",res=0.5, lat=-1, lon=-78)
r

#keep only bio1 and bio12, which are mean annual temp. and annual precip.
r <- r[[c(1,12)]] 
names(r) <- c("Temp","Prec")
r

#plot
samples <- read.csv("Final_observations.csv",  fileEncoding="UTF-8-BOM")


head(samples)
m <- mapview(samples, xcol = "longitude", ycol = "latitude", 
             zcol = "type", crs = 4326, grid = FALSE, alpha = 0, 
             map.types = c("Esri.WorldShadedRelief", "CartoDB.Positron"), color = "grey40")
m

#extract variables
long <- samples$longitude
lats <- samples$latitude
df <- as.data.table(samples)

head(df)
coordinates(df) <- ~ longitude+latitude

head(samples)

df <- SpatialPoints(coords = samples[, 10:11])
summary(df)

coords <- data.frame(x=long, y=lats)

points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84")) #r@crs
head(points)

#now with CHELSA data

valuesc1 <- raster::extract(chelsa_1, points)
valuesc12 <- raster::extract(chelsa_12, points)
values_rsds_mean <- raster::extract(rsds_mean, points)


df_temp <- cbind.data.frame(coordinates(points),valuesc1)
df12 <- cbind.data.frame(coordinates(points),valuesc12)
df1 <- cbind.data.frame(coordinates(points),values_rsds_mean)
df1

#In some cases, problems have been reported using older versions of GDAL or ArcGIS. 
#In this case the scale and offset has to be set manually for the variables.
#This can be done by first multiplying the raster values with the 'scale' value and then adding
#the 'offset' value

df1$valuesc1 <- (df1$valuesc1*0.1)-273.15
df12$valuesc12 <- (df12$valuesc12*0.1)+0

head(df1) #all data
head(df12)

#compare the two datasets
par(mfrow=c(2,2))
hist(df$Temp, main="WorldClim annual temp", xlim=c(120,260))
hist(df$Prec, main="WorldClim annual rainfall", xlim=c(500,6000))
hist(df1$valuesc1, main="CHELSA annual temp", xlim=c(12,26))
hist(df12$valuesc12, main="CHELSA annual rainfall", xlim=c(500,6000))

#plot
r1 <- crop(r[[(1)]], extent(-80,-79, -4.5, -3.5))
r12 <- crop(r[[(2)]], extent(-80,-79, -4.5, -3.5))

#-78.4,-77.5, -1.7, -0.9
#-80,-79, -4.5, -3.5
#-82,-76, -5, 0

plot(r1)
plot(points,add=T)

plot(r12)
plot(points,add=T, main="WorldClim annual rainfall")

chelsa_1 <- crop(chelsa_1, extent(-78.3,-77.3, -1.8, -0.8))
chelsa_12 <- crop(chelsa_12, extent(-78.3,-77.3, -1.8, -0.8))
rsds_mean1 <- crop(rsds_mean, extent(-78.3,-77.3, -1.8, -0.8))
npp1 <- crop(npp, extent(-78.3,-77.3, -1.8, -0.8))

chelsa_1 <- (chelsa_1*0.1)-273.15
chelsa_12<- (chelsa_12*0.1)+0

plot(chelsa_1)
plot(points,add=T, main="CHELSA annual temp")

plot(chelsa_12)
plot(points,add=T, main="CHELSA annual temp")

plot(npp1)
plot(points,add=T, main="CHELSA rsds mean")
points
dev.off()


#Merge env data with original dataframe

df2 <- cbind.data.frame(long, df1$values_rsds_mean, df12$valuesc12, df_temp$valuesc1)

#correct names
names(df2)[names(df2) == 'df1$values_rsds_mean'] <- 'values_rsds_mean'
names(df2)[names(df2) == 'df12$valuesc12'] <- 'annual_rain'
names(df2)[names(df2) == 'df_temp$valuesc1'] <- 'mean_temp'

#make temp in correct units
df2$mean_temp <- df2$mean_temp/100

################## FOREST COVER ##################

#Forest cover raster (in .tif format) was obtained from the Global Forest Watch database.
#I specifically donwloaded the .tif file for forest cover for the appropriate coordinates (00N_080W)
library(raster)
forest_cover <- raster("00N_080W.tif") #THIS IS THE ONE I NEED!

plot(forest_cover)

df2 <- read.csv("eye_data_with_env_variables.csv")

# Convert data frame to SpatialPoints
coordinates(df2) <- ~ longitude + latitude

# Extract forest cover values for these points
forest_cover_values <- extract(forest_cover, df2)
hist(forest_cover_values)

length(forest_cover_values)

# Combine coordinates with forest cover values
results <- data.frame(cbind(df2, forest_cover_values))
colnames(results)[42] <- "forest_cover"

#write.csv(results, "eye_data_with_env_variables_and_forest_cover.csv")


## Analysis of correlations between environmental variables and eye size ####

# import data
df2 <- read.csv("eye_data_with_env_variables_and_forest_cover.csv")
head(df2)

#Model eye size (corneal area) along with each environmental variable ####

library(lme4)
library(lmerTest)  # Extends lme4 to provide p-values for mixed models

# ELEVATION #

# Fit the mixed-effects model
model <- lmer(CA_residuals ~ Elevation + (1 | location), data = df2, REML = TRUE)

plot(model)

# Summary of the model
summary(model)

# Extract the estimate (slope) for the fixed effect (Elevation)
model_summary <- summary(model)
estimate <- model_summary$coefficients[2, "Estimate"]  # Extract the estimate value for Elevation
p_value <- model_summary$coefficients[2, "Pr(>|t|)"]   # Extract the p-value for Elevation
# Generate predictions from the model
df2$predicted <- predict(model)

#theme settings for plots
par(family = "serif")

# Create a scatter plot with a regression line for the mixed model
ggplot(df2, aes(x = Elevation, y = CA_residuals)) +
  geom_point(alpha = 0.2) +  # Scatter plot of individual data points
  geom_line(aes(y = predicted), color = "black", size = 1) +  # Fitted line for the fixed effect
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.2) +  # Global regression line for fixed effect
  labs(x = "Elevation (m)",
       y = "CA residuals") +
  theme_classic() +
  theme(legend.position = "none") +  # Remove legend for site colors if not needed
  # Annotate estimate value and p-value on the plot
  annotate("text", x = Inf, y = Inf, label = paste("Estimate =", round(estimate, 4)), 
           hjust = 1.9, vjust = 2, size = 3.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("p-val =", round(p_value, 6)), 
           hjust = 1.1, vjust = 2, size = 3.5, color = "black")+
  theme(legend.position = "none")  # Remove legend for site colors if not needed


# TEMPERATURE #

# Fit the mixed-effects model
model <- lmer(CA_residuals ~ mean_temp + (1 | location), data = df2)

plot(model)
# Summary of the model
summary(model)

# Extract the estimate (slope) for the fixed effect (Elevation)
model_summary <- summary(model)
estimate <- model_summary$coefficients[2, "Estimate"]  # Extract the estimate value for Elevation
p_value <- model_summary$coefficients[2, "Pr(>|t|)"]   # Extract the p-value for Elevation
# Generate predictions from the model
df2$predicted <- predict(model)

# Create a scatter plot with a regression line for the mixed model
ggplot(df2, aes(x = mean_temp, y = CA_residuals)) +
  geom_point(alpha = 0.2) +  # Scatter plot of individual data points
  geom_line(aes(y = predicted), color = "black", size = 1) +  # Fitted line for the fixed effect
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.2) +  # Global regression line for fixed effect
  labs(x = "Mean temperature (°C)",
       y = "CA residuals") +
  theme_classic() +
  theme(legend.position = "none") +  # Remove legend for site colors if not needed
  # Annotate estimate value and p-value on the plot
  annotate("text", x = Inf, y = Inf, label = paste("Estimate =", round(estimate, 4)), 
           hjust = 1.8, vjust = 2, size = 3.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("p-val =", round(p_value, 6)), 
           hjust = 1.1, vjust = 2, size = 3.5, color = "black")+
  theme(legend.position = "none")  # Remove legend for site colors if not needed

# RSDS #
head(df2)
# Fit the mixed-effects model
model <- lmer(CA_residuals ~ values_rsds_mean + (1 | location), data = df2)

plot(model)
# Summary of the model
summary(model)

# Extract the estimate (slope) for the fixed effect (Elevation)
model_summary <- summary(model)
estimate <- model_summary$coefficients[2, "Estimate"]  # Extract the estimate value for Elevation
p_value <- model_summary$coefficients[2, "Pr(>|t|)"]   # Extract the p-value for Elevation
# Generate predictions from the model
df2$predicted <- predict(model)

# Create a scatter plot with a regression line for the mixed model
ggplot(df2, aes(x = values_rsds_mean, y = CA_residuals)) +
  geom_point(alpha = 0.2) +  # Scatter plot of individual data points
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.2) +  # Global regression line for fixed effect
  labs(x = "rsds (W m-2)",
       y = "CA residuals") +
  theme_classic() +
  theme(legend.position = "none") +  # Remove legend for site colors if not needed
  # Annotate estimate value and p-value on the plot
  annotate("text", x = Inf, y = Inf, label = paste("Estimate =", round(estimate, 6)), 
           hjust = 1.8, vjust = 2, size = 3.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("p-val =", round(p_value, 2)), 
           hjust = 1.1, vjust = 2, size = 3.5, color = "black")+
  theme(legend.position = "none")  # Remove legend for site colors if not needed

# RAINFALL #

# Fit the mixed-effects model
model <- lmer(CA_residuals ~ annual_rain + (1 | location), data = df2)

plot(model)
# Summary of the model
summary(model)

# Extract the estimate (slope) for the fixed effect (Elevation)
model_summary <- summary(model)
estimate <- model_summary$coefficients[2, "Estimate"]  # Extract the estimate value for Elevation
p_value <- model_summary$coefficients[2, "Pr(>|t|)"]   # Extract the p-value for Elevation
# Generate predictions from the model
df2$predicted <- predict(model)

# Create a scatter plot with a regression line for the mixed model
ggplot(df2, aes(x = annual_rain, y = CA_residuals)) +
  geom_point(alpha = 0.2) +  # Scatter plot of individual data points
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.2) +  # Global regression line for fixed effect
  labs(x = "Annual rainfall (mm)",
       y = "CA residuals") +
  theme_classic() +
  theme(legend.position = "none") +  # Remove legend for site colors if not needed
  # Annotate estimate value and p-value on the plot
  annotate("text", x = Inf, y = Inf, label = paste("Estimate =", round(estimate, 8)), 
           hjust = 1.7, vjust = 2, size = 3.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("p-val =", round(p_value, 2)), 
           hjust = 1.1, vjust = 2, size = 3.5, color = "black")+
  theme(legend.position = "none")  # Remove legend for site colors if not needed

# FOREST COVER #

# Fit the mixed-effects model
model <- lmer(CA_residuals ~ forest_cover + (1 | location), data = df2)

plot(model)
# Summary of the model
summary(model)

# Extract the estimate (slope) for the fixed effect (Elevation)
model_summary <- summary(model)
estimate <- model_summary$coefficients[2, "Estimate"]  # Extract the estimate value for Elevation
p_value <- model_summary$coefficients[2, "Pr(>|t|)"]   # Extract the p-value for Elevation
# Generate predictions from the model
df2$predicted <- predict(model)

# Create a scatter plot with a regression line for the mixed model
ggplot(df2, aes(x = forest_cover, y = CA_residuals)) +
  geom_point(alpha = 0.2) +  # Scatter plot of individual data points
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.2) +  # Global regression line for fixed effect
  labs(x = "Forest cover (%)",
       y = "CA residuals") +
  theme_classic() +
  theme(legend.position = "none") +  # Remove legend for site colors if not needed
  # Annotate estimate value and p-value on the plot
  annotate("text", x = Inf, y = Inf, label = paste("Estimate =", round(estimate, 4)), 
           hjust = 1.8, vjust = 2, size = 3.5, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("p-val =", round(p_value, 2)), 
           hjust = 1.1, vjust = 2, size = 3.5, color = "black")+
  theme(legend.position = "none")  # Remove legend for site colors if not needed

#### PCA for environmental variables along sites ####
#using https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples 

# Rearrange dataframe per site
fc <- as.data.frame(df2 %>% group_by(location,Elevation) %>%
                      summarise(rsds = mean(values_rsds_mean),
                                rainfall = mean(annual_rain),
                                temp = mean(mean_temp),
                                forest_cover_mean = mean(forest_cover)))

library(vegan)
#run PCA
PCA <- rda(fc[,-1], scale = TRUE) # the argument scale standardizes the variables
#inspect
head(summary(PCA))

#Get loadings for the first two PCs.
loadings <- scores(PCA, display = 'species', scaling = 0)
loadings

#Which env variables have the most impact on the top two PCs?
sort (abs (loadings[,1]), decreasing = TRUE) #temp and elevation have the highest correlation with first axis
sort (abs (loadings[,2]), decreasing = TRUE) #rsds and rainfall with second axis

#Get the amount of vaiance explained
explained_variance <- summary(PCA)$cont$importance[2,]
explained_variance

#Variance explained by PC1: 41.63085 %
#Variance explained by PC2: 33.75066 %

#Plot PCA
biplot (PCA, display = 'species', scaling = 'species')
dev.off()

#Plot with samples on it
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/cleanplot.pca.R')
cleanplot.pca (PCA, scaling = 2) 

#end

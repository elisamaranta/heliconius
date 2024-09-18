# HZAR R code for cline analysis of morphological data
# Elisa Mofollon Perez, September 2024.
#
# This script is used to conduct cline analyses of the corneal are and facet count residuals obtained
# from linear models producted in the first scipt (1. eye_morph_stats)
# The script can be run with either residuals obtained using abdomen or tibia length as the allometric control. 
# 
# Following code from:
# https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/ClineAnalysis/clines_highFSTsites.sh

setwd("C:/Users/USER/OneDrive - CACTUS/Documents/OneDrive/EES Year 2/IRT3/Measurement data")

#Load packages
library(ggplot2)
library(hzar)
library(dplyr)
library(ggpubr)
library(stats)
library(RColorBrewer)
library(dplyr)
library(scales)
require(grid)

#install.packages("doMC", repos="http://R-Forge.R-project.org")
#devtools::install_version("hzar",version="0.2-5",repos = "http://cran.us.r-project.org")
require(doMC)
require(hzar)

#theme settings for plots
pub<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour = NA,fill=NA),
  panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  axis.title.x=element_text(colour="grey20",size=12,face="plain",hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(colour="grey20",size=12,face="plain",hjust=0.5,vjust=1,angle=90),
  axis.text.x=element_text(colour="grey20",angle=0,size=10),
  axis.text.y=element_text(colour="grey20",angle=0,size=10),
  axis.ticks=element_line(colour="black"),
  text = element_text(family = "Arial"))

# Here, the script is run with the abdomen residuals dataframe produced in 1.Eye_size_analysis script
# But the script was equally run with the tibia residuals

df <- read.csv("eye_data_with_residuals_abdomen_only.csv",  fileEncoding="UTF-8-BOM")
#df <- read.csv("eye_data_with_residuals_tibia_only.csv",  fileEncoding="UTF-8-BOM")


#Initial exploration
as.data.frame(df %>% group_by(type) %>%
                      summarise(nSamples = n()))


ggplot(df, aes(x=elevation, y=FC_residuals, colour=sex)) +
  geom_point() + geom_smooth(method = lm)

#Load admixture values from Wing Shape paper (Montejo-Kovacevich et al. 2021 Mol Ecol)
ngs <- read.csv("wing_shape_ngs_admix.csv")
head(ngs)

#id         sex.cov    area aspect.ratio sex.cov.1 sex ngs.admix  subsp
#1 CAM016019    male 494.735        2.138      male   0 0.5677630 hybrid
#2 CAM016020    male 353.010        2.103      male   0 0.9999998 hybrid

# Merge the 'admixture' and 'aspect.ratio' columns from ngs into my dataset based on ID
df <- merge(df, ngs[, c("id", "aspect.ratio","ngs.admix")], by.x = "ID", by.y = "id", all.x = TRUE)

hist(df$aspect.ratio)
hist(df$ngs.admix)

samples <- df %>% drop_na(latitude, longitude) #from our samples data
counts <- data.frame(df %>% count(elevation, longitude, type))

# Plot of sampling across elevation and longitude
ggplot(counts, aes(longitude, elevation, size = n, color=type)) + 
  geom_point(alpha = 0.9) +
  scale_color_manual(values=c("hybrid" = "#A6CEE3", "lativitta" = "#1F78B4", "notabilis" = "#B2DF8A"),
                     name="")  +
  scale_size_continuous(name = "Samples", range = c(3,12), breaks = c(5,20,40)) +
  labs(x = "", y = "Elevation")  +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "grey", stroke = 1, colour = "black")),
         color = guide_legend(override.aes = list(size = 5))) +
         #color = "none") +  # Remove the category legend
  theme(legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.box = "horizontal") 

par(mfrow=c(1,2))

#How does increase in eye size with elevation vary with subspecies
head(df)
ggplot(df, aes(x=type, y=ngs.admix, colour=type)) +
  geom_boxplot() + geom_point(alpha=0.3, size=2)

#how do residuals correlate?

plot(df$FC_residuals ~ df$CA_residuals, main = "Abdomen model residuals", ylab="FC residuals", xlab="CA residuals")
plot(df1$FC_residuals ~ df1$CA_residuals, main = "Tibia model residuals", ylab="FC residuals", xlab="CA residuals")

cor.test(df$FC_residuals, df$CA_residuals, method = "pearson")
##t = 24.986, df = 443, p-value < 2.2e-16
##  cor = 0.7648144 

# What is the standard deviation in the residuals like across the transect? ####

#Calculate sd across sites
fc <- as.data.frame(df %>% group_by(location,Transect.Position) %>%
                      summarise(fc_obs_mean = mean(FC_residuals),
                                fc_sd = sd(FC_residuals),
                                ca_obs_mean = mean(CA_residuals),
                                ca_sd = sd(CA_residuals),
                                nSamples = n()) %>%
                      rename(locality = location,
                             distance = Transect.Position))

# Order by distance and remove off-transect sites.
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

# Plot
fc %>%
  ggplot(aes(x=distance, y=fc_obs_mean))+
  geom_point(col="red")+
  geom_errorbar(aes(ymin = fc_obs_mean-fc_sd, ymax = fc_obs_mean+fc_sd), 
                col=alpha("grey30", 0.2), size=2) + theme_classic()

fc %>%
  ggplot(aes(x=distance, y=ca_obs_mean))+
  geom_point(col="red")+
  geom_errorbar(aes(ymin = ca_obs_mean-ca_sd, ymax = ca_obs_mean+ca_sd), 
                col=alpha("grey30", 0.2), size=2) + theme_classic()


# plot ID of samples across transect for closer inspection of any potential outliers
ggplot(df, aes(Transect.Position, FC_residuals, label=ID)) +
  geom_point() +
  geom_text(hjust = -0.2, size=2)


# 1. CLINE ANALYSIS FOR FACET COUNT ####

## Format data into correct layout


# Reformint into the following layout:
# locality # distance(km) # obs_mean # obs_var # nSamples #

# Facet count dataframe
fc <- as.data.frame(df %>% group_by(location, Transect.Position) %>%
  summarise(obs_mean = mean(FC_residuals),
            obs_var = var(FC_residuals),
            nSamples = n()) %>%
  rename(locality = location,
         distance = Transect.Position))

fc
# Order by distance and remove missing transect positions (i.e. off-transect sites) 
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

#now plot with mean variable
ggplot() +
  geom_point(data = df, aes(Transect.Position, FC_residuals), color="grey", alpha=0.3, size=2) +
  geom_point(data = fc, aes(distance, obs_mean), color = "red", pch=19, size=2.5) +
  labs(x = "Distance (km)", y = "Facet count residuals") +
  theme(text = element_text(family = "Times New Roman", face = "plain"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))

# Set chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(978,544,99,596,528,124), 
              B=c(544,99,596,528,124,978), 
              C=c(99,596,528,124,978,544), 
              D=c(978,99,596,528,124,544), 
              E=c(124,596,978,99,528,544))

# Prepare an object holding everything
erato <- list()
## Space to hold the observed data
erato$obs <- list();
## Space to hold the models to fit
erato$models <- list();
## Space to hold the compiled fit requests
erato$fitRs <- list();
## Space to hold the output data chains
erato$runs <- list();
## Space to hold the analysed data
erato$analysis <- list();


# Generate a hzar.obsData object 
erato$obs <- hzar.doNormalData1DPops(distance = as.list(fc$distance), muObs = as.list(fc$obs_mean),
                                     varObs = as.list(fc$obs_var), nEff = as.list(fc$nSamples),
                                     siteID=paste("P",1:length(as.list(fc$distance)),sep=""))

## Look at a graph of the observed data
hzar.plot.obsData(erato$obs, ylim = c(-0.07,0.07), xlim = c(0,85));
erato$obs$frame
## Make a helper function to define tails and name of models
erato.loadAdaAmodel <- function(tails,id=paste(tails))
  erato$models[[id]] <<- hzar.makeCline1DNormal(erato$obs, tails)

# Add five possible models
erato.loadAdaAmodel("none","modelI")
erato.loadAdaAmodel("right","modelII")
erato.loadAdaAmodel("left","modelIII")
erato.loadAdaAmodel("mirror","modelIV")
erato.loadAdaAmodel("both","modelV")


## Check the default settings
print(erato$models)

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km

erato$models <- sapply(erato$models,
                       hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting 
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
erato$fitRs$init <- sapply(erato$models,
                               hzar.first.fitRequest.gC,
                               obsData=erato$obs,
                               verbose=FALSE,
                               simplify=FALSE)

#erato$fitRs$init <- sapply(erato$models,
#                           hzar.first.fitRequest.old.ML, obsData = erato$obs, 
#                           simplify=FALSE, verbose = FALSE)

## Check fit request settings
print(erato$fitRs$init)

## Run each model for an initial chain
erato$runs$init <- list()
erato$runs$init$modelI <- hzar.doFit(erato$fitRs$init$modelI)
erato$runs$init$modelII <- hzar.doFit(erato$fitRs$init$modelII)
erato$runs$init$modelIII <- hzar.doFit(erato$fitRs$init$modelIII)
erato$runs$init$modelIV <- hzar.doFit(erato$fitRs$init$modelIV)
erato$runs$init$modelV <- hzar.doFit(erato$fitRs$init$modelV)

## Plot the trace of a model
plot(hzar.mcmc.bindLL(erato$runs$init$modelI))

## Compile a new set of fit requests using the initial chains
erato$fitRs$chains <- lapply(erato$runs$init,hzar.next.fitRequest)


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
erato$fitRs$chains <- hzar.multiFitRequest(erato$fitRs$chains,
                                           each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])


##WARNING THIS TAKES A LONG TIME TO RUN###
# Run a chain of 3 runs for every fit request
erato$runs$chains <- hzar.doChain.multi(erato$fitRs$chains,
                                        doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Did modelII converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]]))))


## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
erato$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
erato$analysis$initDGs$modelI <-
  hzar.dataGroup.add(erato$runs$init$modelI)
erato$analysis$initDGs$modelII <-
  hzar.dataGroup.add(erato$runs$init$modelII)
erato$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(erato$runs$init$modelIII)
erato$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(erato$runs$init$modelIV)
erato$analysis$initDGs$modelV <-
  hzar.dataGroup.add(erato$runs$init$modelV)

 ## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
erato$analysis$oDG <-hzar.make.obsDataGroup(erato$analysis$initDGs)

# Use the same model names
erato$analysis$oDG <- hzar.copyModelLabels(erato$analysis$initDGs,
                                               erato$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
erato$analysis$oDG <-hzar.make.obsDataGroup(lapply(erato$runs$chains,
                                                       hzar.dataGroup.add),erato$analysis$oDG)

## Check to make sure that there are only 5 hzar.dataGroup
## in the hzar.obsDataGroup object.
print(summary(erato$analysis$oDG$data.groups))


## Do model selection based on the AICc scores
print(erato$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(erato$analysis$oDG));

#AICc
#modelI   -1206.003
#modelII  -1200.914
#modelIII -1201.779
#modelIV  -1201.892
#modelV   -1197.488

## Compare the 3 cline models to the null model graphically
hzar.plot.cline(erato$analysis$oDG, ylim = c(-0.08,0.08));

## Print out the model with the minimum AICc score
print(erato$analysis$model.name <-
        rownames(erato$analysis$AICcTable)[[which.min(erato$analysis$AICcTable$AICc)]])
# "modelI"

# Compute delta AICc (should be above 2)
deltaAICc<-sort(erato$analysis$AICcTable$AICc)[2]-sort(erato$analysis$AICcTable$AICc)[1]
deltaAICc
#4.111092

## Extract the hzar.dataGroup object for the selected model
erato$analysis$model.selected <-
  erato$analysis$oDG$data.groups[[erato$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(erato$analysis$model.selected,
                         names(erato$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(erato$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(erato$analysis$model.selected, ylim = c(-0.08,0.08));

# Open a pdf file
#pdf("clineplot.pdf") 

## Plot the 95% credible cline region for the selected model
par(mfrow=c(1,1))
hzar.plot.fzCline(erato$analysis$model.selected, ylim = c(-0.08,0.08));

#Extract cline information into an object for plotting with specific format
fc_abd_cline_all_pos <- erato$analysis$model.selected


#fc_abd <- hzar.get.ML.cline(erato$analysis$model.selected)
#fc_abd_var <- hzar.getLLCutParam(erato$analysis$model.selected,
#                               names(erato$analysis$model.selected$data.param))

# 2. CLINE ANALYSIS FOR CORNEAL AREA #### 

# Facet count dataframe
fc <- as.data.frame(df %>% group_by(location,Transect.Position) %>%
                      summarise(obs_mean = mean(CA_residuals),
                                obs_var = var(CA_residuals),
                                nSamples = n()) %>%
                      rename(locality = location,
                             distance = Transect.Position))

head(fc)
# Order by distance
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

#now plot with mean residual values and all values
ggplot() +
  geom_point(data = df, aes(Transect.Position, CA_residuals), color="grey", alpha=0.3, size=2) +
  geom_point(data = fc, aes(distance, obs_mean), color = "red", pch=19, size=2.5) +
  #geom_smooth(data = fc, aes(distance, obs_mean), color = "grey10", size=1, method="lm")+
  labs(x = "Distance (km)", y = "Corneal area residuals") +
  theme(text = element_text(family = "Times New Roman", face = "plain"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))

# Set chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),
              B=c(528,124,978,544,99,596),
              C=c(124,978,544,99,596,528))

# Prepare an object holding everything
erato <- list()
## Space to hold the observed data
erato$obs <- list();
## Space to hold the models to fit
erato$models <- list();
## Space to hold the compiled fit requests
erato$fitRs <- list();
## Space to hold the output data chains
erato$runs <- list();
## Space to hold the analysed data
erato$analysis <- list();


# Generate a hzar.obsData object (HERE I SHOULD GET MISSING TRANSECT VALUES)
erato$obs <- hzar.doNormalData1DPops(distance = as.list(fc$distance), muObs = as.list(fc$obs_mean),
                                     varObs = as.list(fc$obs_var), nEff = as.list(fc$nSamples),
                                     siteID=paste("P",1:length(as.list(fc$distance)),sep=""))

## Look at a graph of the observed data
hzar.plot.obsData(erato$obs, ylim = c(-0.08,0.08));
erato$obs$frame
## Make a helper function to define tails and name of models
erato.loadAdaAmodel <- function(tails,id=paste(tails))
  erato$models[[id]] <<- hzar.makeCline1DNormal(erato$obs, tails)

# Add five possible models
erato.loadAdaAmodel("none","modelI")
erato.loadAdaAmodel("right","modelII")
erato.loadAdaAmodel("left","modelIII")
erato.loadAdaAmodel("mirror","modelIV")
erato.loadAdaAmodel("both","modelV")


## Check the default settings
print(erato$models)

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km.

erato$models <- sapply(erato$models,
                       hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting 
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
erato$fitRs$init <- sapply(erato$models,
                           hzar.first.fitRequest.gC,
                           obsData=erato$obs,
                           verbose=FALSE,
                           simplify=FALSE)

## Check fit request settings
print(erato$fitRs$init)

## Run each model for an initial chain
erato$runs$init <- list()
erato$runs$init$modelI <- hzar.doFit(erato$fitRs$init$modelI)
erato$runs$init$modelII <- hzar.doFit(erato$fitRs$init$modelII)
erato$runs$init$modelIII <- hzar.doFit(erato$fitRs$init$modelIII)
erato$runs$init$modelIV <- hzar.doFit(erato$fitRs$init$modelIV)
erato$runs$init$modelV <- hzar.doFit(erato$fitRs$init$modelV)

## Plot the trace of a model
plot(hzar.mcmc.bindLL(erato$runs$init$modelI))

## Compile a new set of fit requests using the initial chains
erato$fitRs$chains <- lapply(erato$runs$init,hzar.next.fitRequest)


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
erato$fitRs$chains <- hzar.multiFitRequest(erato$fitRs$chains,
                                           each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])

##WARNING THIS TAKES A LONG TIME TO RUN###
# Run a chain of 3 runs for every fit request
erato$runs$chains <- hzar.doChain.multi(erato$fitRs$chains,
                                        doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Did modelII converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]]))))


## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
erato$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
erato$analysis$initDGs$modelI <-
  hzar.dataGroup.add(erato$runs$init$modelI)
erato$analysis$initDGs$modelII <-
  hzar.dataGroup.add(erato$runs$init$modelII)
erato$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(erato$runs$init$modelIII)
erato$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(erato$runs$init$modelIV)
erato$analysis$initDGs$modelV <-
  hzar.dataGroup.add(erato$runs$init$modelV)

## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
erato$analysis$oDG <-hzar.make.obsDataGroup(erato$analysis$initDGs)

# Use the same model names
erato$analysis$oDG <- hzar.copyModelLabels(erato$analysis$initDGs,
                                           erato$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
erato$analysis$oDG <-hzar.make.obsDataGroup(lapply(erato$runs$chains,
                                                   hzar.dataGroup.add),erato$analysis$oDG)

## Check to make sure that there are only 5 hzar.dataGroup
## in the hzar.obsDataGroup object.
print(summary(erato$analysis$oDG$data.groups))


## Do model selection based on the AICc scores
print(erato$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(erato$analysis$oDG));


## Compare the 3 cline models to the null model graphically
hzar.plot.cline(erato$analysis$oDG, ylim = c(-0.08,0.08));

## Print out the model with the minimum AICc score
print(erato$analysis$model.name <-
        rownames(erato$analysis$AICcTable)[[which.min(erato$analysis$AICcTable$AICc)]])
# "modelI"

# Compute delta AICc (should be above 2)
deltaAICc<-sort(erato$analysis$AICcTable$AICc)[2]-sort(erato$analysis$AICcTable$AICc)[1]
deltaAICc
#4.111092

## Extract the hzar.dataGroup object for the selected model
erato$analysis$model.selected <-
  erato$analysis$oDG$data.groups[[erato$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(erato$analysis$model.selected,
                         names(erato$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(erato$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(erato$analysis$model.selected, ylim = c(-0.08,0.08));

# Open a pdf file
#pdf("clineplot.pdf") 

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(erato$analysis$model.selected, ylim = c(-0.08,0.08));

#Extract cline information for plotting
ca_abd_cline <- erato$analysis$model.selected

##Plot the two clines overlaid

grouped_clines <- list(fc_abd_cline, ca_abd_cline)

class(grouped_clines)

hzar.overPlot.fzCline(grouped_clines, fzClineSet = sapply(grouped_clines,
                                                          hzar.getCredParamRed,
                                                          simplify = FALSE))

#Plot the two clines separately
par(mfrow=c(1,2))
# Set graphical parameters
par(family = "serif", cex = 1, cex.axis = 0.8, cex.lab = 0.8, font = 2, font.axis = 1, font.lab = 1, mgp = c(2, 0.5, 0))

dev.off()
hzar.plot.fzCline(fc_abd_cline_all_pos,pch=3,size=0.6,xlab="Distance (km)", ylab="Residuals", col="black", ylim = c(-0.08, 0.08), fzCol = alpha("black", 0.2))
hzar.plot.fzCline(ca_abd_cline,pch=4,xlab="Distance (km)", ylab="Residuals", col="black", ylim = c(-0.08, 0.08), fzCol = alpha("black", 0.2))
legend("top", legend=c("Corneal area", "Facet count"), inset = c(0, -0.25), col=c("black", "red"), lty=1:1, xpd = TRUE)

hzar.plot.fzCline(ca_abd_cline, ylim = c(-0.08, 0.08), pch=1, main="Corneal area residuals", ylab="Log10(CA) residuals")
hzar.plot.fzCline(fc_abd_cline_all_pos, ylim = c(-0.08, 0.08), main="Facet count residuals", ylab="Log10(FC) residuals")


# 3. CLINE ANALYSIS FOR FACET COUNT (TIBIA RESIDUALS) ####

df <- read.csv("eye_data_with_residuals_tibia_only.csv",  fileEncoding="UTF-8-BOM")

# Initial graphical exploration of data
par(mfrow=c(1,2))

plot(log.facet.count ~ Transect.Position, data = df)

ggplot(df, aes(Transect.Position, df$FC_residuals, label=ID)) +
  geom_point() +
  geom_text(hjust = -0.2, size=2)

plot(FC_residuals ~ Elevation, data = df)
plot(CA_residuals ~ Elevation, data = df)

# Facet count dataframe
fc <- as.data.frame(df %>% group_by(location,Transect.Position) %>%
                      summarise(obs_mean = mean(FC_residuals),
                                obs_var = var(FC_residuals),
                                nSamples = n()) %>%
                      rename(locality = location,
                             distance = Transect.Position))

head(fc)
# Order by distance
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

# Set chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),
              B=c(528,124,978,544,99,596),
              C=c(124,978,544,99,596,528))

# Prepare an object holding everything
erato <- list()
## Space to hold the observed data
erato$obs <- list();
## Space to hold the models to fit
erato$models <- list();
## Space to hold the compiled fit requests
erato$fitRs <- list();
## Space to hold the output data chains
erato$runs <- list();
## Space to hold the analysed data
erato$analysis <- list();


# Generate a hzar.obsData object
erato$obs <- hzar.doNormalData1DPops(distance = as.list(fc$distance), muObs = as.list(fc$obs_mean),
                                     varObs = as.list(fc$obs_var), nEff = as.list(fc$nSamples),
                                     siteID=paste("P",1:length(as.list(fc$distance)),sep=""))

## Look at a graph of the observed data
hzar.plot.obsData(erato$obs, ylim = c(-0.07,0.07));
erato$obs$frame
## Make a helper function to define tails and name of models
erato.loadAdaAmodel <- function(tails,id=paste(tails))
  erato$models[[id]] <<- hzar.makeCline1DNormal(erato$obs, tails)

# Add five possible models
erato.loadAdaAmodel("none","modelI")
erato.loadAdaAmodel("right","modelII")
erato.loadAdaAmodel("left","modelIII")
erato.loadAdaAmodel("mirror","modelIV")
erato.loadAdaAmodel("both","modelV")


## Check the default settings
print(erato$models)

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km.

erato$models <- sapply(erato$models,
                       hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting 
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
erato$fitRs$init <- sapply(erato$models,
                           hzar.first.fitRequest.gC,
                           obsData=erato$obs,
                           verbose=FALSE,
                           simplify=FALSE)

## Check fit request settings
print(erato$fitRs$init)

## Run each model for an initial chain
erato$runs$init <- list()
erato$runs$init$modelI <- hzar.doFit(erato$fitRs$init$modelI)
erato$runs$init$modelII <- hzar.doFit(erato$fitRs$init$modelII)
erato$runs$init$modelIII <- hzar.doFit(erato$fitRs$init$modelIII)
erato$runs$init$modelIV <- hzar.doFit(erato$fitRs$init$modelIV)
erato$runs$init$modelV <- hzar.doFit(erato$fitRs$init$modelV)

## Plot the trace of a model
plot(hzar.mcmc.bindLL(erato$runs$init$modelI))

## Compile a new set of fit requests using the initial chains
erato$fitRs$chains <- lapply(erato$runs$init,hzar.next.fitRequest)


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
erato$fitRs$chains <- hzar.multiFitRequest(erato$fitRs$chains,
                                           each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])


##WARNING THIS TAKES A LONG TIME TO RUN###
# Run a chain of 3 runs for every fit request
erato$runs$chains <- hzar.doChain.multi(erato$fitRs$chains,
                                        doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Did modelII converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]]))))


## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
erato$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
erato$analysis$initDGs$modelI <-
  hzar.dataGroup.add(erato$runs$init$modelI)
erato$analysis$initDGs$modelII <-
  hzar.dataGroup.add(erato$runs$init$modelII)
erato$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(erato$runs$init$modelIII)
erato$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(erato$runs$init$modelIV)
erato$analysis$initDGs$modelV <-
  hzar.dataGroup.add(erato$runs$init$modelV)

## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
erato$analysis$oDG <-hzar.make.obsDataGroup(erato$analysis$initDGs)

# Use the same model names
erato$analysis$oDG <- hzar.copyModelLabels(erato$analysis$initDGs,
                                           erato$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
erato$analysis$oDG <-hzar.make.obsDataGroup(lapply(erato$runs$chains,
                                                   hzar.dataGroup.add),erato$analysis$oDG)

## Check to make sure that there are only 5 hzar.dataGroup
## in the hzar.obsDataGroup object.
print(summary(erato$analysis$oDG$data.groups))


## Do model selection based on the AICc scores
print(erato$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(erato$analysis$oDG));

#AICc
#modelI   -1024.946
#modelII  -1039.630
#modelIII -1020.759
#modelIV  -1021.550
#modelV   -1035.428

## Compare the 3 cline models to the null model graphically
hzar.plot.cline(erato$analysis$oDG, ylim = c(-0.08,0.08));
dev.off()
## Print out the model with the minimum AICc score
print(erato$analysis$model.name <-
        rownames(erato$analysis$AICcTable)[[which.min(erato$analysis$AICcTable$AICc)]])
# "modelI"

# Compute delta AICc (should be above 2)
deltaAICc<-sort(erato$analysis$AICcTable$AICc)[2]-sort(erato$analysis$AICcTable$AICc)[1]
deltaAICc
#4.201981

## Extract the hzar.dataGroup object for the selected model
erato$analysis$model.selected <-
  erato$analysis$oDG$data.groups[[erato$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(erato$analysis$model.selected,
                         names(erato$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(erato$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(erato$analysis$model.selected, ylim = c(-0.08,0.08));

# Open a pdf file
#pdf("clineplot.pdf") 
dev.off()
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(erato$analysis$model.selected, ylim = c(-0.08,0.08));

#Extract cline information into an object
fc_tib_cline01 <- erato$analysis$model.selected


# 4. CLINE ANALYSIS FOR CORNEAL AREA (TIBIA RESIDUALS) ####

# Facet count dataframe
fc <- as.data.frame(df %>% group_by(location,Transect.Position) %>%
                      summarise(obs_mean = mean(CA_residuals),
                                obs_var = var(CA_residuals),
                                nSamples = n()) %>%
                      rename(locality = location,
                             distance = Transect.Position))

head(fc)
# Order by distance
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

# Set chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),
              B=c(528,124,978,544,99,596),
              C=c(124,978,544,99,596,528))

# Prepare an object holding everything
erato <- list()
## Space to hold the observed data
erato$obs <- list();
## Space to hold the models to fit
erato$models <- list();
## Space to hold the compiled fit requests
erato$fitRs <- list();
## Space to hold the output data chains
erato$runs <- list();
## Space to hold the analysed data
erato$analysis <- list();


# Generate a hzar.obsData object (HERE I SHOULD GET MISSING TRANSECT VALUES)
erato$obs <- hzar.doNormalData1DPops(distance = as.list(fc$distance), muObs = as.list(fc$obs_mean),
                                     varObs = as.list(fc$obs_var), nEff = as.list(fc$nSamples),
                                     siteID=paste("P",1:length(as.list(fc$distance)),sep=""))

## Look at a graph of the observed data
hzar.plot.obsData(erato$obs, ylim = c(-0.08,0.08));
erato$obs$frame
## Make a helper function to define tails and name of models
erato.loadAdaAmodel <- function(tails,id=paste(tails))
  erato$models[[id]] <<- hzar.makeCline1DNormal(erato$obs, tails)

# Add five possible models
erato.loadAdaAmodel("none","modelI")
erato.loadAdaAmodel("right","modelII")
erato.loadAdaAmodel("left","modelIII")
erato.loadAdaAmodel("mirror","modelIV")
erato.loadAdaAmodel("both","modelV")


## Check the default settings
print(erato$models)

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km.
#erato$models <- sapply(erato$models,
#                       hzar.model.addBoxReq,300,1400,simplify=FALSE)

erato$models <- sapply(erato$models,
                       hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting 
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
erato$fitRs$init <- sapply(erato$models,
                           hzar.first.fitRequest.gC,
                           obsData=erato$obs,
                           verbose=FALSE,
                           simplify=FALSE)

#erato$fitRs$init <- sapply(erato$models,
#                           hzar.first.fitRequest.old.ML, obsData = erato$obs, 
#                           simplify=FALSE, verbose = FALSE)

## Check fit request settings
print(erato$fitRs$init)

## Run each model for an initial chain
erato$runs$init <- list()
erato$runs$init$modelI <- hzar.doFit(erato$fitRs$init$modelI)
erato$runs$init$modelII <- hzar.doFit(erato$fitRs$init$modelII)
erato$runs$init$modelIII <- hzar.doFit(erato$fitRs$init$modelIII)
erato$runs$init$modelIV <- hzar.doFit(erato$fitRs$init$modelIV)
erato$runs$init$modelV <- hzar.doFit(erato$fitRs$init$modelV)

## Plot the trace of a model
plot(hzar.mcmc.bindLL(erato$runs$init$modelI))

## Compile a new set of fit requests using the initial chains
erato$fitRs$chains <- lapply(erato$runs$init,hzar.next.fitRequest)


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
erato$fitRs$chains <- hzar.multiFitRequest(erato$fitRs$chains,
                                           each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])

##WARNING THIS TAKES A LONG TIME TO RUN###
# Run a chain of 3 runs for every fit request
erato$runs$chains <- hzar.doChain.multi(erato$fitRs$chains,
                                        doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Did modelII converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]]))))


## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
erato$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
erato$analysis$initDGs$modelI <-
  hzar.dataGroup.add(erato$runs$init$modelI)
erato$analysis$initDGs$modelII <-
  hzar.dataGroup.add(erato$runs$init$modelII)
erato$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(erato$runs$init$modelIII)
erato$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(erato$runs$init$modelIV)
erato$analysis$initDGs$modelV <-
  hzar.dataGroup.add(erato$runs$init$modelV)

## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
erato$analysis$oDG <-hzar.make.obsDataGroup(erato$analysis$initDGs)

# Use the same model names
erato$analysis$oDG <- hzar.copyModelLabels(erato$analysis$initDGs,
                                           erato$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
erato$analysis$oDG <-hzar.make.obsDataGroup(lapply(erato$runs$chains,
                                                   hzar.dataGroup.add),erato$analysis$oDG)

## Check to make sure that there are only 5 hzar.dataGroup
## in the hzar.obsDataGroup object.
print(summary(erato$analysis$oDG$data.groups))


## Do model selection based on the AICc scores
print(erato$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(erato$analysis$oDG));


## Compare the 3 cline models to the null model graphically
hzar.plot.cline(erato$analysis$oDG, ylim = c(-0.06,0.06));
dev.off()
## Print out the model with the minimum AICc score
print(erato$analysis$model.name <-
        rownames(erato$analysis$AICcTable)[[which.min(erato$analysis$AICcTable$AICc)]])
# "modelI"

# Compute delta AICc (should be above 2)
deltaAICc<-sort(erato$analysis$AICcTable$AICc)[2]-sort(erato$analysis$AICcTable$AICc)[1]
deltaAICc
#4.111092

## Extract the hzar.dataGroup object for the selected model
erato$analysis$model.selected <-
  erato$analysis$oDG$data.groups[[erato$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(erato$analysis$model.selected,
                         names(erato$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(erato$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(erato$analysis$model.selected, ylim = c(-0.06,0.06));

# Open a pdf file
#pdf("clineplot.pdf") 
dev.off()
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(erato$analysis$model.selected, ylim = c(-0.06,0.06));

#Extract cline information
ca_tib_cline01 <- erato$analysis$model.selected

#plot the two clines separately
hzar.plot.fzCline(fc_tib_cline01,pch=3,xlab="Distance (km)", ylab="Residuals", col="black", ylim = c(-0.08, 0.08), fzCol = alpha("red", 0.1))
hzar.plot.fzCline(ca_tib_cline01,pch=3,xlab="Distance (km)", ylab="Residuals", col = "black",  ylim = c(-0.08, 0.08), fzCol = alpha("black", 0.2))
legend("top", legend=c("Corneal area", "Facet count"), inset = c(0, -0.25), col=c("black", "red"), lty=1:1, xpd = TRUE)
dev.off()




# 5. CLINES FOR WING SHAPE TO COMPARE ####

# aspect ratio dataframe
fc <- as.data.frame(df %>% group_by(location,Transect.Position) %>%
                      summarise(obs_mean = mean(aspect.ratio),
                                obs_var = var(aspect.ratio),
                                nSamples = n()) %>%
                      rename(locality = location,
                             distance = Transect.Position))

head(fc)
# Order by distance
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

# Set chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),
              B=c(528,124,978,544,99,596),
              C=c(124,978,544,99,596,528))

# Prepare an object holding everything
erato <- list()
## Space to hold the observed data
erato$obs <- list();
## Space to hold the models to fit
erato$models <- list();
## Space to hold the compiled fit requests
erato$fitRs <- list();
## Space to hold the output data chains
erato$runs <- list();
## Space to hold the analysed data
erato$analysis <- list();


# Generate a hzar.obsData object (HERE I SHOULD GET MISSING TRANSECT VALUES)
erato$obs <- hzar.doNormalData1DPops(distance = as.list(fc$distance), muObs = as.list(fc$obs_mean),
                                     varObs = as.list(fc$obs_var), nEff = as.list(fc$nSamples),
                                     siteID=paste("P",1:length(as.list(fc$distance)),sep=""))

## Look at a graph of the observed data
hzar.plot.obsData(erato$obs, ylim = c(2.05,2.20));
erato$obs$frame
## Make a helper function to define tails and name of models
erato.loadAdaAmodel <- function(tails,id=paste(tails))
  erato$models[[id]] <<- hzar.makeCline1DNormal(erato$obs, tails)

# Add five possible models
erato.loadAdaAmodel("none","modelI")
erato.loadAdaAmodel("right","modelII")
erato.loadAdaAmodel("left","modelIII")
erato.loadAdaAmodel("mirror","modelIV")
erato.loadAdaAmodel("both","modelV")


## Check the default settings
print(erato$models)

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km.
erato$models <- sapply(erato$models,
                       hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting 
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
erato$fitRs$init <- sapply(erato$models,
                           hzar.first.fitRequest.gC,
                           obsData=erato$obs,
                           verbose=FALSE,
                           simplify=FALSE)

## Check fit request settings
print(erato$fitRs$init)

## Run each model for an initial chain
erato$runs$init <- list()
erato$runs$init$modelI <- hzar.doFit(erato$fitRs$init$modelI)
erato$runs$init$modelII <- hzar.doFit(erato$fitRs$init$modelII)
erato$runs$init$modelIII <- hzar.doFit(erato$fitRs$init$modelIII)
erato$runs$init$modelIV <- hzar.doFit(erato$fitRs$init$modelIV)
erato$runs$init$modelV <- hzar.doFit(erato$fitRs$init$modelV)

## Plot the trace of a model
plot(hzar.mcmc.bindLL(erato$runs$init$modelI))

## Compile a new set of fit requests using the initial chains
erato$fitRs$chains <- lapply(erato$runs$init,hzar.next.fitRequest)


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
erato$fitRs$chains <- hzar.multiFitRequest(erato$fitRs$chains,
                                           each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])

##WARNING THIS TAKES A LONG TIME TO RUN###
# Run a chain of 3 runs for every fit request
erato$runs$chains <- hzar.doChain.multi(erato$fitRs$chains,
                                        doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Did modelII converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]]))))


## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
erato$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
erato$analysis$initDGs$modelI <-
  hzar.dataGroup.add(erato$runs$init$modelI)
erato$analysis$initDGs$modelII <-
  hzar.dataGroup.add(erato$runs$init$modelII)
erato$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(erato$runs$init$modelIII)
erato$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(erato$runs$init$modelIV)
erato$analysis$initDGs$modelV <-
  hzar.dataGroup.add(erato$runs$init$modelV)

## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
erato$analysis$oDG <-hzar.make.obsDataGroup(erato$analysis$initDGs)

# Use the same model names
erato$analysis$oDG <- hzar.copyModelLabels(erato$analysis$initDGs,
                                           erato$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
erato$analysis$oDG <-hzar.make.obsDataGroup(lapply(erato$runs$chains,
                                                   hzar.dataGroup.add),erato$analysis$oDG)

## Check to make sure that there are only 5 hzar.dataGroup
## in the hzar.obsDataGroup object.
print(summary(erato$analysis$oDG$data.groups))


## Do model selection based on the AICc scores
print(erato$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(erato$analysis$oDG));


## Compare the 3 cline models to the null model graphically
hzar.plot.cline(erato$analysis$oDG);
dev.off()
## Print out the model with the minimum AICc score
print(erato$analysis$model.name <-
        rownames(erato$analysis$AICcTable)[[which.min(erato$analysis$AICcTable$AICc)]])
# "modelI"

# Compute delta AICc (should be above 2)
deltaAICc<-sort(erato$analysis$AICcTable$AICc)[2]-sort(erato$analysis$AICcTable$AICc)[1]
deltaAICc
#4.111092

## Extract the hzar.dataGroup object for the selected model
erato$analysis$model.selected <-
  erato$analysis$oDG$data.groups[[erato$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(erato$analysis$model.selected,
                         names(erato$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(erato$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(erato$analysis$model.selected);

# Open a pdf file
#pdf("clineplot.pdf") 
dev.off()
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(erato$analysis$model.selected);

#Extract cline information and plot
aspect_ratio_cline <- erato$analysis$model.selected

hzar.plot.fzCline(aspect_ratio_cline,pch=3,xlab="Distance (km)", ylab="Residuals", col="black", ylim = c(-0.08, 0.08), fzCol = alpha("black", 0.1))

# 6. CLINE FOR ADMIXTURE (FROM NGSADMIX) ####

#First we need to flip the admixture values around so that they match the cline direction of eye size
#This should only change the direction of the cline.
df$ngs.admix <- 1-df$ngs.admix 

# aspect ratio dataframe
fc <- as.data.frame(df %>% group_by(location,Transect.Position) %>%
                      summarise(obs_mean = mean(ngs.admix),
                                obs_var = var(ngs.admix),
                                nSamples = n()) %>%
                      rename(locality = location,
                             distance = Transect.Position))

head(fc)
# Order by distance
fc <- na.omit(fc)
fc <- fc[order(-fc$distance), ]

# Set chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(978,544,99,596,528,124), 
              B=c(544,99,596,528,124,978), 
              C=c(99,596,528,124,978,544), 
              D=c(978,99,596,528,124,544), 
              E=c(124,596,978,99,528,544))

# Prepare an object holding everything
erato <- list()
## Space to hold the observed data
erato$obs <- list();
## Space to hold the models to fit
erato$models <- list();
## Space to hold the compiled fit requests
erato$fitRs <- list();
## Space to hold the output data chains
erato$runs <- list();
## Space to hold the analysed data
erato$analysis <- list();


# Generate a hzar.obsData object (HERE I SHOULD GET MISSING TRANSECT VALUES)
erato$obs <- hzar.doNormalData1DPops(distance = as.list(fc$distance), muObs = as.list(fc$obs_mean),
                                     varObs = as.list(fc$obs_var), nEff = as.list(fc$nSamples),
                                     siteID=paste("P",1:length(as.list(fc$distance)),sep=""))

## Look at a graph of the observed data
hzar.plot.obsData(erato$obs);
erato$obs$frame

## Make a helper function to define tails and name of models
erato.loadAdaAmodel <- function(tails,id=paste(tails))
  erato$models[[id]] <<- hzar.makeCline1DNormal(erato$obs, tails)

# Add five possible models
erato.loadAdaAmodel("none","modelI")
erato.loadAdaAmodel("right","modelII")
erato.loadAdaAmodel("left","modelIII")
erato.loadAdaAmodel("mirror","modelIV")
erato.loadAdaAmodel("both","modelV")


## Check the default settings
print(erato$models)

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km.

erato$models <- sapply(erato$models,
                       hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting 
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
erato$fitRs$init <- sapply(erato$models,
                           hzar.first.fitRequest.gC,
                           obsData=erato$obs,
                           verbose=FALSE,
                           simplify=FALSE)

#erato$fitRs$init <- sapply(erato$models,
#                           hzar.first.fitRequest.old.ML, obsData = erato$obs, 
#                           simplify=FALSE, verbose = FALSE)

## Check fit request settings
print(erato$fitRs$init)

## Run each model for an initial chain
erato$runs$init <- list()
erato$runs$init$modelI <- hzar.doFit(erato$fitRs$init$modelI)
erato$runs$init$modelII <- hzar.doFit(erato$fitRs$init$modelII)
erato$runs$init$modelIII <- hzar.doFit(erato$fitRs$init$modelIII)
erato$runs$init$modelIV <- hzar.doFit(erato$fitRs$init$modelIV)
erato$runs$init$modelV <- hzar.doFit(erato$fitRs$init$modelV)

## Plot the trace of a model
plot(hzar.mcmc.bindLL(erato$runs$init$modelI))

## Compile a new set of fit requests using the initial chains
erato$fitRs$chains <- lapply(erato$runs$init,hzar.next.fitRequest)


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
erato$fitRs$chains <- hzar.multiFitRequest(erato$fitRs$chains,
                                           each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(erato$fitRs$chains), function(x) erato$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])

##WARNING THIS TAKES A LONG TIME TO RUN###
# Run a chain of 3 runs for every fit request
erato$runs$chains <- hzar.doChain.multi(erato$fitRs$chains,
                                        doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Did modelII converge? ---> YES
summary(do.call(mcmc.list,
                lapply(erato$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]]))))


## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
erato$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
erato$analysis$initDGs$modelI <-
  hzar.dataGroup.add(erato$runs$init$modelI)
erato$analysis$initDGs$modelII <-
  hzar.dataGroup.add(erato$runs$init$modelII)
erato$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(erato$runs$init$modelIII)
erato$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(erato$runs$init$modelIV)
erato$analysis$initDGs$modelV <-
  hzar.dataGroup.add(erato$runs$init$modelV)

## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
erato$analysis$oDG <-hzar.make.obsDataGroup(erato$analysis$initDGs)

# Use the same model names
erato$analysis$oDG <- hzar.copyModelLabels(erato$analysis$initDGs,
                                           erato$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
erato$analysis$oDG <-hzar.make.obsDataGroup(lapply(erato$runs$chains,
                                                   hzar.dataGroup.add),erato$analysis$oDG)

## Check to make sure that there are only 5 hzar.dataGroup
## in the hzar.obsDataGroup object.
print(summary(erato$analysis$oDG$data.groups))


## Do model selection based on the AICc scores
print(erato$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(erato$analysis$oDG));


## Compare the 3 cline models to the null model graphically
hzar.plot.cline(erato$analysis$oDG);
dev.off()
## Print out the model with the minimum AICc score
print(erato$analysis$model.name <-
        rownames(erato$analysis$AICcTable)[[which.min(erato$analysis$AICcTable$AICc)]])
# "modelI"

# Compute delta AICc (should be above 2)
deltaAICc<-sort(erato$analysis$AICcTable$AICc)[2]-sort(erato$analysis$AICcTable$AICc)[1]
deltaAICc
#4.111092

## Extract the hzar.dataGroup object for the selected model
erato$analysis$model.selected <-
  erato$analysis$oDG$data.groups[[erato$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(erato$analysis$model.selected,
                         names(erato$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(erato$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(erato$analysis$model.selected);

# Open a pdf file
#pdf("clineplot.pdf") 
dev.off()
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(erato$analysis$model.selected);

#Extract cline information
admix_cline <- erato$analysis$model.selected

# Set graphical parameters
par(family = "serif", cex = 1.1, cex.axis = 1, cex.lab = 1.2, font = 2, font.axis = 1, font.lab = 1, mgp = c(2, 0.5, 0))

#plot
hzar.plot.fzCline(admix_cline,pch=4,xlab="Distance (km)", ylab="Admixture proportion", col = "black", fzCol = alpha("black", 0.2))

#end

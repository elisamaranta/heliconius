# Analysis of eye size across Heliconius erato hybrid zone
# Elisa Mogollon Perez, Septemeber 2024.

setwd("C:/Users/USER/OneDrive - CACTUS/Documents/OneDrive/EES Year 2/IRT3/Measurement data")
getwd()

#load packages
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(rstatix)
library(effects)
library(data.table)
library(tidyverse)
library(lme4)
library(nlme)
library(gridExtra)
library(patchwork)
library(emmeans)
library(scales)

#color selection for consistenty with other plots
library(RColorBrewer)

#theme settings for plots
require(grid)
pub<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour = NA,fill=NA),
  panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  axis.title.x=element_text(size=15,face="bold",hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=15,face="bold",hjust=0.5,vjust=1,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=15),
  axis.text.y=element_text(colour="black",angle=0,size=15),
  axis.ticks=element_line(colour="black"))


#import data
eyes <- read.csv("NEW_Lativ_Notab_eye_morphology_data.csv",  fileEncoding="UTF-8-BOM")
head(eyes)

#Remove last two columns
long <- eyes[c(1:20)]

######create new composite variables of facet count & corneal area####
#when left eye was missing or damaged, right is used
long$cornea.area <- ifelse(is.na(long$L_area), long$R_area, long$L_area)
long$facet.count <- ifelse(is.na(long$L_whole_count), long$R_whole_count, long$L_whole_count)

#visually inspect distribution
par(mfrow=c(3,2))

#log-transformed variables
long$log.facet.count <- log(long$facet.count)
long$log.cornea.area <- log(long$cornea.area)
long$log.tibia.length <- log(long$tibia_length)
long$log.abdomen <- log(long$abdomen_length)
long$log.head <- log(long$inter_eye_width)
long$log.facet.area <- log(long$cornea.area/long$facet.count)

#plot
hist(long$log.facet.count)
hist(long$log.cornea.area)
hist(long$log.tibia.length)
hist(long$log.abdomen)
hist(long$log.head)
hist(long$log.facet.area)

#sanity checks for correlation between left and right side
ggplot(long, aes(x=L_whole_count, y=R_whole_count)) + 
  geom_smooth(method = 'lm') + geom_point()  + theme(legend.position = "top")
cor.test(long$L_whole_count,long$R_whole_count)
#cor = 0.95

ggplot(long, aes(x=L_area, y=R_area)) + 
  geom_smooth(method = 'lm') + geom_point()  + theme(legend.position = "top")
cor.test(long$L_area,long$R_area)
#cor = 0.98

## 1) Initial exploarion analysis ####

print(unique(long$type))
#rename hybrids
long <- long %>% mutate(type = dplyr::recode(type, "notabilisxlativitta" = "hybrid"))
head(long)

#filter single missing data sample
long <- long %>% filter(!is.na(long$sex),
                        type != "checkdata base")

#How many individuals per subspp.?
counts <- long %>% count(type)
counts

#    hybrid   f 130
#    hybrid   m 291
# lativitta   f   7
# lativitta   m  22
# notabilis   f  15
# notabilis   m  19

#Where do we find hybrids?
long %>% group_by(type) %>% summarise(unique_sites = n_distinct(location, na.rm=TRUE))

#type      unique_sites
#<chr>            <int>
#  1 hybrid             33
#2 lativitta            3
#3 notabilis            3
head(long)

#Get sites samples as csv
sites <- long %>% group_by(location, latitude, longitude, elevation) %>% summarise(individuals = n_distinct(ID, na.rm=TRUE))

sites <- sites[order(-sites$elevation), ]
#write.csv(sites, "sites_samples.csv")

hyb <- subset(long, long$type == "hybrid")

max(hyb$elevation)
min(hyb$elevation)

#What are the differences between species?
long %>% group_by(type) %>% summarise(mean_CA = mean(cornea.area, na.rm=TRUE),
                                      sd_CA = sd(cornea.area, na.rm=TRUE),
                                      mean_FC = mean(facet.count, na.rm=TRUE),
                                      sd_FC = sd(facet.count, na.rm=TRUE),
                                      mean_abd = mean(abdomen_length, na.rm=TRUE),
                                      sd_abd = sd(abdomen_length, na.rm=TRUE),
                                      mean_head = mean(inter_eye_width, na.rm=TRUE))

#type      mean_CA sd_CA mean_FC sd_FC mean_abd sd_abd
#<chr>       <dbl> <dbl>   <dbl> <dbl>    <dbl>  <dbl>
#  1 hybrid     7.76 0.891  13845. 1045.     16.5   1.20
#2 lativitta    8.07 0.832  14124. 1171.     16.6   1.06
#3 notabilis    7.25 0.633  13298.  794.     16.6   1.01


#Plot allometric relationships

#log abdomen x log tibia by sex
p1 <- ggplot(long, aes(x=log.abdomen, y=log.tibia.length, color=sex)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = sex)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log hind tibia length") + xlab("log abdomen length") + 
  theme(legend.text=element_text(size=12),)+
  theme(text = element_text(family = "Times New Roman", size = 14))

p1

#log corneal area x log tibia by sex
p1 <- ggplot(long, aes(x=log.tibia.length, y=log.cornea.area, color=type)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = type)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  #scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  scale_colour_manual(values=c("#A6CEE3", "#1F78B4", "#61D04F")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("log hind tibia length") + ylab("log cornea area") + 
  theme(legend.text=element_text(size=12),)+
  theme(text = element_text(family = "Times New Roman", size = 14))

p1


# Corneal area x log tibia by sex

# Custom facet labels
facet_labels <- c(m = "Males", f = "Females")

# Define a lighter color palette using RColorBrewer
light_palette <- rev(brewer.pal(n = 9, name = "RdYlGn"))  
brewer.pal(7, "RdYlGn")

custom_colors <- rev(c("#FFFF00", "#66BD63","#006837")) #"#FFFFBF"

p1 <- ggplot(long, aes(x=log.tibia.length, y=log.cornea.area, color=elevation)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = sex)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  #scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  scale_color_gradientn(colors = custom_colors, name = "Elevation (m)", oob = scales::squish) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log corneal area") + xlab("log hind tibia length") + labs(color="Elevation")+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14))+
  theme(text = element_text(family = "Times New Roman", size = 14)) + facet_wrap(~sex, labeller = as_labeller(facet_labels))

p1

#corneal area x abdomen length
p1 <- ggplot(long, aes(x=log.abdomen, y=log.facet.count, color=elevation)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = sex)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  #scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  scale_color_gradientn(colors = custom_colors, name = "Elevation (m)", oob = scales::squish) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log facet count") + xlab("log abdomen length") + labs(color="Elevation")+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14))+
  theme(text = element_text(family = "Times New Roman", size = 14)) + facet_wrap(~sex, labeller = as_labeller(facet_labels))

p1

#corneal area x head size
p1 <- ggplot(long, aes(x=log.head, y=log.cornea.area, color=type)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = type)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  #scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("Log corneal area") + xlab("Log inter-ocular distance") + labs(color="Elevation")+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14))+
  theme(text = element_text(family = "Times New Roman", size = 14))

p1

#model corneal area and head size * type/elevation interaction

m1 <- lm(log.cornea.area ~ log.head * elevation, data=long)

plot(m1, which = 2)
summary(m1) #no sig. interaction. 

# facet count x log tibia by sex
p1 <- ggplot(long, aes(x=log.tibia.length, y=log.facet.count, color=elevation)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = sex)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  #scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  scale_color_gradientn(colors = custom_colors, name = "Elevation (m)", oob = scales::squish) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log facet count") + xlab("log hind tibia length") + labs(color="Elevation")+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size=14))+
  theme(text = element_text(family = "Times New Roman", size = 14)) + facet_wrap(~sex, labeller = as_labeller(facet_labels))
p1

# facet count x corneal area
p1 <- ggplot(long, aes(x=log.cornea.area, y=log.facet.count, color=sex)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = sex)) + theme_classic()+
  #geom_violin(position=position_dodge(width=1)) +
  #geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  #geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  scale_colour_manual(labels = c("f"="Females","m"="Males"), values=c("#F8766D", "#00BFC4")) +
  #scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log facet count") + xlab("log corneal area") + 
  theme(legend.text=element_text(size=12),)+
  theme(text = element_text(family = "Times New Roman", size = 14))

p1

#what are the correlation?
cor.test(df$log.cornea.area, df$log.facet.count, method = "pearson")

cor.test(df$log.abdomen, df$log.tibia.length, method = "pearson")

cor.test(df$log.cornea.area, df$log.abdomen, method = "pearson")

cor.test(df$log.facet.count, df$log.abdomen, method = "pearson")

cor.test(df$log.cornea.area, df$log.tibia.length, method = "pearson")

cor.test(df$log.facet.count, df$log.tibia.length, method = "pearson")



# abdomen size by type or sex
p1 <- ggplot(long, aes(x=type, y=log.abdomen, fill=type)) +
  geom_violin(position=position_dodge(width=1)) + theme_classic()+
  geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A")) +
  scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log10 abdomen length") + xlab("")+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size = 14))

p1

# abdomen size by type or sex
p1 <- ggplot(long, aes(x=type, y=log.tibia.length, fill=type)) +
  geom_violin(position=position_dodge(width=1)) + theme_classic()+
  geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A")) +
  scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("log10 tibia length") + xlab("")+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size = 14))

p1

## 2) Are there significant differences in body measurements? #####

### Plot raw data
dev.off()
hist(long$log.abdomen)

### Abdomen length model
m1 <- lm(log.abdomen ~ type * sex + elevation, data = long)

#mixed model 
mixed <- lmer(log.abdomen ~  type * sex + elevation + (1 | location), data = long, REML = TRUE)
lmerTest::ranova(mixed) #location does not have a significant effect on the model.

plot(m1, which = 1)
plot(m1, which = 2)
hist(resid(m1))

#Full model
m1 <- lm(log.abdomen ~ type * sex * elevation, data = long)

drop1(m1, test = "Chisq") #type:sex:elevation is ns
m1 <- update(m1,.~.-type:sex:elevation)
drop1(m1, test = "Chisq") #type:elevation is ns
m1 <- update(m1,.~.-type:elevation) 
drop1(m1, test = "Chisq") #sex:elevation is ns
m1 <- update(m1,.~.-sex:elevation) 
drop1(m1, test = "Chisq") #type:sex is ns
m1 <- update(m1,.~.-type:sex) 
drop1(m1, test="Chisq") #elevation is ns
m1 <- update(m1,.~.-elevation) 
drop1(m1, test="Chisq") #type is ns
m1 <- update(m1,.~.-type)
drop1(m1, test="Chisq") #sex IS SIGNIFICANT

plot(allEffects(m1))
plot(m1, which = 2) 
plot(m1, which = 1)
hist(resid(m1))

summary(m1) #sex is the only sig. variable for abdomen length
car::Anova(m1, test.statistic="Chisq", method = "fdr")
# sex   p=3.684e-05 ***

#effect means
model_means <- emmeans(object = m1, specs = ~ sex, type = "response") 
model_means

#p-value comparison table
pwpm(model_means, adjust="fdr",  diffs = F)

#effect size
eff_size(model_means, sigma = sigma(m1), edf = df.residual(m1))

#add letters to each mean to indicate difference and plot
library(multcomp)
library(multcompView)

model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means as a dataframe
facet <- as.data.frame(model_means_cld)

facet$sex <- factor(facet$sex, c("male", "female"))

#plot sex effects
facet %>% ggplot(aes(sex,facet$emmean))+
  geom_point(aes(color=sex), size=4) + geom_errorbar(aes(x=sex, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10 (abdomen length)")+
  scale_x_discrete(labels=c("female", "male"))+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))



## Tibia length model
m1 <- lm(log.tibia.length ~ type * sex * elevation, data = long)

#mixed model 
model <- lmer(log.tibia.length ~  type * sex * elevation + (1 | location), data = long, REML = TRUE)
lmerTest::ranova(model)

plot(m1, which = 1)
plot(m1, which = 2)
hist(resid(m1))

drop1(m1, test = "Chisq") #type:sex:elevation is ns
m1 <- update(m1,.~.-type:sex:elevation)
drop1(m1, test = "Chisq") #type:elevation is ns
m1 <- update(m1,.~.-type:sex) 
drop1(m1, test = "Chisq") #sex:elevation is ns
m1 <- update(m1,.~.-sex:elevation) 
drop1(m1, test = "Chisq") #type:sex is ns
m1 <- update(m1,.~.-type:elevation) 
drop1(m1, test="Chisq") #elevation is ns
m1 <- update(m1,.~.-type) 
drop1(m1, test="Chisq") #type is ns
m1 <- update(m1,.~.-elevation)
drop1(m1, test="Chisq") #sex is ns 

plot(m1, which = 2) 
hist(resid(m1))

summary(m1) #sex is ns, and low R-squared
car::Anova(m1, test.statistic="Chisq")

#effect means
model_means <- emmeans(object = m1, specs = ~ sex, type = "response") 
model_means

#p-value comparison table
pwpm(model_means, adjust="fdr",  diffs = F)

#effect size
eff_size(model_means, sigma = sigma(m1), edf = df.residual(m1))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means at a dataframe
facet <- as.data.frame(model_means_cld)
facet
#species
#color selection for consistenty with other papers
library(RColorBrewer)
display.brewer.pal(12, 'Paired')
brewer.pal(12, 'Paired')

facet$sex <- factor(facet$sex, c("male", "female"))

facet %>% ggplot(aes(sex,facet$emmean))+
  geom_point(aes(color=sex), size=4) + geom_errorbar(aes(x=sex, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10 (tibia length)")+
  scale_x_discrete(labels=c("female", "male"))+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## Inter-ocular distance
m1 <- lm(log.head ~ type * sex * elevation, data = long)

#mixed model 
model <- lmer(log.head ~  type * sex * elevation + (1 | location), data = long, REML = TRUE)
lmerTest::ranova(model) #location effect is not significant.
anova(model, m1)

plot(m1, which = 1)
plot(m1, which = 2)
hist(resid(m1))

drop1(m1, test = "Chisq") #type:sex:elevation is ns
m1 <- update(m1,.~.-type:sex:elevation)
drop1(m1, test = "Chisq") #type:elevation is ns
m1 <- update(m1,.~.-type:sex) 
drop1(m1, test = "Chisq") #sex:elevation is ns
m1 <- update(m1,.~.-type:elevation) 
drop1(m1, test = "Chisq") #type:sex is ns
m1 <- update(m1,.~.-sex:elevation) 
drop1(m1, test="Chisq") #elevation is ns
m1 <- update(m1,.~.-sex) 
drop1(m1, test="Chisq") #type is ns
m1 <- update(m1,.~.-elevation)
drop1(m1, test="Chisq") #type IS SIGNIFICANT

plot(m1, which = 2) 
hist(resid(m1))

summary(m1)
car::Anova(m1, test.statistic="Chisq", method = "fdr")

#pairwise comparisons
model_means <- emmeans(object = m1, specs = ~ type, type = "response") 
model_means

contrast(model_means)
#                  estimate      SE  df t.ratio p.value
#hybrid effect    -0.02014 0.00584 433  -3.445  0.0019
#lativitta effect  0.00469 0.00899 433   0.522  0.6021
#notabilis effect  0.01544 0.00831 433   1.858  0.0958

#p-value comparison table
pwpm(model_means, adjust="fdr",  diffs = F)

#           hybrid lativitta notabilis
#hybrid    [0.327]    0.1531    0.0048 
#lativitta           [0.352]    1.0000
#notabilis                     [0.363]

#effect size
eff_size(model_means, sigma = sigma(m1), edf = df.residual(m1))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="fdr",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means at a dataframe
facet <- as.data.frame(model_means_cld)
facet
#species
#color selection for consistenty with other papers
library(RColorBrewer)
display.brewer.pal(12, 'Paired')
brewer.pal(12, 'Paired')

#plot
facet %>% ggplot(aes(type,facet$emmean))+
  geom_point(aes(color=type), size=4) + geom_errorbar(aes(x=type, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10\n (inter-ocular distance)")+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("subspecies effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## 3) What is the relationship between facet count and corneal area?####

d1 <- lm(log.cornea.area~log.facet.count*type*sex+elevation, data=long)

#mixed model 
model <- lmer(log.cornea.area~log.facet.count*type*sex+elevation+ (1 | location), data = long, REML = TRUE)
lmerTest::ranova(model) #location does not have a significant effect

anova(model, d1)

par(mfrow=c(2,2)) 
plot(d1)
hist(resid(d1))

drop1(d1, test="Chisq")
d2 <- update(d2,.~.-log.facet.count:type:sex)
drop1(d2, test="Chisq")
d3 <- update(d2,.~.-log.facet.count:sex:elevation)
drop1(d3, test="Chisq")
d4 <- update(d3,.~.-log.facet.count:sex)
drop1(d4, test="Chisq")
d5 <- update(d4,.~.-type:sex:elevation)
drop1(d5, test="Chisq")
d6 <- update(d5,.~.-log.facet.count:type:elevation)
drop1(d6, test="Chisq")
d7 <- update(d6,.~.-sex:type)
drop1(d7, test="Chisq") #only log.facet.count and elevation explain log.cornea.area
d8 <- update(d7,.~.-type:elevation)
drop1(d8, test="Chisq")
d9 <- update(d8,.~.-type:sex)
drop1(d9, test="Chisq")
d10 <- update(d9,.~.-log.facet.count:type)
drop1(d10, test="Chisq")
d11 <- update(d10,.~.-log.facet.count:elevation)
drop1(d11, test="Chisq")
d12 <- update(d11,.~.-type)
drop1(d12, test="Chisq")
d13 <- update(d12,.~.-sex:elevation)
drop1(d13, test="Chisq")
d14 <- update(d13,.~.-sex)
drop1(d14, test="Chisq") #log.facet.count and elevation are both significant

plot(d14)
hist(resid(d14))
dev.off()

car::Anova(d14, test="Chisq")
#log.facet.count 4.6301   1 1563.4293 < 2.2e-16 ***
#elevation       0.0295   1    9.9718  0.001696 **

### DO THE SPECIES HAVE DIFFERENT EYE SIZES? ###

#Cornea area
d1 <- lm(log.cornea.area~log.abdomen+type*sex, data=long)

#Here we cannot test site as a random effect because it acts as an 'elevation' effect if we don't include elevation.

par(mfrow=c(2,2)) 
plot(d1)
hist(resid(d1)) #normal

drop1(d1, test="Chisq")
d2 <- update(d1,.~.-type:sex)
drop1(d2, test="Chisq")  #all three sign, no interactions.

plot(d2)
hist(resid(d2))
dev.off()

plot(allEffects(d2))

car::Anova(d2, test="Chisq", method = "fdr")
#Response: log.cornea.area
#Sum Sq  Df F value    Pr(>F)    
#log.abdomen 2.56362   1 490.961 < 2.2e-16 ***
#  type        0.12061   2  11.549 1.293e-05 ***
#  sex         0.58568   1 112.163 < 2.2e-16 ***
#  Residuals   2.29752 440                      

#effect means
model_means <- emmeans(object = d2, specs = ~ type, type = "response") 
model_means

#p-value comparison table
pwpm(model_means, adjust="fdr",  diffs = F)

#lativitta notabilis notabilisxlativitta
#lativitta             [2.053]    <.0001              0.0942
#notabilis                       [1.969]              <.0001
#notabilisxlativitta                                 [2.028]

#effect size
eff_size(model_means, sigma = sigma(m1), edf = df.residual(m1))

## facet count
d1 <- lm(log.facet.count~log.abdomen+type*sex, data=long)
par(mfrow=c(2,2)) 
plot(d1)
hist(resid(d1)) #normal

drop1(d1, test="Chisq") #type*sex interaction is significant.
d2 <- update(d1,.~.-type:sex)
drop1(d2, test="Chisq")  #all three sign, no interactions.

plot(d2)
hist(resid(d2))
dev.off()

plot(allEffects(d2))

car::Anova(d2, test="Chisq", method = "fdr")
#Response: log.cornea.area
#Sum Sq  Df F value    Pr(>F)    
#log.tibia.length 2.1618   1  445.19 < 2.2e-16 ***
#  type             0.1009   2   10.39 4.038e-05 ***
#  sex              0.7544   1  155.36 < 2.2e-16 ***
#  Residuals        1.8549 382                      


#effect means
model_means <- emmeans(object = d2, specs = ~ type, type = "response") 
model_means

#p-value comparison table
pwpm(model_means, adjust="fdr",  diffs = F)

#lativitta notabilis notabilisxlativitta
#lativitta             [9.529]    0.0056              0.4941
#notabilis                       [9.487]              0.0017
#notabilisxlativitta                                 [9.521]

#effect size
eff_size(model_means, sigma = sigma(m1), edf = df.residual(m1))

## 4) How does eye size vary with elevation? USING ABDOMEN####

## Full model for facet count
m1 <- lm(log.facet.count ~ log.abdomen * sex * elevation, data = long)

#mixed model 
model <- lmer(log.facet.count ~ log.abdomen * sex * elevation + (1 | location), data = long, REML = TRUE)
lmerTest::ranova(model) #location has no effect

plot(m1, which=1)
plot(m1, which=2)
hist(resid(m1))

drop1(m1, test="Chisq") #abdomen*sex*elevation is ns
m2 <- update(m1,.~.-log.abdomen:sex:elevation) 
drop1(m2, test="Chisq") #sex*elevation is ns
m2 <- update(m2,.~.-sex:elevation) 
drop1(m2, test="Chisq") #abdomen*sex is ns
m2 <- update(m2,.~.-log.abdomen:sex) 
drop1(m2, test="Chisq") #abdomen*elevation is ns
m2 <- update(m2,.~.-log.abdomen:elevation) 
drop1(m2, test="Chisq") #all significant

library(effects)
plot(allEffects(m2))
plot(m2, which = 2) 
hist(resid(m2))

Anova(m2, test.statistic="Chisq") 
summary(m2)

car::Anova(m2, test.statistic="Chisq")
#Response: log.facet.count
#Sum Sq  Df F value    Pr(>F)    
#log.abdomen 0.70454   1 249.222 < 2.2e-16 ***
#  sex         0.33560   1 118.714 < 2.2e-16 ***
#  elevation   0.06295   1  22.269 3.187e-06 ***
#  Residuals   1.24670 441                      

## Sex effect

#pairwise comparisons
model_means <- emmeans(object = m2, specs = ~ sex, type = "response") 
model_means
#p-value comparison table
pwpm(model_means, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means, sigma = sigma(m2), edf = df.residual(m2))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means at a dataframe
facet <- as.data.frame(model_means_cld)
facet
#species
#color selection for consistenty with other papers
library(RColorBrewer)
display.brewer.pal(12, 'Paired')
brewer.pal(12, 'Paired')

#plot sex effect
facet %>% ggplot(aes(sex,facet$emmean))+
  geom_point(aes(color=sex), size=4) + geom_errorbar(aes(x=sex, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10 (facet count)")+
  scale_x_discrete(labels=c("female", "male"))+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size=18))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## Elevation effect

coef(summary(m2))[ , "Estimate"]

el <- Effect("elevation", m2)
el <- data.frame(Effect("elevation", m2))

ggplot(data = el, aes(x=elevation, y=fit, ymin=lower, ymax=upper))+
  geom_line(col="#1F78B4", size=1) + geom_ribbon(alpha=0.4, fill="#1F78B4") + ylab("e.m.m. log10 (facet count)") + 
  xlab("Elevation") + ggtitle("elevation effect")+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 2)) + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## Full model for CORNEAL AREA

m1 <- lm(log.cornea.area ~ log.abdomen + sex * elevation, data = long)

#mixed model 
model <- lmer(log.cornea.area ~ log.abdomen + sex * elevation + (1 | location), data = long, REML =TRUE)
lmerTest::ranova(model)

plot(m1, which=1)
plot(m1, which=2)
hist(resid(m1))

drop1(m1, test="Chisq") #sex:elevation is ns
m2 <- update(m1,.~.-sex:elevation) 
drop1(m2, test="Chisq") #all significant

plot(allEffects(m2))
plot(m2, which = 2) 
hist(resid(m2))

summary(m2)
car::Anova(m2, test.statistic="Chisq")


## Sex effect
model_means <- emmeans(object = m2, specs = ~ sex, type = "response") 
model_means
#p-value comparison table
pwpm(model_means, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means, sigma = sigma(m2), edf = df.residual(m2))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means at a dataframe
facet <- as.data.frame(model_means_cld)
facet

#plots sex effect
facet %>% ggplot(aes(sex,facet$emmean))+
  geom_point(aes(color=sex), size=4) + geom_errorbar(aes(x=sex, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10 (corneal area)")+
  scale_x_discrete(labels=c("female", "male"))+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size=18))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## Elevation effect
coef(summary(m2))[ , "Estimate"]

el <- Effect("elevation", m2)
el <- data.frame(Effect("elevation", m2))

ggplot(data = el, aes(x=elevation, y=fit, ymin=lower, ymax=upper))+
  geom_line(col="#1F78B4", size=1) + geom_ribbon(alpha=0.4, fill="#1F78B4") + ylab("e.m.m. log10 (corneal area)") + 
  xlab("Elevation") + ggtitle("elevation effect")+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 2)) + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## 5) How does eye size vary with elevation? USING TIBIA ####

#Full model for FACET COUNT

m1 <- lm(log.facet.count ~ log.tibia.length + sex * elevation, data = long)

plot(m1, which=1)
plot(m1, which=2)
hist(resid(m1))

drop1(m1, test="Chisq") #sex:elevation is ns
m2 <- update(m1,.~.-sex:elevation) 
drop1(m2, test="Chisq") #all significant

plot(allEffects(m2))
plot(m2, which = 2) 
hist(resid(m2))

Anova(m2, test.statistic="Chisq") 
summary(m2)

car::Anova(m2, test.statistic="fdr")

## Sex effect

#pairwise comparisons
model_means <- emmeans(object = m2, specs = ~ sex, type = "response") 
model_means
#p-value comparison table
pwpm(model_means, adjust="bonferroni",  diffs = F)
#f       m
#f [9.483]  <.0001
#m         [9.553]

#effect size
eff_size(model_means, sigma = sigma(m2), edf = df.residual(m2))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means at a dataframe
facet <- as.data.frame(model_means_cld)
facet

# Plot sex effect
facet %>% ggplot(aes(sex,facet$emmean))+
  geom_point(aes(color=sex), size=4) + geom_errorbar(aes(x=sex, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10 (facet count)")+
  scale_x_discrete(labels=c("female", "male"))+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size=18))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## Elevation effect

coef(summary(m2))[ , "Estimate"]

el <- Effect("elevation", m2)
el <- data.frame(Effect("elevation", m2))

ggplot(data = el, aes(x=elevation, y=fit, ymin=lower, ymax=upper))+
  geom_line(col="#1F78B4", size=1) + geom_ribbon(alpha=0.4, fill="#1F78B4") + ylab("e.m.m. log10 (facet count)") + 
  xlab("Elevation") + ggtitle("elevation effect")+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 2)) + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))

## Full model for CORNEAL AREA

m1 <- lm(log.cornea.area ~ log.tibia.length + sex * elevation, data = long)

plot(m1, which=1)
plot(m1, which=2)
hist(resid(m1))

drop1(m1, test="Chisq") #sex:elevation is ns
m2 <- update(m1,.~.-sex:elevation) 
drop1(m2, test="Chisq") #all significant

plot(allEffects(m2))
plot(m2, which = 2) 
hist(resid(m2))

summary(m2)
car::Anova(m2, test.statistic="Chisq")

## Sex effect

#pairwise comparisons
model_means <- emmeans(object = m2, specs = ~ sex, type = "response") 
model_means
#p-value comparison table
pwpm(model_means, adjust="bonferroni",  diffs = F)
#f       m
#f [9.483]  <.0001
#m         [9.553]

#effect size
eff_size(model_means, sigma = sigma(m2), edf = df.residual(m2))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

model_means_cld
#set marginal means at a dataframe
facet <- as.data.frame(model_means_cld)
facet

#plot elevation effect
facet %>% ggplot(aes(sex,facet$emmean))+
  geom_point(aes(color=sex), size=4) + geom_errorbar(aes(x=sex, ymin=facet$upper.CL, ymax=facet$lower.CL), width=0.2, color="grey50")+
  geom_text(aes(label = .group, y = facet$emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=5)  +
  labs(x="",y="e.m.m. log10 (corneal area)")+
  scale_x_discrete(labels=c("female", "male"))+
  scale_fill_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size=18))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## Elevation effect
coef(summary(m2))[ , "Estimate"]

el <- Effect("elevation", m2)
el <- data.frame(Effect("elevation", m2))

ggplot(data = el, aes(x=elevation, y=fit, ymin=lower, ymax=upper))+
  geom_line(col="#1F78B4", size=1) + geom_ribbon(alpha=0.4, fill="#1F78B4") + ylab("e.m.m. log10 (corneal area)") + 
  xlab("Elevation") + ggtitle("elevation effect")+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 2)) + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


## 6) Obtaining the residuals of eye size

#######################################
####Obtaining residuals of eye size####
#######################################


#Join data with transect position data
info <- read.csv("Meier2021_supplementary_sample_list.csv")

#remove unecessary columns
info <- info[c(1,7:14)] %>% rename(ID = ID.code)
head(info)

# Merge dataframes by different ID column names
long1 <- merge(long, info, by = "ID")
head(long1)

## 6.1 Residuals modeled by log10 abdomen length
FC_abd <- lm(log.facet.count ~ log.abdomen + sex, data = long1)
CA_abd <- lm(log.cornea.area ~ log.abdomen + sex, data = long1)
FA_abd <- lm(log.facet.area ~ log.abdomen + sex, data = long1)

par(mfrow=c(2,2))
plot(FC_abd)

plot(FC_abd, which = 1)
plot(FC_abd, which = 2)
plot(FC_abd, which = 3)
plot(FC_abd, which = 5)

summary(FC_abd)
summary(CA_abd)
summary(FA_abd)

drop1(FC_abd, test = "Chisq") #all significant
drop1(CA_abd, test = "Chisq") #all significant
drop1(FA_abd, test = "Chisq") #all significant

FC_residuals <- residuals(FC_abd)
length(FC_residuals) #445

CA_residuals <- residuals(CA_abd)
length(CA_residuals) #445

FA_residuals <- residuals(FA_abd)
length(FA_residuals) #445

# Create a new dataframe with residuals and non-missing observations
long2 <- long1[!is.na(long1$log.abdomen), ]
long2 <- long2[!is.na(long2$facet.count), ] # Select only non-missing observations
long2$FC_residuals <- FC_residuals
long2$CA_residuals <- CA_residuals
long2$FA_residuals <- FA_residuals # Assign residuals to non-missing observations

head(long2)
#write.csv(long2, "eye_data_residuals_abdomen_only.csv")

#are the residuals normally distributed?
shapiro.test(long2$FC_residuals) #yes
shapiro.test(long2$CA_residuals) #no

dev.off()

#explore data along transect
par(mfrow=c(1,2))
plot(long2$CA_residuals ~ long2$Transect.Position)
abline(h=0, col="blue")

plot(long2$FC_residuals ~ long2$Transect.Position)
abline(h=0, col="blue")

dev.off()

#Facet count residuals by subspp.
p1 <- ggplot(long2, aes(x=type, y=FC_residuals, fill=type)) +
  geom_violin(position=position_dodge(width=1)) + theme_classic()+
  geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A")) +
  scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("Residuals") + xlab("") + 
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size = 18))

p1

#Corneal area residuals by subspp.
p1 <- ggplot(long2, aes(x=type, y=CA_residuals, fill=type)) +
  geom_violin(position=position_dodge(width=1)) + theme_classic()+
  geom_jitter(position = position_dodge(width=1), alpha = 0.8, size=2, col="grey20") + 
  geom_boxplot(position = position_dodge(width=1), width=0.2, alpha = 0.5, col="white") +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A")) +
  scale_x_discrete(limits=c("lativitta", "hybrid", "notabilis")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  ylab("Residuals") + xlab("")  +
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(text = element_text(family = "Times New Roman", size = 18))

p1

# 6.2 Residuals modeled by log10 tibia length
FC_tib <- lm(log.facet.count ~ log.tibia.length + sex, data = long1)
CA_tib <- lm(log.cornea.area ~ log.tibia.length + sex, data = long1)

summary(FC_tib)
summary(CA_tib)

drop1(FC_tib, test = "Chisq") #all significant
drop1(CA_tib, test = "Chisq") #all significant

FC_residuals <- residuals(FC_tib)
length(FC_residuals) #387 - we loose almost 100 samples!

CA_residuals <- residuals(CA_tib)
length(CA_residuals) #387 - we loose almost 100 samples!

# Create a new dataframe with residuals and non-missing observations
long2 <- long1[!is.na(long1$log.tibia.length), ]
long2 <- long2[!is.na(long2$facet.count), ] # Select only non-missing observations
long2$FC_residuals <- FC_residuals
long2$CA_residuals <- CA_residuals  # Assign residuals to those observations

head(long2)
#write.csv(long2, "eye_data_with_residuals_tibia_only.csv")

## Do subspecies differ in their residuals?
# Perform ANOVA
anova_result <- aov(CA_residuals ~ type, data = long2)

# Post-hoc test using Tukey's HSD
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

library(multcompView)
# Extract the Tukey test results and generate significance letters
tukey_letters <- multcompLetters4(anova_result, tukey_result)
tukey_letters

#lativitta    hybrid notabilis 
#"a"       "a"       "b" 

#Now for CORNEAL AREA
anova_result <- aov(CA_residuals ~ type, data = long2)

# Post-hoc test using Tukey's HSD
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

library(multcompView)

# Extract the Tukey test results and generate significance letters
tukey_letters <- multcompLetters4(anova_result, tukey_result)
tukey_letters

#lativitta    hybrid notabilis 
#"a"       "a"       "b" 


# plot residuals along elevation gradient
ggplot(long2,aes(Elevation, FC_residuals))+geom_point(aes(color = type)) + 
  theme_classic() +geom_smooth(method='lm', color="grey")+ labs(title="log(facet count) ~ log(tibia) + sex", y="Residuals")+
  scale_color_manual(values=c("#56B4E9",  "grey70", "#E69F00")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank())

ggplot(long2,aes(Elevation, CA_residuals))+geom_point(aes(color = type)) + 
  theme_classic() +geom_smooth(method='lm', color="grey")+ labs(title="log(cornea area) ~ log(tibia) + sex", y="Residuals")+
  scale_color_manual(values=c("#56B4E9",  "grey70", "#E69F00")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        legend.title = element_blank())

## 6.3 Test that sex and body size do not have an effect on the residuals

#ca
m1 <- lm(CA_residuals ~ sex + log.abdomen, data = long2)
summary(m1)

#fc
m1 <- lm(FC_residuals ~ sex + log.abdomen, data = long2)
summary(m1)

#endEye

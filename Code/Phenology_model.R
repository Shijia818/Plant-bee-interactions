
library(dplyr)
library(lme4)
library(lmerTest)
library(ape)
library(performance)
library(sf)
library(stars)
library(spaMM)
library(biomod2)
library(gam)
library(stringr)
library(raster)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stars)
library(emmeans)
library(glmmTMB)
library(minpack.lm)
library(maxlike)

data.plant <- read.csv("./Analyses/Model/data.plant.final.csv") 
data.bee <- read.csv("./Analyses/Model/data.bee.final.csv")

############# correlation among predictors #################
cor.test(data.plant$altitude,data.plant$T_anm)
cor.test(data.plant$latitude,data.plant$P_anm)
cor.test(data.plant$T_anm,data.plant$P_anm)
cor.test(data.plant$latitude,data.plant$T_anm)
cor.test(data.plant$altitude,data.plant$P_anm)
cor.test(data.bee$T_anm,data.bee$P_anm)
cor.test(data.bee$altitude,data.bee$T_anm)
cor.test(data.bee$altitude,data.bee$P_anm)
cor.test(data.bee$latitude,data.bee$T_anm)
cor.test(data.bee$latitude,data.bee$P_anm)

data.plant <- read.csv("./Analyses/Model/data.plant.final.csv") 
d.plant <- data.plant %>% filter(flower > 0) %>% filter(flower > (bud + fruit))
quantiles <- d.plant %>% group_by(Species) %>% summarize(quantile_95 = quantile(flower, probs = 0.95),
                                                         quantile_5_date = quantile(julian, probs = 0.05),
                                                         quantile_95_date = quantile(julian, probs = 0.95))

d.plant <- left_join(d.plant, quantiles, by = "Species")
d.plant <- d.plant %>% mutate(threshold = quantile_95 * 0.05)
d.plant.filter <- d.plant %>% filter(flower >= threshold, julian >= quantile_5_date, julian <= quantile_95_date)

model.plant <- lmerTest::lmer(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) +
                                scale(T_anm):scale(P_anm) + 
                                (1 | Species) + 
                                (0 + scale(T_anm) | Species) +
                                (0 + scale(P_anm) | Species), data = d.plant.filter, REML = T, 
                              control = lmerControl(optimizer = "bobyqa",
                                                    optCtrl = list(maxfun = 2e5)))

summary(model.plant)
model_performance(model.plant)


data.bee.filter <- read.csv("./Analyses/Results_bio4/Bee_probability/data.bee.filter.csv")
model.bee <- lmerTest::lmer(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) +
                              scale(T_anm):scale(P_anm) +
                              (1 | Species) + 
                              (0 + scale(T_anm) | Species) +
                              (0 + scale(P_anm) | Species), data = data.bee.filter, REML = T, 
                            control = lmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 2e5)))

summary(model.bee)
model_performance(model.bee)


data.generalist <- data.bee.filter[which(data.bee.filter$Species != "Andrena_violae"), ]
data.specialist <- data.bee.filter[which(data.bee.filter$Species == "Andrena_violae"), ]

model.specialist <- lm(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) + 
                         scale(T_anm):scale(P_anm), data = data.specialist)

summary(model.specialist)
model_performance(model.specialist)


model.generalist <- lmerTest::lmer(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) +
                                     scale(T_anm):scale(P_anm) +
                                     (1 | Species) + 
                                     (0 + scale(T_anm) | Species) +
                                     (0 + scale(P_anm) | Species), data = data.generalist, REML = T, 
                                   control = lmerControl(optimizer = "bobyqa",
                                                         optCtrl = list(maxfun = 2e5)))

summary(model.generalist)
model_performance(model.generalist)
















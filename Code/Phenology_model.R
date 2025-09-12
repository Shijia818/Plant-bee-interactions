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
library(Metrics)
library(dplyr)

setwd("D:/Mutualistic_interaction")

################### Plant Phenology model ##########################

data.plant <- read.csv("./Analyses/Model/data.plant.final.csv")

d.plant <- data.plant %>% filter(flower > 0) %>% filter(flower > (bud + fruit))
quantiles <- d.plant %>% group_by(Species) %>% summarize(quantile_95 = quantile(flower, probs = 0.95),
                                                         quantile_5_date = quantile(julian, probs = 0.05),
                                                         quantile_95_date = quantile(julian, probs = 0.95))

d.plant <- left_join(d.plant, quantiles, by = "Species")
d.plant <- d.plant %>% mutate(threshold = quantile_95 * 0.05)
d.plant.filter <- d.plant %>% filter(flower >= threshold, julian >= quantile_5_date, julian <= quantile_95_date)


model.plant <- fitme(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) +
                             scale(T_anm):scale(P_anm) +
                             (1 | Species) + (0 + scale(T_anm) | Species) + (0 + scale(P_anm) | Species) + AR1(1 | year) + Matern(1 | longitude + latitude),
                           data = d.plant.filter,
                           family = gaussian(),
                           method = "REML")

VarF <- var(predict(model.plant, re.form = NA))
VarR <- sum(model.plant$lambda)
VarE <- sigma(model.plant)^2

R2_marginal <- VarF / (VarF + VarR + VarE)
R2_conditional <- (VarF + VarR) / (VarF + VarR + VarE)

df <- d.plant.filter
n <- nrow(df)
preds <- numeric(n)

failures <- logical(n)

for(i in 1:n){
  train <- df[-i, ]
  test <- df[i, , drop = FALSE]
  
  
  fit <- tryCatch({model.plant}, error = function(e) NULL)
  
  if (!is.null(fit)) {
    preds[i] <- predict(fit, newdata = test, allow.new.levels = TRUE)
  } else {
    failures[i] <- TRUE
    preds[i] <- NA
  }
}

actual <- df$julian[!failures]
preds <- preds[!failures]

rmse_val <- rmse(actual, preds)
mae_val <- mae(actual, preds)
r2_val <- cor(actual, preds)^2


##################### bee phenology model #################

data.bee.filter <- read.csv("./Analyses/Results_bio4/Bee_probability/data.bee.filter.csv")

data.generalist <- data.bee.filter[which(data.bee.filter$Species != "Andrena_violae"), ]
data.specialist <- data.bee.filter[which(data.bee.filter$Species == "Andrena_violae"), ]

model.generalist <- fitme(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) +
                                     scale(T_anm):scale(P_anm) +
                                     (1 | Species) + 
                                     (0 + scale(T_anm) | Species) +
                                     (0 + scale(P_anm) | Species) + AR1(1 | year) + Matern(1 | longitude + latitude), data = data.generalist, family = gaussian(),method = "REML")

model.specialist <- fitme(julian ~ scale(mean.T) + scale(mean.P) + scale(T_anm) + scale(P_anm) + 
                         scale(T_anm):scale(P_anm) + AR1(1 | year) + Matern(1 | longitude + latitude), family = gaussian(), method = "REML", data = data.specialist)


df.bee <- data.generalist
n <- nrow(df.bee)
preds <- numeric(n)

failures <- logical(n)

for(i in 1:n){
  train <- df.bee[-i, ]
  test <- df.bee[i, , drop = FALSE]
  
  fit <- tryCatch({model.generalist}, error = function(e) NULL)
  
  if (!is.null(fit)) {
    preds[i] <- predict(fit, newdata = test, allow.new.levels = TRUE)
  } else {
    failures[i] <- TRUE
    preds[i] <- NA
  }
  print(i)
}

actual <- df.bee$julian[!failures]
preds <- preds[!failures]

rmse_val <- rmse(actual, preds)
mae_val <- mae(actual, preds)
r2_val <- cor(actual, preds)^2

#########################    predict plant phenology at each grid cell ##################################

biomod.spdis.plant <- read.csv("./Analyses/Model/biomod.spdis.plant.csv",row.names = 1)
data.XY <- read.csv("./Analyses/SDM/Climate/grid40km.csv", row.names = 1)
data.grid.T <- read.csv("./Analyses/SDM/Climate/data.grid.T.csv")
data.grid.P <- read.csv("./Analyses/SDM/Climate/data.grid.P.csv")
res <- which(rownames(data.XY) %in% data.grid.P$Gridcode == F)
data.XY <- data.XY[-res,]


T_ano_current_plant <- data.grid.T$bio1 - data.grid.T$Ave
T_ano_GI_plant <- data.grid.T$GI - data.grid.T$Ave
T_ano_HD_plant <- data.grid.T$HD - data.grid.T$Ave
T_ano_IN_plant <- data.grid.T$IN - data.grid.T$Ave

P_ano_current_plant <- data.grid.P$bio12 / data.grid.P$Ave
P_ano_GI_plant <- data.grid.P$GI / data.grid.P$Ave
P_ano_HD_plant <- data.grid.P$HD / data.grid.P$Ave
P_ano_IN_plant <- data.grid.P$IN / data.grid.P$Ave

mean.T <- data.grid.T$Ave
mean.P <- data.grid.P$Ave

data.predict.plant <- cbind(data.XY, mean.T, mean.P, T_ano_current_plant,T_ano_GI_plant,T_ano_HD_plant,T_ano_IN_plant,P_ano_current_plant,P_ano_GI_plant,P_ano_HD_plant,P_ano_IN_plant)

Species <- colnames(biomod.spdis.plant)
test <- rbind()
for (i in 1:length(Species)){
  test <- rbind(test,data.predict.plant)}

Species.rev <- rep(Species,nrow(data.XY)) %>% sort() %>% cbind(test)
colnames(Species.rev) [1] <- "Species"

model.current <- Species.rev[,c("Species", "longitude","latitude","mean.T","mean.P","T_ano_current_plant","P_ano_current_plant")]
colnames(model.current)[6:7] <- c("T_anm","P_anm")
model.current$year <- 1970

predicted.current <- predict(model.plant, newdata = model.current, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.current, "predVar"))

lower.CI <- predicted.current - 1.96 * se
upper.CI <- predicted.current + 1.96 * se

plant.prediction.current <- cbind(predicted.current, lower.CI, upper.CI)
colnames(plant.prediction.current) <- c("prediction","lower.CI","upper.CI")
write.csv(plant.prediction.current, file = "./Revision/Results/plant.prediction.current.csv")


###########################################################

model.future.GI <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_GI_plant","P_ano_GI_plant")]
colnames(model.future.GI)[6:7] <- c("T_anm","P_anm")
model.future.GI$year <- 2070


predicted.future.GI <- predict(model.plant, newdata = model.future.GI, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.GI, "predVar"))

lower.CI <- predicted.future.GI - 1.96 * se
upper.CI <- predicted.future.GI + 1.96 * se

plant.prediction.GI <- cbind(predicted.future.GI, lower.CI, upper.CI)
colnames(plant.prediction.GI) <- c("prediction","lower.CI","upper.CI")

###########################################

model.future.HD <- Species.rev[,c("Species", "latitude","longitude", "mean.T","mean.P","T_ano_HD_plant","P_ano_HD_plant")]
colnames(model.future.HD)[6:7] <- c("T_anm","P_anm")
model.future.HD$year <- 2070
predicted.future.HD <- predict(model.plant, newdata = model.future.HD, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.HD, "predVar"))

lower.CI <- predicted.future.HD - 1.96 * se
upper.CI <- predicted.future.HD + 1.96 * se

plant.prediction.HD <- cbind(predicted.future.HD, lower.CI, upper.CI)
colnames(plant.prediction.HD) <- c("prediction","lower.CI","upper.CI")


################################################

model.future.IN <- Species.rev[,c("Species", "latitude","longitude", "mean.T","mean.P","T_ano_IN_plant","P_ano_IN_plant")]
colnames(model.future.IN)[6:7] <- c("T_anm","P_anm")
model.future.IN$year <- 2070
predicted.future.IN <- predict(model.plant, newdata = model.future.IN, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.IN, "predVar"))

lower.CI <- predicted.future.IN - 1.96 * se
upper.CI <- predicted.future.IN + 1.96 * se

plant.prediction.IN <- cbind(predicted.future.IN, lower.CI, upper.CI)
colnames(plant.prediction.IN) <- c("prediction","lower.CI","upper.CI")


################################ predict bee phenology at each grid cell #############

biomod.spdis.bee <- read.csv("./Analyses/Model/biomod.spdis.bee.csv", row.names =1)
data.predict.bee <- cbind(data.XY, mean.T, mean.P, T_ano_current_bee,T_ano_GI_bee,T_ano_HD_bee,T_ano_IN_bee,P_ano_current_bee,P_ano_GI_bee,P_ano_HD_bee,P_ano_IN_bee)
Species <- colnames(biomod.spdis.bee)[which(colnames(biomod.spdis.bee) != "Andrena_violae")]
test <- rbind()
for (i in 1:length(Species)){
  test <- rbind(test,data.predict.bee)}

Species.rev <- rep(Species,nrow(data.XY)) %>% sort() %>% cbind(test)
colnames(Species.rev) [1] <- "Species"

model.current <- Species.rev[,c("Species","latitude","longitude", "mean.T","mean.P","T_ano_current_bee","P_ano_current_bee")]
colnames(model.current)[6:7] <- c("T_anm","P_anm")
model.current$year <- 1970

predicted.current <- predict(model.generalist, newdata = model.current, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.current, "predVar"))

lower.CI <- predicted.current - 1.96 * se
upper.CI <- predicted.current + 1.96 * se

generalist.prediction.current <- cbind(predicted.current, lower.CI, upper.CI)
colnames(generalist.prediction.current) <- c("prediction","lower.CI","upper.CI")


##########################################

model.future.GI <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_GI_bee","P_ano_GI_bee")]
colnames(model.future.GI)[6:7] <- c("T_anm","P_anm")
model.future.GI$year <- 2070

predicted.future.GI <- predict(model.generalist, newdata = model.future.GI, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.GI, "predVar"))

lower.CI <- predicted.future.GI - 1.96 * se
upper.CI <- predicted.future.GI + 1.96 * se

generalist.prediction.GI <- cbind(predicted.future.GI, lower.CI, upper.CI)
colnames(generalist.prediction.GI) <- c("prediction","lower.CI","upper.CI")

#############################################


model.future.HD <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_HD_bee","P_ano_HD_bee")]
colnames(model.future.HD)[6:7] <- c("T_anm","P_anm")
model.future.HD$year <- 2070

predicted.future.HD <- predict(model.generalist, newdata = model.future.HD, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.HD, "predVar"))

lower.CI <- predicted.future.HD - 1.96 * se
upper.CI <- predicted.future.HD + 1.96 * se

generalist.prediction.HD <- cbind(predicted.future.HD, lower.CI, upper.CI)
colnames(generalist.prediction.HD) <- c("prediction","lower.CI","upper.CI")

######################################

model.future.IN <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_IN_bee","P_ano_IN_bee")]
colnames(model.future.IN)[6:7] <- c("T_anm","P_anm")
model.future.IN$year <- 2070

predicted.future.IN <- predict(model.generalist, newdata = model.future.IN, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.IN, "predVar"))

lower.CI <- predicted.future.IN - 1.96 * se
upper.CI <- predicted.future.IN + 1.96 * se

generalist.prediction.IN <- cbind(predicted.future.IN, lower.CI, upper.CI)
colnames(generalist.prediction.IN) <- c("prediction","lower.CI","upper.CI")



############ specialist bee ###########

data.predict.bee <- cbind(data.XY, mean.T, mean.P, T_ano_current_bee,T_ano_GI_bee,T_ano_HD_bee,T_ano_IN_bee,P_ano_current_bee,P_ano_GI_bee,P_ano_HD_bee,P_ano_IN_bee)
Species <- colnames(biomod.spdis.bee)[which(colnames(biomod.spdis.bee) == "Andrena_violae")]
Species.rev <- cbind(Species, data.predict.bee)

model.current <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_current_bee","P_ano_current_bee")]
colnames(model.current)[6:7] <- c("T_anm","P_anm")
model.current$year <- 1970

predicted.current <- predict(model.specialist, newdata = model.current, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.current, "predVar"))

lower.CI <- predicted.current - 1.96 * se
upper.CI <- predicted.current + 1.96 * se

specialist.prediction.current <- cbind(predicted.current, lower.CI, upper.CI)
colnames(specialist.prediction.current) <- c("prediction","lower.CI","upper.CI")



##########################################

model.future.GI <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_GI_bee","P_ano_GI_bee")]
colnames(model.future.GI)[6:7] <- c("T_anm","P_anm")
model.future.GI$year <- 2070

predicted.future.GI <- predict(model.specialist, newdata = model.future.GI, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.GI, "predVar"))

lower.CI <- predicted.future.GI - 1.96 * se
upper.CI <- predicted.future.GI + 1.96 * se

specialist.prediction.GI <- cbind(predicted.future.GI, lower.CI, upper.CI)
colnames(specialist.prediction.GI) <- c("prediction","lower.CI","upper.CI")


############################################

model.future.HD <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_HD_bee","P_ano_HD_bee")]
colnames(model.future.HD)[6:7] <- c("T_anm","P_anm")
model.future.HD$year <- 2070


predicted.future.HD <- predict(model.specialist, newdata = model.future.HD, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.HD, "predVar"))

lower.CI <- predicted.future.HD - 1.96 * se
upper.CI <- predicted.future.HD + 1.96 * se

specialist.prediction.HD <- cbind(predicted.future.HD, lower.CI, upper.CI)
colnames(specialist.prediction.HD) <- c("prediction","lower.CI","upper.CI")

############################################

model.future.IN <- Species.rev[,c("Species","latitude","longitude","mean.T","mean.P","T_ano_IN_bee","P_ano_IN_bee")]
colnames(model.future.IN)[6:7] <- c("T_anm","P_anm")
model.future.IN$year <- 2070


predicted.future.IN <- predict(model.specialist, newdata = model.future.IN, variances = list(predVar = TRUE))
se <- sqrt(attr(predicted.future.IN, "predVar"))

lower.CI <- predicted.future.IN - 1.96 * se
upper.CI <- predicted.future.IN + 1.96 * se

specialist.prediction.IN <- cbind(predicted.future.IN, lower.CI, upper.CI)
colnames(specialist.prediction.IN) <- c("prediction","lower.CI","upper.CI")













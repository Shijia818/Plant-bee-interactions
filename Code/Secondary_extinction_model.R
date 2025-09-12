
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

data.Rsquare <- read.csv("./Revision/data.Rsquare.power.csv",row.names = 1)
result <- as.data.frame(which(data.Rsquare < 0.5 | is.na(data.Rsquare), arr.ind = TRUE))
result$Plant <- rownames(data.Rsquare)[result$row]
result$name <- colnames(data.Rsquare)[result$col]
result <- result[, c("Plant", "name")]

file <- list.files("./Revision/Results/Probability/Median/")
data.res <- c()
for (i in 1:length(file)){
  data.file <- read.csv(paste("./Revision/Results/Probability/Median/",file[i],sep = ""))
  colnames(data.file)[1] <- "Gridcode"
  sp <- str_remove(file[i], ".csv$")
  data.file$Plant <- sp
  data.file <- data.file[, c("Plant", setdiff(names(data.file), "Plant"))]
  data.res <- rbind(data.res, data.file)
}


data.res_t <- pivot_longer(data.res,col = colnames(data.res)[3:9])

spnames <- colnames(data.res)[3:9]
type <- c(rep("generalist",6),"specialist")
data.type <- as.data.frame(cbind(spnames,type))
data.type.1 <- data.type[match(data.res_t$name,data.type$spnames),]
data.res_t$type <- data.type.1$type

data.res_t_filter <- anti_join(data.res_t,result,by = c("Plant", "name"))
data.res_t_filter <- data.res_t_filter[!is.na(data.res_t_filter$value), ]


file.mismatch <- list.files("./Revision/Results/Phenology_mismatch/Median/")
data.mismatch <- c()

for (i in 1:length(file.mismatch)){
  data.file <- read.csv(paste("./Revision/Results/Phenology_mismatch/Medain/",file.mismatch[i],sep = ""))
  colnames(data.file)[1] <- "Gridcode"
  sp <- str_remove(file.mismatch[i], ".csv$")
  data.file$Plant <- sp
  data.file <- data.file[, c("Plant", setdiff(names(data.file), "Plant"))]
  data.mismatch <- rbind(data.mismatch, data.file)
}


data.mismatch_t <- pivot_longer(data.mismatch,col = colnames(data.mismatch)[3:9])
spnames <- colnames(data.mismatch)[3:9]
type <- c(rep("generalist",6),"specialist")
data.type <- as.data.frame(cbind(spnames,type))
data.type.1 <- data.type[match(data.mismatch_t$name,data.type$spnames),]
data.mismatch_t$type <- data.type.1$type

data.mismatch_t_filter <- anti_join(data.mismatch_t,result,by = c("Plant", "name"))

data.res_t_filter$phenological_mismatch <- data.mismatch_t_filter$value
res <- which(data.res_t_filter$value >= 1)
data.res_t_filter_final <- data.res_t_filter[-res,]


fit.res <- glmmTMB(value ~ scale(latitude) * type + (1 | name) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = data.res_t_filter_final)
summary(fit.res)


data.res_t_filter_ge <- data.res_t_filter_final[which(data.res_t_filter_final$type == "generalist"),]
fit.res.mismatch_ge <- glmmTMB(value ~  scale(phenological_mismatch) * scale(latitude) + (1 | name) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = data.res_t_filter_ge)

data.res_t_filter_spe <- data.res_t_filter_final[which(data.res_t_filter_final$type == "specialist"),]
fit.res.mismatch_spe <- glmmTMB(value ~ scale(phenological_mismatch) * scale(latitude) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = data.res_t_filter_spe)

############################

generalists <- data.res_t_filter_final[data.res_t_filter_final$type == "generalist", ]
min_generalist_rows <- generalists[ave(generalists$phenological_mismatch, generalists$Plant, generalists$Gridcode, FUN = min) == generalists$phenological_mismatch, ]
specialists <- data.res_t_filter_final[data.res_t_filter_final$type == "specialist", ]
new_df <- rbind(min_generalist_rows, specialists)
fit.res <- glmmTMB(value ~ latitude * type + (1 | name) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = new_df)





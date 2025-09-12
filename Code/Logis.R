
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

convert_doy_to_angle <- function(doy){
  angle <- (doy / 365) * 360
  return(angle)
}

circular_overlap_ratio <- function(start1, end1, start2, end2) {
  
  get_segments <- function(start, end) {
    if (start <= end) {
      return(list(c(start, end)))
    } else {
      return(list(c(start, 360), c(0, end)))
    }
  }
  
  segs1 <- get_segments(start1, end1)
  segs2 <- get_segments(start2, end2)
  
  overlap_length <- 0
  for (s1 in segs1) {
    for (s2 in segs2) {
      overlap_start <- max(s1[1], s2[1])
      overlap_end <- min(s1[2], s2[2])
      if (overlap_start < overlap_end) {
        overlap_length <- overlap_length + (overlap_end - overlap_start)
      }
    }
  }
  
  flower_length <- if (start1 <= end1) {
    end1 - start1
  } else {
    (360 - start1) + end1
  }
  
  overlap_ratio <- overlap_length / flower_length
  
  return(overlap_ratio)
} 

plant_prediction <- read.csv("./Revision/Data/plant_prediction.csv")
generalist_prediction <- read.csv("./Revision/Data/generalist_prediction.csv")
specialist_prediction <- read.csv("./Revision/Data/specialist_prediction.csv")

file1 <- list.files("./Revision/Results/Probability/Current")
file2 <- list.files("./Revision/Results/phenology_mismatch/Current")


for(j in 1:length(file1)){
  data.XY <- read.csv("./Analyses/SDM/Climate/grid40km.csv", row.names = 1)
  spnames <- str_remove(file1[j],".csv$")
  data.prob.summary <- read.csv(paste("./Revision/Results/Probability/Current/",file1[j],sep=""),row.names = 1)
  data.mismatch <- read.csv(paste("./Revision/Results/Phenology_mismatch/Current/",file2[j],sep = ""),row.names =1)
  
  
  data.sp.plant <- plant_prediction[which(plant_prediction$Species %in% spnames),]
  
  data.sp.plant.rev <- rbind()
  for (k in 1:6){
    data.sp.plant.rev <- rbind(data.sp.plant.rev,data.sp.plant)}
  
  
  ###########################################################
  mismatch_generalist_GI <- mapply(circular_overlap_ratio, data.sp.plant.rev$lower_CI_GI, data.sp.plant.rev$upper_CI_GI, generalist_prediction$lower_CI_GI, generalist_prediction$upper_CI_GI)
  mismatch_specialist_GI <- mapply(circular_overlap_ratio, data.sp.plant$lower_CI_GI, data.sp.plant$upper_CI_GI, specialist_prediction$lower_CI_GI, specialist_prediction$upper_CI_GI)
  
  mismatch_generalist_HD <- mapply(circular_overlap_ratio, data.sp.plant.rev$lower_CI_HD, data.sp.plant.rev$upper_CI_HD, generalist_prediction$lower_CI_HD, generalist_prediction$upper_CI_HD)
  mismatch_specialist_HD <- mapply(circular_overlap_ratio, data.sp.plant$lower_CI_HD, data.sp.plant$upper_CI_HD, specialist_prediction$lower_CI_HD, specialist_prediction$upper_CI_HD)
  
  mismatch_generalist_IN <- mapply(circular_overlap_ratio, data.sp.plant.rev$lower_CI_IN, data.sp.plant.rev$upper_CI_IN, generalist_prediction$lower_CI_IN, generalist_prediction$upper_CI_IN)
  mismatch_specialist_IN <- mapply(circular_overlap_ratio, data.sp.plant$lower_CI_IN, data.sp.plant$upper_CI_IN, specialist_prediction$lower_CI_IN, specialist_prediction$upper_CI_IN)
  
  data.sp.plant.rev$mismatch_generalist_GI <- mismatch_generalist_GI
  data.sp.plant.rev$mismatch_generalist_HD <- mismatch_generalist_HD
  data.sp.plant.rev$mismatch_generalist_IN <- mismatch_generalist_IN
  data.sp.plant.rev$mismatch_generalist_median <- apply(data.sp.plant.rev[,c("mismatch_generalist_GI","mismatch_generalist_HD","mismatch_generalist_IN")],1,median)
  
  data.sp.plant$mismatch_specialist_GI <- mismatch_specialist_GI
  data.sp.plant$mismatch_specialist_HD <- mismatch_specialist_HD
  data.sp.plant$mismatch_specialist_IN <- mismatch_specialist_IN
  data.sp.plant$mismatch_specialist_median <- apply(data.sp.plant[,c("mismatch_specialist_GI","mismatch_specialist_HD","mismatch_specialist_IN")],1,median)
  
  
  data.pro.future <- as.data.frame(matrix(NA,1157,6))
  rownames(data.pro.future) <- data.sp.plant$Gridcode
  colnames(data.pro.future) <- colnames(data.prob.summary)[1:6]
  
  for(i in 1:ncol(data.pro.future)){
    y <- data.prob.summary[,i]
    x <- data.mismatch[,i]
    
    fit.glm <- try(glm(y ~ x, family = binomial(link = "logit")))
    
    
    if (inherits(fit.glm, "try-error")) {
      warning(paste("nlsLM failed for column", i, "- skipping."))
      next
    }
    
    mismatch_gen <- data.sp.plant.rev$mismatch_generalist_median[((i-1)*nrow(data.pro.future)+1):(i*nrow(data.pro.future))]
    mismatch_gen <- 1- mismatch_gen
    predicted_value <- predict(fit.glm, newdata = list(x = mismatch_gen), type = "response")
    data.pro.future[,i] <- predicted_value
  }
  
  mismatch_spe <- data.sp.plant$mismatch_specialist_median
  mismatch_spe <- 1- mismatch_spe
  y <- data.prob.summary[,7]
  x <- data.mismatch[,7]
  

  fit.glm <- try(glm(y ~ x, family = binomial(link = "logit")))
  
  if (inherits(fit.glm, "try-error")) {
    warning(paste("nlsLM failed for column", i, "- skipping."))
    next
  }
  
  predicted_spe <- predict(fit.glm, newdata = list(x = mismatch_spe), type = "response")
  
  data.pro.future$Andrena_violae <- predicted_spe
  
  res <- which(rownames(data.XY) %in% data.sp.plant$Gridcode == F)
  data.XY <- data.XY[-res,]
  data.pro.future$latitude <- data.XY$latitude
  
  write.csv(data.pro.future, file = paste("./Revision/Results/Probability/glm/",spnames,".csv",sep = ""))
  
  
  data.mismatch.future <- as.data.frame(matrix(NA,1157,6))
  rownames(data.mismatch.future) <- data.sp.plant$Gridcode
  colnames(data.mismatch.future) <- c("Andrena_carlini", "Andrena_miserabilis","Osmia_atriventris","Osmia_bucephala","Osmia_lignaria","Osmia_pumila")
  
  for(n in 1:ncol(data.mismatch.future)){
    step1 <- data.sp.plant.rev$mismatch_generalist_median[((n-1)*nrow(data.mismatch.future)+1):(n*nrow(data.mismatch.future))]
    step2 <- 1-step1
    data.mismatch.future[,n] <- step2
  }
  
  step1.1 <- data.sp.plant$mismatch_specialist_median
  step2.1 <- 1-step1.1
  data.mismatch.future$Andrena_violae <- step2.1
  
  write.csv(data.mismatch.future, file = paste("./Revision/Results/Phenology_mismatch/glm/",spnames,".csv",sep = ""))
}


##########################################


data.Rsquare <- read.csv("./Revision/data.Rsquare.glm.csv",row.names = 1)
result <- as.data.frame(which(data.Rsquare < 0.2 | is.na(data.Rsquare), arr.ind = TRUE))
result$Plant <- rownames(data.Rsquare)[result$row]
result$name <- colnames(data.Rsquare)[result$col]
result <- result[, c("Plant", "name")]

file <- list.files("./Revision/Results/Probability/glm/")
data.res <- c()
for (i in 1:length(file)){
  data.file <- read.csv(paste("./Revision/Results/Probability/glm/",file[i],sep = ""))
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


file.mismatch <- list.files("./Revision/Results/Phenology_mismatch/glm/")
data.mismatch <- c()

for (i in 1:length(file.mismatch)){
  data.file <- read.csv(paste("./Revision/Results/Phenology_mismatch/glm/",file.mismatch[i],sep = ""))
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











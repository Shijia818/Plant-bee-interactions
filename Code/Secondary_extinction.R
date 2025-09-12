
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

################ calculate phenological overlap between plants and bees and secondary extinction risk ######

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

biomod.prob.plant <- read.csv("./Analyses/SDM/Results/Plant/MAXNET/pseudo_absence/probability/proj_current_MAXNET.csv", row.names =1)
plant.prediction.current <- read.csv("./Revision/Results/Plant/plant.prediction.current.csv")
generalist.prediction.current <- read.csv("./Revision/Results/Generalist/generalist.prediction.current.csv")
specialist.prediction.current <- read.csv("./Revision/Results/Specialist/specialist.prediction.current.csv")

plant.prediction.current$plant_start_angle <- convert_doy_to_angle(plant.prediction.current$lower_CI_plant)
plant.prediction.current$plant_end_angle <- convert_doy_to_angle(plant.prediction.current$upper_CI_plant)

generalist.prediction.current$generalist_start_angle <- convert_doy_to_angle(generalist.prediction.current$lower_CI_bee)
generalist.prediction.current$generalist_end_angle <- convert_doy_to_angle(generalist.prediction.current$upper_CI_bee)

specialist.prediction.current$specialist_start_angle <- convert_doy_to_angle(specialist.prediction.current$lower_CI_bee)
specialist.prediction.current$specialist_end_angle <- convert_doy_to_angle(specialist.prediction.current$upper_CI_bee)

sp.plant <- unique(plant.prediction.current$Species)

for(i in 1:length(sp.plant)){
  
  spname <- sp.plant[i]
  data.sp.plant <- plant.prediction.current[which(plant.prediction.current$Species %in% spname),]
  
  data.sp.plant.rev <- rbind()
  for (k in 1:6){
    data.sp.plant.rev <- rbind(data.sp.plant.rev,data.sp.plant)}
  
  mismatch_generalist <- mapply(circular_overlap_ratio, data.sp.plant.rev$plant_start_angle, data.sp.plant.rev$plant_end_angle, generalist.prediction.current$generalist_start_angle, generalist.prediction.current$generalist_end_angle)
  data.sp.plant.rev$mismatch_generalist <- mismatch_generalist
  
  mismatch_specialist <- mapply(circular_overlap_ratio, data.sp.plant$plant_start_angle, data.sp.plant$plant_end_angle, specialist.prediction.current$specialist_start_angle, specialist.prediction.current$specialist_end_angle)
  data.sp.plant$mismatch_specialist <- mismatch_specialist
  
  biomod.prob.plant.e <- biomod.prob.plant[,which(colnames(biomod.prob.plant) %in% spname),drop = FALSE]
  biomod.prob.plant.e <- biomod.prob.plant.e[which(rownames(biomod.prob.plant.e) %in% data.sp.plant$Gridcode),]
  
  data.prob.summary <- as.data.frame(matrix(NA,1157,6))
  rownames(data.prob.summary) <- data.sp.plant$Gridcode
  colnames(data.prob.summary) <- unique(generalist.prediction.current$Species)
  
  for(j in 1:ncol(data.prob.summary)){
    step1 <- data.sp.plant.rev$mismatch_generalist[((j-1)*nrow(data.prob.summary)+1):(j*nrow(data.prob.summary))]
    step2 <- step1
    step3 <- biomod.prob.plant.e * step2
    step4 <- biomod.prob.plant.e - step3   ### secondary extinction risk
    data.prob.summary[,j] <- step4
  }
  
  step1.1 <- data.sp.plant$mismatch_specialist
  step2.1 <- step1.1
  step3.1 <- biomod.prob.plant.e * step2.1
  step4.1 <- biomod.prob.plant.e - step3.1
  data.prob.summary$Andrena_violae <- step4.1
  
  write.csv(data.prob.summary, file = paste("./Revision/Results/Probability/Current/",spname,".csv",sep = ""))
  
  data.mismatch <- as.data.frame(matrix(NA,1157,6))
  rownames(data.mismatch) <- data.sp.plant$Gridcode
  colnames(data.mismatch) <- unique(generalist.prediction.current$Species)
  
  for(n in 1:ncol(data.mismatch)){
    step1 <- data.sp.plant.rev$mismatch_generalist[((n-1)*nrow(data.mismatch)+1):(n*nrow(data.mismatch))]
    step2 <- 1-step1
    data.mismatch[,n] <- step2
  }
  
  step1.1 <- data.sp.plant$mismatch_specialist
  step2.1 <- 1-step1.1
  data.mismatch$Andrena_violae <- step2.1
  
  
  write.csv(data.mismatch, file = paste("./Revision/Results/Phenology_mismatch/Current/",spname,".csv",sep = ""))
}

  
#################### future projection of  ############################

plant_prediction <- read.csv("./Revision/Results/Plant/plant.prediction.future.csv")
generalist_prediction <- read.csv("./Revision/Results/Generalist/generalist.prediction.future.csv")
specialist_prediction <- read.csv("./Revision/Results/Specialist/specialist.prediction.future.csv")

plant_prediction$plant_start_angle_GI <- convert_doy_to_angle(plant_prediction$lower_CI_GI)
plant_prediction$plant_end_angle_GI <- convert_doy_to_angle(plant_prediction$upper_CI_GI)

plant_prediction$plant_start_angle_HD <- convert_doy_to_angle(plant_prediction$lower_CI_HD)
plant_prediction$plant_end_angle_HD <- convert_doy_to_angle(plant_prediction$upper_CI_HD)

plant_prediction$plant_start_angle_IN <- convert_doy_to_angle(plant_prediction$lower_CI_IN)
plant_prediction$plant_end_angle_IN <- convert_doy_to_angle(plant_prediction$upper_CI_IN)



generalist_prediction$generalist_start_angle_GI <- convert_doy_to_angle(generalist_prediction$lower_CI_GI)
generalist_prediction$generalist_end_angle_GI <- convert_doy_to_angle(generalist_prediction$upper_CI_GI)

generalist_prediction$generalist_start_angle_HD <- convert_doy_to_angle(generalist_prediction$lower_CI_HD)
generalist_prediction$generalist_end_angle_HD <- convert_doy_to_angle(generalist_prediction$upper_CI_HD)

generalist_prediction$generalist_start_angle_IN <- convert_doy_to_angle(generalist_prediction$lower_CI_IN)
generalist_prediction$generalist_end_angle_IN <- convert_doy_to_angle(generalist_prediction$upper_CI_IN)


specialist_prediction$specialist_start_angle_GI <- convert_doy_to_angle(specialist_prediction$lower_CI_GI)
specialist_prediction$specialist_end_angle_GI <- convert_doy_to_angle(specialist_prediction$upper_CI_GI)

specialist_prediction$specialist_start_angle_HD <- convert_doy_to_angle(specialist_prediction$lower_CI_HD)
specialist_prediction$specialist_end_angle_HD <- convert_doy_to_angle(specialist_prediction$upper_CI_HD)

specialist_prediction$specialist_start_angle_IN <- convert_doy_to_angle(specialist_prediction$lower_CI_IN)
specialist_prediction$specialist_end_angle_IN <- convert_doy_to_angle(specialist_prediction$upper_CI_IN)


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
  
  data.sp.plant.rev$mismatch_generalist_GI <- abs(mismatch_generalist_GI)
  data.sp.plant.rev$mismatch_generalist_HD <- abs(mismatch_generalist_HD)
  data.sp.plant.rev$mismatch_generalist_IN <- abs(mismatch_generalist_IN)
  data.sp.plant.rev$mismatch_generalist_median <- apply(data.sp.plant.rev[,c("mismatch_generalist_GI","mismatch_generalist_HD","mismatch_generalist_IN")],1,median)
  
  data.sp.plant$mismatch_specialist_GI <- abs(mismatch_specialist_GI)
  data.sp.plant$mismatch_specialist_HD <- abs(mismatch_specialist_HD)
  data.sp.plant$mismatch_specialist_IN <- abs(mismatch_specialist_IN)
  data.sp.plant$mismatch_specialist_median <- apply(data.sp.plant[,c("mismatch_specialist_GI","mismatch_specialist_HD","mismatch_specialist_IN")],1,median)
  
  
  data.pro.future <- as.data.frame(matrix(NA,1157,6))
  rownames(data.pro.future) <- data.sp.plant$Gridcode
  colnames(data.pro.future) <- colnames(data.prob.summary)[1:6]
  
  for(i in 1:ncol(data.pro.future)){
    y <- data.prob.summary[,i]
    x <- data.mismatch[,i]
   #fit.logis <- try(nlsLM(y ~ 1 / (1 + exp(-(a + b * x))), start = list(a = 0, b = 1)))
    fit.power <- try(nlsLM(y ~ a * x^b, start = list(a = 2, b = 0.5)))
    
    
    if (inherits(fit.power, "try-error")) {
      warning(paste("nlsLM failed for column", i, "- skipping."))
      next
    }
    
    mismatch_gen <- data.sp.plant.rev$mismatch_generalist_median[((i-1)*nrow(data.pro.future)+1):(i*nrow(data.pro.future))]
    mismatch_gen <- 1- mismatch_gen
    predicted_value <- predict(fit.power, newdata = list(x = mismatch_gen))
    data.pro.future[,i] <- predicted_value
  }
  
  mismatch_spe <- data.sp.plant$mismatch_specialist_median
  mismatch_spe <- 1- mismatch_spe
  y <- data.prob.summary[,7]
  x <- data.mismatch[,7]
  
  fit.power <- try(nlsLM(y ~ a * x^b, start = list(a = 2, b = 0.5)))
  
  
  if (inherits(fit.power, "try-error")) {
    warning(paste("nlsLM failed for column", i, "- skipping."))
    next
  }
  
  predicted_spe <- predict(fit.power, newdata = list(x = mismatch_spe))
  
  data.pro.future$Andrena_violae <- predicted_spe
  
  res <- which(rownames(data.XY) %in% data.sp.plant$Gridcode == F)
  data.XY <- data.XY[-res,]
  data.pro.future$latitude <- data.XY$latitude
  
  write.csv(data.pro.future, file = paste("./Revision/Results/Probability/power/",spnames,".csv",sep = ""))
  
  
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
  
  write.csv(data.mismatch.future, file = paste("./Revision/Results/Phenology_mismatch/power/",spnames,".csv",sep = ""))
}


 ###################################################

spnames_plants <- str_remove(file1,".csv$")

data.Rsquare <- as.data.frame(matrix(NA,length(spnames_plants),7))
rownames(data.Rsquare) <- spnames_plants
colnames(data.Rsquare) <- c("Andrena_carlini", "Andrena_miserabilis","Osmia_atriventris","Osmia_bucephala","Osmia_lignaria","Osmia_pumila","Andrena_violae")

for (j in 1:nrow(data.Rsquare)) {
  data.prob.summary <- read.csv(paste0("./Revision/Results/Probability/Current/", file1[j]), row.names = 1)
  data.mismatch <- read.csv(paste0("./Revision/Results/Phenology_mismatch/Current/", file2[j]), row.names = 1)
  
  for (i in 1:ncol(data.Rsquare)) {
    y <- data.prob.summary[, i]
    x <- data.mismatch[, i]
    
    valid <- which(!is.na(x) & !is.na(y) & x > 0 & y >= 0)
    x_valid <- x[valid]
    y_valid <- y[valid]
    
    if (length(x_valid) >= 5 && sd(x_valid) > 0 && sd(y_valid) > 0) {
      
      fit.power <- tryCatch({
        nlsLM(y_valid ~ a * x_valid^b, start = list(a = 2, b = 0.5))
      }, error = function(e) NULL)
      
      if (!is.null(fit.power)) {
        y_pred <- predict(fit.power)
        RSS <- sum((y_valid - y_pred)^2)
        TSS <- sum((y_valid - mean(y_valid))^2)
        R2 <- 1 - (RSS / TSS)
        data.Rsquare[j, i] <- R2
      }
    }
  }
}


write.csv(data.Rsquare, file = "./Revision/data.Rsquare.power.csv")








biomod.prob.plant <- read.csv("./Analyses/SDM/Results/Plant/MAXNET/pseudo_absence/probability/proj_current_MAXNET.csv", row.names =1) ## primary extinction risk
plant_prediction <- read.csv("./Analyses/Results/prediction/plant.prediction.final.csv")
generalist_prediction <- read.csv("./Analyses/Results/prediction/generalist.prediction.final.csv")
specialist_prediction <- read.csv("./Analyses/Results/prediction/specialist.prediction.final.csv")

################################ transform day of year to angle######################

convert_doy_to_angle <- function(doy){
  angle <- (doy / 365) * 360
  return(angle)
}

plant_prediction$Current_angle <- sapply(plant_prediction$Current, convert_doy_to_angle)
plant_prediction$GI_angle <- sapply(plant_prediction$GI, convert_doy_to_angle)
plant_prediction$HD_angle <- sapply(plant_prediction$HD, convert_doy_to_angle)
plant_prediction$IN_angle <- sapply(plant_prediction$IN, convert_doy_to_angle)

generalist_prediction$Current_angle <- sapply(generalist_prediction$Current, convert_doy_to_angle)
generalist_prediction$GI_angle <- sapply(generalist_prediction$GI, convert_doy_to_angle)
generalist_prediction$HD_angle <- sapply(generalist_prediction$HD, convert_doy_to_angle)
generalist_prediction$IN_angle <- sapply(generalist_prediction$IN, convert_doy_to_angle)

specialist_prediction$Current_angle <- sapply(specialist_prediction$Current, convert_doy_to_angle)
specialist_prediction$GI_angle <- sapply(specialist_prediction$GI, convert_doy_to_angle)
specialist_prediction$HD_angle <- sapply(specialist_prediction$HD, convert_doy_to_angle)
specialist_prediction$IN_angle <- sapply(specialist_prediction$IN, convert_doy_to_angle)


########################
# angle_trans <- function(angle){
#if(angle > 360) return(angle - 360) else
#return(angle)}

###############################################################

min_angle_diff <- function(angle1, angle2){
  diff <- abs(angle1 - angle2)
  min_diff <- min(diff, 360 - diff)
  if(angle1 - angle2 < 0) return(-min_diff) else
    return(min_diff)
}

sp.plant <- unique(plant_prediction$Species)

for(i in 1:length(sp.plant)){
  
  spname <- sp.plant[i]
  data.sp.plant <- plant_prediction[which(plant_prediction$Species %in% spname),]
  
  data.sp.plant.rev <- rbind()
  for (k in 1:6){
    data.sp.plant.rev <- rbind(data.sp.plant.rev,data.sp.plant)}
  
  
  mismatch_generalist <- mapply(min_angle_diff, data.sp.plant.rev$Current_angle, generalist_prediction$Current_angle)
  data.sp.plant.rev$mismatch_generalist <- mismatch_generalist
  
  mismatch_specialist <- mapply(min_angle_diff, data.sp.plant$Current_angle, specialist_prediction$Current_angle)
  data.sp.plant$mismatch_specialist <- mismatch_specialist
  
  ###########################
  
  biomod.prob.plant.e <- biomod.prob.plant[,which(colnames(biomod.prob.plant) %in% spname),drop = FALSE]
  biomod.prob.plant.e <- biomod.prob.plant.e[which(rownames(biomod.prob.plant.e) %in% data.sp.plant$Gridcode),]
  
  
  
  data.prob.summary <- as.data.frame(matrix(NA,1157,6))
  rownames(data.prob.summary) <- data.sp.plant$Gridcode
  colnames(data.prob.summary) <- unique(generalist_prediction$Species)
  
  for(j in 1:ncol(data.prob.summary)){
    step1 <- data.sp.plant.rev$mismatch_generalist[((j-1)*nrow(data.prob.summary)+1):(j*nrow(data.prob.summary))]
    step2 <- 1- abs(step1 / 180)
    step3 <- biomod.prob.plant.e * step2
    step4 <- biomod.prob.plant.e - step3   ### secondary extinction risk
    data.prob.summary[,j] <- step4
  }
  
  ################### specialist #####################
  step1.1 <- data.sp.plant$mismatch_specialist
  step2.1 <- 1 - abs(step1.1 / 180)
  step3.1 <- biomod.prob.plant.e * step2.1
  step4.1 <- biomod.prob.plant.e - step3.1
  data.prob.summary$Andrena_violae <- step4.1
  
  write.csv(data.prob.summary, file = paste("./Analyses/Results/probability/MAXNET/Current/",spname,".csv",sep = ""))
  
  
  ##########################################
  data.mismatch <- as.data.frame(matrix(NA,1157,6))
  rownames(data.mismatch) <- data.sp.plant$Gridcode
  colnames(data.mismatch) <- unique(generalist_prediction$Species)
  
  for(n in 1:ncol(data.mismatch)){
    step1 <- data.sp.plant.rev$mismatch_generalist[((n-1)*nrow(data.mismatch)+1):(n*nrow(data.mismatch))]
    step2 <- abs(step1)
    data.mismatch[,n] <- step2
  }
  
  step1.1 <- data.sp.plant$mismatch_specialist
  step2.1 <- abs(step1.1)
  data.mismatch$Andrena_violae <- step2.1
  
  write.csv(data.mismatch, file = paste("./Analyses/Results/phenology_mismatch/MAXNET/Current/",spname,".csv",sep = ""))
}

###################################################

file1 <- list.files("./Analyses/Results/probability/MAXNET/Current")
file2 <- list.files("./Analyses/Results/phenology_mismatch/MAXNET/Current")


for(j in 1:length(file1)){
  data.XY <- read.csv("./Analyses/SDM/Climate/grid40km.csv", row.names = 1)
  spnames <- str_remove(file1[j],".csv$")
  data.prob.summary <- read.csv(paste("./Analyses/Results/probability/MAXNET/Current/",file1[j],sep=""),row.names = 1)
  data.mismatch <- read.csv(paste("./Analyses/Results/phenology_mismatch/MAXNET/Current/",file2[j],sep = ""),row.names =1)
  
  #fit.lm <- lm(y ~ x)  
  #fit.power <- nls(y ~ a*x^b, start = list(a=2, b= 1.5))
  #fit.loess <- loess(y ~ x)
  
  data.sp.plant <- plant_prediction[which(plant_prediction$Species %in% spnames),]
  
  data.sp.plant.rev <- rbind()
  for (k in 1:6){
    data.sp.plant.rev <- rbind(data.sp.plant.rev,data.sp.plant)}
  
  
  ###########################################################
  mismatch_generalist_GI <- mapply(min_angle_diff, data.sp.plant.rev$GI_angle, generalist_prediction$GI_angle)
  mismatch_specialist_GI <- mapply(min_angle_diff, data.sp.plant$GI_angle, specialist_prediction$GI_angle)
  
  mismatch_generalist_HD <- mapply(min_angle_diff, data.sp.plant.rev$HD_angle, generalist_prediction$HD_angle)
  mismatch_specialist_HD <- mapply(min_angle_diff, data.sp.plant$HD_angle, specialist_prediction$HD_angle)
  
  mismatch_generalist_IN <- mapply(min_angle_diff, data.sp.plant.rev$IN_angle, generalist_prediction$IN_angle)
  mismatch_specialist_IN <- mapply(min_angle_diff, data.sp.plant$IN_angle, specialist_prediction$IN_angle)
  
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
    fit.power <- nlsLM(y ~ a * x^b, start = list(a = 2, b = 0.5))
    #fit.spline <- lm(y ~ ns(x, df = 4))
    #fit.lm <- lm(y ~x)
    mismatch_gen <- data.sp.plant.rev$mismatch_generalist_median[((i-1)*nrow(data.pro.future)+1):(i*nrow(data.pro.future))]
    mismatch_gen <- abs(mismatch_gen)
    predicted_value <- predict(fit.power, newdata = list(x = mismatch_gen))
    data.pro.future[,i] <- predicted_value
  }
  
  mismatch_spe <- data.sp.plant$mismatch_specialist_median
  mismatch_spe <- abs(mismatch_spe)
  y <- data.prob.summary[,7]
  x <- data.mismatch[,7]
  fit.power <- nlsLM(y ~ a * x^b, start = list(a = 2, b = 0.5))
  #fit.spline <- lm(y ~ ns(x, df = 4))
  #fit.lm <- lm(y ~x)
  predicted_spe <- predict(fit.power, newdata = list(x = mismatch_spe))
  
  data.pro.future$Andrena_violae <- predicted_spe
  
  res <- which(rownames(data.XY) %in% data.sp.plant$Gridcode == F)
  data.XY <- data.XY[-res,]
  data.pro.future$latitude <- data.XY$latitude
  write.csv(data.pro.future, file = paste("./Analyses/Results/probability/MAXNET/Median/",spnames,".csv",sep = ""))
}

############################ calculate R square of both linear and power function ######

spnames_plants <- str_remove(file1,".csv$")

data.Rsquare <- as.data.frame(matrix(NA,length(spnames_plants),7))
rownames(data.Rsquare) <- spnames_plants
colnames(data.Rsquare) <- c("Andrena_carlini", "Andrena_miserabilis","Osmia_atriventris","Osmia_bucephala","Osmia_lignaria","Osmia_pumila","Andrena_violae")

for(j in 1:nrow(data.Rsquare)){
  data.prob.summary <- read.csv(paste("./Analyses/Results/probability/MAXNET/Current/",file1[j],sep=""),row.names = 1)
  data.mismatch <- read.csv(paste("./Analyses/Results/phenology_mismatch/MAXNET/Current/",file2[j],sep = ""),row.names =1)
  for(i in 1:ncol(data.Rsquare)){
    y <- data.prob.summary[,i]
    x <- data.mismatch[,i]
    fit.power <- nlsLM(y ~ a * x^b, start = list(a = 2, b = 0.5))
    y_pred <- predict(fit.power)
    RSS <- sum((y - y_pred)^2)
    TSS <- sum((y - mean(y))^2)
    R2 <- 1 - (RSS/TSS)
    data.Rsquare[j, i] <- R2
  }
}
write.csv(data.Rsquare, file = "./Analyses/Results/data.Rsquare.power.MAXNET.csv")

for(j in 1:nrow(data.Rsquare)){
  data.prob.summary <- read.csv(paste("./Analyses/Results/probability/MAXNET/Current/",file1[j],sep=""),row.names = 1)
  data.mismatch <- read.csv(paste("./Analyses/Results/phenology_mismatch/MAXNET/Current/",file2[j],sep = ""),row.names =1)
  for(i in 1:ncol(data.Rsquare)){
    y <- data.prob.summary[,i]
    x <- data.mismatch[,i]
    fit.lm <- lm(y ~ x)
    y_pred <- predict(fit.lm)
    model_summary <- summary(fit.lm)
    r_squared <- model_summary$r.squared
    data.Rsquare[j, i] <- r_squared
  }
}

write.csv(data.Rsquare, file = "./Analyses/Results/data.Rsquare.lm.MAXNET.csv")


############## We used the results from power function #########

################ First, we excluded plant-pollinator combinations with low R2 ########

data.Rsquare <- read.csv("./Analyses/Results/data.Rsquare.power.MAXNET.csv",row.names = 1)
result <- as.data.frame(which(data.Rsquare  < 0.5, arr.ind = TRUE))
result$Plant <- rownames(data.Rsquare)[result$row]
result$name <- colnames(data.Rsquare)[result$col]
result <- result[, c("Plant", "name")]

file <- list.files("./Analyses/Results/probability/MAXNET/Median/")
data.res <- c()
for (i in 1:length(file)){
  data.file <- read.csv(paste("./Analyses/Results/probability/MAXNET/Median/",file[i],sep = ""))
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


file.mismatch <- list.files("./Analyses/Results/phenology_mismatch/MAXNET/Median/")
data.mismatch <- c()
for (i in 1:length(file.mismatch)){
  data.file <- read.csv(paste("./Analyses/Results/phenology_mismatch/MAXNET/Median/",file.mismatch[i],sep = ""))
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


fit.res <- glmmTMB(value ~ scale(latitude) * type + (1 | name) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = data.res_t_filter)
summary(fit.res)
fit.res_no_interaction <- glmmTMB(value ~ latitude + type + (1 | name) + (1 | Plant), ziformula = ~1,
                                  family = beta_family(link = "logit"), 
                                  data = data.res_t_filter)
anova(fit.res, fit.res_no_interaction)

fit.res_no_type <- glmmTMB(value ~ latitude + (1 | name) + (1 | Plant), 
                           family = beta_family(link = "logit"), 
                           data = data.res_t_filter)
anova(fit.res, fit.res_no_type)

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}
result <- inv_logit(-0.155707)
print(result)


######### how the relationship between secondary extinction risk and phenological mismatch varies with latitude ######## 

data.res_t_filter_ge <- data.res_t_filter[which(data.res_t_filter$type == "generalist"),]
fit.res.mismatch_ge <- glmmTMB(value ~  scale(phenological_mismatch) * scale(latitude) + (1 | name) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = data.res_t_filter_ge)

data.res_t_filter_spe <- data.res_t_filter[which(data.res_t_filter$type == "specialist"),]
fit.res.mismatch_spe <- glmmTMB(value ~ scale(phenological_mismatch) * scale(latitude) + (1 | Plant), ziformula = ~1, family = beta_family(link = "logit"), data = data.res_t_filter_spe)


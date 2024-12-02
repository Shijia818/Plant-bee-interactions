

############## species distribution model ###########

future.climate.GI <- read.csv("D:/phenology/Climate/Future Climate/csv/ssp585/GISS-E2-1-G/GISS-E2-1-G.csv", row.names = 1)
future.climate.HD <- read.csv("D:/phenology/Climate/Future Climate/csv/ssp585/HadGEM3/HadGEM3.csv", row.names = 1)
future.climate.IN <- read.csv("D:/phenology/Climate/Future Climate/csv/ssp585/INM-CM4-8/INM-CM4-8.csv", row.names = 1)

######################  make species distribution maps and run SDM ######################

biomod.spdis.plant <- read.csv("./Analyses/Model/biomod.spdis.plant.csv", row.names = 1)
current.climate <- read.csv("./Analyses/Model/current_climate.csv", row.names = 1)
data.soil <- read.csv("./Analyses/Model/data.soil.csv", row.names = 1)
current.climate.plant <- cbind(current.climate, data.soil)

######################### PCA analyses #################
################# current plant climate ################
results.plant.PCA <- prcomp(current.climate.plant, scale = TRUE)
summary(results.plant.PCA)

grid.current.climate.plant <- as.data.frame(matrix(NA,nrow(current.climate.plant),6))
rownames(grid.current.climate.plant) <- rownames(current.climate.plant)
colnames(grid.current.climate.plant) <- c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6")
grid.current.climate.plant[,1:6] <- (-1*results.plant.PCA$x)[,1:6]
write.csv(grid.current.climate.plant,file = "./Analyses/SDM/Climate/Plant/current_climate.csv")



###################### future plant climate ###################
future.climate.HD.1 <- cbind(future.climate.HD, data.soil)
results.plant.future.PCA <- prcomp(future.climate.HD.1, scale = TRUE)
summary(results.plant.future.PCA)

grid.future.climate.plant <- as.data.frame(matrix(NA,nrow(future.climate.HD.1),6))
rownames(grid.future.climate.plant) <- rownames(future.climate.HD.1)
colnames(grid.future.climate.plant) <- c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6")
grid.future.climate.plant[,1:6] <- (-1*results.plant.future.PCA$x)[,1:6]
write.csv(grid.future.climate.plant,file = "./Analyses/SDM/Climate/Plant/future_climate_HD.csv")



######################  build species distribution models ##############
###################### MaxEnt, GLM, GAM ##############

biomod.spdis.plant <- read.csv("./Analyses/Model/biomod.spdis.plant.csv",row.names = 1)
data.XY <- read.csv("D:/phenology/Data/datasets/grid40km.csv", row.names = 1)
data.climate <- read.csv("./Analyses/SDM/Climate/Plant/current_climate.csv", row.names = 1)
GI85 <- read.csv("./Analyses/SDM/Climate/Plant/future_climate_GI.csv", row.names = 1)
HD85 <- read.csv("./Analyses/SDM/Climate/Plant/future_climate_HD.csv", row.names = 1)
IN85 <- read.csv("./Analyses/SDM/Climate/Plant/future_climate_IN.csv", row.names = 1)

data.sp <- as.data.frame(biomod.spdis.plant)
setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/true_absence")

for(i in 1:ncol(data.sp)){
  myRespName<- colnames(data.sp)[i]
  myResp<-as.numeric(data.sp[,myRespName])
  myRespxy<-data.XY
  myExpl<-data.climate
  # system("icacls .\\Viola.hirsutula\\ /grant Everyone:(OI)(CI)F /T")
  
  myBiomodData<-BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespxy,
                                     resp.name = myRespName)
  
  
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                      modeling.id = paste(myRespName,"Firstmodeling",sep=""),
                                      models = c("GAM", "GLM", "MAXNET"),
                                      CV.strategy = "random",
                                      CV.nb.rep = 30,
                                      CV.perc = 0.8,
                                      OPT.strategy = "bigboss",
                                      metric.eval  = c('TSS','ROC'),
                                      var.import = 3,
                                      CV.do.full.models = FALSE)
  
  myBiomodProjection <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          models.chosen = 'all',
                                          compress = T,
                                          build.clamping.mask = TRUE)
  
  Proj_GI <-BIOMOD_Projection(bm.mod = myBiomodModelOut, models.chosen = 'all', new.env = GI85,
                              proj.name = "2070_ssp85_GI", metric.binary = "TSS", build.clamping.mask = TRUE, compress = TRUE)
  Proj_HD <-BIOMOD_Projection(bm.mod = myBiomodModelOut, models.chosen = 'all', new.env = HD85,
                              proj.name = "2070_ssp85_HD", metric.binary = "TSS", build.clamping.mask = TRUE, compress = TRUE)
  Proj_IN <-BIOMOD_Projection(bm.mod = myBiomodModelOut, models.chosen = 'all', new.env = IN85,
                              proj.name = "2070_ssp85_IN", metric.binary = "TSS", build.clamping.mask = TRUE, compress = TRUE)
  
}

################ transform true absences into potential pseudo-absences ##### 
setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/pseudo_absence")
for(i in 1:ncol(data.sp)){
  myRespName <- colnames(data.sp)[i]
  myResp <- as.numeric(data.sp[,myRespName])
  myRespxy <- data.XY
  myExpl <- data.climate
  
  myResp.PA <- ifelse(myResp == 1, 1, NA)
  
  res <- length(which(myResp == 1))
  
  myBiomodData<-BIOMOD_FormatingData(resp.var = myResp.PA,
                                     expl.var = myExpl,
                                     resp.xy = myRespxy,
                                     resp.name = myRespName,
                                     PA.nb.rep = 10,
                                     PA.nb.absences = 2 * res,
                                     PA.strategy = 'random') 
  
  
  myBiomodModelOut <- BIOMOD_Modeling(bm.format  = myBiomodData,
                                      modeling.id = paste(myRespName,"Firstmodeling",sep=""),
                                      models = c("GAM", "GLM", "MAXNET"),
                                      CV.strategy = "random",
                                      CV.nb.rep = 10,
                                      CV.perc = 0.8,
                                      OPT.strategy = "bigboss",
                                      metric.eval  = c('TSS','ROC'),
                                      var.import = 3,
                                      CV.do.full.models = FALSE)
  
  myBiomodProjection <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = 'all',
                                          compress = T,
                                          build.clamping.mask = TRUE)
  
  Proj_GI <- BIOMOD_Projection(bm.mod = myBiomodModelOut, models.chosen = 'all', new.env = GI85,
                               proj.name = "2070_ssp85_GI", metric.binary = "TSS", build.clamping.mask = TRUE, compress = TRUE)
  Proj_HD <- BIOMOD_Projection(bm.mod = myBiomodModelOut, models.chosen = 'all', new.env = HD85,
                               proj.name = "2070_ssp85_HD", metric.binary = "TSS", build.clamping.mask = TRUE, compress = TRUE)
  Proj_IN <- BIOMOD_Projection(bm.mod = myBiomodModelOut, models.chosen = 'all', new.env = IN85,
                               proj.name = "2070_ssp85_IN", metric.binary = "TSS", build.clamping.mask = TRUE, compress = TRUE)
  
} 

#################################### deal with distribution results #################

setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/true_absence")
spnames <- list.files()
biomod.spdis.plant <- read.csv("D:/Mutualistic_interaction/Analyses/Model/biomod.spdis.plant.csv", row.names =1)

for(i in 1:length(spnames)){
  sp.names <- spnames[i]
  dat.current <- matrix(NA,1158,1)
  rownames(dat.current) <- rownames(biomod.spdis.plant)
  colnames(dat.current) <- sp.names
  model.out <- get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  model.value <- get_evaluations(model.out)
  TSS <- model.value[which(model.value$metric.eval == "TSS"),] 
  prediction <- get(load(paste(sp.names,"/","proj_current","/","proj_current_",sp.names,".RData",sep="")))
  
  ########################### deal with GAM ################
  
  data.gam <- TSS[which(TSS$algo == "GAM"),]
  res.gam <-which(data.gam$validation < 0.5)
  if(length(res.gam)!= 0){TSS.gam <- data.gam[-res.gam,]}else{TSS.gam <- data.gam}
  
  run.names <- TSS.gam$run
  prediction.gam <- prediction[which(prediction$algo == "GAM"),]
  prediction.gam.1 <- prediction.gam[which(prediction.gam$run %in% run.names),][,c("run","algo","points","pred")]
  prediction.gam.final <- pivot_wider(prediction.gam.1, names_from = run, values_from = pred)
  pred.median <- apply(prediction.gam.final[,3:ncol(prediction.gam.final), drop = FALSE],1,median)
  
  threshold <- bm_FindOptimStat(metric.eval = "TSS", obs = biomod.spdis.plant[,i], fit = pred.median, nb.thresh = 200, mpa.perc = 0.9)
  presence <- bm_BinaryTransformation(pred.median, threshold$cutoff, do.filtering = FALSE)
  
  dat.current[,1] <- pred.median/1000
  save(dat.current, file = paste("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/GAM/true_absence/probability/","proj_current_",sp.names,".RData",sep=""))
}

###########################################

setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/GAM/true_absence/probability")
sp <- list.files()
dis.sp <- c()
for(i in 1:length(sp)){dis.sp.t <- get(load(sp[i]))
dis.sp<-cbind(dis.sp,dis.sp.t)}
data.sp<-dis.sp[,match(colnames(biomod.spdis.plant),colnames(dis.sp))]
write.csv(data.sp,file = "proj_current_GAM.csv")


############################## get future distribution #################

setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/true_absence")
spnames <- list.files()

for(i in 1:length(spnames)){
  sp.names <- spnames[i]
  dat.future <- matrix(NA,1158,1)
  rownames(dat.future) <- rownames(biomod.spdis.plant)
  colnames(dat.future) <- sp.names
  model.out <- get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  model.value <- get_evaluations(model.out)
  TSS <- model.value[which(model.value$metric.eval == "TSS"),] 
  prediction <- get(load(paste(sp.names,"/","proj_current","/","proj_current_",sp.names,".RData",sep="")))
  prediction.IN <- get(load(paste(sp.names,"/","proj_2070_ssp85_IN","/","proj_2070_ssp85_IN_",sp.names,".RData",sep="")))
  
  
  data.gam <- TSS[which(TSS$algo == "GAM"),]
  res.gam <-which(data.gam$validation < 0.5)
  if(length(res.gam)!= 0){TSS.gam <- data.gam[-res.gam,]}else{TSS.gam <- data.gam}
  
  run.names <- TSS.gam$run
  prediction.gam <- prediction[which(prediction$algo == "GAM"),]
  prediction.gam.1 <- prediction.gam[which(prediction.gam$run %in% run.names),][,c("run","algo","points","pred")]
  prediction.gam.final <- pivot_wider(prediction.gam.1, names_from = run, values_from = pred)
  pred.median <- apply(prediction.gam.final[,3:ncol(prediction.gam.final), drop = FALSE],1,median)
  threshold <- bm_FindOptimStat(metric.eval = "TSS", obs = biomod.spdis.plant[,i], fit = pred.median, nb.thresh = 200, mpa.perc = 0.9)
  
  
  prediction.future.gam <- prediction.IN[which(prediction.IN$algo == "GAM"),]
  prediction.future.gam.1 <- prediction.future.gam[which(prediction.future.gam$run %in% run.names),][,c("run","algo","points","pred")]
  prediction.future.gam.final <- pivot_wider(prediction.future.gam.1, names_from = run, values_from = pred)
  pred.future.median <- apply(prediction.future.gam.final[,3:ncol(prediction.future.gam.final), drop = FALSE],1,median)
  presence <- bm_BinaryTransformation(pred.future.median, threshold$cutoff, do.filtering = FALSE)
  
  dat.future[,1] <- pred.future.median/1000
  save(dat.future, file = paste("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/GAM/true_absence/probability/","proj_future_IN_",sp.names,".RData",sep=""))
  
}

setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/GAM/true_absence/probability")
sp <- list.files()
dis.sp <- c()
for(i in 4:length(sp)){dis.sp.t <- get(load(sp[i]))
dis.sp<-cbind(dis.sp,dis.sp.t)}
data.sp<-dis.sp[,match(colnames(biomod.spdis.plant),colnames(dis.sp))]
write.csv(data.sp,file = "proj_future_IN_GAM.csv")


#################################### extract critical parameters of SDM ###########
################# TSS, true_absence ########
setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/true_absence")
spnames <- list.files()
data.model <- as.data.frame(matrix(NA, 90, 25))
colnames(data.model) <- c("RUNS","algo",spnames)
runs <- paste("RUN", 1:30, sep = "")
repeated_runs <- rep(runs, each = 3)
data.model$RUNS <- repeated_runs
algo <- c("GAM", "GLM", "MAXNET")
data.model$algo <- rep(algo, times = 10)

for(i in 1:length(spnames)){
  sp.names <- spnames[i]
  model.out <- get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  model.value <- get_evaluations(model.out)
  TSS <- model.value[which(model.value$metric.eval == "TSS"),]
  data.model[, i+2] <- TSS$validation
}

#write.csv(data.model, file = "D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/TSS_true.csv")

data.summary <- as.data.frame(matrix(NA, 23, 4))
colnames(data.summary) <- c("Species", "GAM", "GLM", "MAXNET")
data.summary$Species <- spnames

for(i in 1:length(spnames)){
  dat <- data.model[,c(1:2, i+2), drop = FALSE]
  res.gam <- which(dat$algo == "GAM" & dat[,3] >= 0.5)
  res.glm <- which(dat$algo == "GLM" & dat[,3] >= 0.5)
  res.maxnet <- which(dat$algo == "MAXNET" & dat[,3] >= 0.5)
  
  data.summary[i,2] <- length(res.gam)
  data.summary[i,3] <- length(res.glm)
  data.summary[i,4] <- length(res.maxnet)
  
} 
write.csv(data.summary, file = "D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/TSS_true_summary.csv")  

################################## TSS pseudo_absence ############
setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/pseudo_absence")
spnames <- list.files()

data.model <- as.data.frame(matrix(NA, 300, 25))
colnames(data.model) <- c("RUNS","algo",spnames)
PA <- paste("PA", 1:10, sep = "")
RUN <- paste("RUN", 1:10, sep = "")

combinations <- outer(PA, RUN, paste, sep = "_")
combinations <- as.vector(t(combinations))
repeated_combinations <- rep(combinations, each = 3)
data.model$RUNS <- repeated_combinations

algo <- c("GAM", "GLM", "MAXNET")
data.model$algo <- rep(algo, times = 100)

for(i in 1:length(spnames)){
  sp.names <- spnames[i]
  model.out <- get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  model.value <- get_evaluations(model.out)
  TSS <- model.value[which(model.value$metric.eval == "TSS"),]
  data.model[, i+2] <- TSS$validation
}

write.csv(data.model, file = "D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/TSS_pseudo.csv")

data.summary <- as.data.frame(matrix(NA, 23, 4))
colnames(data.summary) <- c("Species", "GAM", "GLM", "MAXNET")
data.summary$Species <- spnames

for(i in 1:length(spnames)){
  dat <- data.model[,c(1:2, i+2), drop = FALSE]
  res.gam <- which(dat$algo == "GAM" & dat[,3] >= 0.5)
  res.glm <- which(dat$algo == "GLM" & dat[,3] >= 0.5)
  res.maxnet <- which(dat$algo == "MAXNET" & dat[,3] >= 0.5)
  
  data.summary[i,2] <- length(res.gam)
  data.summary[i,3] <- length(res.glm)
  data.summary[i,4] <- length(res.maxnet)
  
} 
write.csv(data.summary, file = "D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/TSS_pseudo_summary.csv")  

#################### get TSS values ##############
setwd("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/pseudo_absence")
spnames <- list.files()

for(i in 1:length(spnames)){
  sp.names <- spnames[i]
  model.out <- get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  model.value <- get_evaluations(model.out)
  TSS <- model.value[which(model.value$metric.eval == "TSS"),][ ,c(1,4,9,10)]
  TSS$species <- sp.names
  TSS <- TSS[, c("species", setdiff(names(TSS), "species"))]
  write.csv(TSS, file = paste("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/TSS_pseudo/",sp.names,".csv", sep = ""))
}

proj_future_GLM <- read.csv("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/GLM/pseudo_absence/probability/proj_future_HD_GLM.csv", row.names =1)
proj_future_GAM <- read.csv("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/GAM/pseudo_absence/probability/proj_future_HD_GAM.csv",row.names =1)
proj_future_MAXNET <- read.csv("D:/Mutualistic_interaction/Analyses/SDM/Results/Plant/MAXNET/pseudo_absence/probability/proj_future_HD_MAXNET.csv",row.names =1)
cor_results <- mapply(cor, proj_future_GAM, proj_future_MAXNET)








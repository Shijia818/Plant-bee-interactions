

library(ecospat)
setwd("D:/Mutualistic_interaction/Climate/total")

data1 <- read.csv("./Andrena_carlini/Tem_ave.csv")
data2 <- read.csv("./Andrena_miserabilis/Tem_ave.csv")
data3 <- read.csv("./Andrena_violae/Tem_ave.csv")
data4 <- read.csv("./Osmia_atriventris/Tem_ave.csv")
data5 <- read.csv("./Osmia_bucephala/Tem_ave.csv")
data6 <- read.csv("./Osmia_lignaria/Tem_ave.csv")
data7 <- read.csv("./Osmia_pumila/Tem_ave.csv")

data.total.ave <- rbind(data1,data2,data3,data4,data5,data6,data7)
bio1 <- apply(data.total.ave[,15:26],1,mean)
bio4 <- apply(data.total.ave[,15:26],1,sd)

T.max.month <- apply(data.total.ave[,15:26],1,which.max)
T.min.month <- apply(data.total.ave[,15:26],1,which.min)


data1.1 <- read.csv("./Andrena_carlini/Tem_max.csv")
data2.1 <- read.csv("./Andrena_miserabilis/Tem_max.csv")
data3.1 <- read.csv("./Andrena_violae/Tem_max.csv")
data4.1 <- read.csv("./Osmia_atriventris/Tem_max.csv")
data5.1 <- read.csv("./Osmia_bucephala/Tem_max.csv")
data6.1 <- read.csv("./Osmia_lignaria/Tem_max.csv")
data7.1 <- read.csv("./Osmia_pumila/Tem_max.csv")

data.total.max <- rbind(data1.1,data2.1,data3.1,data4.1,data5.1,data6.1,data7.1)
df.max <- data.total.max[,15:26]
values.max <- mapply(function(row, col)df.max[row, col], 1:nrow(df.max),T.max.month)


data1.2 <- read.csv("./Andrena_carlini/Tem_min.csv")
data2.2 <- read.csv("./Andrena_miserabilis/Tem_min.csv")
data3.2 <- read.csv("./Andrena_violae/Tem_min.csv")
data4.2 <- read.csv("./Osmia_atriventris/Tem_min.csv")
data5.2 <- read.csv("./Osmia_bucephala/Tem_min.csv")
data6.2 <- read.csv("./Osmia_lignaria/Tem_min.csv")
data7.2 <- read.csv("./Osmia_pumila/Tem_min.csv")


data.total.min <- rbind(data1.2,data2.2,data3.2,data4.2,data5.2,data6.2,data7.2)
df.min <- data.total.min[,15:26]
values.min <- mapply(function(row, col)df.min[row, col], 1:nrow(df.min),T.min.month)

diff <- values.max - values.min


################## determine bee probability of occurrence #######

data.bee.probability <- read.csv("./Analyses/Bee_probability/data.bee.final.csv")
unique_species <- unique(data.bee.probability$Species)
df <- data.bee.probability[,c("Species","longitude","latitude","BIO4","BIO10","BIO15","BIO16","BIO17")]
#partition_data <- function(df, n_splits = 100, calibration_ratio = 0.8)

n_splits = 100
calibration_ratio = 0.8

  for(i in 1:length(unique_species)){
  #mpa <- c()
  species_splits <- list()
  spnames <- unique_species[i]
  sp_data <- df[df$Species == spnames, ]
  
  xmn_1 <- min(sp_data$longitude)
  xmx_1 <- max(sp_data$longitude)
  ymn_1 <- min(sp_data$latitude)
  ymx_1 <- max(sp_data$latitude)
  r_template_1 <- raster(nrows = 50,ncols = 50 , xmn = xmn_1, xmx = xmx_1, ymn = ymn_1, ymx = ymx_1)
  df_env_presence <- sp_data[,c("BIO4","BIO10","BIO15","BIO16","BIO17")]
  rl_presence <- lapply(df_env_presence, function(x) {
    m1 <- matrix(x, nrow = nrow(r_template_1), ncol = ncol(r_template_1), byrow = TRUE)
    r1 <- raster(m1, template = r_template_1)
    return(r1)
  })
  rs1 <- stack(rl_presence)
  
    for(j in 1:n_splits){
      index <- sample(1:nrow(sp_data), size = round(calibration_ratio * nrow(sp_data)))
      calibration_data <- sp_data[index, ]
      evaluation_data <- sp_data[-index, ]
      
      species_splits[[j]] <- list(
        calibration = calibration_data,
        evaluation = evaluation_data
      )
    }
      
  for(k in 1:n_splits){
    calibration_data <- species_splits[[k]]$calibration
    
    xmn <- min(calibration_data$longitude)
    xmx <- max(calibration_data$longitude)
    ymn <- min(calibration_data$latitude)
    ymx <- max(calibration_data$latitude)
    r_template <- raster(nrows = 50,ncols = 50, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
    df_env <- calibration_data[,c("BIO4","BIO10","BIO15","BIO16","BIO17")]
    rl <- lapply(df_env, function(x) {
      m <- matrix(x, nrow = nrow(r_template), ncol = ncol(r_template), byrow = TRUE)
      r <- raster(m, template = r_template)
      return(r)
    })
    rs <- stack(rl)
    
    model  <- maxlike(formula = ~ BIO4 + BIO10 + BIO15 + BIO16 + BIO17,  # 模型公式
      points = cbind(calibration_data$longitude, calibration_data$latitude),  # 存在点数据
      rasters = rs,  # 栅格栈
      link = "logit"
    )
   
    psi.hat <- predict(model, newdata = rs1)
    coordinates(sp_data) <- ~longitude+latitude
    proj4string(sp_data) <- proj4string(psi.hat)
    extracted_values <- raster::extract(psi.hat, sp_data)
    
    sp_data <- data.frame(sp_data)
    
    sp_data$predicted_probabilities <- extracted_values
    
    obs <- sp_data$predicted_probabilities
    mpa.value <- ecospat.mpa(obs,perc = 0.9)
    if(mpa.value < 0.7){
      write.csv(sp_data, paste("./Analyses/Bee_probability/Results/",spnames,"/","file_",k,".csv",sep = ""))}
    #mpa[k] <- mpa.value
  }
  #write.csv(data.frame(mpa), file = paste("./Analyses/Bee_probability/MPA/",spnames,".csv",sep = ""))
}
  #for(r in 1:length(mpa)){
    #if(mpa[r] <= quantile(mpa,0.25)){
      #write.csv(sp_data, paste("./Analyses/Bee_probability/Results/",spnames,"/","file_", r,".csv",sep = ""))}
    #}
  #}    

files <- list.files("./Analyses/Bee_probability/Results")
spnames <- files
for(i in 1:length(spnames)){
  res <- list.files(paste("./Analyses/Bee_probability/Results/",spnames[i],sep = ""))
  file1 <- read.csv(paste("./Analyses/Bee_probability/Results/",spnames[i],"/", res[1],sep = ""), row.names = 1)
  prediction <- file1$predicted_probabilities
  for(j in 2:length(res)){
    data <- read.csv(paste("./Analyses/Bee_probability/Results/",spnames[i],"/", res[j],sep = ""), row.names = 1)
    prediction.1 <- data$predicted_probabilities
    prediction <- cbind(prediction, prediction.1)
  }
 data.final <- as.data.frame(apply(prediction,1,median))
 colnames(data.final) <- spnames[i]
 write.csv(data.final, file = paste("./Analyses/Bee_probability/Results/",spnames[i],".csv",sep = ""))
}
      
data.bee.final <- read.csv("./Analyses/Results_bio4/Bee_probability/data.bee.final.csv")
spnames <- unique(data.bee.final$Species)
for(i in 1:length(spnames)){
  data <- data.bee.final[which(data.bee.final$Species == spnames[i]),]
  res <- quantile(data$Probability, c(0.05,0.95))
  pos <- which(data$Probability >= res[1] & data$Probability <= res[2])
  data.final <- data[pos, ]
  write.csv(data.final, file = paste("./Analyses/Results_bio4/Bee_probability/Results/",spnames[i],".csv",sep = "")) 
}

file <- list.files("./Analyses/Results_bio4/Bee_probability/Results")
data.bee.filter <- c()
for(j in 1:length(file)){
  dat <- read.csv(paste("./Analyses/Results_bio4/Bee_probability/Results/",file[j],sep = ""))
  data.bee.filter <- rbind(data.bee.filter,dat)
}










## Author: Karthik Srinivasan
## Submission date: Dec 13,2016.
## Code for course project - STAT 574 
## Title: "FINDING THE RIGHT CUMULATIVE LAG FOR ENVIRONMENT EFFECTS ON INDIVIDUAL WELLBEING
#          - An application of lasso-based feature extraction for linear mixed effect models"
#     

  library(nlme)
  library(mgcv)
  library(lmerTest)
  library(tseries)
  library(plm)
  library(foreign)
  library(gplots)
  library(data.table)
  library(Hmisc)
  library(gtools)
  library(randomForest)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(psych)
  library(lme4)
  library(MASS)
  library(glmnet)
  library(lars)
  library(psychometric)
  library(Metrics) 
  library(MuMIn)

  data_f2 <- read.csv("STAT574datset_KarthikS.csv")
  data_f2$P_ID <- as.factor(data_f2$P_ID)

## Four different outcomes in dataset - Imputing missing values

  data_f2$RMSSD_5[is.na(data_f2$RMSSD_5)] <- data_f2$RMSSD[is.na(data_f2$RMSSD_5)]
  data_f2$SDNN_5[is.na(data_f2$SDNN_5)] <- data_f2$SDNN[is.na(data_f2$SDNN_5)]
  data_f2$n_HF_5[is.na(data_f2$n_HF_5)] <- data_f2$n_HF[is.na(data_f2$n_HF_5)]
  data_f2$n_LF_5[is.na(data_f2$n_LF_5)] <- data_f2$n_LF[is.na(data_f2$n_LF_5)]
  data_f2$LF.HF_5[is.na(data_f2$LF.HF_5)] <- data_f2$LF.HF[is.na(data_f2$LF.HF_5)]

  data_f <- data_f2[,names(data_f2)%in% c("Date","P_ID","Sound","Temperature","Pressure","CO2","Relative_humidity","Light","Absolute_humidity","SDNN_5","RMSSD_5","LF.HF_5","n_HF_5","ToD","DoW","BMI","Age","Gender")]

  # Check the number of missing values in dataset
  sapply(data_f, function(y) sum(length(which(is.na(y)))))
    ## Light and Absolute humidity have missing values. We will handle these variables later.
  
############# LAGS - Create 5 minute lags as new columns, to create the accumulated lag values (as sum of cells in a row)

  create_lags <- function(data_in,n)
  {
    data_in <- data.table(data_in)
    data_in$Date <- as.Date(data_in$Date)
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"CO2",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(CO2,i), by=list(P_ID,Date)] 
    }
    
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"Sound",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(Sound,i), by=list(P_ID,Date)] 
    }
    
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"Relative_humidity",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(Relative_humidity,i), by=list(P_ID,Date)] 
    }
    
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"Pressure",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(Pressure,i), by=list(P_ID,Date)] 
    }
    
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"Temperature",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(Temperature,i), by=list(P_ID,Date)] 
    }
    
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"Light",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(Light,i), by=list(P_ID,Date)] 
    }
    
    for (i in 1:n)
    {
      IEQcol_temp <- paste("lag",i,"Absolute_humidity",sep=".")
      data_in <- data_in[, (IEQcol_temp):=shift(Absolute_humidity,i), by=list(P_ID,Date)] 
    }
    
    return(data_in) 
  }

## Date pre-processing ##
    data_p_count <- aggregate(RMSSD~P_ID,data=data_f2,NROW)
    dt_l <- create_lags(data_f,12)
    ## Remove participants with less than 10 points (or less than 50 minutes of data)
    dt_p_min <- dt_l[dt_l$P_ID %in% data_p_count$P_ID[data_p_count$RMSSD > 12],]
    sapply(dt_p_min, function(y) sum(length(which(is.na(y)))))
    ## Since creating Lags creates missing values (initial), we assume value hasn't changed for initial points (less than 5%)
    dt_p_min2 <- as.data.frame(dt_p_min) 

  col_names <- c("Pressure","Relative_humidity","Temperature", "CO2","Sound","Light","Absolute_humidity")
  for(name in col_names)
  {
    for(col in grep(name, names(dt_p_min)))
    {
      dt_p_min2[is.na(dt_p_min2[,col]),col] <- dt_p_min2[is.na(dt_p_min2[,col]),name]
    }
  }

  sapply(dt_p_min2, function(y) sum(length(which(is.na(y)))))  

  data_lag <- dt_p_min2
  data_lag$rown <- 1:nrow(data_lag)
  
  data_lag$ToD <- relevel(data_lag$ToD,ref= "Morning")
  data_lag$DoW <- relevel(data_lag$DoW,ref= "Monday")
  data_lag$Gender <- relevel(data_lag$Gender,ref= "Male")
  
  data_lag$co2org <- data_lag$CO2   ## To compare later
  data_lag$CO2 <- 0.01*data_lag$CO2 

## Data scaling for the model ##

  ## Step I: Scale the inputs ## Even though it doesnt matter as we are considering only one input type and its lags
  data_lag_del <- data_lag[!is.na(data_lag$lag.12.Relative_humidity ),] 
  
  ieq_list <- list()
  
  ## Take only instantaneous, 30 min. and 60 minutes cumulative lags for comparison 
  col_names <- c("Temperature", "CO2","Sound","Relative_humidity","Pressure","Light","Absolute_humidity")
  for(var_list in 1:length(col_names))
  {
    ieq_list[[var_list]] <- data_lag_del[,grep(col_names[var_list], names(data_lag_del))][,c(1,7,13)]
  }
  
  ## Average lag (from 1 to 10)
  
  acc_means <- function(dat1)
  {
    acc_num <- ncol(dat1)
    for(i in 2:acc_num)
    {
      k <- dat1[,c(1:i)]   
      dat1[,paste(names(dat1)[1],".alag.",i-1,sep="")] <- rowMeans(k,na.rm=T)
    }
    dat1 <- dat1[,-c(2:acc_num)]
    return(dat1)
  }
  
  acc_means_list <- lapply(ieq_list,function(y) acc_means(y))
  acc_means_col <- do.call(cbind,acc_means_list)
  
  data_acc_var <- cbind(data_lag_del[,c("P_ID","rown")],acc_means_col) 

  ## Group/Person mean centering for the multilevel model
  data_n <- split(data_acc_var,data_acc_var$P_ID)
  data_m <- lapply(data_n,function(y) {
    z <- sapply(y[,-c(1,2)],function(x) scale(x,center=T,scale=F))
    return(data.frame(P_ID = y[,1],rown=y[,2],z))
  }
  )
  data_lc <- do.call(rbind,data_m)
  
  names(data_lc) <- paste(names(data_acc_var),"_lc",sep="")
  data_lag_s <- left_join(data_lag_del,data_lc,by=c("rown"="rown_lc"))

  ## We remove the original variables from the dataset for simplicity
  data_lag_s_only <- data_lag_s[,-c(12:ncol(data_lag))]

####  
## Training-test data: Divide to training and test dataset (75:25) : Stratified sampling to ensure each participant is present
####
  
  data_p_count <- aggregate(RMSSD_5~P_ID,data=data_lag_s_only,NROW)
  data_lag_s_del2 <- data_lag_s_only[data_lag_s_only$P_ID %in% 
                                       data_p_count$P_ID[data_p_count$RMSSD_5 > 15],]
  data_lag_s_del2$P_ID <- factor(data_lag_s_del2$P_ID)
  
  data_subs <- split(data_lag_s_del2,data_lag_s_del2$P_ID)
  len_clusters <- sapply(data_subs, function(x) return(nrow(x)))
  
  ds_tr <- list()
  ds_te <- list()
  
  set.seed(12345)
  for(list in 1:length(data_subs))
  {
    x <- data_subs[[list]]
    sample_tr <- sample(nrow(x),size=round(0.75*nrow(x)),replace=F)
    ds_tr[[list]] <- x[sample_tr,]
    ds_te[[list]] <- x[-sample_tr,]
  }
  ds_train <- do.call(rbind,ds_tr)
  ds_test <- do.call(rbind,ds_te)

  ## Split input variables to lag groups
  
  X_list <- list()
  
  col_names <- c("Temperature", "CO2","Sound","Relative_humidity","Pressure","Light","Absolute_humidity")
  for(var_list in 1:length(col_names))
  {
    X_list[[var_list]] <- ds_train[,grep(col_names[var_list], names(ds_train))]
  }
  
  Y_list <- ds_train[,c("SDNN_5","RMSSD_5","n_HF_5","LF.HF_5")]

####  
####
## Function for feature selection in multilevel data (MAIN idea of project) 
####
####  
  find_opt_lag <- function(X,Y,cov,split_var)
    {
      ds_train_sub_X <- split(X,split_var) 
      ds_train_sub_Y <- split(Y,split_var)
      cov_flatten <- as.data.frame(model.matrix(~.,data=cov)[,-1])
      ds_train_sub_cov <- split(cov_flatten,split_var)
      coef_list <- list()
      weight <- list()
      mse <- list()
      set.seed(12345)
      for(i in 1: length(ds_train_sub_Y))
      {
        x <- as.matrix(cbind(ds_train_sub_X[[i]],ds_train_sub_cov[[i]]))
        y <- as.vector(ds_train_sub_Y[[i]])
        try(cvfit <- cv.glmnet(x, y,family="gaussian",alpha=1,nfolds=10,grouped=F))
        coef_list[[i]] <- coef(cvfit, s = "lambda.min")
        weight[[i]] <- length(y)
        mse[[i]] <- min(cvfit$cvm)
      }
      mse_v <- as.matrix(unlist(mse))
      coef_list_df <- do.call(cbind,coef_list)
      wt <- as.matrix(unlist(weight))
      adj_wt <- wt/(mse_v)^0.5
      coef_list_df_w <- (coef_list_df %*% adj_wt)/length(adj_wt) 
      return(round(coef_list_df_w,6))
    }  

## Run the feature selection for different input sets (each of these, we consider the three versions and select optimal one)
## Time of day (ToD) and Day of week (DoW) are important covariates and hence included in the models

  temperature_list_SDNN <- find_opt_lag(X_list[[1]],Y_list[[1]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  temperature_list_RMSSD <- find_opt_lag(X_list[[1]],Y_list[[2]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  temperature_list_LFHF <- find_opt_lag(X_list[[1]],Y_list[[3]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  temperature_list_n_HF <- find_opt_lag(X_list[[1]],Y_list[[4]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  
  CO2_list_SDNN <- find_opt_lag(X_list[[2]],Y_list[[1]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  CO2_list_RMSSD <- find_opt_lag(X_list[[2]],Y_list[[2]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  CO2_list_LFHF <- find_opt_lag(X_list[[2]],Y_list[[3]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  CO2_list_n_HF <- find_opt_lag(X_list[[2]],Y_list[[4]],ds_train[,c("ToD","DoW")],ds_train$P_ID)
  
  Sound_list_SDNN <- find_opt_lag(X_list[[3]],Y_list[[1]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Sound_list_RMSSD <- find_opt_lag(X_list[[3]],Y_list[[2]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Sound_list_LFHF <- find_opt_lag(X_list[[3]],Y_list[[3]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Sound_list_n_HF <- find_opt_lag(X_list[[3]],Y_list[[4]],ds_train[,c("ToD","DoW")],ds_train$P_ID)
  
  Relative_humidity_list_SDNN <- find_opt_lag(X_list[[4]],Y_list[[1]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Relative_humidity_list_RMSSD <- find_opt_lag(X_list[[4]],Y_list[[2]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Relative_humidity_list_LFHF <- find_opt_lag(X_list[[4]],Y_list[[3]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Relative_humidity_list_n_HF <- find_opt_lag(X_list[[4]],Y_list[[4]],ds_train[,c("ToD","DoW")],ds_train$P_ID)
  
  Pressure_list_SDNN <- find_opt_lag(X_list[[5]],Y_list[[1]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Pressure_list_RMSSD <- find_opt_lag(X_list[[5]],Y_list[[2]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Pressure_list_LFHF <- find_opt_lag(X_list[[5]],Y_list[[3]],ds_train[,c("ToD","DoW")],ds_train$P_ID) 
  Pressure_list_n_HF <- find_opt_lag(X_list[[5]],Y_list[[4]],ds_train[,c("ToD","DoW")],ds_train$P_ID)

## Since light and absolute humidity have missing values, we consider only the subset of data that has populated values for them

  light_in <- X_list[[6]][!is.na(X_list[[6]][,1]),]
  light_out <- lapply(Y_list,function(y) return(y[!is.na(X_list[[6]][,1])]))
  light_PID <- factor(ds_train$P_ID[!is.na(X_list[[6]][,1])])
  light_cov <- ds_train[,c("ToD","DoW")][!is.na(X_list[[6]][,1]),]
  
  Light_list_SDNN <- find_opt_lag(light_in,light_out[[1]],light_cov,light_PID) 
  Light_list_RMSSD <- find_opt_lag(light_in,light_out[[2]],light_cov,light_PID) 
  Light_list_LFHF <- find_opt_lag(light_in,light_out[[3]],light_cov,light_PID) 
  Light_list_n_HF <- find_opt_lag(light_in,light_out[[4]],light_cov,light_PID)
  
  abs_humidity_in <- X_list[[7]][!is.na(X_list[[7]][,1]),]
  abs_humidity_out <- lapply(Y_list,function(y) return(y[!is.na(X_list[[7]][,1])]))
  abs_humidity_PID <- factor(ds_train$P_ID[!is.na(X_list[[7]][,1])])
  abs_humidity_cov <- ds_train[,c("ToD","DoW")][!is.na(X_list[[7]][,1]),]
  
  Absolute_humidity_list_SDNN <- find_opt_lag(abs_humidity_in,abs_humidity_out[[1]],abs_humidity_cov,abs_humidity_PID) 
  Absolute_humidity_list_RMSSD <- find_opt_lag(abs_humidity_in,abs_humidity_out[[2]],abs_humidity_cov,abs_humidity_PID) 
  Absolute_humidity_list_LFHF <- find_opt_lag(abs_humidity_in,abs_humidity_out[[3]],abs_humidity_cov,abs_humidity_PID) 
  Absolute_humidity_list_n_HF <- find_opt_lag(abs_humidity_in,abs_humidity_out[[4]],abs_humidity_cov,abs_humidity_PID)
  
  data_out <- cbind(temperature_list_SDNN,temperature_list_RMSSD,temperature_list_n_HF,temperature_list_LFHF,
                    Sound_list_SDNN,Sound_list_RMSSD,Sound_list_n_HF,Sound_list_LFHF,
                    CO2_list_SDNN,CO2_list_RMSSD,CO2_list_n_HF,CO2_list_LFHF,
                    Pressure_list_SDNN,Pressure_list_RMSSD,Pressure_list_n_HF,Pressure_list_LFHF,
                    Relative_humidity_list_SDNN,Relative_humidity_list_RMSSD,Relative_humidity_list_n_HF,Relative_humidity_list_LFHF,
                    Absolute_humidity_list_SDNN,Absolute_humidity_list_RMSSD,Absolute_humidity_list_n_HF,Absolute_humidity_list_LFHF,
                    Light_list_SDNN,Light_list_RMSSD,Light_list_n_HF,Light_list_LFHF)

#### NOTE:                                                                #####
### This method gives us the best versions of inputs as follows:          #####
### Sound - Instantaneous (for all outcomes)                              #####
### Temperature - 60 minutes cumulative lagged version (for all outcomes) #####
### Relative humidity - 60 minutes cumulative lagged version (for all outcomes) 
### Pressure - 60 minutes cumulative lagged version (for all outcomes)    #####
### Light - 60 minutes cumulative lagged version (for all outcomes)       #####
### Absolute humidity - 60 minutes cumulative lagged version (for all outcomes)

######
### Verification - Comparing model performance to models derived from 
### other feature selection methods or pre-determined feature sets
######


## Test the performance of HLM models with variables

    error_compute <- function(Yhat,Y)
    {
      c(rmse(Yhat,Y),mae(Yhat,Y)) 
    }
    

## INSTANTANEOUS ##
  mod_RMSSD_inst <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_inst <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_inst <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_inst <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_inst),AIC(mod_SDNN_inst),AIC(mod_n_HF_inst),AIC(mod_LF.HF_inst))
  c(r.squaredGLMM(mod_RMSSD_inst)[2],r.squaredGLMM(mod_SDNN_inst)[2],r.squaredGLMM(mod_n_HF_inst)[2],r.squaredGLMM(mod_LF.HF_inst)[2])
  
  inst_error <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_inst,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_inst,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_inst,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_inst,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(inst_error) <- c("RMSE","MAE")
#write.csv(inst_error,"Instantaneous_input_error.csv",row.names=F)

## Considering all possible variables as second feature set for comparison ##
  
  mod_RMSSD_all <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc
                       , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_all <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc
                      , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_all <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc
                      , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_all <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc
                       , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_all),AIC(mod_SDNN_all),AIC(mod_n_HF_all),AIC(mod_LF.HF_all))
  c(r.squaredGLMM(mod_RMSSD_all)[2],r.squaredGLMM(mod_SDNN_all)[2],r.squaredGLMM(mod_n_HF_all)[2],r.squaredGLMM(mod_LF.HF_all)[2])
  
  
  allvar_error <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_all,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_all,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_all,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_all,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(allvar_error) <- c("RMSE","MAE")

#write.csv(allvar_error,"Allinputs_error.csv",row.names=F)

####                                               #####  
### Variables identified using our PROPOSED APPROACH ### 
####                                               #####

  mod_RMSSD_ML <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_ML <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_ML <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_ML <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_ML),AIC(mod_SDNN_ML),AIC(mod_n_HF_ML),AIC(mod_LF.HF_ML))
  c(r.squaredGLMM(mod_RMSSD_ML)[2],r.squaredGLMM(mod_SDNN_ML)[2],r.squaredGLMM(mod_n_HF_ML)[2],r.squaredGLMM(mod_LF.HF_ML)[2])
  
  
  ML_error <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_ML,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_ML,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_ML,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_ML,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(ML_error) <- c("RMSE","MAE")

#write.csv(ML_error,"Mixedlags_input_error.csv",row.names=F)

## Instataneous - only fixed effects, BASELINE!!

  mod_RMSSD_instf <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ 1 )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_instf <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ 1 )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_instf <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ 1 )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_instf <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ 1 )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_instf),AIC(mod_SDNN_instf),AIC(mod_n_HF_instf),AIC(mod_LF.HF_instf))
  c(r.squaredGLMM(mod_RMSSD_instf)[2],r.squaredGLMM(mod_SDNN_instf)[2],r.squaredGLMM(mod_n_HF_instf)[2],r.squaredGLMM(mod_LF.HF_instf)[2])
  
  inst_errorf <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_instf,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_instf,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_instf,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_instf,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(inst_errorf) <- c("RMSE","MAE")
#write.csv(inst_errorf,"Instantaneous_only_fixed_error.csv",row.names=F)

##### All 30 minutes acc. lag variables (manual feature selection) 

  mod_RMSSD_30 <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc  , random = list(P_ID=pdDiag(~ Sound.alag.1_lc +   Temperature.alag.1_lc + Relative_humidity.alag.1_lc + Pressure.alag.1_lc + CO2.alag.1_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_30 <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc  , random = list(P_ID=pdDiag(~ Sound.alag.1_lc +   Temperature.alag.1_lc + Relative_humidity.alag.1_lc + Pressure.alag.1_lc + CO2.alag.1_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_30 <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc  , random = list(P_ID=pdDiag(~ Sound.alag.1_lc +   Temperature.alag.1_lc + Relative_humidity.alag.1_lc + Pressure.alag.1_lc + CO2.alag.1_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_30 <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.1_lc +  Temperature.alag.1_lc + Relative_humidity.alag.1_lc + CO2.alag.1_lc + Pressure.alag.1_lc  , random = list(P_ID=pdDiag(~ Sound.alag.1_lc +   Temperature.alag.1_lc + Relative_humidity.alag.1_lc + Pressure.alag.1_lc + CO2.alag.1_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_30),AIC(mod_SDNN_30),AIC(mod_n_HF_30),AIC(mod_LF.HF_30))
  c(r.squaredGLMM(mod_RMSSD_30)[2],r.squaredGLMM(mod_SDNN_30)[2],r.squaredGLMM(mod_n_HF_30)[2],r.squaredGLMM(mod_LF.HF_30)[2])
  
  
  I30_error <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_30,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_30,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_30,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_30,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(I30_error) <- c("RMSE","MAE")
#write.csv(I30_error,"Instantaneous_30area_error.csv",row.names=F)

##### All 60 minutes acc. lag variables (manual feature selection) 

  mod_RMSSD_60 <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc  , random = list(P_ID=pdDiag(~ Sound.alag.2_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure.alag.2_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_60 <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc  , random = list(P_ID=pdDiag(~ Sound.alag.2_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure.alag.2_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_60 <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc  , random = list(P_ID=pdDiag(~ Sound.alag.2_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure.alag.2_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_60 <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.2_lc +  Temperature.alag.2_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure.alag.2_lc  , random = list(P_ID=pdDiag(~ Sound.alag.2_lc +   Temperature.alag.2_lc + Relative_humidity.alag.2_lc + Pressure.alag.2_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_60),AIC(mod_SDNN_60),AIC(mod_n_HF_60),AIC(mod_LF.HF_60))
  c(r.squaredGLMM(mod_RMSSD_60)[2],r.squaredGLMM(mod_SDNN_60)[2],r.squaredGLMM(mod_n_HF_60)[2],r.squaredGLMM(mod_LF.HF_60)[2])
  
  I60_error <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_60,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_60,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_60,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_60,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(I60_error) <- c("RMSE","MAE")
#write.csv(I60_error,"Instantaneous_60area_error.csv",row.names=F)

## Competing feature selection method: Stepwise approach - Code for stepwise approach not shown as a trivial plug and check model fit manual approach was used
#### Idetified following cumulative lags for Temperature, Sound, CO2, Relative humidity: 
###  2,3,3,3;1,1,1,1; 2,3,3,3; 1,1,3,2 for RMSSD, SDNN, n_HF, LF/HF

  mod_RMSSD_SW <- lme(fixed = RMSSD_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.2_lc +  Temperature.alag.1_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound.alag.2_lc +   Temperature.alag.1_lc + Relative_humidity.alag.2_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_SDNN_SW <- lme(fixed = SDNN_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity_lc + CO2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity_lc + Pressure_lc + CO2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_n_HF_SW <- lme(fixed = n_HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound.alag.2_lc +  Temperature.alag.1_lc + Relative_humidity.alag.2_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound.alag.2_lc +   Temperature.alag.1_lc + Relative_humidity.alag.2_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  mod_LF.HF_SW <- lme(fixed = LF.HF_5 ~ 1+  ToD + DoW + Age + BMI + Gender  + Sound_lc +  Temperature_lc + Relative_humidity.alag.1_lc + CO2.alag.2_lc + Pressure_lc  , random = list(P_ID=pdDiag(~ Sound_lc +   Temperature_lc + Relative_humidity.alag.1_lc + Pressure_lc + CO2.alag.2_lc )),data=ds_train,control=lmeControl(opt="optim"),method="ML",na.action=na.exclude)  
  
  c(AIC(mod_RMSSD_SW),AIC(mod_SDNN_SW),AIC(mod_n_HF_SW),AIC(mod_LF.HF_SW))
  c(r.squaredGLMM(mod_RMSSD_SW)[2],r.squaredGLMM(mod_SDNN_SW)[2],r.squaredGLMM(mod_n_HF_SW)[2],r.squaredGLMM(mod_LF.HF_SW)[2])
  
  
  ISW_error <- data.frame(
    rbind(
      error_compute(predict(mod_RMSSD_SW,newdata = ds_test),ds_test$RMSSD_5),
      error_compute(predict(mod_SDNN_SW,newdata = ds_test),ds_test$SDNN_5),
      error_compute(predict(mod_n_HF_SW,newdata = ds_test),ds_test$n_HF_5),
      error_compute(predict(mod_LF.HF_SW,newdata = ds_test),ds_test$LF.HF_5)
    )
  )
  names(ISW_error) <- c("RMSE","MAE")
#write.csv(ISW_error,"Stepwise_error.csv",row.names=F)



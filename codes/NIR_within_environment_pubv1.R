####################### Within-Environment Prediction ##########################
################################################################################

library(BGLR)
library(dplyr)
library(stringr)
library(doBy)
library(caret)
library(tidyverse)
library(splitstackshape)
library(dplyr)
library(BGGE)
library(prospectr)
library(progress)

traitnames <- c("yield", "pht", "dy")

pheno_geno_NIR = read.csv("blues_pheno_NIR.csv", check.names = FALSE)
names(pheno_geno_NIR) = tolower(names(pheno_geno_NIR))

envnames = c("18COL", "18CS", "18GC", "18VC", "19COL", "19CS", "19TA", "19VC")

# env =4
for (env in 1:length(envnames)) {
  
  if (env == 1) {
    pheno <- subset(pheno_geno_NIR, env == "18COL")
  } else if (env ==2) {
    pheno <- subset(pheno_geno_NIR, env == "18CS")
  } else if (env ==3) {
    pheno <- subset(pheno_geno_NIR, env == "18GC")
  }else if (env ==4) {
    pheno <- subset(pheno_geno_NIR, env == "18VC")
  }else if (env ==5) {
    pheno <- subset(pheno_geno_NIR, env == "19COL")
  }else if (env ==6) {
    pheno <- subset(pheno_geno_NIR, env == "19CS")
  }else if (env ==7) {
    pheno <- subset(pheno_geno_NIR, env == "19TA")
  }else {
    pheno <- subset(pheno_geno_NIR, env == "19VC")
  }
  
  rownames(pheno) = pheno$pedigree
  
  NIR.d1a<-scale(savitzkyGolay(pheno[,-c(1:20)],m=1,p=1,w=11)) #First derivative
  ZN1 = tcrossprod(as.matrix(NIR.d1a)/ncol(as.matrix(NIR.d1a)))
  dim(ZN1)
  
  NIR.d2a<-scale(savitzkyGolay(pheno[,-c(1:20)],m=1,p=2,w=11)) #Second derivative
  ZN2 = tcrossprod(as.matrix(NIR.d2a)/ncol(as.matrix(NIR.d2a)))
  dim(ZN2)
  
  col = colnames(pheno)
  factors = col[1:8] # independent variables 
  numerics = col[9: length(col)] # dependent var- agronomic data, composition data and NIR wavelength
  
  
  # Constructing marker data
  
  tmp<-strsplit(pheno$pedigree,"/")
  
  P1<-rep(NA,length(tmp))
  P2<-rep(NA,length(tmp))
  
  for(i in 1:length(tmp))
  {
    P1[i]=tmp[[i]][1]
    P2[i]=tmp[[i]][2]
  }
  
  P1<-as.factor(P1)
  P2<-as.factor(P2)
  cross<-as.factor(pheno$pedigree)
  
  Z1=model.matrix(~P1-1)
  dim(Z1)
  
  Z2=model.matrix(~P2-1)
  dim(Z2)
  
  
  Z3=model.matrix(~cross-1)
  dim(Z3)
  
  load("X_G_Jales.RData")
  K1=G[levels(P1), levels(P1)]
  dim(K1)
  rownames(K1)
  K2=G[levels(P2), levels(P2)]
  dim(K2)
  rownames(K2)
  
  K3=kronecker(K1,K2)
  dim(K3)
  
  K1star=Z1%*%K1%*%t(Z1) #female
  K2star=Z2%*%K2%*%t(Z2) #male
  K3star=Z3%*%K3%*%t(Z3) #female * male
  
  
  ###cleaning data
  pheno[,factors] = lapply(pheno[,factors], as.factor)
  pheno[,numerics] = lapply(pheno[,numerics], as.numeric)
  
  ##set ETAs predictors for models
  
  
  ##### NIR1 (first derivative)
  Eta1<-list(list(K=ZN1,model="RKHS"))    #NIR1
  
  
  ##### NIR2 (second derivative)
  Eta2<-list(list(K=ZN2,model="RKHS"))     #NIR2
  
  #### only G #######
  Eta3 =   list(list(K=K1star,model="RKHS"),     #Female
                list(K=K2star,model="RKHS"),
                list(K=K3star,model="RKHS"))
  
  ### G + NIR1 ###########
  Eta4<-list(list(K=K1star,model="RKHS"),     #Female
             list(K=K2star, model = "RKHS"),  #Male
             list(K=K3star, model = "RKHS"),  #Female x Male
             list(K=ZN1,    model="RKHS"))        #NIR1
  
  ### G NIR2 ###########
  Eta5<-list(list(K=K1star,model="RKHS"),     #Female
             list(K=K2star,model = "RKHS"),  #Male
             list(K=K3star, model = "RKHS"),  #Female x Male
             list(K=ZN2, model="RKHS"))     #NIR2 x E
  
  Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5)
  
  
  # tr=3; env = env
  set.seed(1234)
  
  for (tr in 1:length(traitnames)) {
    if (tr == 1) {
      pheno_data <- read.csv("pheno_yield.csv")
    } else if (tr == 2) {
      pheno_data <- read.csv("pheno_pht.csv")
    } else {
      pheno_data <- read.csv("pheno_dy.csv")
    }
      
      for (env in 1:length(envnames)) {
        Phenotype_data1 = pheno_data[pheno_data$env == envnames[env],]
      }
  }
  
    # cross-validation #
    hybrid = as.character(unique(Phenotype_data1$pedigree))
    
    set.seed(123)
    cycles = 10
    CV2 = list()
    CV1 = list()
    
    #MODEL =3; rep_num=1
    for (MODEL in 1:length(Models)) {  
      
      for (rep_num in 1:50) {
        
        test_geno <- sample(hybrid, 30)  
        train_geno <-setdiff(hybrid, test_geno)
        
        CV_Data_1_2<-Phenotype_data1
        CV_Data_1_2$Y<-NA
        CV_Data_1_2$Y[CV_Data_1_2$pedigree%in%train_geno]<-CV_Data_1_2[,1][CV_Data_1_2$pedigree%in%train_geno] 
        
        y_t<-as.numeric(CV_Data_1_2$Y)
        fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10) #nIter=5000,burnIn=1000, thin =10
        
        CV_Data_1_2$yhat <- fit$yHat
      
        # CV1
        df_test <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% test_geno)
        CV1[[(rep_num)]] <- as.data.frame(df_test %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
        
        
        #CV2
        df_train <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% train_geno)
        CV2[[(rep_num)]] <- as.data.frame(df_train %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      }
      
      if (rep_num == 50) {
        CV1out <- plyr::ldply(CV1, data.frame)
        CV2out <- plyr::ldply(CV2, data.frame)
        write.csv(CV1out, file = paste("ACC_", envnames[env], "_", traitnames[tr],"_CV1_model", MODEL, ".csv", sep=""), row.names = F)
        write.csv(CV2out, file = paste("ACC_", envnames[env], "_", traitnames[tr],"_CV2_model", MODEL, ".csv", sep=""), row.names = F)
      }
    }
  }
####################################################################################################

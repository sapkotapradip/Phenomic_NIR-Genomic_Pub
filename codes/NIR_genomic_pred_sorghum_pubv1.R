################################### Sapkota, Pradip ####################
################################# NIRBLUP Pradip 2024 ##################
################################# Agronomic traits #####################
################################# multi- environment ###################

####################### Across-Environment Prediction ##########################
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

a = subset(pheno_geno_NIR, env =="18COL")
b = subset(pheno_geno_NIR, env =="18CS")
c = subset(pheno_geno_NIR, env =="18GC")
d = subset(pheno_geno_NIR, env =="18VC")
e = subset(pheno_geno_NIR, env =="19COL")
f = subset(pheno_geno_NIR, env =="19CS")
g = subset(pheno_geno_NIR, env =="19TA")
h = subset(pheno_geno_NIR, env =="19VC")


rownames(a) = a$pedigree
rownames(b) = b$pedigree
rownames(c) = c$pedigree
rownames(d) = d$pedigree
rownames(e) = e$pedigree
rownames(f) = f$pedigree
rownames(g) = g$pedigree
rownames(h) = h$pedigree


NIR.d1a<-scale(savitzkyGolay(a[,-c(1:20)],m=1,p=1,w=11)) #First derivative
NIR.d1b<-scale(savitzkyGolay(b[,-c(1:20)],m=1,p=1,w=11))
NIR.d1c<-scale(savitzkyGolay(c[,-c(1:20)],m=1,p=1,w=11))
NIR.d1d<-scale(savitzkyGolay(d[,-c(1:20)],m=1,p=1,w=11))
NIR.d1e<-scale(savitzkyGolay(e[,-c(1:20)],m=1,p=1,w=11))
NIR.d1f<-scale(savitzkyGolay(f[,-c(1:20)],m=1,p=1,w=11))
NIR.d1g<-scale(savitzkyGolay(g[,-c(1:20)],m=1,p=1,w=11))
NIR.d1h<-scale(savitzkyGolay(h[,-c(1:20)],m=1,p=1,w=11))

NIR.d1 = rbind(NIR.d1a, NIR.d1b, NIR.d1c, NIR.d1d, NIR.d1e, NIR.d1f, NIR.d1g, NIR.d1h)
ZN1 = tcrossprod(as.matrix(NIR.d1)/ncol(as.matrix(NIR.d1)))
dim(ZN1)

NIR.d2a<-scale(savitzkyGolay(a[,-c(1:20)],m=1,p=2,w=11)) #Second derivative
NIR.d2b<-scale(savitzkyGolay(b[,-c(1:20)],m=1,p=2,w=11))
NIR.d2c<-scale(savitzkyGolay(c[,-c(1:20)],m=1,p=2,w=11))
NIR.d2d<-scale(savitzkyGolay(d[,-c(1:20)],m=1,p=2,w=11))
NIR.d2e<-scale(savitzkyGolay(e[,-c(1:20)],m=1,p=2,w=11))
NIR.d2f<-scale(savitzkyGolay(f[,-c(1:20)],m=1,p=2,w=11))
NIR.d2g<-scale(savitzkyGolay(g[,-c(1:20)],m=1,p=2,w=11))
NIR.d2h<-scale(savitzkyGolay(h[,-c(1:20)],m=1,p=2,w=11))

NIR.d2 = rbind(NIR.d2a, NIR.d2b, NIR.d2c, NIR.d2d, NIR.d2e, NIR.d2f, NIR.d2g, NIR.d2h)
ZN2 = tcrossprod(as.matrix(NIR.d2)/ncol(as.matrix(NIR.d2)))
dim(ZN2)

col = colnames(pheno_geno_NIR)
factors = col[1:8] # independent variables 
numerics = col[9: length(col)] # dependent var- agronomic data, composition data and NIR wavelength


# Constructing marker data
tmp<-strsplit(pheno_geno_NIR$pedigree,"/")

P1<-rep(NA,length(tmp))
P2<-rep(NA,length(tmp))

for(i in 1:length(tmp))
{
  P1[i]=tmp[[i]][1]
  P2[i]=tmp[[i]][2]
}

P1<-as.factor(P1)
P2<-as.factor(P2)
cross<-as.factor(pheno_geno_NIR$pedigree)

Site=as.factor(pheno_geno_NIR$env)

ZE=model.matrix(~Site-1)

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

ZEZEt=tcrossprod(ZE) #env

K4=K1star*ZEZEt #female *env
K5=K2star*ZEZEt #male*env
K6=K3star*ZEZEt #female*male*env

###cleaning data
pheno_geno_NIR[,factors] = lapply(pheno_geno_NIR[,factors], as.factor)
pheno_geno_NIR[,numerics] = lapply(pheno_geno_NIR[,numerics], as.numeric)

######################### interaction kernel for phenomic kernels #############
ZNZE1 = ZN1 * ZEZEt #interaction of 1st der with env
ZNZE2 = ZN2 * ZEZEt #interaction of 1st der with env


##set ETAs predictors for models


##### NIR1 (first derivative)
Eta1<-list(E=list(X=ZE,model="BRR"),      #Env
           P=list(K=ZN1,model="RKHS"),    #NIR1
           PE=list(K=ZNZE1,model="RKHS")) #NIR1 x Env


##### NIR2 (second derivative)
Eta2<-list(E=list(X=ZE,model="BRR"),       #Env
           P=list(K=ZN2,model="RKHS"),     #NIR2
           PE=list(K=ZNZE2,model="RKHS"))  #NIR2 x Env


#### only G #######
Eta3 =   list(list(X=ZE,model="BRR"),     #Env
         list(K=K1star,model="RKHS"),     #Female
         list(K=K2star,model="RKHS"),
         list(K=K3star,model="RKHS"),
         list(K=K4,model="RKHS"),
         list(K=K5,model="RKHS"),
         list(K=K6,model="RKHS"))

### G + NIR1 ###########
Eta4<-list(list(X=ZE, model = "BRR"),       #Env
           list(K=K1star,model="RKHS"),     #Female
           list(K=K2star, model = "RKHS"),  #Male
           list(K=K3star, model = "RKHS"),  #Female x Male
           list(K=K4,model="RKHS"),         #Female x Env
           list(K=K5,model="RKHS"),         #Male x Env
           list(K=K6,model="RKHS"),         #Female x Male x Env
           list(K=ZN1,model="RKHS"),       #NIR1
           list(K=ZNZE1,model="RKHS"))     #NIR1 x E

### G NIR2 ###########
Eta5<-list(list(X=ZE, model = "BRR"),       #E
           list(K=K1star,model="RKHS"),     #Female
           list(K=K2star,model = "RKHS"),  #Male
           list(K=K3star, model = "RKHS"),  #Female x Male
           list(K=K4,model="RKHS"),         #Male x Env
           list(K=K5,model="RKHS"),         #Female x Env
           list(K=K6,model="RKHS"),         #Female x Male x Env
           list(K=ZN2,model="RKHS"),       #NIR2
           list(K=ZNZE2,model="RKHS"))     #NIR2 x E

Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5)

# tr=1
set.seed(1234)
for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==2) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  # cross-validation #
  hybrid = as.character(unique(pheno$pedigree))
  Phenotype_data1 = pheno
  
  cycles = 10
  CV2 = list()
  CV1 = list()
  CV0_18CS = list()
  CV00_18CS = list()
  CV0_19COL = list()
  CV00_19COL = list()
  
  #MODEL =3; rep_num=2
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
      CV1[[(rep_num)]] <- as.data.frame(df_test %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      
      
      #CV2
      df_train <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% train_geno)
      CV2[[(rep_num)]] <- as.data.frame(df_train %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      
      # Untested Environment as "18CS"###################
      ################################################
      
      CV_Data_00_0<-Phenotype_data1
      CV_Data_00_0$Y0<-NA
      
      CV_Data_00_0$Y0[CV_Data_00_0$pedigree%in%train_geno] <- CV_Data_00_0[,1][CV_Data_00_0$pedigree%in%train_geno]
      CV_Data_00_0$Y0[CV_Data_00_0$env ==  '18CS' ] <- NA 
      y_t1<-as.numeric(CV_Data_00_0$Y0)
      
      fit1<-BGLR(y=y_t1,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      CV_Data_00_0$yhat1 <- fit1$yHat
      
      # CV0 #
      df_train1 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% train_geno)
      CV0_18CS[[(rep_num)]] <- as.data.frame(df_train1 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat1,use = "complete.obs")))
      
      #CV_00 #
      df_test1 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% test_geno)
      CV00_18CS[[(rep_num)]] <- as.data.frame(df_test1 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat1,use = "complete.obs")))
      
      # Untested Environment as "19COL"##################
      ###############################################
      
      CV_Data_00_0<-Phenotype_data1
      CV_Data_00_0$Y0<-NA
      
      CV_Data_00_0$Y0[CV_Data_00_0$pedigree%in%train_geno] <- CV_Data_00_0[,1][CV_Data_00_0$pedigree%in%train_geno]
      CV_Data_00_0$Y0[CV_Data_00_0$env ==  '19COL' ] <- NA 
      y_t2<-as.numeric(CV_Data_00_0$Y0)
      
      fit2<-BGLR(y=y_t2,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      CV_Data_00_0$yhat2 <- fit2$yHat
      
      # CV0 #
      df_train2 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% train_geno)
      CV0_19COL[[(rep_num)]] <- as.data.frame(df_train2 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat2,use = "complete.obs")))
      
      #CV_00 #
      df_test2 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% test_geno)
      CV00_19COL[[(rep_num)]] <- as.data.frame(df_test2 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat2,use = "complete.obs")))
      
      }
    
    if (rep_num == 50) {
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      write.csv(CV1out, file = paste("ACC_", traitnames[tr],"_CV1_", MODEL, ".csv", sep=""), row.names = F)
      write.csv(CV2out, file = paste("ACC_", traitnames[tr],"_CV2_", MODEL, ".csv", sep=""), row.names = F)
      
      
      CV00out <- plyr::ldply(CV00_18CS, data.frame)
      CV0out  <- plyr::ldply(CV0_18CS, data.frame)
      write.csv(CV00out, file = paste("ACC_", traitnames[tr],"_CV00_18CS_model", MODEL, ".csv", sep=""), row.names = F)
      write.csv(CV0out, file = paste("ACC_", traitnames[tr],"_CV0_18CS_model", MODEL, ".csv", sep=""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_19COL, data.frame)
      CV0out  <- plyr::ldply(CV0_19COL, data.frame)
      write.csv(CV00out, file = paste("ACC_", traitnames[tr],"_CV00_19COL_model", MODEL, ".csv", sep=""), row.names = F)
      write.csv(CV0out, file = paste("ACC_", traitnames[tr],"_CV0_19COL_model", MODEL, ".csv", sep=""), row.names = F)
    }
  }
}

### Plot the result for agronomic traits
df = read.csv("Pred.ability.CV1.CV2.CV0.CV00.csv")
head(df)

df = df[df$select == "1",]

df = df %>%
  filter(trait %in% c("yield", "dy", "pht")) %>%
  mutate(
    trait = case_when(
      trait == "yield" ~ "GY",
      trait == "dy" ~ "DA",
      trait == "pht" ~ "PH",
      TRUE ~ trait
    ),
    model = case_when(
      model == "Genomic" ~ "G",
      model == "Phenomic" ~ "P",
      model == "Genomic + Phenomic" ~ "G + P",
      TRUE ~ model
    )
  )

df <- as.data.frame(df %>%  dplyr::group_by(model,CV,trait) %>% 
                      dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                       SD = sd(cor, na.rm=TRUE)))


df$CV <- factor(df$CV, levels =  c("CV2", "CV1", "CV0", "CV00"))
df$model = factor(df$model, levels = c("G", "P", "G + P"))
df$trait = factor(df$trait, levels = c("GY", "DA", "PH"))


head(df)
library(ggplot2)

p = ggplot(df, aes(trait, y=M, fill=model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(M,2)  ), hjust=3, color="white",
            position = position_dodge(0.9), angle = 90,size=3.5)+
  geom_errorbar(aes(ymin=M, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~CV)+
  scale_y_continuous("Prediction ability")+
  xlab("Agronomic traits") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),    
    axis.title.y = element_text(size=12, face="bold", colour = "black"),    
    axis.text.x = element_text(size=7, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=7, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 8, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 8, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #axis.text.x.bottom = element_blank()
  )



p = p + guides(fill = guide_legend(title = "Models"))

jpeg("Pred.ability_agronomic_traits.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()



######################################################################################################
######################################################################################################
### This script is used only for grain characterstics traits; which is only available for four out
## of eight environments 18CS, 18VC, 19COL, 19VC

################################### Sapkota, Pradip ###########################
################################# NIRBLUP Pradip 2023 #########################
################################# Kernel_physical factors #####################
################################# multi- environment ##########################

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

traitnames <- c("kernel_hardness" , "kernel_diameter", "kernel_weight")

pheno_geno_NIR = read.csv("blues_pheno_NIR_physical_factors.csv", check.names = FALSE)
names(pheno_geno_NIR) = tolower(names(pheno_geno_NIR))

a = subset(pheno_geno_NIR, env =="18CS")
b = subset(pheno_geno_NIR, env =="18VC")
c = subset(pheno_geno_NIR, env =="19COL")
d = subset(pheno_geno_NIR, env =="19VC")

rownames(a) = a$pedigree
rownames(b) = b$pedigree
rownames(c) = c$pedigree
rownames(d) = d$pedigree


NIR.d1a<-scale(savitzkyGolay(a[,-c(1:22)],m=1,p=1,w=11)) #First derivative
NIR.d1b<-scale(savitzkyGolay(b[,-c(1:22)],m=1,p=1,w=11))
NIR.d1c<-scale(savitzkyGolay(c[,-c(1:22)],m=1,p=1,w=11))
NIR.d1d<-scale(savitzkyGolay(d[,-c(1:22)],m=1,p=1,w=11))

NIR.d1 = rbind(NIR.d1a, NIR.d1b, NIR.d1c, NIR.d1d)
ZN1 = tcrossprod(as.matrix(NIR.d1)/ncol(as.matrix(NIR.d1)))
dim(ZN1)

NIR.d2a<-scale(savitzkyGolay(a[,-c(1:23)],m=1,p=2,w=11)) #Second derivative
NIR.d2b<-scale(savitzkyGolay(b[,-c(1:23)],m=1,p=2,w=11))
NIR.d2c<-scale(savitzkyGolay(c[,-c(1:23)],m=1,p=2,w=11))
NIR.d2d<-scale(savitzkyGolay(d[,-c(1:23)],m=1,p=2,w=11))

NIR.d2 = rbind(NIR.d2a, NIR.d2b, NIR.d2c, NIR.d2d)
ZN2 = tcrossprod(as.matrix(NIR.d2)/ncol(as.matrix(NIR.d2)))
dim(ZN2)

col = colnames(pheno_geno_NIR)
factors = col[1:7] # independent variables 
numerics = col[8: length(col)] # dependent var- agronomic data, composition data and NIR wavelength


# Constructing marker data

tmp<-strsplit(pheno_geno_NIR$pedigree,"/")

P1<-rep(NA,length(tmp))
P2<-rep(NA,length(tmp))

for(i in 1:length(tmp))
{
  P1[i]=tmp[[i]][1]
  P2[i]=tmp[[i]][2]
}

P1<-as.factor(P1)
P2<-as.factor(P2)
cross<-as.factor(pheno_geno_NIR$pedigree)

Site=as.factor(pheno_geno_NIR$env)

ZE=model.matrix(~Site-1)

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

ZEZEt=tcrossprod(ZE) #env

K4=K1star*ZEZEt #female *env
K5=K2star*ZEZEt #male*env
K6=K3star*ZEZEt #female*male*env

###cleaning data
pheno_geno_NIR[,factors] = lapply(pheno_geno_NIR[,factors], as.factor)
pheno_geno_NIR[,numerics] = lapply(pheno_geno_NIR[,numerics], as.numeric)

######################### interaction kernel for phenomic kernels #############
ZNZE1 = ZN1 * ZEZEt #interaction of 1st der with env
ZNZE2 = ZN2 * ZEZEt #interaction of 1st der with env


##set ETAs predictors for models


##### NIR1 (first derivative)
Eta1<-list(E=list(X=ZE,model="BRR"),      #Env
           P=list(K=ZN1,model="RKHS"),    #NIR1
           PE=list(K=ZNZE1,model="RKHS")) #NIR1 x Env


##### NIR2 (second derivative)
Eta2<-list(E=list(X=ZE,model="BRR"),       #Env
           P=list(K=ZN2,model="RKHS"),     #NIR2
           PE=list(K=ZNZE2,model="RKHS"))  #NIR2 x Env


#### only G #######
Eta3 =   list(list(X=ZE,model="BRR"),     #Env
              list(K=K1star,model="RKHS"),     #Female
              list(K=K2star,model="RKHS"),
              list(K=K3star,model="RKHS"),
              list(K=K4,model="RKHS"),
              list(K=K5,model="RKHS"),
              list(K=K6,model="RKHS"))

### G + NIR1 ###########
Eta4<-list(list(X=ZE, model = "BRR"),       #Env
           list(K=K1star,model="RKHS"),     #Female
           list(K=K2star, model = "RKHS"),  #Male
           list(K=K3star, model = "RKHS"),  #Female x Male
           list(K=K4,model="RKHS"),         #Female x Env
           list(K=K5,model="RKHS"),         #Male x Env
           list(K=K6,model="RKHS"),         #Female x Male x Env
           list(K=ZN1,model="RKHS"),       #NIR1
           list(K=ZNZE1,model="RKHS"))     #NIR1 x E

### G NIR2 ###########
Eta5<-list(list(X=ZE, model = "BRR"),       #E
           list(K=K1star,model="RKHS"),     #Female
           list(K=K2star,model = "RKHS"),  #Male
           list(K=K3star, model = "RKHS"),  #Female x Male
           list(K=K4,model="RKHS"),         #Male x Env
           list(K=K5,model="RKHS"),         #Female x Env
           list(K=K6,model="RKHS"),         #Female x Male x Env
           list(K=ZN2,model="RKHS"),       #NIR2
           list(K=ZNZE2,model="RKHS"))     #NIR2 x E

Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5)

# tr=3
set.seed(1234)
for (tr in 1:length(traitnames)) {
  
  if (tr == 1) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else if (tr ==2) {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  } else {
    pheno <- read.csv(file = paste("pheno_", traitnames[tr], ".csv", sep = ""))
  }
  
  # cross-validation #
  hybrid = as.character(unique(pheno$pedigree))
  Phenotype_data1 = pheno
  
  set.seed(123)
  cycles = 10
  CV2 = list()
  CV1 = list()
  CV0_18CS = list()
  CV0_19COL = list()
  CV00_18CS = list()
  CV00_19COL = list()
  
  
  #MODEL =1; rep_num=1
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
      CV1[[(rep_num)]] <- as.data.frame(df_test %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      
      
      #CV2
      df_train <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% train_geno)
      CV2[[(rep_num)]] <- as.data.frame(df_train %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat,use = "complete.obs")))
      
      # Untested Environment as "18CS"###################
      ################################################
      
      CV_Data_00_0<-Phenotype_data1
      CV_Data_00_0$Y0<-NA
      
      CV_Data_00_0$Y0[CV_Data_00_0$pedigree%in%train_geno]<-CV_Data_00_0[,1][CV_Data_00_0$pedigree%in%train_geno]
      CV_Data_00_0$Y0[grepl("18CS",CV_Data_00_0$env )]<- NA # Train-test split based environment 18CS, 
      y_t1<-as.numeric(CV_Data_00_0$Y0)
      
      fit1<-BGLR(y=y_t1,ETA = Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      CV_Data_00_0$yhat1 <- fit1$yHat
      
      # CV0 #
      df_train1 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% train_geno)
      CV0_18CS[[(rep_num)]] <- as.data.frame(df_train1 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat1,use = "complete.obs")))
      
      #CV_00 #
      df_test1 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% test_geno)
      CV00_18CS[[(rep_num)]] <- as.data.frame(df_test1 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat1,use = "complete.obs")))
      
      # Untested Environment as "19COL"##################
      ###################################################
      
      CV_Data_00_0<-Phenotype_data1
      CV_Data_00_0$Y0<-NA
      
      CV_Data_00_0$Y0[CV_Data_00_0$pedigree%in%train_geno]<-CV_Data_00_0[,1][CV_Data_00_0$pedigree%in%train_geno]
      CV_Data_00_0$Y0[grepl("19COL",CV_Data_00_0$env )]<- NA # Train-test split based environment 18CS, 
      y_t2<-as.numeric(CV_Data_00_0$Y0)
      
      fit2<-BGLR(y=y_t2,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      CV_Data_00_0$yhat2 <- fit2$yHat
      
      # CV0 #
      df_train2 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% train_geno)
      CV0_19COL[[(rep_num)]] <- as.data.frame(df_train2 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat2,use = "complete.obs")))
      
      #CV_00 #
      df_test2 <- subset(CV_Data_00_0, CV_Data_00_0$pedigree %in% test_geno)
      CV00_19COL[[(rep_num)]] <- as.data.frame(df_test2 %>% group_by(env) %>% dplyr::summarize(cor=cor(blue, yhat2,use = "complete.obs")))
    }
    
    if (rep_num == 50) {
      CV1out <- plyr::ldply(CV1, data.frame)
      CV2out <- plyr::ldply(CV2, data.frame)
      write.csv(CV1out, file = paste("ACC_", traitnames[tr],"_CV1_", MODEL, ".csv", sep=""), row.names = F)
      write.csv(CV2out, file = paste("ACC_", traitnames[tr],"_CV2_", MODEL, ".csv", sep=""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_18CS, data.frame)
      CV0out  <- plyr::ldply(CV0_18CS, data.frame)
      write.csv(CV00out, file = paste("ACC_", traitnames[tr],"_CV00_18CS_model", MODEL, ".csv", sep=""), row.names = F)
      write.csv(CV0out, file = paste("ACC_", traitnames[tr],"_CV0_18CS_model", MODEL, ".csv", sep=""), row.names = F)
      
      CV00out <- plyr::ldply(CV00_19COL, data.frame)
      CV0out  <- plyr::ldply(CV0_19COL, data.frame)
      write.csv(CV00out, file = paste("ACC_", traitnames[tr],"_CV00_19COL_model", MODEL, ".csv", sep=""), row.names = F)
      write.csv(CV0out, file = paste("ACC_", traitnames[tr],"_CV0_19COL_model", MODEL, ".csv", sep=""), row.names = F)
    }
  }
}


#################################################NIR within env ################
############ Sapkota Pradip 2024 ###############################################
##################v2############################################################
#######################Agronomic traits#########################################

#loading libraries
library(BGLR)
library(dplyr)
library(stringr)
library(doBy)
library(caret)
library(tidyverse)
library(splitstackshape)
library(BGGE)
library(prospectr)
library(progress)

traitnames <- c("yield", "pht", "dy")
envnames = c("18COL", "18CS", "18GC", "18VC", "19COL", "19CS", "19TA", "19VC")

pheno_geno_NIR = read.csv("blues_pheno_NIR.csv", check.names = FALSE)
names(pheno_geno_NIR) = tolower(names(pheno_geno_NIR))


for (env in 1:length(envnames)){
  for (tr in 1:length(traitnames)) {
    tryCatch({
      cat("Processing environment:", envnames[env], "and trait:", traitnames[tr], "\n")
      
      if (env == 1) {
        pheno <- subset(pheno_geno_NIR, env == "18COL")
      } else if (env == 2) {
        pheno <- subset(pheno_geno_NIR, env == "18CS")
      } else if (env == 3) {
        pheno <- subset(pheno_geno_NIR, env == "18GC")
      } else if (env == 4) {
        pheno <- subset(pheno_geno_NIR, env == "18VC")
      } else if (env == 5) {
        pheno <- subset(pheno_geno_NIR, env == "19COL")
      } else if (env == 6) {
        pheno <- subset(pheno_geno_NIR, env == "19CS")
      } else if (env == 7) {
        pheno <- subset(pheno_geno_NIR, env == "19TA")
      } else {
        pheno <- subset(pheno_geno_NIR, env == "19VC")
      }
      
      cat("Phenotype data dimensions for env:", envnames[env], "are", dim(pheno), "\n")
      
      if (tr == 1) {
        pheno_data <- read.csv("pheno_yield.csv")
      } else if (tr == 2) {
        pheno_data <- read.csv("pheno_pht.csv")
      } else {
        pheno_data <- read.csv("pheno_dy.csv")
      }
      
      cat("Trait data dimensions for trait:", traitnames[tr], "are", dim(pheno_data), "\n")
      
      Phenotype_data1 = pheno_data[pheno_data$env == envnames[env],]
      cat("Filtered phenotype data dimensions for env:", envnames[env], "are", dim(Phenotype_data1), "\n")
      
      rownames(pheno) = pheno$pedigree
      
      NIR.d1a <- scale(savitzkyGolay(pheno[,-c(1:20)], m = 1, p = 1, w = 11))
      ZN1 = tcrossprod(as.matrix(NIR.d1a) / ncol(as.matrix(NIR.d1a)))
      
      NIR.d2a <- scale(savitzkyGolay(pheno[,-c(1:20)], m = 1, p = 2, w = 11))
      ZN2 = tcrossprod(as.matrix(NIR.d2a) / ncol(as.matrix(NIR.d2a)))
      
      col = colnames(pheno)
      factors = col[1:8]
      numerics = col[9:length(col)]
      
      tmp <- strsplit(pheno$pedigree, "/")
      
      P1 <- rep(NA, length(tmp))
      P2 <- rep(NA, length(tmp))
      
      for (i in 1:length(tmp)) {
        P1[i] = tmp[[i]][1]
        P2[i] = tmp[[i]][2]
      }
      
      P1 <- as.factor(P1)
      P2 <- as.factor(P2)
      cross <- as.factor(pheno$pedigree)
      
      Z1 = model.matrix(~P1 - 1)
      Z2 = model.matrix(~P2 - 1)
      Z3 = model.matrix(~cross - 1)
      
      load("X_G_Jales.RData")
      K1 = G[levels(P1), levels(P1)]
      K2 = G[levels(P2), levels(P2)]
      
      K3 = kronecker(K1, K2)
      
      K1star = Z1 %*% K1 %*% t(Z1)
      K2star = Z2 %*% K2 %*% t(Z2)
      K3star = Z3 %*% K3 %*% t(Z3)
      
      pheno[,factors] = lapply(pheno[,factors], as.factor)
      pheno[,numerics] = lapply(pheno[,numerics], as.numeric)
      
      Eta1 = list(list(K = ZN1, model = "RKHS"))
      Eta2 = list(list(K = ZN2, model = "RKHS"))
      Eta3 = list(list(K = K1star, model = "RKHS"), list(K = K2star, model = "RKHS"), list(K = K3star, model = "RKHS"))
      Eta4 = list(list(K = K1star, model = "RKHS"), list(K = K2star, model = "RKHS"), list(K = K3star, model = "RKHS"), list(K = ZN1, model = "RKHS"))
      Eta5 = list(list(K = K1star, model = "RKHS"), list(K = K2star, model = "RKHS"), list(K = K3star, model = "RKHS"), list(K = ZN2, model = "RKHS"))
      
      Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5)
      
      hybrid = as.character(unique(Phenotype_data1$pedigree))
      
      set.seed(1234)
      cycles = 10
      CV2 = list()
      CV1 = list()
      
      set.seed(1234)
      for (MODEL in 1:length(Models)) {
        for (rep_num in 1:50) {
          test_geno <- sample(hybrid, 30)
          train_geno <- setdiff(hybrid, test_geno)
          
          CV_Data_1_2 <- Phenotype_data1
          CV_Data_1_2$Y <- NA
          CV_Data_1_2$Y[CV_Data_1_2$pedigree %in% train_geno] <- CV_Data_1_2[,1][CV_Data_1_2$pedigree %in% train_geno]
          
          y_t <- as.numeric(CV_Data_1_2$Y)
          fit <- BGLR(y = y_t, ETA = Models[[MODEL]], nIter = 5000, burnIn = 1000, thin = 10)
          
          CV_Data_1_2$yhat <- fit$yHat
          
          df_test <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% test_geno)
          CV1[[(rep_num)]] <- as.data.frame(df_test %>% dplyr::summarize(cor = cor(blue, yhat, use = "complete.obs")))
          
          df_train <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% train_geno)
          CV2[[(rep_num)]] <- as.data.frame(df_train %>% dplyr::summarize(cor = cor(blue, yhat, use = "complete.obs")))
        }
        
        if (rep_num == 50) {
          CV1out <- plyr::ldply(CV1, data.frame)
          CV2out <- plyr::ldply(CV2, data.frame)
          write.csv(CV1out, file = paste("ACC_", envnames[env], "_", traitnames[tr], "_CV1_model", MODEL, ".csv", sep = ""), row.names = F)
          write.csv(CV2out, file = paste("ACC_", envnames[env], "_", traitnames[tr], "_CV2_model", MODEL, ".csv", sep = ""), row.names = F)
        }
      }
    }, error = function(e) {
      cat("Error in processing environment:", envnames[env], "and trait:", traitnames[tr], "\n")
      cat("Error message:", e$message, "\n")
    })
  }
}


############# End of the code ######################################################################
#################################################NIR within env ####################################
############ Sapkota Pradip 2024 ###################################################################
##################v2################################################################################
######################Single Environment########################################
##################### Kernel Physical Factors ############################################

#loading libraries
library(BGLR)
library(dplyr)
library(stringr)
library(doBy)
library(caret)
library(tidyverse)
library(splitstackshape)
library(BGGE)
library(prospectr)
library(progress)

traitnames <- c("kernel_hardness" , "kernel_diameter", "kernel_weight")
envnames = c("18CS", "18VC", "19COL", "19VC")

pheno_geno_NIR = read.csv("blues_pheno_NIR_physical_factors.csv", check.names = FALSE)
names(pheno_geno_NIR) = tolower(names(pheno_geno_NIR))

#env =1; tr =1
for (env in 1:length(envnames)){
  for (tr in 1:length(traitnames)) {
    tryCatch({
      cat("Processing environment:", envnames[env], "and trait:", traitnames[tr], "\n")
      
      if (env == 1) {
        pheno <- subset(pheno_geno_NIR, env == "18CS")
      } else if (env == 2) {
        pheno <- subset(pheno_geno_NIR, env == "18VC")
      } else if (env == 3) {
        pheno <- subset(pheno_geno_NIR, env == "19COL")
      } else {
        pheno <- subset(pheno_geno_NIR, env == "19VC")
      }
      
      cat("Phenotype data dimensions for env:", envnames[env], "are", dim(pheno), "\n")
      
      if (tr == 1) {
        pheno_data <- read.csv("pheno_yield.csv")
      } else if (tr == 2) {
        pheno_data <- read.csv("pheno_pht.csv")
      } else {
        pheno_data <- read.csv("pheno_dy.csv")
      }
      
      cat("Trait data dimensions for trait:", traitnames[tr], "are", dim(pheno_data), "\n")
      
      Phenotype_data1 = pheno_data[pheno_data$env == envnames[env],]
      cat("Filtered phenotype data dimensions for env:", envnames[env], "are", dim(Phenotype_data1), "\n")
      
      rownames(pheno) = pheno$pedigree
      
      NIR.d1a <- scale(savitzkyGolay(pheno[,-c(1:22)], m = 1, p = 1, w = 11))
      ZN1 = tcrossprod(as.matrix(NIR.d1a) / ncol(as.matrix(NIR.d1a)))
      dim(ZN1)
      
      NIR.d2a <- scale(savitzkyGolay(pheno[,-c(1:22)], m = 1, p = 2, w = 11))
      ZN2 = tcrossprod(as.matrix(NIR.d2a) / ncol(as.matrix(NIR.d2a)))
      dim(ZN2)
      
      
      col = colnames(pheno)
      factors = col[1:7]
      numerics = col[8:length(col)]
      
      tmp <- strsplit(pheno$pedigree, "/")
      
      P1 <- rep(NA, length(tmp))
      P2 <- rep(NA, length(tmp))
      
      for (i in 1:length(tmp)) {
        P1[i] = tmp[[i]][1]
        P2[i] = tmp[[i]][2]
      }
      
      P1 <- as.factor(P1)
      P2 <- as.factor(P2)
      cross <- as.factor(pheno$pedigree)
      
      Z1 = model.matrix(~P1 - 1); dim(Z1)
      Z2 = model.matrix(~P2 - 1); dim(Z2)
      Z3 = model.matrix(~cross - 1); dim(Z3)
      
      load("X_G_Jales.RData")
      K1 = G[levels(P1), levels(P1)]; dim(K1)
      K2 = G[levels(P2), levels(P2)]; dim(K1)
      
      K3 = kronecker(K1, K2); dim(K3)
      
      K1star = Z1 %*% K1 %*% t(Z1); dim(K1star)
      K2star = Z2 %*% K2 %*% t(Z2); dim(K2star)
      K3star = Z3 %*% K3 %*% t(Z3); dim(K3star)
      
      pheno[,factors] = lapply(pheno[,factors], as.factor)
      pheno[,numerics] = lapply(pheno[,numerics], as.numeric)
      
      Eta1 = list(list(K = ZN1, model = "RKHS"))
      Eta2 = list(list(K = ZN2, model = "RKHS"))
      Eta3 = list(list(K = K1star, model = "RKHS"), list(K = K2star, model = "RKHS"), list(K = K3star, model = "RKHS"))
      Eta4 = list(list(K = K1star, model = "RKHS"), list(K = K2star, model = "RKHS"), list(K = K3star, model = "RKHS"), list(K = ZN1, model = "RKHS"))
      Eta5 = list(list(K = K1star, model = "RKHS"), list(K = K2star, model = "RKHS"), list(K = K3star, model = "RKHS"), list(K = ZN2, model = "RKHS"))
      
      Models <- list(Eta1, Eta2, Eta3, Eta4, Eta5)
      
      hybrid = as.character(unique(Phenotype_data1$pedigree))
      
      set.seed(1234)
      cycles = 10
      CV2 = list()
      CV1 = list()
      
      #MODEL =1; rep_num =1
      set.seed(1234)
      for (MODEL in 1:length(Models)) {
        for (rep_num in 1:50) {
          test_geno <- sample(hybrid, 30)
          train_geno <- setdiff(hybrid, test_geno)
          
          CV_Data_1_2 <- Phenotype_data1
          CV_Data_1_2$Y <- NA
          CV_Data_1_2$Y[CV_Data_1_2$pedigree %in% train_geno] <- CV_Data_1_2[,1][CV_Data_1_2$pedigree %in% train_geno]
          
          y_t <- as.numeric(CV_Data_1_2$Y)
          fit <- BGLR(y = y_t, ETA = Models[[MODEL]], nIter = 5000, burnIn = 1000, thin = 10)
          
          CV_Data_1_2$yhat <- fit$yHat
          
          df_test <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% test_geno)
          CV1[[(rep_num)]] <- as.data.frame(df_test %>% dplyr::summarize(cor = cor(blue, yhat, use = "complete.obs")))
          
          df_train <- subset(CV_Data_1_2, CV_Data_1_2$pedigree %in% train_geno)
          CV2[[(rep_num)]] <- as.data.frame(df_train %>% dplyr::summarize(cor = cor(blue, yhat, use = "complete.obs")))
        }
        
        if (rep_num == 50) {
          CV1out <- plyr::ldply(CV1, data.frame)
          CV2out <- plyr::ldply(CV2, data.frame)
          write.csv(CV1out, file = paste("ACC_", envnames[env], "_", traitnames[tr], "_CV1_model", MODEL, ".csv", sep = ""), row.names = F)
          write.csv(CV2out, file = paste("ACC_", envnames[env], "_", traitnames[tr], "_CV2_model", MODEL, ".csv", sep = ""), row.names = F)
        }
      }
    }, error = function(e) {
      cat("Error in processing environment:", envnames[env], "and trait:", traitnames[tr], "\n")
      cat("Error message:", e$message, "\n")
    })
  }
}

############### End of coding ##################################################################

######################## Plotting Prediction abilities ############# 
################within environments ################################
library(plyr)
library(readr)
library(dplyr)
library(stringr)

list_csv_files <- list.files("path_to_single_env_all/")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2
write.csv(df2, "Pred.ability.CV1.CV2.single.env.csv")

# Load necessary library
library(dplyr)

# Read the CSV file
df = read.csv("Pred.ability.CV1.CV2.single.env.csv")
df
df$file_name = gsub("kernel_hardness", "khi", df$file_name)
df$file_name = gsub("kernel_diameter", "kd", df$file_name)
df$file_name = gsub("kernel_weight", "kw", df$file_name)


df1 = df %>%
  mutate(split_file_name = str_split(file_name, "_", simplify = TRUE)) %>%
  mutate(location = split_file_name[, 2],
         traits = split_file_name[, 3],
         CV = split_file_name[, 4],
         model = str_remove(split_file_name[, 5], ".csv")) %>%
  select(-split_file_name)

# Subset data to include only model1, model3, and model4, and rename models
df <- df1 %>%
  filter(model %in% c("model1", "model3", "model4")) %>%
  mutate(
    model = case_when(
      model == "model1" ~ "P",
      model == "model3" ~ "G",
      model == "model4" ~ "G + P",
      TRUE ~ model
    ),
    traits = case_when(
      traits == "dy" ~ "DA",
      traits == "pht" ~ "PH",
      traits == "yield" ~ "GY",
      traits == "kd" ~ "KD",
      traits == "khi" ~ "KHI",
      traits == "kw" ~ "KW",
      TRUE ~ traits)
  )

df

df <- as.data.frame(df %>%  dplyr::group_by(model,CV,traits, location) %>% 
                      dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                       SD = sd(cor, na.rm=TRUE)))


df$CV <- factor(df$CV, levels =  c("CV2", "CV1"))
df$model = factor(df$model, levels = c("G", "P", "G + P"))
df$traits = factor(df$traits, levels = c("GY", "PH", "DA", "KHI", "KD", "KW"))
df$location = factor(df$location, levels = c("18COL", "18CS", "18GC", "18VC",
                                             "19COL", "19CS", "19TA", "19VC"))
df= df[df$CV =="CV1",]

head(df)
library(ggplot2)


# Create the plot
p =
ggplot(df, aes(x = model, y = M, color = M, ymin = M - SD, ymax = M + SD)) +
  geom_point(size = 4) +
  geom_errorbar(width = 0.2) +
  scale_color_gradient(low = "red", high = "blue", name = "Prediction Ability") +
  facet_grid(location ~ traits) +
  theme_bw()+
  xlab("Prediction Models") +
  ylab("Prediction Accuracy") +
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=15, face="bold", colour = "black"),    
    axis.title.y = element_text(size=15, face="bold", colour = "black"),    
    axis.text.x = element_text(size=11, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=11, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 11, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 11, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )


p = p + guides(fill = guide_legend(title = "Models"))

jpeg("Pred.single.env.jpeg",width = 12,height =9,units = "in", res=600)
p
dev.off()

head(df)


g = df[df$model == "G",]

write.csv(df, "pred.acc.single.csv")


ggplot(df, aes(models, pa, fill = models)) +
  geom_boxplot(colour = "blue")+
  labs(y = "Grain Yield", x= "")+
  scale_y_continuous(limits = c(0.2,1))+
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=20, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=25, face="bold", colour = "black"),    
    axis.title.y = element_text(size=25, face="bold", colour = "black"),    
    axis.text.x = element_text(size=22, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=22, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 15, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 15, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.4),
    axis.line.y = element_line(color="black", size = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, size=0.4),
    axis.text.x.bottom = element_blank()
  )









# Create the plot
p = ggplot(df, aes(model, y = M, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(M, 2)), hjust = -0.2, color = "black",
            position = position_dodge(0.9), angle = 90, size = 4) +
  geom_errorbar(aes(ymin = M, ymax = M + SD), width = 0.2,
                position = position_dodge(0.9)) +
  theme_bw() +
  facet_grid(location ~ traits) +
  scale_y_continuous("Prediction ability", limits = c(0, 1)) +
  xlab("Prediction Models") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 10),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 10, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "black"),
    strip.text.y = element_text(size = 10, face = "bold", colour = "black"),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  ) +
  guides(fill = guide_legend(title = "Models"))

# Save the plot as a JPEG file
jpeg("Pred.single.env.jpeg", width = 9, height = 6, units = "in", res = 600)
p
dev.off()

############## for kernel characterstics

df = read.csv("Pred.ability.CV1.CV2.CV0.CV00.csv")
head(df)

df = df[df$select == "1",]

df = df %>%
  filter(trait %in% c("hardness", "diameter", "weight")) %>%
  mutate(
    trait = case_when(
      trait == "hardness" ~ "KHI",
      trait == "diameter" ~ "KD",
      trait == "weight" ~ "KW",
      TRUE ~ trait
    ),
    model = case_when(
      model == "Genomic" ~ "G",
      model == "Phenomic" ~ "P",
      model == "Genomic + Phenomic" ~ "G + P",
      TRUE ~ model
    )
  )

df <- as.data.frame(df %>%  dplyr::group_by(model,CV,trait) %>% 
                      dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                       SD = sd(cor, na.rm=TRUE)))


df$CV <- factor(df$CV, levels =  c("CV2", "CV1", "CV0", "CV00"))
df$model = factor(df$model, levels = c("G", "P", "G + P"))
df$trait = factor(df$trait, levels = c("KHI", "KD", "KW"))


head(df)
library(ggplot2)

p = ggplot(df, aes(trait, y=M, fill=model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(M,2)  ), hjust=3, color="white",
            position = position_dodge(0.9), angle = 90,size=3.5)+
  geom_errorbar(aes(ymin=M, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  facet_grid(~CV)+
  scale_y_continuous("Prediction ability")+
  xlab("Kernel Physical Characterstics") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  
  
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=5, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=12, face="bold", colour = "black"),    
    axis.title.y = element_text(size=12, face="bold", colour = "black"),    
    axis.text.x = element_text(size=7, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=7, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 8, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 8, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    #axis.text.x.bottom = element_blank()
  )



p = p + guides(fill = guide_legend(title = "Models"))

jpeg("Pred.ability_kernel_traits.jpeg",width = 9,height =4,units = "in", res=600)
p
dev.off()

################ End of the coding ##################################################################### Plotting #######################
library(readr)
library(tidyverse)
list_csv_files = list.files(path = "path_to_CV1_CV2")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2

library(readr)
library(tidyverse)
list_csv_files = list.files(path = "path_to_CV0_CV00")
df2 <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
df2


############################# Variance Components for agronomic traits #################
########################################################################################
library(lme4)
library(lmerTest)
pheno = read.csv("pheno_NIR_data_8env.csv", check.names = FALSE)
colSums(is.na(pheno))
names(pheno) = tolower(names(pheno))
#pheno= pheno[,-c(22:24)]
#pheno = pheno[,-19]
#pheno = na.omit(pheno)


#### converting variables in factors

pheno$pedigree <- as.factor(pheno$pedigree)
pheno$env <- as.factor(pheno$env)
pheno$rep <- as.factor(pheno$rep)
pheno$female <- as.factor(pheno$female)
pheno$male <- as.factor(pheno$male)

# random effects model to extract variance components

model <- lmer(y1 ~ (1|female) + (1|male) + (1|pedigree) + (1|rep/env) +
                (1|env) + (1|female:env) + (1|male:env) + (1|female:male:env), 
              pheno)


ranova(model)
summary(model)
GCAmales <- ranef(model)$male
GCAfemales <- ranef(model)$female
SCA <- ranef(model)$"pedigree"
GXE = ranef(model)$"female:male:env"
GCAmales
GCAfemales
SCA

### calculation of variance heritability and CV component
variancehard <- print(VarCorr(model), com = c("Variance", "Std.Dev."))
variances = as.data.frame(variancehard)[,c(1,4,5)]
hybridvar = variances[6,2] + variances[7,2] + variances[2,2]
hybridenvvar = variances[1,2] + variances[3,2] +variances[4,2]
errorvar = variances[10,2]

addvar = variances[6,2] + variances[7,2]

cve = sqrt(errorvar)/mean(pheno$dy)

repeatability = hybridvar/(hybridvar+(hybridenvvar/8)+errorvar/17) #BSH
NSH = addvar/(hybridvar+(hybridenvvar/8)+errorvar/17) #NSH
repeatability
NSH


########################Summary statistics ####################################
pheno = read.csv("blues_pheno_NIR.csv", check.names = FALSE)



########################## visualizing varcomp ################################
##############################visualizing variance component##########
library(ggplot2)
library(plyr)


data = read.csv("variance_compo_8envagro4kernel.csv", check.names = FALSE)

data$value = as.numeric(data$value)


data$traits <- factor(data$traits, levels = c("Grain Yield", "Days to Anthesis",
                                              "Plant Height", "Kernel Hardness Index",
                                              "Kernel Diameter", "Kernel Weight"))

data$source <- factor(data$source, levels = c("male", "female", "pedigree",
                                              "env", "male:env", "female:env",
                                              "female:male:env", "rep",
                                              "Residual"))
data$source <- revalue(data$source, c("male" = "Male", "female" = "Female",
                                  "pedigree" = "Pedigree", "env" = "Environment",
                                  "male:env" = "Male X Environment", "female:env"=
                                    "Female X Environment", 
                                  "female:male:env" = "Female X Male X Environment",
                                  "rep" = "Replicaton", 
                                  "Residual" = "Residual"))




a <- ggplot(data, aes(fill=source, y=value, x=traits)) +
  geom_bar(position="stack", stat="identity")+
  labs(y= "Variance Components", x = "Traits")+
  # Add the point indicating the single value
  #geom_point(aes(x = category, y = single_value), color = "red", size = 3) +
  # Add labels to the points
  #geom_text(aes(label = single_value), vjust = -0.5)+
  guides(fill=guide_legend(title="Source of Variation")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  theme(
    # LABELS APPEARANCE
    plot.title = element_text(size=10, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=8, face="bold", colour = "black"),    
    axis.title.y = element_text(size=15, face="bold", colour = "black"),    
    axis.text.x = element_text(size=7, face="bold", colour = "black"), 
    # axis.text.y = element_text(size=22,  colour = "black"), # unbold
    axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 5, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 5, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.4),
    axis.line.y = element_line(color="black", size = 0.4),
    panel.border = element_rect(colour = "black", fill=NA, size=0.4),
    #axis.text.x.bottom = element_blank()
  )

a

jpeg("Varcomp.jpeg",width = 9,height =4,units = "in", res=600)
a
dev.off()

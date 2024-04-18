#!/usr/bin/env Rscript

## Julia Carrasco Zanini ---- 10/06/2021

## This script run feature selection of proteome determinants, 
## gets a performance metric in training and testing,
## and finally fits one joint model for selected proteome determinants to calculate the explained variance

## Set working directory
setwd("~/rds/rds-rjh234-mrc-epid/Studies/People/Julia/proteome_determinants/bin/")

options(stringsAsFactors = F)

rm(list=ls())

## load packages needed 
library(readstata13)
library(RNOmni)
library(caret)
library(glmnet)
library(doMC)
library(tidyr)
library(dplyr)
library(data.table)

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## define somamer
soma <- args[1]

## read in phenotype input file
pheno <- as.data.frame(fread(file="../data_input/Phenos_prots_scocres_input_file_proteome_determinats.txt", sep = "\t"))

## read protein labels
p.label <- read.delim("../../data/Fenland/phase1/SOMAscan_Assay_v4_Annotations_version3.3.2_tab.txt", sep = "\t", he=T)
p.label$SeqId <- as.character(p.label$SeqId)

## generate MRCseqid
p.label$MRC_seqid <- paste0("SeqId_", gsub("-","_", p.label$SeqId))

## subset to human and protein - 4979 proteins
p.label <- subset(p.label, Organism=="human"&Type=="Protein")

prots <- p.label$MRC_seqid

## read in phenotype labels
t.labels <- read.delim("../data_input/proteome_determinants_labels.txt", sep = "\t", he=T)

## split into training and testing
## set seed to get the same split
set.seed(3456)
trainIndex <- createDataPartition(pheno$age, p=0.75, list = F, times = 1)

## Training set
pheno.Train <- pheno[trainIndex,]
## Test set 
pheno.Test <- pheno[-trainIndex,]

## read in function to run feature selection 
source("functions/feature_rank_LASSO_gaussian_cis_trans_scores.R")

## parallelize
registerDoMC(32)
## run feature selection over 100 boots 
fs.rank.f <- lasso.fs.rank(dat = pheno.Train, feat = c(t.labels$pheno.labels,"noise.1","noise.2","noise.3"), outc = soma, n=100)
write.table(fs.rank.f[[soma]], file = paste0("../data_output/feature_selection/",soma,"_FS_determinants.txt"), sep = "\t", row.names = T)

####----    Filtering variables from feature selection and comupting LASSO R2 in test set   ----#### 

## subset labels into continuous or binary
t.labels.cont <- subset(t.labels, class=="continuous")
t.labels.cat <- subset(t.labels, class=="binary")

## define list of categorical and continuos variables
cont.vars <- lapply(c(soma), function(x){
  ## define the number of boots the highest ranked noise variable was selected in
  noise.rank <- grep("noise",row.names(fs.rank.f[[x]]))
  thresh <- max(fs.rank.f[[x]][noise.rank, "select"])
  ## if thresh == 100 ; the compute variance explained for other variables selected in 100 % of boots
  if(thresh ==100 ){
    tmp <- subset(fs.rank.f[[x]], select==100)
    tmp$vars <- row.names(tmp)
    tmp <- subset(tmp, !vars %in%t.labels.cat$pheno.labels& (!vars%in%c("noise.1","noise.2","noise.3")))
    return(tmp$vars)
  } else{
    tmp <- fs.rank.f[[x]][which(fs.rank.f[[x]][,"select"]>thresh),]
    tmp$vars <- row.names(tmp)
    tmp <- subset(tmp, !vars %in%t.labels.cat$pheno.labels& (!vars%in%c("noise.1","noise.2","noise.3")))
    return(tmp$vars) 
  }
})
names(cont.vars) <- c(soma)

## create alternate version for cis trans score names
cont.vars.2 <- cont.vars
cont.vars.2 <- lapply(c(soma), function(x){
  cont.vars.2[[x]] <- gsub(".SeqId_[[:digit:]]+_[[:digit:]]+","",cont.vars.2[[x]])
})
names(cont.vars.2) <- c(soma)

## binary
cat.vars <- lapply(c(soma), function(x){
  ## define the number of boots the highest ranked noise variable was selected in
  noise.rank <- grep("noise",row.names(fs.rank.f[[x]]))
  thresh <- max(fs.rank.f[[x]][noise.rank, "select"])
  ## if thresh == 100 ; the compute variance explained for other variables selected in 100 % of boots
  if(thresh ==100 ){
    tmp <- subset(fs.rank.f[[x]], select==100)
    tmp$vars <- row.names(tmp)
    tmp <- subset(tmp, vars %in%t.labels.cat$pheno.labels& (!vars%in%c("noise.1","noise.2","noise.3")))
    return(tmp$vars)
  } else{
    tmp <- fs.rank.f[[x]][which(fs.rank.f[[x]][,"select"]>thresh),]
    tmp$vars <- row.names(tmp)
    tmp <- subset(tmp, vars %in%t.labels.cat$pheno.labels& (!vars%in%c("noise.1","noise.2","noise.3")))
    return(tmp$vars) 
  }
})
names(cat.vars) <- c(soma)

## Training and Test R2
if(length(cont.vars[[soma]])>0&length(cont.vars[[soma]])>0){
  las.train <- train(pheno.Train[, c(cont.vars[[soma]],cont.vars[[soma]])], pheno.Train[ , soma],
                     family="gaussian",
                     method = "glmnet",
                     tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                     trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                     ))
  r2.train <- R2(pred = predict(las.train$finalModel, s= las.train$finalModel$lambdaOpt, newx = as.matrix(pheno.Train[,c(cont.vars[[soma]],cont.vars[[soma]])])),
                 obs = pheno.Train[,soma])
  r2.Test <- R2(pred = predict(las.train$finalModel, s= las.train$finalModel$lambdaOpt, newx = as.matrix(pheno.Test[,c(cont.vars[[soma]],cont.vars[[soma]])])),
                obs = pheno.Test[,soma])
} else{print("Next step")}

## second step
if(length(cont.vars[[soma]])>0&length(cont.vars[[soma]])==0){
  las.train <- train(pheno.Train[, c(cont.vars[[soma]])], pheno.Train[ , soma],
                     family="gaussian",
                     method = "glmnet",
                     tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                     trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                     ))
  r2.train <- R2(pred = predict(las.train$finalModel, s= las.train$finalModel$lambdaOpt, newx = as.matrix(pheno.Train[,c(cont.vars[[soma]])])),
                 obs = pheno.Train[,soma])
  r2.Test <- R2(pred = predict(las.train$finalModel, s= las.train$finalModel$lambdaOpt, newx = as.matrix(pheno.Test[,c(cont.vars[[soma]])])),
                obs = pheno.Test[,soma])
}else{print("Next step")}

## third step
if(length(cont.vars[[soma]])==0&length(cont.vars[[soma]])>0) {
  las.train <- train(pheno.Train[, c(cont.vars[[soma]])], pheno.Train[ , soma],
                     family="gaussian",
                     method = "glmnet",
                     tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                     trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                     ))
  r2.train <- R2(pred = predict(las.train$finalModel, s= las.train$finalModel$lambdaOpt, newx = as.matrix(pheno.Train[,c(cont.vars[[soma]])])),
                 obs = pheno.Train[,soma])
  r2.Test <- R2(pred = predict(las.train$finalModel, s= las.train$finalModel$lambdaOpt, newx = as.matrix(pheno.Test[,c(cont.vars[[soma]])])),
                obs = pheno.Test[,soma])
} else{print("Next step")}

lasso.perf <- data.frame(prot=soma,Train=r2.train, Test=r2.Test)
colnames(lasso.perf) <- c("prot","R2.Train","R2.Test")
## write train test R2
write.table(lasso.perf, file = paste0("../data_output/train_test_R2/",soma,"_performance.txt"), sep = "\t", row.names = F)

## turn categorical variables into factors
for (i in t.labels.cat$pheno.labels) {
  pheno[,i] <- as.factor(pheno[,i])
}

## create t.labels for cis and trans scores
t.labels[62,] <- c("cis.score","","0","grs","continuous") 
t.labels[63,] <- c("trans.score","","0","grs","continuous") 

####----    Variance Explained in all using selected features   ----#### 
##call variance partition function
source("functions/variance_partition.R")
source("functions/variance_partition_cont.R")

library(variancePartition)
require(doSNOW)

cl <- makeCluster(1)
registerDoSNOW(cl)

## run variance partition one joint model 
var.res <- lapply(c(soma), function(x){
  ## if there are binary variables run this version
  if(length(cat.vars[[x]])>0){
    tmp <- expl.var.part(pheno,x,cont.vars[[x]],cat.vars[[x]])
  } else{
    tmp <- expl.var.part.cont(pheno,x,cont.vars[[x]])
  }
  ## remove seqid from cis and trans score names
  colnames(tmp) <- gsub(".SeqId_[[:digit:]]+_[[:digit:]]+","",colnames(tmp))
  ## generate empty data frame for variables not included in the mode 
  ex.vars <- t.labels$pheno.labels[which(!t.labels$pheno.labels%in%cont.vars.2[[x]]&!t.labels$pheno.labels%in%cat.vars[[x]])]
  vec <- rep(NA,length(ex.vars))
  names(vec) <- ex.vars
  vec <- t(as.data.frame(vec))
  row.names(vec) <- x
  ## merge with results 
  tmp <- cbind(tmp, vec)
  tmp <- as.data.frame(tmp[,c("n","Residuals",t.labels$pheno.labels)])
  colnames(tmp) <- x
  return(tmp)
})

var.res.all <- data.frame(row.names(var.res[[1]]),var.res[[1]])
colnames(var.res.all) <- c("var",soma)
##save results
write.table(var.res.all, file = paste0("../data_output/variance_explained/",soma,"_VarianceExplained.txt"), sep = "\t", row.names = F)
##############################################
#### Julia Carrasco-Zanini    07/05/2021  ####
##############################################

## This function runs feature selection by LASSO regression by SUBSAMPLING (but can be modified to run by bootstrapping)
## The data it takes as input is a data drame containing all variables needed, assumes cleaning, transformations and formatting has been done previously,
## Also assumes splitting of the training sets has been done (i.e. the input should be the training set only)
## Functions enables parallelization, but number of cores to be used needs to be registered prior running the function


lasso.fs.rank <- function(dat, feat, outc, n){
  ## dat - dataframe with all variables
  ## feat - vector of features used as predictors
  ## outc - vector of variables used as outcome
  ## n - # of iterations to be done for the lasso loop 
  
  library(caret)
  library(glmnet)
  
  ## draw subsamples of the data 
  jj <- lapply(1:n, function(x) sample(nrow(dat), round(nrow(dat)*0.80), replace = F))
  
  ## LASSO loop over the subsamples
  cf <- lapply(outc, function(x){
    ## lasso loop 
    res <- lapply(1:n, function(i){
      print(i)
      ## run lasso 
      if(paste0("cis.score.",x)%in%colnames(dat)&paste0("trans.score.",x)%in%colnames(dat)){
        las.morb <- tryCatch({train(dat[jj[[i]], c(feat, paste0("cis.score.",x),paste0("trans.score.",x))], dat[jj[[i]], x],
                                    family="gaussian",
                                    method = "glmnet",
                                    tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                                    trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                                    ))
        })
      }
      if(paste0("cis.score.",x)%in%colnames(dat)&!paste0("trans.score.",x)%in%colnames(dat)){
        las.morb <- tryCatch({train(dat[jj[[i]], c(feat, paste0("cis.score.",x))], dat[jj[[i]], x],
                                    family="gaussian",
                                    method = "glmnet",
                                    tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                                    trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                                    ))
        })
      }
      if(!paste0("cis.score.",x)%in%colnames(dat)&paste0("trans.score.",x)%in%colnames(dat)){
        las.morb <- tryCatch({train(dat[jj[[i]], c(feat, paste0("trans.score.",x))], dat[jj[[i]], x],
                                    family="gaussian",
                                    method = "glmnet",
                                    tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                                    trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                                    ))
        })
        }
      if(!paste0("cis.score.",x)%in%colnames(dat)&!paste0("trans.score.",x)%in%colnames(dat)){
        las.morb <- tryCatch({train(dat[jj[[i]], feat], dat[jj[[i]], x],
                                    family="gaussian",
                                    method = "glmnet",
                                    tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                                    trControl = trainControl(method="repeatedcv", number=10, repeats = 10, allowParallel=T
                                    ))
        })
      }      
      
      ## collect results for feature selection in the final model
      tmp <- as.matrix(coef(las.morb$finalModel, s = las.morb$finalModel$lambdaOpt))
      return(tmp)
    })
    
    ## generate a feature selection 0/1 matrix
    fs <- lapply(1:n, function(x){
      tmp <- cbind(ifelse(res[[x]][,1]!=0,1,0))
      return(tmp)
    })
    ## make intp data.frame
    fs <- do.call(cbind,fs)
    
    ## order by number of times a feature was selected and remove the intercept
    fs <- fs[-1,]
    fs <- as.data.frame(fs)
    fs$select <- rowSums(fs[,1:n])
    fs <- fs[order(fs$select, decreasing = T),]   
    
    return(fs)
  })
  names(cf) <- outc
  return(cf)
}
##########################################
## function to compute exlpained variance
## by different factors at once

expl.var.part <- function(dat, out, c.cont, c.cat){
  
  ## 'dat'    -- data frame containing all needed data
  ## 'out'    -- vector of outcome variables to be tested
  ## 'c.cont' -- continous traits model as fixed effects
  ## 'c.cat'  -- categorical traits model as random effect
  
  ## load package to do so
  library(variancePartition)
  
  ## code categorical variables
  # dat[,c.cat] <- apply(dat[,c.cat,drop=F], 2, as.factor)
  
  ## create formula object
  form <- paste0("~ ", paste(c.cont, collapse = " + "), " + ", paste(paste0("(1|", c.cat, ")"), collapse = " + "))
  cat("using the following formula for variance partitioning:\n", form,"\n")
  
  ## prepare storage
  res           <- array(data=NA, dim=c(length(out), length(c(c.cat, c.cont))+2))
  rownames(res) <- out
  
  print(dim(na.omit(dat[, c(out[1], c.cat, c.cont)])))

  ## run simple loop over all variables
  for(j in 1:length(out)){
    if(j%%10 == 0){
      cat(".")
    }
    if(j%%100==0){
      cat(j,"\n")
    }

    ## run within tryCatch to avoid douplings
    tmp <- tryCatch({
      ## run the analyses; cave: omit NA first!!
      tmp       <- na.omit(dat[, c(out[j], c.cat, c.cont)])
      ## store number of observations used
      n         <- nrow(tmp)
      ## cave transpose the outcome matrix
      tmp       <- fitExtractVarPartModel(t(tmp[,1,drop=F]), form, tmp[, -1], showWarnings = F)
      tmp       <- c(n,unlist(tmp[1,]))
    }, error=function(e){
      cat("no compuation on ", out[j], "\n")
      return(rep(NA, length(c(c.cat, c.cont))+2))
    }
    )

    ## store the results
    res[j,] <- tmp
  }

  ## store the names
  colnames(res) <- c("n", names(tmp)[-1])
  
  cat("\n")
  
  return(res)  
}
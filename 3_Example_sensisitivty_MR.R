## This script runs sensitivity MR analyses for BMD to proteins
## (only proteins for which bone mass was the major determinant)

rm(list = ls())

setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Julia/proteome_determinants/bin/")

options(stringsAsFactors = F)

library(RadialMR)
library(MendelianRandomization)
library(TwoSampleMR)

## read in determinants file 
params <- read.delim("mrc_seqid_param_file.txt", sep = "\ ", he=F) 
head(params)
var.expl <- read.delim("../data_output/variance_explained_all_proteins.txt", sep ="\t", header = T)

var.res.all <- subset(var.expl, var!="n"&var!="Residuals"&var!="TestSite.2"&var!="TestSite.3"&var!="AppDate_Attended_ROUNDED"&var!="AppDate.spring"&
                        var!="AppDate.summer"&var!="AppDate.winter"&var!="PC1.20"&var!="PC1.5"&var!="PC1.005")

u.labels.s <- sapply(params$V1, function(x){
  ii <- which.max(var.res.all[,x])
  deter <- var.res.all$var[ii]
  return(deter)
})

mod.deter <- as.data.frame(cbind(u.labels.s))
mod.deter$protein <- row.names(mod.deter)

## filter bone mass as a determinants
# length(grep("AD04i_iDEXA_Total_bone_mass", mod.deter$u.labels.s))
exp.prots <- c(mod.deter$protein[grep("AD04i_iDEXA_Total_bone_mass", mod.deter$u.labels.s)])
gwas.exp.p <- gsub("SeqId_","res_invn_X", exp.prots)

## BMD instruments
library(MRInstruments)
data(gwas_catalog)

exp.stats <- subset(gwas_catalog, STUDY.ACCESSION=="GCST001482")
exp.stats <- subset(exp.stats,pval<5e-8)
exp.stats <- format_data(exp.stats)
exp.stats <- clump_data(exp.stats, clump_r2 = 0.01)

mr.res <- read.delim("../data_output/MR_outliers_summary_nominal_results.txt", sep = "\t")
mr.res <- subset(mr.res,u.labels.s=="AD04i_iDEXA_Total_bone_mass")
mr.res$gwas.prot <- gsub("SeqId_","res_invn_X", mr.res$prot) 
gwas.exp.p <- gwas.exp.p[which(gwas.exp.p%in%mr.res$gwas.prot)]

## filter in the instruments for the required proteins
radial.res <- lapply(gwas.exp.p, function(x){
  ## read in some info
  res.soma <- read_outcome_data(
    snps = exp.stats$SNP,
    filename = paste0("zcat ~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Meta-Analysis/Fenland_auto_chrX_filtered/output/",
                      x,"_Fenland_MA_auto_chrX_filtered.txt.gz"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "Pvalue",
    samplesize_col = " TotalSampleSize")
  ## subset to variants founs in soma gwas
  exp.stats <- subset(exp.stats, SNP%in%res.soma$SNP)
  
  ## harmonize
  tsmrdat <- TwoSampleMR::harmonise_data(exposure_dat = exp.stats,
                                         outcome_dat = res.soma )
  ##run radial MR
  res.mr <- ivw_radial(r_input = tsmrdat, alpha = 0.05,
                       weights = 1, tol = 0.0001)
  ## sensitivity analyses
  if(res.mr$outliers!="No significant outliers"){
    out.dat <- subset(tsmrdat,!SNP%in%res.mr$outliers )
  } else{
    out.dat <- tsmrdat}
  
  ## egger, normal IVW, median
  sens.mr <- mr(out.dat)
  
  egger.int <- mr_pleiotropy_test(out.dat)
  egger.int$nsnp <- NA
  egger.int$method <- "egger_intercept"
  
  egger.int <- egger.int[,c("id.exposure","id.outcome","outcome","exposure","method","nsnp","egger_intercept","se","pval")]
  colnames(egger.int)[7] <- "b"
  
  sens.mr <- rbind(sens.mr, egger.int)
  sens.mr$id.exposure <- x
  sens.mr$id.outcome <- "Bone Mineral Density"
  return(sens.mr)
})
radial.res <- do.call(rbind, radial.res)
## save results for the sensitivity analyses
write.table(radial.res, file="../data_output/outlier_reomval/sensitivity_MR_BMD_res.txt", sep = "\t", row.names = F)

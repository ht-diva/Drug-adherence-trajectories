#### create covariates file ####
rm(list=ls())
gc()

setwd("/media/volume/mferro/scripts")
.libPaths("/shared-directory/sd-tools/apps/R/lib/")
library(feather)
library(data.table)
library(dplyr)

char = fread("/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_SUMMARY.csv.gz",
             select = c("PATIENT_ID","AGE"))
#read clusters
groups_MPR = fread("/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_CLUSTERS.csv") %>% 
  dplyr::select(FINREGISTRYID = PATIENT_ID, T_GROUPS = DESC) %>% as.data.table()

#read minimal
mp = fread("/home/mattferr/Projects/SD-Connect/project_2007099/processed_data/FinRegistry_v01/minimal_phenotype/minimal_phenotype_2023-08-14.csv.gz",
           select = c("FINREGISTRYID","DATE_OF_BIRTH","SEX"))[
             FINREGISTRYID %in% groups_MPR$FINREGISTRYID]

cov = merge(groups_MPR,mp,all.x = T) %>%
  merge(char,all.x = T, by.x = "FINREGISTRYID", by.y = "PATIENT_ID") 

cov = cov[T_GROUPS == "High" | T_GROUPS == "Dec"]

invisible(gc(full = T))


# read matrixes
mat = fread("/media/volume/mferro/data/ADHERENCE/socioecomic_covariates.csv")
cov = merge(cov,mat, by = "FINREGISTRYID")

#read ICC
cci = fread("/media/volume/mferro/data/ADHERENCE/STATINS/CCI.csv")
colnames(cci) = c("FINREGISTRYID","CCI_score")
cov = merge(cov,cci, by = "FINREGISTRYID")

# read endpoints
end = fread("/media/volume/mferro/data/ADHERENCE/diseases_covariates.csv")
cov = merge(cov,end, by = "FINREGISTRYID")

# read meds
med = fread("/media/volume/mferro/data/ADHERENCE/medication_covariates.csv") 
cov = merge(cov,med, by = "FINREGISTRYID")

uniqueN(cov)
cov = cov[AGE > 18]
uniqueN(cov)

cov[,YOB := year(DATE_OF_BIRTH)]
cov[,DATE_OF_BIRTH := NULL]

cov = cov %>%
  mutate_if(is.logical, as.numeric)

fwrite(cov,"/media/volume/mferro/data/ADHERENCE/ALL_COVARIATES_MPR.tsv",
       sep = "\t", col.names = T, row.names = F, quote = F)

rm(list=ls())
gc(full = T)
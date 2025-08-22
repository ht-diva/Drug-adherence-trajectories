setwd("/media/volume/mferro/scripts")
rm(list=ls())
gc()


install.packages("/media/volume/mferro/scripts/CCI/ICCI/ICCI_2.3.0.tar.gz",
                 "/shared-directory/sd-tools/apps/R/lib/",
                 repos = NULL, type="source")


.libPaths("/shared-directory/sd-tools/apps/R/lib/")
library(feather)
library(ICCI)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(lubridate)
library(comorbidity)
setDTthreads(0)

long <- fread("/media/volume/mferro/data/CCI/longitudinal_adherence_statins_all.csv.gz")
colnames(long) = c("FINREGISTRYID", "SOURCE", "EVENT_AGE", "EVENT_DATE", "CODE1", "ICDVER")

# Adherence
adh <-  fread("/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_SUMMARY.csv.gz", select = c("PATIENT_ID", "first"))


long <- merge(long, adh, by.x = "FINREGISTRYID" , by.y = "PATIENT_ID")[,ICDVER := fifelse(ICDVER == "O3", "10", ICDVER)] %>% 
  select(ID = FINREGISTRYID, PRIMARY_ICD = CODE1, ICD_VERSION = ICDVER, EVENT_AGE = EVENT_AGE, startDATE = first) 

setDT(long)
long = long[!is.na(EVENT_AGE) & EVENT_AGE < 150 & ICD_VERSION %in% c("9","10")]
  
  
score_data <- ICCI::calc_cci(long, exp_end = long$startDATE)
  
# Some indivuduals won't have records in the longitudinal files, hence no comorbidities and CCI == 0
d <- merge(adh, score_data, by.x = "PATIENT_ID" , by.y = "ID", all.x = T)[,CCI_score := fifelse(is.na(CCI), 0, CCI)] %>% 
  select(PATIENT_ID, CCI_score)
  
fwrite(d,"/media/volume/mferro/data/ADHERENCE/STATINS/CCI.feather")
rm(list=ls())
gc()
  
  
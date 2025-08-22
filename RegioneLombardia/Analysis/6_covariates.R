#### setting ####
rm(list = ls())
invisible(gc(full = T))

library(dplyr)
library(lubridate)
library(data.table)
library(haven)
library(ggplot2)
library(naniar)

setwd("N:/output/scripts/acorbetta/")

#### CREATE COVARIATES FILE ####

# outcome
groups = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_CLUSTERS.csv",
               select = c("PATIENT_ID","DESC"),
               col.names = c("COD_SOGGETTO","T_GROUPS"))

groups = groups[T_GROUPS == "High" | T_GROUPS == "Dec"]
  
scores = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_FPCSCORES.csv")
colnames(scores)[1] = "COD_SOGGETTO"

# get medications
meds = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/medication_covariates_all.csv")
meds$C10 = NULL
# get diseases 
dis = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/diseases_covariates_all.csv")

# get visits
vis = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/visits_covariates_all.csv")
    
# get socioeco
mat = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/socioeco_covariates_all.csv")


#merge data
data = merge(scores, groups, by = "COD_SOGGETTO", all.x = T) %>%
  merge(meds, by = "COD_SOGGETTO", all.x = T) %>%
  merge(dis, by = "COD_SOGGETTO", all.x = T) %>%
  merge(mat, by = "COD_SOGGETTO", all.x = T) %>%
  merge(vis, by = "COD_SOGGETTO", all.x = T)

data[, ATS := fifelse(ATS == "",NA_character_,ATS)]
data[,CONTINENT := relevel(as.factor(CONTINENT), ref = "European")]
data[,ATS := relevel(as.factor(ATS),ref = "ATS DELLA CITTA' METROPOLITANA DI MILANO")]
data[,SEX := fifelse(SEX == "M",1,0)]

colnames(data)
colnames(data) = make.names(colnames(data))

library(mice)
library(future.apply)
future::plan("multisession")
cov_imp = futuremice(data[,-c("COD_SOGGETTO", "T_GROUPS", "PC1",
                              "PC2", "PC3", "PC4", "PC5", "PC6")], method = "pmm", n.core = 5, 
                     future.plan = "multisession", m = 5,
                     parallelseed  = 1)


data2 = cbind(data[,c("COD_SOGGETTO", "T_GROUPS", "PC1",
                      "PC2", "PC3", "PC4", "PC5", "PC6")], complete(cov_imp))


fwrite(data2,"N:/output/data/acorbetta/STATINS/ADHERENCE/ALL_COVARIATES_MPR.tsv",
       sep = "\t", col.names = T, row.names = F, quote = F)

rm(list=ls())
gc(full = T)


#### SET UP ####
rm(list=ls())
gc(full = T)


setwd("/media/volume/mferro/scripts/")

#.libPaths("/shared-directory/sd-tools/apps/R/lib/")
library(data.table)
library(dplyr)
library(ggplot2)
library(mice)


groups_MPR = fread("/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_CLUSTERS.csv") %>% 
  dplyr::select(FINREGISTRYID = PATIENT_ID, T_GROUPS = DESC) %>% as.data.table()
groups_MPR = groups_MPR[T_GROUPS == "High" | T_GROUPS == "Dec"]


#### read data ####
cov = fread("/media/volume/mferro/data/ADHERENCE/covariates_adherence_statins.csv-matrix.csv")[FINREGISTRYID %in% groups_MPR$FINREGISTRYID] %>%
  dplyr::select(-SEX,-date_of_birth, - start_of_followup, - end_of_followup) 


fwrite(list(colnames(cov)[-1]),"/media/volume/mferro/summary_stats/covariates_socioeco.csv")


### PREPROCESS DATA ####

#correct mr variables
#mr = (colnames(cov)[grepl("mr_", colnames(cov))])

cov = cov %>%
  mutate_at(vars(starts_with("mr_")), function(x){ifelse(x == "0" | x == "" | x=="False", FALSE, TRUE)}) %>%
  #keep over 15
  # keep only big zones
  mutate(zip_code = substr(as.character(zip_code),1,1),
         children = if_else(children > 0 ,1,0)) %>%
  #factorize 
  mutate_if(~length(unique(.)) < 4, factor)

#change var names
cov = cov %>%
  rename(residential_care_housing_under_65yo = `247_residential_care_housing_under_65yo`,
         psychiatric_residential_care = `247_psychiatric_residential_care`) %>%
  as.data.frame()

# get single-value variables
uni = sapply(cov, function(x) nlevels(x))
singleval = names(uni[uni ==1])
#remove those
cov = cov[, !names(cov) %in% singleval]

#re-code variables
cov$urban_rural_class_code = as.factor(ifelse(cov$urban_rural_class_code == "", NA, as.numeric(substr(cov$urban_rural_class_code,2,3))))
cov$zip_code = as.factor(ifelse(cov$zip_code == "", NA, (cov$zip_code)))
cov$mothertongue_other = if_else(cov$mothertongue_other == T, cov$mothertongue_other, cov$mothertongue_rus)
gc(full = T)

#### missing patterns ####
library(ggplot2)
library(UpSetR)
library(naniar)

svg("/media/volume/mferro/plots/socioeco_miss.svg", width = 10, height = 10)
gg_miss_upset(cov)
dev.off()

#remove obs with more than 80% miss
cov = cov[which(rowMeans(!is.na(cov)) > 0.8), which(colMeans(!is.na(cov)) > 0.8)]

# get numeric vars
num = cov %>%
  select_if(is.numeric)
scaled_num = scale(num)

r2 = cor(num, use  = "pairwise.complete.obs")
fwrite(r2,"/media/volume/mferro/summary_stats/socioeco_numeric_cormat.csv")


library(ggcorrplot)

svg("/media/volume/mferro/plots/socioeco_corrplot.svg", width = 10, height = 10)
ggcorrplot(r2, type = "lower", lab = F )
dev.off()

# Find pairs with high correlation
high_cor <- which(abs(r2) > 0.7 & r2 != 1, arr.ind = TRUE)

# Extract variable names for high correlations
high_cor_pairs <- data.frame(
  Var1 = rownames(r2)[high_cor[, 1]],
  Var2 = colnames(r2)[high_cor[, 2]],
  Correlation = r2[high_cor]
)

# Remove duplicates (since correlation matrix is symmetric)
high_cor_pairs <- high_cor_pairs[high_cor_pairs$Var1 < high_cor_pairs$Var2, ]

# Print the results
fwrite(high_cor_pairs,"/media/volume/mferro/summary_stats/socioeco_high_cor_pairs.csv")



scaled_num = scaled_num %>%
  data.table()%>%
  dplyr::select(-average_income_of_inhabitants, -total_labor_income, 
                -economic_dependency_ratio, general_at_risk_of_poverty_rate_for_the_municipality,
                - miscarriages, -terminated_pregnancies)

fac = cov %>%
  select_if(is.factor)

gmcc = function(factor_variable) {
  table_factor = table(factor_variable)
  most_common = names(table_factor)[which.max(table_factor)]
  rel_freq = table_factor[most_common] /sum(table_factor)
  return(rel_freq)
}


rel_freq_a = sapply(fac, gmcc)

vk = names(rel_freq_a[rel_freq_a < 0.99])


vk = sapply(vk, function(x) sub("\\..*$","",x))
vk = c(unique(vk) , "mothertongue_other")
fac = fac[,names(fac) %in% vk]

rm_fac = names(fac)[!names(fac) %in% vk]
fwrite(list(rm_fac),"/media/volume/mferro/summary_stats/socioeco_removed_factors.csv")

#### continue ####
cov = cbind(scaled_num,fac) %>% as.data.frame()
rm(fac)
rm(num)
rm(scaled_num)

gc(full = T)
library(mice)

svg("/media/volume/mferro/plots/socioeco_miss.svg", width = 10, height = 10)
gg_miss_upset(cov, nsets =5)
dev.off()

#### impute missing ####

library(future.apply)
future::plan("multisession")
cov_imp = futuremice(cov, method = "pmm", n.core = 5, future.plan = "multisession", m = 5)


cov_imputed = cbind(FINREGISTRYID = groups_MPR$FINREGISTRYID,complete(cov_imp)) #%>%
#  rename(FINREGISTRYID = `cov$FINREGISTRYID`)


fwrite(cov_imputed,"/media/volume/mferro/data/ADHERENCE/socioecomic_covariates.csv")
cov_imputed = fread("/media/volume/mferro/data/ADHERENCE/socioecomic_covariates.csv")
fwrite(list(colnames(cov_imputed)[-1]),"/media/volume/mferro/summary_stats/covariates_socioeco_final_selection.csv")


### CLEAN UP #
rm(list=ls())
gc(full = T)
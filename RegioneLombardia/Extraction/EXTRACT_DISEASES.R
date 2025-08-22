#### setting ####
rm(list = ls())
invisible(gc(full = T))

library(dplyr)
library(lubridate)
library(data.table)
library(parallel)
library(feather)


setwd("N:/output/scripts/acorbetta/")

parse_sas_date = function(data) {
  return(as.Date(strptime(data, "%d%b%Y:%H:%M:%S.%OS")))
}


#### set up parameters ####
#filename wants a csv file without header with id and date of end of evaluation window (END)

#change to csv
filename = "N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv"

#output file
fileout = "N:/output/data/acorbetta/STATINS/ADHERENCE/diseases_data_all.csv"

print("reading anagrafe")
tmp <- fread("N:/output/data/minimum_data.csv",header = F, nrows = 1) %>% as.character()
pt =  tmp[grepl("COD_SOGGETTO|^FIRST", tmp)]
anagrafe <- fread("N:/output/data/minimum_data.csv", select = pt)
anagrafe$FIRST_COVID_DIAG = NULL
anagrafe$FIRST_HYPERC_DIAG = NULL
anagrafe$FIRST_HYPERT_DIAG = NULL

sum_stat = fread(filename)

#remove this
colnames(sum_stat) = c("COD_SOGGETTO", "END")

gc(full=T)

anagrafe = anagrafe[COD_SOGGETTO %in% sum_stat$COD_SOGGETTO]
anagrafe = merge(anagrafe,sum_stat,  by = "COD_SOGGETTO", all.x = T)

cols_to_check <- setdiff(names(anagrafe), c("COD_SOGGETTO", "END"))


for (col in cols_to_check) {
  new_col_name <- sub(".*_(.*)_.*", "\\1", col)
  anagrafe[, (new_col_name) := fifelse(get(col) < END,1 ,0 , na =0)]
}

library(dplyr)
anagrafe = anagrafe %>%
  dplyr::select(-starts_with("FIRST"), - END)


fwrite(anagrafe, fileout)

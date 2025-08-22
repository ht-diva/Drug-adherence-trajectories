
rm(list = ls())
invisible(gc(full = T))

library(data.table)
library(dplyr)
setDTthreads(0)
library(haven)
library(feather)
library(parallel)
library(tidyr)

setwd("N:/output/scripts/acorbetta/")

parse_sas_date = function(data) {
  return(as.Date(strptime(data, "%d%b%Y:%H:%M:%S.%OS")))
}

# set how long the evaluation should be in days
time_of_eval = 3*365.25

# read file
sum_stat = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv")
colnames(sum_stat) = c("COD_SOGGETTO", "first")

#set start and stop of evaluation
sum_stat[,data_iniziale := first - time_of_eval]
sum_stat[,data_finale := first ]
sum_stat$first = NULL



#### read map ####
map = read_sas('N:/input/DATI/htcovid_codifiche.sas7bdat') 
setDT(map)
disc = map[TIPO == "DISCIPLINE", .(CODICE, DESCRIZIONE)]
colnames(disc) = c("DISCIPLINA", "DESCRIZIONE_DISCIPLINA")



#### read and filter ####
# set list of files
string = paste0("AMBULATORIALE_",(year(min(sum_stat$data_iniziale))):year(max(sum_stat$data_finale)), collapse = "|")
amb_list =  list.files("N:/input/DATI/", full.names = T) 
amb_list = amb_list[grep(string,amb_list)]

load_and_filter <- function(file_path, sum_stat) {
  dt <- fread(file_path, select = c("COD_SOGGETTO","DISCIPLINA","QUANTITA","DT_EROGAZIONE"))
  dt = merge(dt,sum_stat, by = "COD_SOGGETTO")
  dt[, DT_EROGAZIONE := parse_sas_date(DT_EROGAZIONE)]
  dt = dt[DT_EROGAZIONE >= data_iniziale]
  dt = dt[DT_EROGAZIONE <= data_finale]
  gc(full=T)
  return((dt[, .(COD_SOGGETTO, DISCIPLINA,QUANTITA)]))
}

#  cluster con 2 core ( noth enough memory for more)
cluster <- makeCluster(7)

# Esportare le variabili e le funzioni necessarie nel cluster
clusterExport(cluster, c("load_and_filter", 'sum_stat', "amb_list","parse_sas_date"))
clusterEvalQ(cluster, {
  library(data.table)
})

# Parallelizzazione del caricamento e filtraggio dei file
# go have some coffe, it's gonna take a while
results <- parLapply(cluster, 1:length(amb_list), function(i) {
  file_path <- amb_list[[i]]
  load_and_filter(file_path, sum_stat)
})

stopCluster(cluster)
rm(cluster)
gc(full=T)


amb = rbindlist(results, use.names = TRUE, fill = TRUE)

rm(results)

amb = merge(amb,disc, by = "DISCIPLINA", all.x = T)

amb[,DISCIPLINA := NULL]

amb2 = amb %>%
  group_by(COD_SOGGETTO) %>%
  pivot_wider(names_from = DESCRIZIONE_DISCIPLINA, values_from = QUANTITA, 
              values_fn = sum, values_fill = 0)

setDT(amb2)
amb2[,`DATO MANCANTE` := NULL]


amb3 =  merge(sum_stat[,1], amb2, by = "COD_SOGGETTO", all.x = T)
amb3[is.na(amb3)] <- 0

fwrite(amb3, "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_ambulatoriale_data_all.csv.gz")

rm(list=ls())
gc(full=T)

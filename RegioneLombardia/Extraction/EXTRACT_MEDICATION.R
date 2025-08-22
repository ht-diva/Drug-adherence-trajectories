#### SET UP ####
rm(list=ls())
invisible(gc(full = T))


setwd("N:/output/scripts/acorbetta/")
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(feather)
library(purrr)
library(parallel)
library(doParallel)
library(tidyr)

parse_sas_date = function(data) {
  return(as.Date(strptime(data, "%d%b%Y:%H:%M:%S.%OS")))
}
print("Setting parameters and starting process ...")

#### SET PARAMS
filename = "N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv"

file_path = "N:/input/DATI"
file_pattern = "^HT_COVID_FAR_TERR"
filenameout = "N:/output/data/acorbetta/STATINS/ADHERENCE/drugs_data.csv_all.gz"

time_of_eval = 365.25 *3


sink("drug_matrix_log.txt")
to = Sys.time()
print("Starting extraction ...")


sum_stat = fread(filename)
colnames(sum_stat) = c("COD_SOGGETTO", "END")
sum_stat[,DATA_START_EVAL := ymd(END) - time_of_eval]


load_and_filter <- function(file_path, lapse) {
  dt = fread(file_path, select = c("COD_SOGGETTO","DT_EROGAZIONE","ATC"))
  dt = merge(dt, lapse, by ="COD_SOGGETTO")
  dt[, DT_EROGAZIONE := parse_sas_date(DT_EROGAZIONE)]
  dt = dt[DT_EROGAZIONE <= END]
  dt = dt[DT_EROGAZIONE >= DATA_START_EVAL]
  dt = dt[ATC != "?" & ATC != "-" & ATC != "#" ]
  dt[, ATC := substr(ATC, 1, 3)]
  return(unique(dt[, .(COD_SOGGETTO, ATC)]))
}


ymax = year(max(sum_stat$END))
ymin = year(min(sum_stat$DATA_START_EVAL))

pattern = paste0("HTCOVID_FAR_TERR_",ymin:ymax,".csv",collapse = "|")

file_paths <- list.files("N:/input/DATI", full=T)
file_paths = file_paths[grepl(pattern, file_paths)]
file_paths = c(file_paths,  "N:/input/DATI/HTCOVID_FILEF.csv")


# Creazione di un cluster con 7 core
cluster <- makeCluster(7)


# Esportare le variabili e le funzioni necessarie nel cluster
clusterExport(cluster, c("load_and_filter", "sum_stat", "file_paths","parse_sas_date"))
clusterEvalQ(cluster, {
  library(data.table)
  library(lubridate)
  library(dplyr)
})

# Parallelizzazione del caricamento e filtraggio dei file
results <- parLapply(cluster, 1:length(file_paths), function(i) {
  file_path <- file_paths[[i]]
  load_and_filter(file_path, sum_stat)
})

# Unire i risultati
farmaceutica <- rbindlist(results, use.names = TRUE, fill = TRUE)[
  , unique(.SD)]

# Arrestare il cluster
stopCluster(cluster)
rm(cluster)

# Pulizia della memoria
rm(results)
gc(full = TRUE)


farmaceutica[,D:=1]
farmaceutica = farmaceutica %>% 
  pivot_wider(id_cols = COD_SOGGETTO, names_from = ATC, values_from = D,
              values_fill = list(D = 0))  



farma = merge(sum_stat, farmaceutica, by = "COD_SOGGETTO", all.x = T)
farma = farma[,-c("END","DATA_START_EVAL")]
farma[is.na(farma)] <- 0


# write out
print(paste0("Writing file: ",filenameout))
fwrite(farma, filenameout, quote = F)


print(paste("Time taken: ", Sys.time() - to))
sink()

rm(list = ls())
invisible(gc(full = T))


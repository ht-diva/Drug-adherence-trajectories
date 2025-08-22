#### setting ####
start_time = Sys.time()
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
fileout = "N:/output/data/acorbetta/STATINS/ADHERENCE/diseases_data_all.csv.gz"

# how much time from END do you want to evaluate MCS retrospectively
# NB: if they don't have enought time they will be discarded
time_of_eval = 365.25 *3

#change to csv
sum_stat = fread(filename)

#remove this
colnames(sum_stat) = c("COD_SOGGETTO", "END")

gc(full=T)



###################### ANAGRAFE ###################### 

print("reading anagrafe")
tmp <- fread("N:/output/data/minimum_data.csv",header = F, nrows = 1) %>% as.character()
pt =  tmp[grepl("COD_SOGGETTO|SESSO|ANNO_NASCITA|DEATH_DATE|^RL_", tmp)]
anagrafe <- fread("N:/output/data/minimum_data.csv", select = pt)


gc(full = T)


#get start and end of assistance
base = fread("N:/input/DATI/HTCOVID_ANAGRAFICA.csv", 
             select = c("COD_SOGGETTO", "DATA_FINE_ASSISTENZA", "DATA_INIZIO_ASSISTENZA"))
gc(full = T)

setDT(anagrafe)
anagrafe = anagrafe[COD_SOGGETTO %in% sum_stat$COD_SOGGETTO]
anagrafe = merge(anagrafe,base,  by = "COD_SOGGETTO", all.x = T)
anagrafe = merge(anagrafe,sum_stat,  by = "COD_SOGGETTO", all.x = T)


rm(base)
rm(sum_stat)
gc(full = T)


#### set up ####
anagrafe[,DATA_INIZIALE := ymd(END)]

# date of evaluation start
anagrafe[,DATA_START_EVAL := ymd(END) - time_of_eval]

# Imposta la data di ingresso come Date
anagrafe[, dt_ingresso := as.Date(DATA_INIZIALE)]

# Imposta l'anno di nascita come Date
anagrafe[, dt_nascita := as.Date(paste0(ANNO_NASCITA, "-01-01"))]

# Imposta la data di inizio assistenza come Date
anagrafe[, dt_inizio_assistenza := as.Date(DATA_INIZIO_ASSISTENZA, format='%d/%m/%Y')]

# Imposta la data di fine assistenza come Date
anagrafe[, dt_fine_assistenza := as.Date(DATA_FINE_ASSISTENZA, format='%d/%m/%Y')]

# Imposta la data di decesso come Date
anagrafe[, dt_decesso := as.Date(DEATH_DATE, format='%d/%m/%Y')]

# Calcola l'eta all'ingresso
anagrafe[, eta := year(dt_ingresso) - year(dt_nascita)]

# Calculate differences in time
anagrafe[, diff_ing_iniass := as.numeric(difftime(dt_ingresso, dt_inizio_assistenza, units = "days"))]
anagrafe[, diff_ing_finass := as.numeric(difftime(dt_fine_assistenza, dt_ingresso, units = "days"))]
anagrafe[, diff_ing_dec := as.numeric(difftime(dt_decesso, dt_ingresso, units = "days"))]


###################### SDO ###################### 
print("reading SDO")

flusso_sdo <- fread("N:/input/DATI/HTCOVID_RICOVERI.csv",
                    select = c("COD_SOGGETTO", "DT_DIMISSIONE",
                               "DIAGP_ID", "DIAG1_ID", "DIAG2_ID", "DIAG3_ID", 
                               "DIAG4_ID", "DIAG5_ID"))[COD_SOGGETTO %in% anagrafe$COD_SOGGETTO]

#set dimissione 
flusso_sdo[,DT_DIMISSIONE := as.Date(parse_sas_date(DT_DIMISSIONE))]

# set the beginning of the follow-up 
flusso_sdo = merge(flusso_sdo, anagrafe[,c("COD_SOGGETTO","dt_ingresso","DATA_START_EVAL")], by = "COD_SOGGETTO")

#keep only info about eval period
flusso_sdo <- flusso_sdo[DATA_START_EVAL < DT_DIMISSIONE & DT_DIMISSIONE <= dt_ingresso]

flusso_sdo = flusso_sdo[, list(COD_SOGGETTO, DT_DIMISSIONE, dt_ingresso,
                               DIAGP_ID, DIAG1_ID, DIAG2_ID, DIAG3_ID, 
                               DIAG4_ID, DIAG5_ID)] 

gc(full = T)


###################### FARMACEUTICA ###################### 

print("reading pharma")

##parallel

# Funzione per caricare e filtrare i dati
load_and_filter <- function(file_path, lapse) {
  dt = fread(file_path, select = c("COD_SOGGETTO","DT_EROGAZIONE","ATC"))
  dt = merge(dt, lapse, by ="COD_SOGGETTO")
  dt[, DT_EROGAZIONE := parse_sas_date(DT_EROGAZIONE)]
  dt = dt[DT_EROGAZIONE <= DATA_INIZIALE]
  dt = dt[DT_EROGAZIONE >= DATA_START_EVAL]
  return(unique(dt[, .(COD_SOGGETTO, ATC)]))
}

# Definizione dei percorsi e delle condizioni

#get maximum

ymax = year(max(anagrafe$dt_ingresso))
ymin = year(min(anagrafe$dt_ingresso) - time_of_eval)


pattern = paste0("HTCOVID_FAR_TERR_",ymin:ymax,".csv",collapse = "|")

file_paths <- list.files("N:/input/DATI", full=T)
file_paths = file_paths[grepl(pattern, file_paths)]
file_paths = c(file_paths,  "N:/input/DATI/HTCOVID_FILEF.csv")


# Creazione di un cluster con 2 core
cluster <- makeCluster(7)

# create lapse dataset
lapse = anagrafe[,c("COD_SOGGETTO","DATA_INIZIALE","DATA_START_EVAL")]

# Esportare le variabili e le funzioni necessarie nel cluster
clusterExport(cluster, c("load_and_filter", "lapse", "file_paths","parse_sas_date"))
clusterEvalQ(cluster, {
  library(data.table)
  library(lubridate)
})

# Parallelizzazione del caricamento e filtraggio dei file
results <- parLapply(cluster, 1:length(file_paths), function(i) {
  file_path <- file_paths[[i]]
  load_and_filter(file_path, lapse)
})

# Unire i risultati
farmaceutica <- rbindlist(results, use.names = TRUE, fill = TRUE)[
  , unique(.SD)]

# Pulizia della memoria
rm(results)
gc(full = TRUE)

# Arrestare il cluster
stopCluster(cluster)
rm(cluster)




###################### MACROPATOLOGIE ######################## 

print("calculating")

convert_IDate_to_Date <- function(dt) {
  # Iterate over each column
  for (col in names(dt)) {
    # Check if the column contains elements of class IDate
    if (all(sapply(dt[[col]], function(x) class(x) == "IDate"))) {
      # Convert IDate to Date
      dt[[col]] <- as.Date(dt[[col]])
    }
  }
  return(dt)
}


# Apply the function to the dataset
#anagrafe <- convert_IDate_to_Date(anagrafe)

# Check the class of each variable after conversion
sapply(anagrafe, class)
anagrafe$DEATH_DATE = as.Date(anagrafe$DEATH_DATE)
sapply(anagrafe, class)
anagrafe = data.table(anagrafe)

sapply(flusso_sdo, class)
flusso_sdo = data.table(flusso_sdo)
farmaceutica = data.table(farmaceutica)
gc(full = T)

#### HIV ####

substr_conditions <- c(
  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 042 & as.numeric(substr(DIAGP_ID, 1, 3)) <= 044)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 042 & as.numeric(substr(DIAG1_ID, 1, 3)) <= 044)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 042 & as.numeric(substr(DIAG2_ID, 1, 3)) <= 044)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 042 & as.numeric(substr(DIAG3_ID, 1, 3)) <= 044)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 042 & as.numeric(substr(DIAG4_ID, 1, 3)) <= 044)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 042 & as.numeric(substr(DIAG5_ID, 1, 3)) <= 044)"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "305AF01", "J05AR01", "J05AR04", "J05AR05", "J05AR06",
  "J05AF03", "P01CX01", "301FA09", "J04AB04", "P01AX06"
)

# data.table syntax                                                      
anagrafe <- anagrafe[, HIV_PRE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
          COD_SOGGETTO %in% unique(farmaceutica[substr(ATC, 1, 7) %in% atc_conditions, COD_SOGGETTO]),
          HIV_PRE := 1]


gc(full = T)
#### TUBERCOLOSI ####
substr_conditions <- c(
  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 010 & as.numeric(substr(DIAGP_ID, 1, 3)) <= 018)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 010 & as.numeric(substr(DIAG1_ID, 1, 3)) <= 018)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 010 & as.numeric(substr(DIAG2_ID, 1, 3)) <= 018)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 010 & as.numeric(substr(DIAG3_ID, 1, 3)) <= 018)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 010 & as.numeric(substr(DIAG4_ID, 1, 3)) <= 018)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 010 & as.numeric(substr(DIAG5_ID, 1, 3)) <= 018)"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("J04AB")

# data.table syntax                                                      
anagrafe <- anagrafe[, TUBERCOLOSI := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
          COD_SOGGETTO %in% unique(farmaceutica[substr(ATC, 1, 5) %in% atc_conditions, COD_SOGGETTO]),
          TUBERCOLOSI := 1]

gc(full = T)



#### LINFOMA  ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('2038','2386','2733')",
  "substr(DIAG1_ID,1,4) %in% c('2038','2386','2733')",
  "substr(DIAG2_ID,1,4) %in% c('2038','2386','2733')",
  "substr(DIAG3_ID,1,4) %in% c('2038','2386','2733')",
  "substr(DIAG4_ID,1,4) %in% c('2038','2386','2733')",
  "substr(DIAG5_ID,1,4) %in% c('2038','2386','2733')",
  "substr(DIAGP_ID,1,5) %in% c('V1071','V1072','V1079')",
  "substr(DIAG1_ID,1,5) %in% c('V1071','V1072','V1079')",
  "substr(DIAG2_ID,1,5) %in% c('V1071','V1072','V1079')",
  "substr(DIAG3_ID,1,5) %in% c('V1071','V1072','V1079')",
  "substr(DIAG4_ID,1,5) %in% c('V1071','V1072','V1079')",
  "substr(DIAG5_ID,1,5) %in% c('V1071','V1072','V1079')",
  "200 <= as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 201",
  "200 <= as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAG1_ID,1,3)) <= 201",
  "200 <= as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAG2_ID,1,3)) <= 201",
  "200 <= as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAG3_ID,1,3)) <= 201",
  "200 <= as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAG4_ID,1,3)) <= 201",
  "200 <= as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAG5_ID,1,3)) <= 201",
  "20200 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 20238",
  "20200 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 20238",
  "20200 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 20238",
  "20200 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 20238",
  "20200 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 20238",
  "20200 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 20238",
  "20250 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 20301",
  "20250 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 20301",
  "20250 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 20301",
  "20250 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 20301",
  "20250 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 20301",
  "20250 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 20301"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, LINFOMA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          LINFOMA := 1]

gc(full = T)

#### CANCRO METASTATICO ####

substr_conditions <- c(
  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 196 & as.numeric(substr(DIAGP_ID, 1, 3)) <= 198)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 196 & as.numeric(substr(DIAG1_ID, 1, 3)) <= 198)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 196 & as.numeric(substr(DIAG2_ID, 1, 3)) <= 198)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 196 & as.numeric(substr(DIAG3_ID, 1, 3)) <= 198)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 196 & as.numeric(substr(DIAG4_ID, 1, 3)) <= 198)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 196 & as.numeric(substr(DIAG5_ID, 1, 3)) <= 198)",
  
  "(as.numeric(substr(DIAGP_ID, 1, 4)) >= 1990 & as.numeric(substr(DIAGP_ID, 1, 4)) <= 1991)",
  "(as.numeric(substr(DIAG1_ID, 1, 4)) >= 1990 & as.numeric(substr(DIAG1_ID, 1, 4)) <= 1991)",
  "(as.numeric(substr(DIAG2_ID, 1, 4)) >= 1990 & as.numeric(substr(DIAG2_ID, 1, 4)) <= 1991)",
  "(as.numeric(substr(DIAG3_ID, 1, 4)) >= 1990 & as.numeric(substr(DIAG3_ID, 1, 4)) <= 1991)",
  "(as.numeric(substr(DIAG4_ID, 1, 4)) >= 1990 & as.numeric(substr(DIAG4_ID, 1, 4)) <= 1991)",
  "(as.numeric(substr(DIAG5_ID, 1, 4)) >= 1990 & as.numeric(substr(DIAG5_ID, 1, 4)) <= 1991)"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, CANCRO_MET := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          CANCRO_MET := 1]

gc(full = T)

#### CANCRO NON METASTATICO ####
substr_conditions <- c(
  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 140 & as.numeric(substr(DIAGP_ID, 1, 3)) <= 172)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 140 & as.numeric(substr(DIAG1_ID, 1, 3)) <= 172)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 140 & as.numeric(substr(DIAG2_ID, 1, 3)) <= 172)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 140 & as.numeric(substr(DIAG3_ID, 1, 3)) <= 172)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 140 & as.numeric(substr(DIAG4_ID, 1, 3)) <= 172)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 140 & as.numeric(substr(DIAG5_ID, 1, 3)) <= 172)",
  
  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 174 & as.numeric(substr(DIAGP_ID, 1, 4)) <= 175)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 174 & as.numeric(substr(DIAG1_ID, 1, 4)) <= 175)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 174 & as.numeric(substr(DIAG2_ID, 1, 4)) <= 175)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 174 & as.numeric(substr(DIAG3_ID, 1, 4)) <= 175)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 174 & as.numeric(substr(DIAG4_ID, 1, 4)) <= 175)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 174 & as.numeric(substr(DIAG5_ID, 1, 4)) <= 175)",

  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 179 & as.numeric(substr(DIAGP_ID, 1, 4)) <= 195)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 179 & as.numeric(substr(DIAG1_ID, 1, 4)) <= 195)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 179 & as.numeric(substr(DIAG2_ID, 1, 4)) <= 195)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 179 & as.numeric(substr(DIAG3_ID, 1, 4)) <= 195)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 179 & as.numeric(substr(DIAG4_ID, 1, 4)) <= 195)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 179 & as.numeric(substr(DIAG5_ID, 1, 4)) <= 195)",
  
  "(substr(DIAGP_ID,1,1) %in% c('V') & substr(DIAGP_ID,2,2) %in% c('10'))",
  "(substr(DIAG1_ID,1,1) %in% c('V') & substr(DIAG1_ID,2,2) %in% c('10'))",
  "(substr(DIAG2_ID,1,1) %in% c('V') & substr(DIAG2_ID,2,2) %in% c('10'))",
  "(substr(DIAG3_ID,1,1) %in% c('V') & substr(DIAG3_ID,2,2) %in% c('10'))",
  "(substr(DIAG4_ID,1,1) %in% c('V') & substr(DIAG4_ID,2,2) %in% c('10'))",
  "(substr(DIAG5_ID,1,1) %in% c('V') & substr(DIAG5_ID,2,2) %in% c('10'))"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, CANCRO_NOMET := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          CANCRO_NOMET := 1]


gc(full = T)
#### DIABETE ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('250') | substr(DIAGP_ID,1,4) %in% c('3620','3572')",
  "substr(DIAG1_ID,1,3) %in% c('250') | substr(DIAG1_ID,1,4) %in% c('3620','3572')",
  "substr(DIAG2_ID,1,3) %in% c('250') | substr(DIAG2_ID,1,4) %in% c('3620','3572')",
  "substr(DIAG3_ID,1,3) %in% c('250') | substr(DIAG3_ID,1,4) %in% c('3620','3572')",
  "substr(DIAG4_ID,1,3) %in% c('250') | substr(DIAG4_ID,1,4) %in% c('3620','3572')",
  "substr(DIAG5_ID,1,3) %in% c('250') | substr(DIAG5_ID,1,4) %in% c('3620','3572')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("A10")

# data.table syntax                                                      
anagrafe <- anagrafe[, DIABETE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
          COD_SOGGETTO %in% unique(farmaceutica[substr(ATC, 1, 3) %in% atc_conditions, COD_SOGGETTO]),
          DIABETE := 1]

gc(full = T)

#### IPOTIROIDISMO  ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('243')",
  "substr(DIAG1_ID,1,3) %in% c('243')",
  "substr(DIAG2_ID,1,3) %in% c('243')",
  "substr(DIAG3_ID,1,3) %in% c('243')",
  "substr(DIAG4_ID,1,3) %in% c('243')",
  "substr(DIAG5_ID,1,3) %in% c('243')",

  "substr(DIAGP_ID,1,4) %in% c('2440','2441','2442','2448','2449')",
  "substr(DIAG1_ID,1,4) %in% c('2440','2441','2442','2448','2449')",
  "substr(DIAG2_ID,1,4) %in% c('2440','2441','2442','2448','2449')",
  "substr(DIAG3_ID,1,4) %in% c('2440','2441','2442','2448','2449')",
  "substr(DIAG4_ID,1,4) %in% c('2440','2441','2442','2448','2449')",
  "substr(DIAG5_ID,1,4) %in% c('2440','2441','2442','2448','2449')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- "substr(ATC, 1, 4) %in% c('H03A','H03B')"

# data.table syntax                                                      
anagrafe <- anagrafe[, IPOTIROIDISMO := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          IPOTIROIDISMO := 1]

gc(full = T)

##### OBESITA ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,5) %in% c('27800','27801')",
  "substr(DIAG1_ID,1,5) %in% c('27800','27801')",
  "substr(DIAG2_ID,1,5) %in% c('27800','27801')",
  "substr(DIAG3_ID,1,5) %in% c('27800','27801')",
  "substr(DIAG4_ID,1,5) %in% c('27800','27801')",
  "substr(DIAG5_ID,1,5) %in% c('27800','27801')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, OBESITA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          OBESITA := 1]

gc(full = T)

#### PERDITA PESO ####

substr_conditions <- c(
  "(as.numeric(substr(DIAGP_ID, 1, 3)) >= 260 & as.numeric(substr(DIAGP_ID, 1, 3)) <= 263)",
  "(as.numeric(substr(DIAG1_ID, 1, 3)) >= 260 & as.numeric(substr(DIAG1_ID, 1, 3)) <= 263)",
  "(as.numeric(substr(DIAG2_ID, 1, 3)) >= 260 & as.numeric(substr(DIAG2_ID, 1, 3)) <= 263)",
  "(as.numeric(substr(DIAG3_ID, 1, 3)) >= 260 & as.numeric(substr(DIAG3_ID, 1, 3)) <= 263)",
  "(as.numeric(substr(DIAG4_ID, 1, 3)) >= 260 & as.numeric(substr(DIAG4_ID, 1, 3)) <= 263)",
  "(as.numeric(substr(DIAG5_ID, 1, 3)) >= 260 & as.numeric(substr(DIAG5_ID, 1, 3)) <= 263)"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, PERDITA_PESO := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          PERDITA_PESO := 1]

gc(full = T)

#### ALTERAZIONE FLUIDI ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('276')",
  "substr(DIAG1_ID,1,3) %in% c('276')",
  "substr(DIAG2_ID,1,3) %in% c('276')",
  "substr(DIAG3_ID,1,3) %in% c('276')",
  "substr(DIAG4_ID,1,3) %in% c('276')",
  "substr(DIAG5_ID,1,3) %in% c('276')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, ALT_FLUIDI := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          ALT_FLUIDI := 1]


gc(full = T)

#### GOTTA ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('274')",
  "substr(DIAG1_ID,1,3) %in% c('274')",
  "substr(DIAG2_ID,1,3) %in% c('274')",
  "substr(DIAG3_ID,1,3) %in% c('274')",
  "substr(DIAG4_ID,1,3) %in% c('274')",
  "substr(DIAG5_ID,1,3) %in% c('274')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('M04AA','M04AB')",
	                  "substr(ATC,1,7) %in% c('M04AC01')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, GOTTA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
          COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          GOTTA := 1]


gc(full = T)

#### COAGULOPATIA #####


substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('286')",
  "substr(DIAG1_ID,1,3) %in% c('286')",
  "substr(DIAG2_ID,1,3) %in% c('286')",
  "substr(DIAG3_ID,1,3) %in% c('286')",
  "substr(DIAG4_ID,1,3) %in% c('286')",
  "substr(DIAG5_ID,1,3) %in% c('286')",
  
  "substr(DIAGP_ID,1,4) %in% c('2871')",
  "substr(DIAG1_ID,1,4) %in% c('2871')",
  "substr(DIAG2_ID,1,4) %in% c('2871')",
  "substr(DIAG3_ID,1,4) %in% c('2871')",
  "substr(DIAG4_ID,1,4) %in% c('2871')",
  "substr(DIAG5_ID,1,4) %in% c('2871')",
  
  "2873 <= as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 2875",
  "2873 <= as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAG1_ID,1,4)) <= 2875",
  "2873 <= as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAG2_ID,1,4)) <= 2875",
  "2873 <= as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAG3_ID,1,4)) <= 2875",
  "2873 <= as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAG4_ID,1,4)) <= 2875",
  "2873 <= as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAG5_ID,1,4)) <= 2875"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, COAGULOPATIA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          COAGULOPATIA := 1]

gc(full = T)

#### ANEMIA ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('2859')",
  "substr(DIAG1_ID,1,4) %in% c('2859')",
  "substr(DIAG2_ID,1,4) %in% c('2859')",
  "substr(DIAG3_ID,1,4) %in% c('2859')",
  "substr(DIAG4_ID,1,4) %in% c('2859')",
  "substr(DIAG5_ID,1,4) %in% c('2859')",
  
  
  "280 <= as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 281",
  "280 <= as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAG1_ID,1,3)) <= 281",
  "280 <= as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAG2_ID,1,3)) <= 281",
  "280 <= as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAG3_ID,1,3)) <= 281",
  "280 <= as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAG4_ID,1,3)) <= 281",
  "280 <= as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAG5_ID,1,3)) <= 281"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('L03AA')",
                    "substr(ATC,1,7) %in% c('B03XA01')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, ANEMIA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          ANEMIA := 1]

gc(full = T)


##### DEMENZA ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('290')",
  "substr(DIAG1_ID,1,3) %in% c('290')",
  "substr(DIAG2_ID,1,3) %in% c('290')",
  "substr(DIAG3_ID,1,3) %in% c('290')",
  "substr(DIAG4_ID,1,3) %in% c('290')",
  "substr(DIAG5_ID,1,3) %in% c('290')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, DEMENZA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          DEMENZA := 1]

gc(full = T)

#### PSICOSI ####

substr_conditions <- c(
  "29910 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 29911",
  "29910 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 29911",
  "29910 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 29911",
  "29910 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 29911",
  "29910 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 29911",
  "29910 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 29911",
  
  "295 <= as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 298",
  "295 <= as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAG1_ID,1,3)) <= 298",
  "295 <= as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAG2_ID,1,3)) <= 298",
  "295 <= as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAG3_ID,1,3)) <= 298",
  "295 <= as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAG4_ID,1,3)) <= 298",
  "295 <= as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAG5_ID,1,3)) <= 298"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('N05AD','N05AA','N05AB','N05AC','N05AX')",
                    "substr(ATC,1,7) %in% c('N05AF04')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, PSICOSI := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          PSICOSI := 1]


gc(full = T)

#### DEPRESSIONE ####


substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('311')",
  "substr(DIAG1_ID,1,3) %in% c('311')",
  "substr(DIAG2_ID,1,3) %in% c('311')",
  "substr(DIAG3_ID,1,3) %in% c('311')",
  "substr(DIAG4_ID,1,3) %in% c('311')",
  "substr(DIAG5_ID,1,3) %in% c('311')",
  
  "substr(DIAGP_ID,1,4) %in% c('3004','3090','3091')",
  "substr(DIAG1_ID,1,4) %in% c('3004','3090','3091')",
  "substr(DIAG2_ID,1,4) %in% c('3004','3090','3091')",
  "substr(DIAG3_ID,1,4) %in% c('3004','3090','3091')",
  "substr(DIAG4_ID,1,4) %in% c('3004','3090','3091')",
  "substr(DIAG5_ID,1,4) %in% c('3004','3090','3091')",
  
  "substr(DIAGP_ID,1,5) %in% c('30112')",
  "substr(DIAG1_ID,1,5) %in% c('30112')",
  "substr(DIAG2_ID,1,5) %in% c('30112')",
  "substr(DIAG3_ID,1,5) %in% c('30112')",
  "substr(DIAG4_ID,1,5) %in% c('30112')",
  "substr(DIAG5_ID,1,5) %in% c('30112')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('N06AA','N04BD','N06AF','N06AG')",
                    "substr(ATC,1,7) %in% c('N06AB03','N06CA03')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, DEPRESSIONE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
          COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          DEPRESSIONE := 1]


gc(full = T)

#### ABUSO ALCOL ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('2911','2912','2915','2918','2919', 'V113')",
  "substr(DIAG1_ID,1,4) %in% c('2911','2912','2915','2918','2919', 'V113')",
  "substr(DIAG2_ID,1,4) %in% c('2911','2912','2915','2918','2919', 'V113')",
  "substr(DIAG3_ID,1,4) %in% c('2911','2912','2915','2918','2919', 'V113')",
  "substr(DIAG4_ID,1,4) %in% c('2911','2912','2915','2918','2919', 'V113')",
  "substr(DIAG5_ID,1,4) %in% c('2911','2912','2915','2918','2919', 'V113')",
  
  "30390 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 30393", 
  "30390 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 30393", 
  "30390 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 30393", 
  "30390 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 30393", 
  "30390 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 30393", 
  "30390 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 30393",
  
  "30500 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 30503", 
  "30500 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 30503", 
  "30500 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 30503", 
  "30500 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 30503", 
  "30500 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 30503", 
  "30500 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 30503"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, AB_ALCOL := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          AB_ALCOL := 1]


gc(full = T)




#### DISTURBO BIPOLARE ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('2960')",
  "substr(DIAG1_ID,1,4) %in% c('2960')",
  "substr(DIAG2_ID,1,4) %in% c('2960')",
  "substr(DIAG3_ID,1,4) %in% c('2960')",
  "substr(DIAG4_ID,1,4) %in% c('2960')",
  "substr(DIAG5_ID,1,4) %in% c('2960')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('N05AN')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, DIST_BIPOLARE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          DIST_BIPOLARE := 1]


gc(full = T)



#### PARALISI ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('3420','3429')",
  "substr(DIAG1_ID,1,4) %in% c('3420','3429')",
  "substr(DIAG2_ID,1,4) %in% c('3420','3429')",
  "substr(DIAG3_ID,1,4) %in% c('3420','3429')",
  "substr(DIAG4_ID,1,4) %in% c('3420','3429')",
  "substr(DIAG5_ID,1,4) %in% c('3420','3429')",
  
  "343 <= as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 344", 
  "343 <= as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAG1_ID,1,3)) <= 344", 
  "343 <= as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAG2_ID,1,3)) <= 344", 
  "343 <= as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAG3_ID,1,3)) <= 344", 
  "343 <= as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAG4_ID,1,3)) <= 344", 
  "343 <= as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAG5_ID,1,3)) <= 344", 
  
  "34210 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 34212", 
  "34210 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 34212", 
  "34210 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 34212", 
  "34210 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 34212", 
  "34210 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 34212", 
  "34210 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 34212"
  
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, PARALISI := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          PARALISI := 1]


gc(full = T)
#### ALTRE MALATTIE NEUROLOGICHE ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('340')",
  "substr(DIAG1_ID,1,3) %in% c('340')",
  "substr(DIAG2_ID,1,3) %in% c('340')",
  "substr(DIAG3_ID,1,3) %in% c('340')",
  "substr(DIAG4_ID,1,3) %in% c('340')",
  "substr(DIAG5_ID,1,3) %in% c('340')",
  
  "substr(DIAGP_ID,1,4) %in% c('3319','3320','3334','3335','3450','3454','3481','3483','3458','7803','7843')",
  "substr(DIAG1_ID,1,4) %in% c('3319','3320','3334','3335','3450','3454','3481','3483','3458','7803','7843')",
  "substr(DIAG2_ID,1,4) %in% c('3319','3320','3334','3335','3450','3454','3481','3483','3458','7803','7843')",
  "substr(DIAG3_ID,1,4) %in% c('3319','3320','3334','3335','3450','3454','3481','3483','3458','7803','7843')",
  "substr(DIAG4_ID,1,4) %in% c('3319','3320','3334','3335','3450','3454','3481','3483','3458','7803','7843')",
  "substr(DIAG5_ID,1,4) %in% c('3319','3320','3334','3335','3450','3454','3481','3483','3458','7803','7843')",
  
  "334 <= as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 335",
  "334 <= as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAG1_ID,1,3)) <= 335",
  "334 <= as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAG2_ID,1,3)) <= 335", 
  "334 <= as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAG3_ID,1,3)) <= 335", 
  "334 <= as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAG4_ID,1,3)) <= 335", 
  "334 <= as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAG5_ID,1,3)) <= 335", 
  
  "3411 <= as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3419", 
  "3411 <= as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAG1_ID,1,4)) <= 3419", 
  "3411 <= as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAG2_ID,1,4)) <= 3419", 
  "3411 <= as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAG3_ID,1,4)) <= 3419", 
  "3411 <= as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAG4_ID,1,4)) <= 3419", 
  "3411 <= as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAG5_ID,1,4)) <= 3419", 
  
  "34510 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 34511", 
  "34510 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 34511", 
  "34510 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 34511", 
  "34510 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 34511", 
  "34510 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 34511", 
  "34510 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 34511",

  "34550 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 34551", 
  "34550 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 34551", 
  "34550 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 34551", 
  "34550 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 34551", 
  "34550 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 34551", 
  "34550 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 34551", 
  
  "34590 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 34591", 
  "34590 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAG1_ID,1,5)) <= 34591", 
  "34590 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAG2_ID,1,5)) <= 34591", 
  "34590 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAG3_ID,1,5)) <= 34591", 
  "34590 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAG4_ID,1,5)) <= 34591", 
  "34590 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAG5_ID,1,5)) <= 34591"
  
)

substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, ALTRE_MAL_NEURO := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          ALTRE_MAL_NEURO := 1]


gc(full = T)

#### GLAUCOMA ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('365')",
  "substr(DIAG1_ID,1,3) %in% c('365')",
  "substr(DIAG2_ID,1,3) %in% c('365')",
  "substr(DIAG3_ID,1,3) %in% c('365')",
  "substr(DIAG4_ID,1,3) %in% c('365')",
  "substr(DIAG5_ID,1,3) %in% c('365')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,4) %in% c('S01E')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, GLAUCOMA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          GLAUCOMA := 1]


gc(full = T)

#### EPILESSIA ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('345')",
  "substr(DIAG1_ID,1,3) %in% c('345')",
  "substr(DIAG2_ID,1,3) %in% c('345')",
  "substr(DIAG3_ID,1,3) %in% c('345')",
  "substr(DIAG4_ID,1,3) %in% c('345')",
  "substr(DIAG5_ID,1,3) %in% c('345')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('N03AA','N03AX')",
                    "substr(ATC,1,7) %in% c('N03AB02','N03AB05','N03AB52')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, EPILESSIA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          EPILESSIA := 1]


gc(full = T)
#### PARKINSON ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('332')",
  "substr(DIAG1_ID,1,3) %in% c('332')",
  "substr(DIAG2_ID,1,3) %in% c('332')",
  "substr(DIAG3_ID,1,3) %in% c('332')",
  "substr(DIAG4_ID,1,3) %in% c('332')",
  "substr(DIAG5_ID,1,3) %in% c('332')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,4) %in% c('N04B')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, PARKINSON := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          PARKINSON := 1]


gc(full = T)
#### AMI  ####

substr_conditions <- c(
  
  "410 <= as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 412", 
  "410 <= as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAG1_ID,1,3)) <= 412", 
  "410 <= as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAG2_ID,1,3)) <= 412", 
  "410 <= as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAG3_ID,1,3)) <= 412", 
  "410 <= as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAG4_ID,1,3)) <= 412", 
  "410 <= as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAG5_ID,1,3)) <= 412"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, AMI := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          AMI := 1]

gc(full = T)
#### HEART FAILURE ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,5) %in% c('39891','40211','40291','40411','40413','40491','40493')",
  "substr(DIAG1_ID,1,5) %in% c('39891','40211','40291','40411','40413','40491','40493')",
  "substr(DIAG2_ID,1,5) %in% c('39891','40211','40291','40411','40413','40491','40493')",
  "substr(DIAG3_ID,1,5) %in% c('39891','40211','40291','40411','40413','40491','40493')",
  "substr(DIAG4_ID,1,5) %in% c('39891','40211','40291','40411','40413','40491','40493')",
  "substr(DIAG5_ID,1,5) %in% c('39891','40211','40291','40411','40413','40491','40493')",
  
  "substr(DIAGP_ID,1,3) %in% c('428')",
  "substr(DIAG1_ID,1,3) %in% c('428')",
  "substr(DIAG2_ID,1,3) %in% c('428')",
  "substr(DIAG3_ID,1,3) %in% c('428')",
  "substr(DIAG4_ID,1,3) %in% c('428')",
  "substr(DIAG5_ID,1,3) %in% c('428')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,4) %in% c('C03C')",
                    "substr(ATC,1,5) %in% c('C01DA','C01AA')",
                    "substr(ATC,1,7) %in% c('C01BA93','C01BA02','C01BA01','C01BA51','C01BA71')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, HF := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          HF := 1]

gc(full = T)


#### ARITMIA ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,5) %in% c('42610','42611','42613','42731','42760')",
  "substr(DIAG1_ID,1,5) %in% c('42610','42611','42613','42731','42760')",
  "substr(DIAG2_ID,1,5) %in% c('42610','42611','42613','42731','42760')",
  "substr(DIAG3_ID,1,5) %in% c('42610','42611','42613','42731','42760')",
  "substr(DIAG4_ID,1,5) %in% c('42610','42611','42613','42731','42760')",
  "substr(DIAG5_ID,1,5) %in% c('42610','42611','42613','42731','42760')",
  
  "substr(DIAGP_ID,1,4) %in% c('4270','4272','4279','7850','V450','V533')",  
  "substr(DIAG1_ID,1,4) %in% c('4270','4272','4279','7850','V450','V533')",  
  "substr(DIAG2_ID,1,4) %in% c('4270','4272','4279','7850','V450','V533')",  
  "substr(DIAG3_ID,1,4) %in% c('4270','4272','4279','7850','V450','V533')",  
  "substr(DIAG4_ID,1,4) %in% c('4270','4272','4279','7850','V450','V533')",  
  "substr(DIAG5_ID,1,4) %in% c('4270','4272','4279','7850','V450','V533')",  
  
  "4262 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4264",
  "4262 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4264",
  "4262 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4264",
  "4262 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4264",
  "4262 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4264",
  "4262 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4264",
  
  "42650 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42653",
  "42650 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42653",
  "42650 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42653",
  "42650 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42653",
  "42650 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42653",
  "42650 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42653",
  
  "4266 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4268",
  "4266 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4268",
  "4266 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4268",
  "4266 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4268",
  "4266 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4268",
  "4266 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4268"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c("substr(ATC,1,5) %in% c('C01BA','C01BC','C01BD')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, ARITMIA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          ARITMIA := 1]

gc(full = T)

#### MALATTIA VALVOLARE ####


substr_conditions <- c(
  "substr(DIAGP_ID,1,5) %in% c('09320','09321','09322','09323','09324')",
  "substr(DIAG1_ID,1,5) %in% c('09320','09321','09322','09323','09324')",
  "substr(DIAG2_ID,1,5) %in% c('09320','09321','09322','09323','09324')",
  "substr(DIAG3_ID,1,5) %in% c('09320','09321','09322','09323','09324')",
  "substr(DIAG4_ID,1,5) %in% c('09320','09321','09322','09323','09324')",
  "substr(DIAG5_ID,1,5) %in% c('09320','09321','09322','09323','09324')",
  
  "substr(DIAGP_ID,1,4) %in% c('V422','V433')",  
  "substr(DIAG1_ID,1,4) %in% c('V422','V433')",  
  "substr(DIAG2_ID,1,4) %in% c('V422','V433')",  
  "substr(DIAG3_ID,1,4) %in% c('V422','V433')",  
  "substr(DIAG4_ID,1,4) %in% c('V422','V433')",  
  "substr(DIAG5_ID,1,4) %in% c('V422','V433')",  
  
  "4240 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4248",
  "4240 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4248",
  "4240 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4248",
  "4240 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4248",
  "4240 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4248",
  "4240 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4248",
  
  "42490 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42491",
  "42490 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42491",
  "42490 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42491",
  "42490 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42491",
  "42490 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42491",
  "42490 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 42491",
  
  "3970 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3971",
  "3970 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3971",
  "3970 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3971",
  "3970 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3971",
  "3970 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3971",
  "3970 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 3971",
  
  "7463 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 7466",
  "7463 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 7466",
  "7463 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 7466",
  "7463 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 7466",
  "7463 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 7466",
  "7463 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 7466",
  
  "394 <=  as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 396",
  "394 <=  as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 396",
  "394 <=  as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 396",
  "394 <=  as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 396",
  "394 <=  as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 396",
  "394 <=  as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 396"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, MAL_VALV := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          MAL_VALV := 1]


gc(full = T)



#### MALATTIA VASCOLARE ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c(440)",
  "substr(DIAG1_ID,1,3) %in% c(440)",
  "substr(DIAG2_ID,1,3) %in% c(440)",
  "substr(DIAG3_ID,1,3) %in% c(440)",
  "substr(DIAG4_ID,1,3) %in% c(440)",
  "substr(DIAG5_ID,1,3) %in% c(440)",
  
  "substr(DIAGP_ID,1,4) %in% c('4412','4414','4417','4419','4471','5571','5579','7854','V434')",  
  "substr(DIAG1_ID,1,4) %in% c('4412','4414','4417','4419','4471','5571','5579','7854','V434')",  
  "substr(DIAG2_ID,1,4) %in% c('4412','4414','4417','4419','4471','5571','5579','7854','V434')",  
  "substr(DIAG3_ID,1,4) %in% c('4412','4414','4417','4419','4471','5571','5579','7854','V434')",  
  "substr(DIAG4_ID,1,4) %in% c('4412','4414','4417','4419','4471','5571','5579','7854','V434')",  
  "substr(DIAG5_ID,1,4) %in% c('4412','4414','4417','4419','4471','5571','5579','7854','V434')",  
  

  "4431 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4431",
  "4431 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4431",
  "4431 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4431",
  "4431 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4431",
  "4431 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4431",
  "4431 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4431"
  

)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, MAL_VASCOL := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          MAL_VASCOL := 1]

gc(full = T)

#### MALATTIA CEREBROVASCOLARE ####
substr_conditions <- c(
  "430 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 438",
  "430 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 438",
  "430 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 438",
  "430 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 438",
  "430 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 438",
  "430 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 438"
  
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, CEREBRO := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          CEREBRO := 1]


gc(full = T)
#### MALATTIA IPERTENSIVA ####
substr_conditions <- c(
  "401 <=  as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 405",
  "401 <=  as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 405",
  "401 <=  as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 405",
  "401 <=  as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 405",
  "401 <=  as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 405",
  "401 <=  as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 405"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,3) %in% c('C04','C07','C08')",
  "substr(ATC,1,4) %in% c('C09A','C09B','C01D','C07E','C03A','C07B','C07D','C02B')",
  "substr(ATC,1,5) %in% c('C02AB','C02LB','C02AA','C02LA','C02CA','C02LE')",
  "substr(ATC,1,7) %in% c('C04AX02','C04AB01','C04AB02','C02AC01','C02LC01','C02LC51','C02CC02','C02LF01')"
)
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, MAL_IPERTENSIVA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          MAL_IPERTENSIVA := 1]

gc(full = T)



#### MALATTIA POLMONARE ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('494')",
  "substr(DIAG1_ID,1,3) %in% c('494')",
  "substr(DIAG2_ID,1,3) %in% c('494')",
  "substr(DIAG3_ID,1,3) %in% c('494')",
  "substr(DIAG4_ID,1,3) %in% c('494')",
  "substr(DIAG5_ID,1,3) %in% c('494')",
  
  "substr(DIAGP_ID,1,4) %in% c('5064')",  
  "substr(DIAG1_ID,1,4) %in% c('5064')",  
  "substr(DIAG2_ID,1,4) %in% c('5064')",  
  "substr(DIAG3_ID,1,4) %in% c('5064')",  
  "substr(DIAG4_ID,1,4) %in% c('5064')",  
  "substr(DIAG5_ID,1,4) %in% c('5064')",  
  
  "490 <=  as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 492",
  "490 <=  as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 492",
  "490 <=  as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 492",
  "490 <=  as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 492",
  "490 <=  as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 492",
  "490 <=  as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 492",
  
  "495 <=  as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 505",
  "495 <=  as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 505",
  "495 <=  as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 505",
  "495 <=  as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 505",
  "495 <=  as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 505",
  "495 <=  as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 505",
  
  "49390 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 49391",
  "49390 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 49391",
  "49390 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 49391",
  "49390 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 49391",
  "49390 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 49391",
  "49390 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 49391",
  
  "4930 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('R03AA','R03AB','R03AC','R03DA','R03DB','R03BA')",
  "substr(ATC,1,7) %in% c('R03DA20','R01AC01','R03BC01','R01AC51','S01GX01','S01GX51')"
)
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, MAL_POLMONARE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          MAL_POLMONARE := 1]

gc(full = T)



#### FIBROSI CISTICA ####
substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('2770')",  
  "substr(DIAG1_ID,1,4) %in% c('2770')",  
  "substr(DIAG2_ID,1,4) %in% c('2770')",  
  "substr(DIAG3_ID,1,4) %in% c('2770')",  
  "substr(DIAG4_ID,1,4) %in% c('2770')",  
  "substr(DIAG5_ID,1,4) %in% c('2770')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('R05CB')",
  "substr(ATC,1,7) %in% c('R05FB01','R05FA01','A09AA02')"
)
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FIBROSI_CISTICA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          FIBROSI_CISTICA := 1]

gc(full = T)

#### ULCERA ####
substr_conditions <- c(
  "531 <=  as.numeric(substr(DIAGP_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 534",
  "531 <=  as.numeric(substr(DIAG1_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 534",
  "531 <=  as.numeric(substr(DIAG2_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 534",
  "531 <=  as.numeric(substr(DIAG3_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 534",
  "531 <=  as.numeric(substr(DIAG4_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 534",
  "531 <=  as.numeric(substr(DIAG5_ID,1,3)) & as.numeric(substr(DIAGP_ID,1,3)) <= 534"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('A02BB','A02BA')",
  "substr(ATC,1,7) %in% c('A02BC01','A02BD05','A02BD01')"
)
atc_conditions <- paste(atc_conditions, collapse = " | ")


# data.table syntax                                                      
anagrafe <- anagrafe[, ULCERA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          ULCERA := 1]

gc(full = T)
#### MALATTIA FEGATO ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('4560','4561','5710','5712','5713','5715','5716','5718','5719','5723','5728','V427')",  
  "substr(DIAG1_ID,1,4) %in% c('4560','4561','5710','5712','5713','5715','5716','5718','5719','5723','5728','V427')",  
  "substr(DIAG2_ID,1,4) %in% c('4560','4561','5710','5712','5713','5715','5716','5718','5719','5723','5728','V427')",  
  "substr(DIAG3_ID,1,4) %in% c('4560','4561','5710','5712','5713','5715','5716','5718','5719','5723','5728','V427')",  
  "substr(DIAG4_ID,1,4) %in% c('4560','4561','5710','5712','5713','5715','5716','5718','5719','5723','5728','V427')",  
  "substr(DIAG5_ID,1,4) %in% c('4560','4561','5710','5712','5713','5715','5716','5718','5719','5723','5728','V427')",
  
  "substr(DIAGP_ID,1,5) %in% c('07032','07033','07054','45620','45621')",  
  "substr(DIAG1_ID,1,5) %in% c('07032','07033','07054','45620','45621')",  
  "substr(DIAG2_ID,1,5) %in% c('07032','07033','07054','45620','45621')",  
  "substr(DIAG3_ID,1,5) %in% c('07032','07033','07054','45620','45621')",  
  "substr(DIAG4_ID,1,5) %in% c('07032','07033','07054','45620','45621')",  
  "substr(DIAG5_ID,1,5) %in% c('07032','07033','07054','45620','45621')",  
  

  "57140 <= as.numeric(substr(DIAGP_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 57149",
  "57140 <= as.numeric(substr(DIAG1_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 57149",
  "57140 <= as.numeric(substr(DIAG2_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 57149",
  "57140 <= as.numeric(substr(DIAG3_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 57149",
  "57140 <= as.numeric(substr(DIAG4_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 57149",
  "57140 <= as.numeric(substr(DIAG5_ID,1,5)) & as.numeric(substr(DIAGP_ID,1,5)) <= 57149",
  
  "4930 <=  as.numeric(substr(DIAGP_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG1_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG2_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG3_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG4_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938",
  "4930 <=  as.numeric(substr(DIAG5_ID,1,4)) & as.numeric(substr(DIAGP_ID,1,4)) <= 4938"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('A06AD')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, MAL_FEGATO := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          MAL_FEGATO := 1]


gc(full = T)
#### CROHN ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('555','556')",  
  "substr(DIAG1_ID,1,3) %in% c('555','556')",  
  "substr(DIAG2_ID,1,3) %in% c('555','556')",  
  "substr(DIAG3_ID,1,3) %in% c('555','556')",  
  "substr(DIAG4_ID,1,3) %in% c('555','556')",  
  "substr(DIAG5_ID,1,3) %in% c('555','556')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,7) %in% c('A07EC01','A07EC03','A07EC02')")

atc_conditions <- paste(atc_conditions, collapse = " | ")

anagrafe <- anagrafe[, CROHN := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          CROHN := 1]

gc(full = T)
#### MALATTIA RENALE ####


substr_conditions <- c(
  "substr(DIAGP_ID,1,3) %in% c('582','585','586','588')",  
  "substr(DIAG1_ID,1,3) %in% c('582','585','586','588')",  
  "substr(DIAG2_ID,1,3) %in% c('582','585','586','588')",  
  "substr(DIAG3_ID,1,3) %in% c('582','585','586','588')",  
  "substr(DIAG4_ID,1,3) %in% c('582','585','586','588')",  
  "substr(DIAG5_ID,1,3) %in% c('582','585','586','588')",
  
  "substr(DIAGP_ID,1,4) %in% c('5830','5831','5834','5837','5838','5846')",  
  "substr(DIAG1_ID,1,4) %in% c('5830','5831','5834','5837','5838','5846')",  
  "substr(DIAG2_ID,1,4) %in% c('5830','5831','5834','5837','5838','5846')",  
  "substr(DIAG3_ID,1,4) %in% c('5830','5831','5834','5837','5838','5846')",  
  "substr(DIAG4_ID,1,4) %in% c('5830','5831','5834','5837','5838','5846')",  
  "substr(DIAG5_ID,1,4) %in% c('5830','5831','5834','5837','5838','5846')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,7) %in% c('V03AE01')")

atc_conditions <- paste(atc_conditions, collapse = " | ")

anagrafe <- anagrafe[, MAL_RENALE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          MAL_RENALE := 1]


gc(full = T)


#### DIALISI RENALE ####


substr_conditions <- c(
  
  "substr(DIAGP_ID,1,4) %in% c('V560')",  
  "substr(DIAG1_ID,1,4) %in% c('V560')",  
  "substr(DIAG2_ID,1,4) %in% c('V560')",  
  "substr(DIAG3_ID,1,4) %in% c('V560')",  
  "substr(DIAG4_ID,1,4) %in% c('V560')",  
  "substr(DIAG5_ID,1,4) %in% c('V560')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")


anagrafe <- anagrafe[, DIALISI_RENALE := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]),
          DIALISI_RENALE := 1]


gc(full = T)
#### MALATTIA REUMATICA ####

substr_conditions <- c(
  "substr(DIAGP_ID,1,4) %in% c('7140','7141','7143','7149')",  
  "substr(DIAG1_ID,1,4) %in% c('7140','7141','7143','7149')",  
  "substr(DIAG2_ID,1,4) %in% c('7140','7141','7143','7149')",  
  "substr(DIAG3_ID,1,4) %in% c('7140','7141','7143','7149')",  
  "substr(DIAG4_ID,1,4) %in% c('7140','7141','7143','7149')",  
  "substr(DIAG5_ID,1,4) %in% c('7140','7141','7143','7149')",
  
  "substr(DIAGP_ID,1,3) %in% c('390','391','725','720')",  
  "substr(DIAG1_ID,1,3) %in% c('390','391','725','720')",  
  "substr(DIAG2_ID,1,3) %in% c('390','391','725','720')",  
  "substr(DIAG3_ID,1,3) %in% c('390','391','725','720')",  
  "substr(DIAG4_ID,1,3) %in% c('390','391','725','720')",  
  "substr(DIAG5_ID,1,3) %in% c('390','391','725','720')"
)
substr_conditions <- paste(substr_conditions, collapse = " | ")

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('M01BA','M01CB')",
  "substr(ATC,1,7) %in% c('P01BA02')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, MAL_REUMATICA := 0]
anagrafe[ COD_SOGGETTO %in% unique(flusso_sdo[eval(parse(text = substr_conditions)), COD_SOGGETTO]) |
            COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
          MAL_REUMATICA := 1]

gc(full = T)
#### FAR MALIGNANCIES #####

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,3) %in% c('L01','A04')",
  "substr(ATC,1,5) %in% c('L03AC','L03AA')",
  "substr(ATC,1,7) %in% c('C07AB05')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FAR_MALIGNANCIES := 0]
anagrafe[COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
         FAR_MALIGNANCIES := 1]


gc(full = T)

#### FAR ANXIETY ####

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('N05BA','N05CD','N05CF','N05BX','N06BX')",
  "substr(ATC,1,7) %in% c('N05BC01','N05BC51','N05CX01')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FAR_ANXIETY := 0]
anagrafe[COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
         FAR_ANXIETY := 1]

gc(full = T)
#### FAR CORONARY ####

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,5) %in% c('B01AA','B01AB','B01AF','B01AE')",
  "substr(ATC,1,7) %in% c('B01AX01','B01AD10','B01AD12','C04AD03','B01AC05')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FAR_CORONARY := 0]
anagrafe[COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
         FAR_CORONARY := 1]

gc(full = T)

#### FAR HYPERLIPIDEMIA ####

# Define the substrings for the ATC condition
atc_conditions <- c(
  "substr(ATC,1,3) %in% c('C10')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FAR_HYPERLIPIDEMIA := 0]
anagrafe[COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
         FAR_HYPERLIPIDEMIA := 1]

gc(full = T)
#### FAR TRANSPLANTATION ####
atc_conditions <- c(
  "substr(ATC,1,7) %in% c('L04AD01','L04AX01')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FAR_TRANSPLANTATION := 0]
anagrafe[COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
         FAR_TRANSPLANTATION := 1]

gc(full = T)
#### FAR PAIN ####

atc_conditions <- c(
  "substr(ATC,1,3) %in% c('N02')",
  "substr(ATC,1,4) %in% c('M01A')")
atc_conditions <- paste(atc_conditions, collapse = " | ")

# data.table syntax                                                      
anagrafe <- anagrafe[, FAR_PAIN := 0]
anagrafe[COD_SOGGETTO %in% unique(farmaceutica[eval(parse(text = atc_conditions)), COD_SOGGETTO]),
         FAR_PAIN := 1]

gc(full = T)


#### calculate ####
anagrafe$MCS <- 
  18*anagrafe$CANCRO_MET +
  11*anagrafe$AB_ALCOL +
  10*anagrafe$CANCRO_NOMET +
  10*anagrafe$TUBERCOLOSI +
  8*anagrafe$PSICOSI +
  8*anagrafe$MAL_FEGATO +
  6*anagrafe$FAR_ANXIETY +
  6*anagrafe$PERDITA_PESO +
  6*anagrafe$DEMENZA +
  5*anagrafe$FAR_MALIGNANCIES +
  5*anagrafe$PARKINSON +
  5*anagrafe$LINFOMA +
  5*anagrafe$PARALISI +
  5*anagrafe$COAGULOPATIA +
  4*anagrafe$ALT_FLUIDI +
  4*anagrafe$MAL_RENALE +
  4*anagrafe$DIALISI_RENALE +
  4*anagrafe$HF +
  3*anagrafe$ALTRE_MAL_NEURO +
  3*anagrafe$MAL_REUMATICA +
  3*anagrafe$ANEMIA +
  3*anagrafe$CEREBRO +
  2*anagrafe$DIABETE +
  2*anagrafe$MAL_VASCOL +
  2*anagrafe$GOTTA +
  2*anagrafe$EPILESSIA +
  2*anagrafe$MAL_POLMONARE +
  2*anagrafe$ULCERA +
  1*anagrafe$AMI +
  1*anagrafe$FAR_CORONARY +
  1*anagrafe$MAL_VALV +
  1*anagrafe$ARITMIA +
  1*anagrafe$OBESITA +
  1*anagrafe$IPOTIROIDISMO



anagrafe = anagrafe %>%
  mutate(MCS_CLA = case_when(MCS >= 20 ~ 5,
                             MCS >= 15 & MCS <20 ~ 4,
                             MCS >= 10 & MCS <15 ~ 3,
                             MCS >= 5 & MCS <10 ~ 2,
                             MCS >= 1 & MCS <5 ~ 1,
                             MCS == 0 ~ 0,
                             TRUE ~ 0,
                             .default = 0))
print("calculate")

print("writing")

anagrafe_s = anagrafe[,c(-4:-30)]
fwrite(anagrafe_s,fileout)
anagrafe_s = fread(fileout)


time_taken = Sys.time() - start_time
print(paste0("Time:",time_taken))
print("finished")


rm(list = ls())
invisible(gc(full = T))

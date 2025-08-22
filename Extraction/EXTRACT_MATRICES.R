#### setting ####
rm(list = ls())
invisible(gc(full = T))
Sys.setlocale("LC_TIME", "English")

library(dplyr)
library(lubridate)
library(data.table)
setDTthreads(0)
library(parallel)
library(feather)
library(haven)


setwd("N:/output/scripts/acorbetta/")

parse_sas_date = function(data) {
  return(as.Date(strptime(data, "%d%b%Y:%H:%M:%S.%OS")))
}

#### set up parameters ####
#filename wants a csv file without header with id and date of end of evaluation window (END)

#change to csv
filename = "N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv"

#output file
fileout = "N:/output/data/acorbetta/STATINS/ADHERENCE/matrices_data_all.csv.gz"

# how much time from END do you want to evaluate matrices retrospectively
# NB: if they don't have enought time they will be discarded
time_of_eval = 365.25 *3

#change to csv
sum_stat = fread(filename)
colnames(sum_stat) = c("COD_SOGGETTO", "END")


gc(full=T)

print("reading anagrafe")
anagrafe <- fread("N:/output/data/minimum_data.csv", select = c("COD_SOGGETTO", "SESSO", "ANNO_NASCITA", "DEATH_DATE")) 

gc(full = T)

#get start and end of assistance
cols = colnames(fread("N:/input/DATI/HTCOVID_ANAGRAFICA.csv", nrows = 1))
base = fread("N:/input/DATI/HTCOVID_ANAGRAFICA.csv", select =  which(cols %in% c("COD_SOGGETTO", "DATA_FINE_ASSISTENZA", "DATA_INIZIO_ASSISTENZA") | grepl("^MMG_ASS", cols)))
gc(full = T)

anagrafe = anagrafe[COD_SOGGETTO %in% sum_stat$COD_SOGGETTO]
anagrafe = merge(anagrafe,base,  by = "COD_SOGGETTO", all.x = T)
anagrafe = merge(anagrafe,sum_stat,  by = "COD_SOGGETTO", all.x = T)


rm(base)
rm(sum_stat)
gc(full = T)

cod = read_sas("N:/input/DATI/htcovid_codifiche.sas7bdat")
setDT(cod)

#### set up ####
anagrafe[,DATA_INIZIALE := ymd(END)]

# date of evaluation start
anagrafe[,DATA_START_EVAL := ymd(END) - time_of_eval]

lapse = anagrafe[,c("COD_SOGGETTO","DATA_INIZIALE","DATA_START_EVAL")]
lapse = lapse[,YEAR := year(DATA_INIZIALE)]

mmg = anagrafe %>%
  select(COD_SOGGETTO, starts_with("MMG_ASS")) %>%
  melt(id.vars = c("COD_SOGGETTO"), measure.vars = patterns("^MMG_ASS"), 
                  variable.name = "YEAR", value.name = "MEDICO")
mmg =  mmg[,YEAR := as.numeric(substr(YEAR,9,12))]

lapse = merge(lapse,mmg, by= c("COD_SOGGETTO", "YEAR"),all.x = T)
rm(mmg)
gc(full=T)

#### ADI #### 
print("reading ADI")

adi <- fread("N:/input/DATI/HTCOVID_ADI_ACCESSI.csv",
             select = c("COD_SOGGETTO", "DATA_ACCESSO", "NUM_ACCESSI","QTA_VISITE_DOMICILIARI"))[COD_SOGGETTO %in% lapse$COD_SOGGETTO]
adi[, DATA_ACCESSO:= parse_sas_date(DATA_ACCESSO)]

adi = merge(adi,lapse, by="COD_SOGGETTO")
adi = adi[DATA_ACCESSO <= DATA_INIZIALE]
adi = adi[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR",'MEDICO')]

adi = adi %>%
  group_by(COD_SOGGETTO) %>%
  summarise(ADI_ACCESSI = sum(NUM_ACCESSI),
            ADI_VISITE = sum(QTA_VISITE_DOMICILIARI))

adi = merge(lapse,adi, by="COD_SOGGETTO",all.x = T, )
adi = adi[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]
colnames(adi)
adi[is.na(adi)] <- 0

#### SOSIA ###################### 
print("reading SOSIA")

sos <- fread("N:/input/DATI/HTCOVID_SOSIA_AMMISSIONI.csv",
             select = c("COD_SOGGETTO", "DATA_AMMISSIONE","DATA_DIMISSIONE"))[COD_SOGGETTO %in% lapse$COD_SOGGETTO]
sos[, DATA_AMMISSIONE:= parse_sas_date(DATA_AMMISSIONE)]

sos = merge(sos,lapse, by="COD_SOGGETTO")
sos = sos[DATA_AMMISSIONE <= DATA_INIZIALE]
sos[, DATA_DIMISSIONE:= parse_sas_date(DATA_DIMISSIONE)]
sos[, SOSIA_TIME :=as.numeric(difftime(pmin(DATA_DIMISSIONE,DATA_INIZIALE, na.rm =T), DATA_AMMISSIONE, units = "days"))]
sos = sos[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]
colnames(sos)

sos = sos %>%
  group_by(COD_SOGGETTO) %>%
  summarise(SOS_TIME = sum(SOSIA_TIME))


sos = merge(lapse,sos, by="COD_SOGGETTO",all.x = T, )
sos = sos[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]
colnames(sos)
sos[is.na(sos)] <- 0


#### MEDICI ####

print("reading MEDICI")

med <- fread("N:/input/DATI/HTCOVID_MEDICI.csv",
             select = c("COD_REG","SESSO","ANNO_NASCITA","ATS","ANNO_LAUREA","ANNO_SPECIALIZZAZIONE","DISTRETTO","NUMERO_ASSISTITI","TIPO_PRESCRITTORE"))

med = merge(lapse,med, by.x="MEDICO", by.y = "COD_REG",all.x = T)
med = med[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]
table(med$SESSO)

#### VACCINI ####
library(tidyr)

cod_vac = cod[TIPO == "VACCINI"]

gc(full=T)
vac = fread("N:/input/DATI/HTCOVID_VACCINAZIONI.csv",
            select = c(1,4,10))[COD_SOGGETTO %in% lapse$COD_SOGGETTO]
vac[, DATA_VACCINO:= parse_sas_date(DATA_VACCINO)]


vac = merge(vac,lapse, by="COD_SOGGETTO")
vac = vac[DATA_VACCINO <= DATA_INIZIALE]
vac = vac[DATA_VACCINO >= DATA_START_EVAL]
vac = vac[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]

vac2 = vac %>%
  group_by(COD_SOGGETTO, VACCINO_ID) %>%
  summarise(N_VACC = n())
vac2 = merge(vac2,cod_vac[,-c("TIPO")], by.x = "VACCINO_ID", by.y = "CODICE", all.x = T)

setDT(vac2)
vac2[, DESCRIZIONE := gsub(" ", "_", DESCRIZIONE)]
vac2[, DESCRIZIONE := gsub(",", "", DESCRIZIONE)]

vac2 = vac2 %>%
  pivot_wider(names_from = DESCRIZIONE, names_prefix = "VACC_", id_cols = "COD_SOGGETTO",
              values_from = N_VACC, values_fill = 0)

vac2$N_VACC = apply(vac2[,-c(1:2)],FUN = sum, MARGIN = 1 )

vac2 = merge(lapse,vac2, by="COD_SOGGETTO",all.x = T )
vac2 = vac2[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO","VACCINO_ID")]
vac2[is.na(vac2)] <- 0



#### ESENZIONI #####

ez = fread("N:/input/DATI/HTCOVID_ESENZIONI.csv",select = c(1,3,4))[COD_SOGGETTO %in% lapse$COD_SOGGETTO]
ez[, DT_INIZIO:= parse_sas_date(DT_INIZIO)]

ez = merge(ez,lapse, by="COD_SOGGETTO")
ez = ez[DT_INIZIO <= DATA_INIZIALE]
ez = ez[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]
ez = ez[,D :=1]
ez = ez[,ESE_ID := sub("\\..*$", "", ESE_ID)]

ez = distinct(ez[,c("COD_SOGGETTO","ESE_ID","D")])


ez2 = ez %>%
  group_by(COD_SOGGETTO) %>%
  summarise(
    # finanancialy expemption
    ES_FIN = fifelse(any(grepl("^E",ESE_ID)),1,0),
    # disability expemption
    ES_INV = fifelse(any(grepl("^I",ESE_ID)),1,0),
    # disease expempition
    ES_DIS = fifelse(any(grepl("^0",ESE_ID)),1,0))

ez = ez %>%
  mutate(ESE_ID = paste0("ES_",ESE_ID))

ez = ez %>%
  dcast(COD_SOGGETTO ~ ESE_ID,
        value.var = "D")

ez = merge(ez2, ez, by="COD_SOGGETTO")
rm(ez2)
gc(full=T)

ez = merge(lapse,ez, by="COD_SOGGETTO",all.x = T, )
ez = ez[,-c("DATA_INIZIALE","DATA_START_EVAL","YEAR","MEDICO")]
ez[is.na(ez)] <- 0

#### put togheter ####

data = merge(adi, sos, by = "COD_SOGGETTO") %>%
  merge(vac2, by = "COD_SOGGETTO") %>%
  merge(ez, by = "COD_SOGGETTO") %>%
  merge(med, by = "COD_SOGGETTO")


fwrite(data, fileout)

rm(list = ls())
invisible(gc(full = T))



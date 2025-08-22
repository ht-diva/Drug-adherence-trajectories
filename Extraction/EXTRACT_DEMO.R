#### create endpoints file ####
rm(list=ls())
gc(full = T)

setwd("N:/output/scripts/acorbetta/")

library(data.table)
setDTthreads(0)
library(dplyr)
library(stringr)
library(tidyr)


parse_sas_date = function(data) {
  return(as.Date(strptime(data, "%d%b%Y:%H:%M:%S.%OS")))
}



#### set up parameters ####
#filename wants a csv file without header with id and date of end of evaluation window (END)

#change to csv
filename = "N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv"

#output file
fileout = "N:/output/data/acorbetta/STATINS/ADHERENCE/demo_data_all.csv.gz"


# read cohort
sum_stat = fread(filename)

#remove this
colnames(sum_stat) = c("COD_SOGGETTO", "END")


# read anagrafe
tmp <- fread("N:/input/DATI/HTCOVID_ANAGRAFICA.csv",header = F, nrows = 1) %>% as.character()

pt =  tmp[!grepl("^STATO|^MMG|RES|DECESSO", tmp)]
anagrafe <- fread("N:/input/DATI/HTCOVID_ANAGRAFICA.csv", select = pt)[COD_SOGGETTO %in% sum_stat$COD_SOGGETTO]
gc(full=T)


# Melt the data into long format
long_data <- melt(
  anagrafe,
  id.vars = "COD_SOGGETTO",                                  # Columns to keep as is
  measure.vars = patterns("^ISTAT_DOM_", "^ATS_ASSISTENZA_"),     # Group columns by patterns
  variable.name = "YEAR",                          # Name for the long variable
  value.name = c("DOM", "ATS")                   # Names for values in long format
)

# Extract the year from column names (e.g., ISTAT_2012 -> 2012)
long_data[, YEAR := as.numeric(gsub(".*_", "", names(anagrafe)[grep("^ISTAT_DOM_", names(anagrafe))])[YEAR])]

sum_stat[,YEAR := year(END)]
long_data = merge(sum_stat, long_data, all.x = T)
summary(long_data)

an2 = anagrafe[,c("COD_SOGGETTO","SESSO","ANNO_NASCITA","NAZIONE_NASCITA")]

df = merge(an2,long_data)
df[,AGE := YEAR - ANNO_NASCITA]
df[,YEAR := NULL]
df[,END := NULL]

cou = fread("N:/output/data/ISTAT_COUNTRIES.tsv")[,c(1,2,4)]
colnames(cou) = c("REGION", "NAZIONE_NASCITA", "NATION")
cou[, NAZIONE_NASCITA := as.numeric(NAZIONE_NASCITA)]
df[,NAZIONE_NASCITA := as.numeric(NAZIONE_NASCITA)]
df[,NAZIONE_NASCITA :=  fifelse(NAZIONE_NASCITA %in% c("-","?","#") | is.na(NAZIONE_NASCITA), 0, NAZIONE_NASCITA)]
df = df %>%
  merge(cou, by = "NAZIONE_NASCITA", all.x = T)

df = df %>%
  mutate(
  REGION = fifelse(is.na(REGION), "Unknown", REGION ),
  NATION = fifelse(is.na(REGION), "Unknown", NATION ),
  CONT = relevel(factor(fcase(grepl("euro", tolower(REGION)) , "European",
                            grepl("ocean|unkown", tolower(REGION)) , "Other",
                            grepl("asia", tolower(REGION)) , "Asian",
                            grepl("africa", tolower(REGION)) , "African",
                            grepl("america", tolower(REGION)) , "American",
                            default = "Other")),ref = "European"))

#### get codifiche ####
map = read_sas('N:/input/DATI/htcovid_codifiche.sas7bdat') %>%
  filter(TIPO %in% c("ENTE", "CENTRO_VACCINALE", "VACCINI",
                     "POLIVALENTE", "FARMACI", "ATC"))

vm = (map[map$TIPO == "VACCINI",])
ats = map[grep("^ATS",map$DESCRIZIONE),] %>%
  rename(ATS_ASSISTENZA_2019 = CODICE,
         ATS = DESCRIZIONE) %>%
  select(-TIPO) %>%
  filter(ATS != "ATS CITTA' METROPOLITANA DI MILANO")

colnames(ats)= c("ATS","ATS_NAME")

df = df %>%
  merge(ats, by = "ATS",all.x = T)
table(df$ATS_NAME, useNA = "ifany")

gc(full = T)


#### define urban_rural
ur = fread("N:/output/data/Urban_rural.txt", col.names = c("DOM","URBAN"))
ur[,DOM := as.numeric(fifelse(DOM == "-", NA, DOM))]
df[,DOM := as.numeric(fifelse(DOM == "-", NA, DOM))]

df = merge(ungroup(df), ur, all.x = T, by = c("DOM") )


# write out
fwrite(df, fileout)

rm(list=ls())
gc(full = T)


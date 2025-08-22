rm(list=ls())
invisible(gc(full = T))

setwd("N:/output/scripts/acorbetta/")
library(data.table)
setDTthreads()
library (dplyr)
library (lubridate)
library(stringr)
library(dtplyr)
#.libPaths ("/shared-directory/sd-tools/apps/R/lib/")
library(feather)
library(haven)
sink(paste0("N:/output/data/acorbetta/STATINS/output_process_purchases",Sys.Date(),".txt"))
to = Sys.time()
Sys.setlocale("LC_TIME", "English")
#
#### DESCRIPTION ####
#this is a list of everything that will be performed on the dataset before calculating the adherence
#FILTERS:
# 1. remove ANJA patient
# 2. keep only patient with more (or equal) than 5 purchases
# 3. keep only trajectories with complete pack size information
# 4. aggregate 0 purchase lags to previous purchase
# 5. remove all patients with daily dosage different than 1
# 6. if change of drug, remove the extra pills based on the last purchase
# 



ti = Sys.time()
#  print (paste0 ('starting file ',i+1, " of ",10))
df = read_feather('N:/output/data/drug_extraction/2024-05-10_STATINS.feather', 
                  columns = c("COD_SOGGETTO","AIC", "ATC", "DT_EROGAZIONE", "QTA")) %>% 
    as.data.table()
colnames(df) =  c("PATIENT_ID","AIC","ATC_CODE","PURCH_DATE","N_PACKS")
df[, PURCH_DATE := as.Date(PURCH_DATE, format = "%d%b%Y")]


#### STEP CODINGS ####
# get pack size from AIC codes
CODES = read_sas("N:/input/DATI/htcovid_codifiche.sas7bdat")

#get AIC for statins
AIC = unique(df$AIC)
codiciAIC = CODES[CODES$CODICE %in% AIC,c(2,3)]
colnames(codiciAIC) = c("AIC", "DESCRIZIONE")

codiciAIC$PACK_SIZE_AIC =  as.numeric(str_extract_all(codiciAIC$DESCRIZIONE,"\\b(\\d+)\\s*(?=CPR|CPS)" ) )
codiciAIC$DOSAGE =  as.numeric(str_extract(codiciAIC$DESCRIZIONE,"\\d+(?=\\s*MG)" ) )

codiciAIC = codiciAIC %>% dplyr::select(-DESCRIZIONE) %>%
    as.data.table()

df = merge(df,codiciAIC, by =  "AIC", all.x = T)

codiciAIC = merge(codiciAIC,distinct(df[,c("AIC","ATC_CODE")]), by = "AIC", all.x = T)  

invisible(gc(full = T))

atc_names = fread ("N:/output/data/ATC_CODES.txt")
colnames (atc_names) = c("ATC_CODE", "ATC_NAME")

df = merge(df, atc_names, all.x = T, by = "ATC_CODE")
df[, STAT_NAME := sub (".*","",ATC_NAME)]
df[, STAT_NAME := sub ( ",","", ATC_NAME) ]


#### STEP 1  ####
print('starting step 1')



pairs_dt <- df[, .(X = pmin(ATC_NAME, as.character (DOSAGE)),
                   Y = pmax(ATC_NAME, as.character (DOSAGE) ) )]
pairs_dt <- as.data.table(table(pairs_dt))
pairs_dt <- pairs_dt[N != 0]

# Remove duplicate pairs
unique_pairs_dt <- distinct(pairs_dt)
unique_pairs_dt


fwrite(unique_pairs_dt, "N:/output/data/acorbetta/STATINS/ADHERENCE/statin_dosage_summary.csv")

invisible(gc(full =T))

#### STEP 2 ####

print('starting step 2')


df = df[order(PATIENT_ID, PURCH_DATE)][
  !is.na(PACK_SIZE_AIC), 
  if (all(!is.na(PACK_SIZE_AIC))) .SD, 
  by = PATIENT_ID
]

invisible(gc(full = T))


#### STEP 3 ####


print("starting step 3")


# Step 3.1: Compute number of pills

df[, N_PILLS := ifelse(N_PACKS == 0, PACK_SIZE_AIC, PACK_SIZE_AIC * N_PACKS)]

# Step 3.2: Merge observations within 7 days and create group_id
df[, group_id := cumsum(c(0, as.numeric(diff(PURCH_DATE)) > 7)), by = PATIENT_ID]

# Step 3.3: Summarise by PATIENT_ID, group_id, and ATC_CODE
df <- df[, .(N_PILLS_AG = sum(N_PILLS),
             PURCH_DATE = first(PURCH_DATE),
             DOSAGE = last(DOSAGE)),
         by = .(PATIENT_ID, group_id, ATC_CODE)]

# Step 3.4: Rename and arrange the data
df <- df[, .(PATIENT_ID, group_id, ATC = ATC_CODE, N_PILLS_AG, PURCH_DATE, DOSAGE)][order(PATIENT_ID, PURCH_DATE)]

# Step 3.5: Select and remove group_id
df <- df[, .(PATIENT_ID, ATC, N_PILLS_AG, PURCH_DATE,DOSAGE)]

# Clean up
print(paste0(uniqueN(df$PATIENT_ID), " indviduals at step 3"))
invisible(gc(full = T))


#### STEP FINAL ####
  
  #as.numeric(difftime(as.IDate("2012-01-10") , as.IDate("2012-11-30"), units = "days"))
  
  # cumulative adherence = tot pills / tot days
  # DAYS_PREV = number of days from purchase in T to purchase in T-1
  # DAYS_NEXT = number of days from purchase in T to purchase in T+1

# Set days between current and previous purchase, difference otherwise
df[, DAYS_PREV := fifelse(seq_len(.N) == 1, 0, 
                          as.numeric(difftime(PURCH_DATE, shift(PURCH_DATE), units = "days"))),
   by = PATIENT_ID]

# Set flags for change in medication and discontinuation
df[, CHANGE := fifelse(ATC == shift(ATC, fill = 0, type = "lag"), 0, 1),
   by = PATIENT_ID]
df[, DISCON := fifelse((DAYS_PREV > 365 + N_PILLS_AG) | seq_len(.N) == 1, 1, 0),
   by = PATIENT_ID]

# Adjust the previous discontinuation count to not count the days in the cumsum
df[, DAYS_PREV := fifelse(DISCON == 1, 0, DAYS_PREV),
   by = PATIENT_ID]

# Adjust the cumsum to not count the pills in the cumsum
df[, N_PILLS_ADJ := fifelse(shift(DISCON, fill = 0, type = "lead") == 1, 0,
                            fifelse(shift(CHANGE, fill = 0, type = "lead") == 0, N_PILLS_AG,
                                    fifelse(N_PILLS_AG < shift(DAYS_PREV, fill = 365, type = "lead"),
                                            N_PILLS_AG,
                                            shift(DAYS_PREV, fill = 0, type = "lead")))),
   by = PATIENT_ID]

# Compute cumulative pills adjusted for discontinuation and change of meds
df[, CUM_PILLS := fifelse(DISCON == 1, NA_real_,
                          cumsum(N_PILLS_ADJ) - N_PILLS_ADJ),
   by = PATIENT_ID]

# Compute cumulative days and adherence rate at each point
df[, CUM_DAYS := cumsum(DAYS_PREV),
   by = PATIENT_ID]

df[, CUM_MPR := fifelse(seq_len(.N) == 1 | DISCON == 1, NA_real_,
                        CUM_PILLS / CUM_DAYS),
   by = PATIENT_ID]

df[, POINT_MPR := fifelse(seq_len(.N) == 1 | DISCON == 1, NA_real_,
                          shift(N_PILLS_ADJ, fill = NA_real_, type = "lag") / DAYS_PREV),
   by = PATIENT_ID]

# Cap MPR values at 1 and flag high dose rates
df[, CUM_MPR_cap := fifelse(CUM_MPR > 1, 1, CUM_MPR),
   by = PATIENT_ID]
df[, POINT_MPR_cap := fifelse(POINT_MPR > 1, 1, POINT_MPR),
   by = PATIENT_ID]
df[, DOS2 := fifelse(median(CUM_MPR, na.rm = T) > 1.2, 1, 0),
   by = PATIENT_ID]


# Adjust previous dosage by shifting forward and filling NA with previous values
df[, DOSAGE_PREV := shift(DOSAGE, fill = NA_real_, type = "lag"), by = PATIENT_ID]




med = df[, .(median_adh = median(POINT_MPR, na.rm = T),
             tot = sum(!is.na(POINT_MPR))),
         by = .(ATC,DOSAGE_PREV)]

fwrite(med,"N:/output/data/acorbetta/STATINS/ADHERENCE/median_dosages_statins.csv")


m1 = merge(data.frame (ATC = unique(df$ATC)), atc_names, all.x = T, by.x = "ATC", by.y = "ATC_CODE")
conversion = merge(m1, unique_pairs_dt, by.x = "ATC_NAME", by.y = "Y") %>%
  arrange (ATC, X)



conversion = data.frame (conversion,
                          INTENSITY = c(1,2,2, # simvastatin
                                        1,2, #lovastatin
                                        1,1,2, #pravastatin
                                        1,1,2, # fluvastatin
                                        2,2,2,3,3,3, #atrovastatin
                                        2,2,3,3,3,2,#Rosuvastatin 
                                        2,3,3, #Simvastatin and ezetimibe
                                        2, # Simvastatin and fenofibrate
                                        3, #Atorvastatin and ezetimibe
                                        3,3,3,3, #Rosuvastatin and ezetimibe
                                        2,3,2, #Rosuvastatin and ezetimibe
                                        2,3, #Rosuvastatin and amlodipine
                                        2,1, #Atorvastatin, amlodipine and perindopril
                                        2,2 #Rosuvastatin and ramipril
                                        )
                          )

conversion$N = NULL                          
colnames (conversion) = c("ATC_NAME", "ATC", "DOSAGE", "STAT_CLASS" )
setDT(conversion)
conversion[, STAT_CLASS := fcase(STAT_CLASS == 1, "Low",
                                STAT_CLASS == 2, "Mid",
                                STAT_CLASS == 3, "High")]

conversion[, DOSAGE := as.numeric(DOSAGE)]

# Merge conversion data back to the main dataframe
df <- merge(df, conversion, all.x = TRUE, by = c("ATC", "DOSAGE"))

# Adjust STAT_CLASS by shifting forward and filling NA with previous values
df[, STAT_CLASS_PREV := shift(STAT_CLASS, fill = NA, type = "lag")]

# Order the data by PATIENT_ID and PURCH_DATE
df <- df[order(PATIENT_ID, PURCH_DATE)]


print("writing...")
filename = paste0("ADH_STATINS_LONG_ALL_MPR")
fwrite(df, paste0("N:/output/data/acorbetta/STATINS/ADHERENCE/",filename,".csv.gz"))
rm(df)
invisible(gc(full = T))

print(paste0("time taken: ", Sys.time() - ti))

print(paste0("time total: ", Sys.time() - t0))


sink()

rm(list=ls())
invisible(gc(full = T))


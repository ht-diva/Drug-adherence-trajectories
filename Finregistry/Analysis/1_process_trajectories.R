rm(list=ls())
invisible(gc(full = T))


library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
#.libPaths("/shared-directory/sd-tools/apps/R/lib/")
#library(feather)

sink("/media/volume/mferro/data/output_process_trajectories.txt")
t0 = Sys.time()
#______________________
####  DESCRIPTION  ####

#this is a list of everything that will be performed on the dataset before calculating the adherence

#FILTERS:

# remove ANJA patient  
# keep only patient with more (or equal) than 5 purchases 
# keep only trajectories with complete pack size information
# aggregate 0 purchase lags to previous purchase
# remove all patients with daily dosage different than 1 
# if change of drug, remove the extra pills based on the last purchase

# what is missing?
# handle big gaps -HANDLED
#text_to_add = readRDS('/media/volume/users/acorbett/migrate/data1/joined_dosage.rds')%>%
#  select(PATIENT_ID, DOSAGE_INSTRUCTION)
#gc(verbose=F)  
#extract information about the dosage and keep only those with 1 tablet a day
#NB: every patient need to have one dosage instruction, pick the first one

#TECNIQUE TO USE: regex 
#str_to_grepl = "1|(yksi)"
#str_to_grepl2 = "(puoli)|½|1/2|1-2|0.5"

#text_to_add = text_to_add %>% 
#  group_by(PATIENT_ID) %>%
#  summarise( DOSE_TXT = first(na.omit(DOSAGE_INSTRUCTION)) ) %>% 
  #extract dosage
#  filter( grepl(str_to_grepl,DOSE_TXT), !grepl(str_to_grepl2,DOSE_TXT)) %>%
#  ungroup()
#gc(verbose=F)



ti = Sys.time()
print(paste0('starting'))
df = fread('/media/volume/mferro/data/DRUGS/STATINS.csv.gz', 
           select = c(1,4,5,7,8), 
           col.names = c("PATIENT_ID","PURCH_DATE","ATC_CODE","VNRO","NPACK"))[,
      PURCH_DATE := ymd(PURCH_DATE)]

df = df[ATC_CODE != "C10AA06"]
df = df[year(PURCH_DATE) >= 1998]


invisible(gc(full = T))

#### STEP CODINGS ####
vnr = read.csv("/media/volume/mferro/data/vnr_dictionary_statins.csv")
setDT(vnr)
vnr = vnr[,c("vnr","ATC","valmiste","vahvuus","pkoko_num")]

vnr2 = read.csv("/media/volume/mferro/data/vnr_dictionary_statins.csv")
setDT(vnr2)
vnr2 = vnr2[,c("vnr","ATC","valmiste","vahvuus","pkoko_num")]
vnr2 = vnr2[vnr != 9101]
vnr2[,vnr := paste0(vnr,0)]

vnr3 = read.csv("/media/volume/mferro/data/vnr_dictionary_statins.csv")
setDT(vnr3)
vnr3 = vnr3[,c("vnr","ATC","valmiste","vahvuus","pkoko_num")]
vnr3 = vnr3[vnr != 9101]
vnr3[,vnr := paste0(vnr,"00")]

vnr4 = distinct(rbind(vnr,vnr2,vnr3))
colnames(vnr4) = c("VNRO","ATC_CODE","STAT_NAME","DOSAGE_STR","PACK_SIZE")


#tmp2 = merge(tmp,vnr4, all.x = T, by.x = "Var1", by.y = "VNRO")
#tofill = tmp2[is.na(tmp2$STAT_NAME),]

#vnr_all[grepl(paste0(tofill$Var1,collapse = "|"),vnr) ,]

df$VNRO = as.character(df$VNRO)

df= merge(df,vnr4,all.x = T, by = "VNRO")
df$ATC_CODE.y = NULL
df = df %>% rename(ATC_CODE = ATC_CODE.x)

atc_names = fread ("/media/volume/mferro/data/ATC_CODES.csv")
colnames (atc_names) = c("ATC_CODE", "ATC_NAME")

df = merge(df, atc_names, all.x = T, by = "ATC_CODE")
df[, STAT_NAME := sub (".*","",ATC_NAME)]
df[, STAT_NAME := sub ( ",","", ATC_NAME) ]


####   STEP 1    ####
print('starting step 1')


# Define a function to extract numbers from a string
extract_numbers <- function(x) {
   as.numeric(sub("(^\\d+).*", "\\1", x))
}


# Read in another data table presumably with statistics
df[, DOSAGE := as.numeric(extract_numbers(DOSAGE_STR))]

# Adjust DOSAGE based on a condition
df[, DOSAGE := ifelse(DOSAGE > 99, DOSAGE / 10, DOSAGE)]

# Create combinations of statin and dosage
pairs_dt <- df[, .(X = pmin(ATC_NAME, as.character(DOSAGE)), 
                     Y = pmax(ATC_NAME, as.character(DOSAGE)))]
pairs_dt <- as.data.table(table(pairs_dt))
pairs_dt <- pairs_dt[N != 0]

# Remove duplicate pairs
unique_pairs_dt <- unique(pairs_dt)

# Print the resulting unique pairs
print(unique_pairs_dt)
fwrite(unique_pairs_dt, "/media/volume/mferro/data/ADHERENCE/statin_dosage_summary.csv")

invisible(gc(full =T))


 ####   STEP 2    ####
print('starting step 2')

# remove also trajectories with incomplete information about the pack size

df = df[order(PATIENT_ID,PURCH_DATE)][
  !is.na(PACK_SIZE),
  if(all(!is.na(PACK_SIZE))) .SD,
  by = PATIENT_ID
]



invisible(gc(full =T))

#_____________________
####   STEP 3    ####
print('starting step 3')
  


#for (i in c(1:length(ids))) {
  #print(paste0("Starting step ",i, " of ", length(ids) ))
# merge same day purchases into one
# Step 3.1: Compute number of pills

df[, N_PILLS := ifelse(NPACK == 0, PACK_SIZE, PACK_SIZE * NPACK)]

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

# Compute median adherence, total non-NA adherence counts by ATC and DOSAGE_PREV
med <- df[, .(
   median_adh = median(POINT_MPR, na.rm = TRUE),
   tot = sum(!is.na(POINT_MPR))
), by = .(ATC, DOSAGE_PREV)] %>%
   na.omit()

fwrite(med, "/media/volume/mferro/data/ADHERENCE/empirical_indications.csv")



# Set weight based on median adherence greater than 1.2
med[, weight := fifelse(median_adh > 1.2, 2, 1)]

# Create a conversion table for medication intensity classification
conversion <- merge(merge(data.frame(ATC = unique(df$ATC)),
                          atc_names, all.x = T, by.x = "ATC", by.y = "ATC_CODE"),
                    unique_pairs_dt[,-c('N')], 
                    by.x = "ATC_NAME", by.y = "Y")

int = fread("/media/volume/mferro/data/ADHERENCE/conversion_statins.csv")

conversion = merge(conversion, int, all.x=T)



colnames(conversion) <- c("ATC_NAME", "ATC", "DOSAGE", "STAT_CLASS")
conversion$DOSAGE <- as.numeric(conversion$DOSAGE)

# Merge conversion data back to the main dataframe
df <- merge(df, conversion, all.x = TRUE, by = c("ATC", "DOSAGE"))

# Adjust STAT_CLASS by shifting forward and filling NA with previous values
df[, STAT_CLASS_PREV := shift(STAT_CLASS, fill = NA, type = "lag")]

# Order the data by PATIENT_ID and PURCH_DATE
df <- df[order(PATIENT_ID, PURCH_DATE)]




print("writing...")
filename = paste0("ADH_STATINS_LONG_ALL_MPR")
fwrite(df, paste0("/media/volume/mferro/data/ADHERENCE/STATINS/",filename,".csv.gz"))
rm(df)
invisible(gc(full = T))
  
print(paste0("time taken: ", Sys.time() - ti))

print(paste0("time total: ", Sys.time() - t0))


sink()

rm(list=ls())
invisible(gc(full = T))



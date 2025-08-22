rm (list=ls())
invisible(gc(full = T))

setwd("N:/output/scripts/acorbetta/")
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
setDTthreads(0)

df = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_LONG_ALL_MPR.csv.gz")

n_start = uniqueN(df$PATIENT_ID)

df = df %>%
  group_by(PATIENT_ID) %>%
  filter(row_number() != 1) %>%
  filter(!any(DISCON==1)) %>%
  ungroup()

setDT(df)
df <- df[
  PURCH_DATE >= as.IDate("2015-01-01") & PURCH_DATE < as.IDate("2023-06-01")
]

n_after = uniqueN(df$PATIENT_ID)

df = df[!is.na(POINT_MPR_cap)]
df = df[,DATE_FIRST := min(PURCH_DATE), by = PATIENT_ID]
df = df[, DAYS := as.numeric(difftime(PURCH_DATE, DATE_FIRST, units = "days" ))]
df = df[DAYS <= 365.25*6]
df[, EW := max(DAYS), by = PATIENT_ID]
df = df[EW > 365.25*5]

dim(df)
uniqueN(df$PATIENT_ID)

df = df %>%
  group_by(PATIENT_ID) %>%
  filter(n() > 5)

n_m5 = uniqueN(df$PATIENT_ID)

fwrite(data.frame(N_START = n_start,
                  N_DISCON = n_after,
                  N_5OBS = n_m5
), "N:/output/data/acorbetta/STATINS/ADHERENCE/selection_cohort.cvs" )



fwrite(df, "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_LONG_FILTERED.csv.gz")


invisible(gc(full = T))


dft = df %>%
  group_by (PATIENT_ID) %>%
  summarise(n_obs = n(),
            first = first(DATE_FIRST) ,
            last = last(PURCH_DATE),
            mean = mean (POINT_MPR_cap, na.rm = T),
            median = median (POINT_MPR_cap,na.rm=T) ,
            sddev = sd (POINT_MPR_cap, na.rm = T),
            time = as.numeric(difftime(last, first, units = "weeks")) ,
            change = sum (CHANGE) ,
            lag = mean (DAYS_PREV) ,
            sd_lag = sd (DAYS_PREV) ,
            seq_dos = paste0(rle(STAT_CLASS[!is.na(STAT_CLASS_PREV)])$values, collapse = "="),
            seqlen = length (rle(STAT_CLASS[!is.na(STAT_CLASS_PREV)])$values))


fwrite(dft, "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_SUMMARY.csv.gz" )

#### ADHERENCE PLOTS ####
library(ggplot2)
library(ggpubr)

start_plot = ggplot(dft, aes(x = first)) +
  geom_density(fill = "dodgerblue") +
  scale_x_date(date_breaks = "1 year") +
  xlab("Date initiation") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

mean_plot = ggplot(dft, aes(x = mean)) +
  geom_density(fill = "dodgerblue2", bounds = c(0,1)) +
  scale_x_continuous() +
  xlab("Adherence (MPR)") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

lag_plot = ggplot(dft, aes(x = lag)) +
  geom_density(fill = "dodgerblue3", bounds = c(0,+Inf), trim = T) +
  scale_x_continuous() +
  xlab("Lag days") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

obs_plot = ggplot(dft, aes(x = n_obs *1.0 )) +
  geom_density(fill = "dodgerblue4", trim = T, bounds = c(0,+Inf), bw = 2) +
  scale_x_continuous() +
  xlab("Number of purchases") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

library(ggpubr)
densities = ggarrange(start_plot, mean_plot, lag_plot,obs_plot)

trace(grDevices:::svg, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/summary_densities.svg", width = 10, height = 10)
densities
dev.off()


# add age and sex distributions
minimal = fread("N:/output/data/minimum_data.csv", select = c("COD_SOGGETTO", "SESSO", "ANNO_NASCITA")) 
minimal = minimal[COD_SOGGETTO %in% dft$PATIENT_ID]

dft = merge(dft, minimal, by.x = "PATIENT_ID", by.y = "COD_SOGGETTO", all.x = T)
setDT(dft)
dft[, AGE := as.numeric(difftime(first, as.Date(paste0(ANNO_NASCITA,"/06/01")), units = "days"))/365.25]


fwrite(dft, "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_SUMMARY.csv.gz" )


library(tidyverse)

summary_stats <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    Q3 = quantile(x, 0.75, na.rm = TRUE)
  )
}



# Apply the function to each column in the data frame
sumstat <- dft %>%
  dplyr::select(-PATIENT_ID, -first, -last, -median, -change, -seq_dos,-time, -seqlen,
                - SESSO, -ANNO_NASCITA) %>%
  reframe(across(everything(), summary_stats, .names = "{.col}_{.fn}")) %>%
  t() %>%
  as.data.frame()

colnames(sumstat) = c("Mean","SD", "1q","Median","3q")
rownames(sumstat) = c("N_OBS","MEAN_ADH", "SD_ADH","MEAN_LAG","SD_LAG","AGE")

fwrite(sumstat, "N:/output/data/acorbetta/STATINS/ADHERENCE/trajectories_descriptives_numeric.csv" , row.names = T)

# summary sequences 

# changes

t_seq = table(dft$seq_dos)
t_seq = t_seq[order(t_seq, decreasing = T)][1:10] %>% data.frame()
colnames(t_seq)  = c("SEQUENCE", "N")
fwrite(t_seq, "N:/output/data/acorbetta/STATINS/ADHERENCE/sequence_table.csv" )
t_len = table(fifelse(dft$seqlen>6,6-1,dft$seqlen-1)) %>% data.frame()
colnames(t_len)  = c("CHANGES", "N")
fwrite(t_len, "N:/output/data/acorbetta/STATINS/ADHERENCE/changes_table.csv" )
t_sex = table(dft$SESSO) %>% data.frame()
t_sex$desc = c("MALE","FEMALE")
colnames(t_sex)  = c("SEX", "N","DESC")
fwrite(t_sex, "N:/output/data/acorbetta/STATINS/ADHERENCE/sex_table.csv" )


age_plot = ggplot(dft[AGE>18], aes(x = AGE, fill = as.factor(SESSO))) +
  geom_density(bw = 2, alpha = 0.8) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = c("Male","Female"), name = "Sex")+
  scale_x_continuous(breaks = seq(20,100,10)) +
  xlab("Age") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/age_sex_densities.svg", width = 12, height = 6)
age_plot
dev.off()


rm(list=ls())
invisible(gc(full = T))



rm (list=ls())
invisible(gc(full = T))

setwd("N:/output/scripts/acorbetta/")
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(dtplyr)
library(haven)
#.libPaths ("/shared-directory/sd-tools/apps/R/lib/")
library(feather)

sink(paste0("N:/output/data/acorbetta/STATINS/output_select_summarise_covid",Sys.Date(),".txt"))
to = Sys.time()

df = read_feather("N:/output/data/acorbetta/STATINS/ADH_STATINS_LONG_ALL.feather")


#fwrite(df, "N:/output/data/acorbetta/ANTIHYPCOL/ADH_ANTIHYPCOL_LONG_ALL_test.feather")
#
#sum = df %>% 
#  group_by(PATIENT_ID) %>%
#  filter(n() > 5) %>%
#  summarise(DISCON = min(DISCON), DOS2 = min(DOS2))

#t1 = Sys.time()
df = df %>%
  lazy_dt() %>%
  filter(PURCH_DATE >= ymd("2014-01-01") & PURCH_DATE <= ymd("2023-06-01") ) %>%
  group_by(PATIENT_ID) %>%
  # filter discontinuing, dosage = 2*day, few observations
  filter(n() > 5, !any(DISCON == 1), !any(DOS2 == 1)) %>%
  # remove the first obs for MPR calculation
  filter(!row_number() == 1) %>%
  # remove WINDOWS excessive then 6 years from 2020
  #filter( as.numeric(difftime(PURCH_DATE, first(PURCH_DATE), units = "d")) < 365 * 6) %>%
  mutate(
    # compute the number of days
    DAYS = as.numeric(difftime(PURCH_DATE, first(PURCH_DATE), units = "d")),
    # compute total number of days
    EW = last(DAYS)) %>%
  # remove people with short EW 
  filter(EW >365*5, last(PURCH_DATE) > ymd("2021-03-01"), first(PURCH_DATE) <  ymd("2018-03-01") ) %>%
  as.data.table()

write_feather(df, "N:/output/data/acorbetta/STATINS/ADH_STATINS_LONG_COVID_FILTERED_2020.feather")
print("file written")
print(paste("time taken for filtering: ", Sys.time() - to))


df = read_feather("N:/output/data/acorbetta/STATINS/ADH_STATINS_LONG_COVID_FILTERED_2020.feather")
invisible(gc(full = T))
df = df %>% 
  group_by(PATIENT_ID) %>%
  summarise(MED_LAG = median(DAYS_PREV), 
            SD_LAG = sd(DAYS_PREV),
            CUM_MPR = median(CUM_MPR_cap),
            POINT_MPR = median(POINT_MPR_cap),
            SD_POINT_MPR = sd(POINT_MPR_cap),
            SD_CUM_MPR = sd(CUM_MPR_cap),
            CHANGE = max(CHANGE),
            NdiffATC = length(unique(ATC)),
            FIRST_PURCH = first(PURCH_DATE),
            LAST_PURCH  = last(PURCH_DATE),
            TOT_PURCH = n(),
            median_PILLS = median(N_PILLS_AG),
            tot_PILLS = sum(N_PILLS_AG)) %>%
  rename(COD_SOGGETTO = PATIENT_ID)

invisible(gc(full = T)) 


mp = fread("N:/pgmsas/minimum_data.csv")  %>%
  filter(COD_SOGGETTO %in% df$COD_SOGGETTO) %>%
  select(COD_SOGGETTO, SESSO, ANNO_NASCITA, NAZIONE_NASCITA)

invisible(gc(full = T))

df = df %>%
  left_join(mp) %>%
  mutate(AGE = year(FIRST_PURCH) - ANNO_NASCITA)
df$FIRST_PURCH = ymd(df$FIRST_PURCH)

rm(mp)
invisible(gc(full = T))

write_feather(df, "N:/output/data/acorbetta/STATINS/ADH_STATINS_SUMMARY_COVID_1DOS_NODIS_ACROSS.feather")
#df = read_feather("N:/output/data/acorbetta/ANTIHYPCOL/ADH_ANTIHYPCOL_SUMMARY_1DOS_NODIS.feather") %>% na.omit()
print(paste0("TOTAL NUMBER OF PATIENTS W/ 5 PURCH, NO DISCON, NO 2X DOSE, AND 5Y EW IS: ", dim(df)[1]))

#library(ggplot2)
#library(gridExtra)

#lag1 = ggplot(df, aes(x = MED_LAG, y = MPR, fill = SESSO, color =SESSO)) +
#  ylim(c(0.4,1))+
#  geom_density_2d(size = 1.1, alpha = 0.7)

#lag2 = ggplot(df, aes(x = SD_LAG, y = MPR, fill = SESSO, color =SESSO)) +
#  ylim(c(0.4,1))+
#  geom_density_2d(size = 1.1, alpha = 0.7)

#age = ggplot(df, aes(y = MPR, x = AGE, color =SESSO, fill = SESSO)) +
#  ylim(c(0.4,1))+
#  geom_density_2d(size = 1.1, alpha = 0.7)

#MPR = ggplot(df, aes(y = MPR, x = SD_MPR, color =SESSO, fill = SESSO)) +
#  ylim(c(0.4,1))+
#  geom_density_2d(size = 1.1, alpha = 0.7)

#change = ggplot(df, aes(x = MPR, fill = factor(CHANGE))) +
#  geom_density(alpha = 0.7) +
#  scale_fill_manual(values = c("blue","yellow"), name = "CHANGE") +
#  facet_wrap(vars(SESSO))

#NAZ = ggplot(df, aes(x = MPR, fill = factor(NAZIONE_NASCITA == 100))) +
#  geom_density(alpha = 0.7) +
#  labs(fill = 'Born in ITA') +
#  scale_fill_manual(values = c("purple","cyan"), labels = c("0","1")) +
#  facet_wrap(vars(SESSO))


#rm(df)
#invisible(gc())
#p_no_legend = list(lag1 + theme(legend.position = "none"),
#                   lag2 + theme(legend.position = "none"),
#                   age + theme(legend.position = "none"),
#                   MPR + theme(legend.position = "none"))

#legend = cowplot::get_legend(p_no_legend[[1]] +
#                             theme(legend.position = "bottom") +
#                             labs(color = "SEX"))


#p_grid = cowplot::plot_grid(plotlist = p_no_legend, ncol = 2)

#biv = cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = c(1,0.1))
#png("N:/output/data/acorbetta/ANTIHYPCOL/ADH_ANTIHYPCOL_PLOTBIV_1DOS_NODIS.png")
#biv
#dev.off()
#print("BIV plot written")

#cat = cowplot::plot_grid(change + theme(legend.position = "bottom"),
#                   NAZ + theme(legend.position = "bottom"),
#                   ncol = 1)

#png("N:/output/data/acorbetta/ANTIHYPCOL/ADH_ANTIHYPCOL_PLOTCAT_1DOS_NODIS.png",)
#cat
#dev.off()
#print("BIV plot written")

print("end")
print(paste("time taken: ", Sys.time() - to))
sink()

rm(list=ls())
invisible(gc(full=T))
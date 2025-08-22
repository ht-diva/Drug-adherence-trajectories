rm (list=ls())
invisible(gc(full=T))

setwd("N:/output/scripts/acorbetta/")
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(dtplyr)
library(haven)
#.libPaths ("/shared-directory/sd-tools/apps/R/lib/")
library(feather)
library(ggplot2)
library(purrr)

df = read_feather("N:/output/data/acorbetta/STATINS/ADH_STATINS_SUMMARY_1DOS_NODIS.feather")

mp = fread("N:/output/data/minimum_data.csv", select = c("COD_SOGGETTO", "SESSO", "ANNO_NASCITA", "NAZIONE_NASCITA")) 

df = df %>%
  select(-ANNO_NASCITA, - SESSO, - NAZIONE_NASCITA, -AGE,)

df = df %>%
  left_join(mp, by = "COD_SOGGETTO")

df = df %>%
  mutate(AGE = year(FIRST_PURCH) - ANNO_NASCITA,
         ITA = fifelse(NAZIONE_NASCITA == "100",1,0))
  
num = df %>%
  select_if(is.numeric)

a = apply(num,2,median,na.rm =T) %>% as.data.frame() %>% 
  cbind(apply(num,2,mean,na.rm =T) %>% as.data.frame()) %>%
  cbind(apply(num,2,sd,na.rm =T) %>% as.data.frame()) %>%
  round(3)
a

fac = df %>%
  select(ITA,SESSO)
b = apply(fac, 2, table) / nrow(df)
rownames(b) = c("0","1")
b

### plot age-pdc by nationality ###
x11()
ggplot(df, aes(x = POINT_MPR, y = AGE, group=factor(ITA), color = factor(ITA))) +
  stat_density_2d() +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_bw()

ggplot(df, aes(x = POINT_MPR, y = AGE, group=factor(SESSO), color = factor(SESSO))) +
  stat_density_2d_() +
  theme_bw()

ggplot(df, aes(POINT_MPR, group=factor(CHANGE), color = factor(CHANGE))) +
  stat_density()

ggplot(df, aes(POINT_MPR)) +
  stat_density()

ggplot(na.omit(df), aes(POINT_MPR, group=factor(CHANGE), color = factor(CHANGE), fill = factor(CHANGE))) +
  stat_density()

ggplot(na.omit(df),aes(POINT_MPR, group=factor(SESSO), color = factor(SESSO), fill = factor(SESSO))) +
  stat_density()


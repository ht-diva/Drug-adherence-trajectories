#### create socioeco file ####
rm(list=ls())
gc(full = T)

setwd("N:/output/scripts/acorbetta/")


library(feather)
library(data.table)
setDTthreads(0)
library(dplyr)
library(stringr)
library(tidyr)

fe = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/matrices_data_all.csv.gz")



med = fe %>%
  select(SESSO, ANNO_NASCITA, ATS, ANNO_LAUREA, ANNO_SPECIALIZZAZIONE, DISTRETTO, NUMERO_ASSISTITI, TIPO_PRESCRITTORE)

fe2 = fe %>%
  select(-SESSO, -ANNO_NASCITA, -ATS, -ANNO_LAUREA, -ANNO_SPECIALIZZAZIONE, -DISTRETTO, -NUMERO_ASSISTITI, -TIPO_PRESCRITTORE)

ch = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv", 
           col.names= c("COD_SOGGETTO","x"))[,1]

fe2 = merge(ch,fe2, all.x = T)
fe2[is.na(fe2)] <- 0

ID = fe2$COD_SOGGETTO

fwrite(list(colnames(fe)[-1]),"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_socioeco_all.csv")

# remove low frequancy diseases
fe_wide = fe2 %>%
  select(-COD_SOGGETTO) %>%
  select_if(function(x){sum(x)/nrow(fe) >= 0.005}) %>%
  as.data.frame()

#isolate total vaccines
nvacc = fe_wide$N_VACC
# mutate vaccines

fe_wide = fe_wide %>%
  select(-ADI_ACCESSI, - SOS_TIME, -N_VACC) %>%
  mutate(VACC_DIFTERITE = fifelse(VACC_DIFTERITE>0,1,0),
         VACC_INFLUENZA = fifelse(VACC_INFLUENZA>0,1,0),
         VACC_TETANO = fifelse(VACC_TETANO>0,1,0),
         VACC_M.I._PNEUMOCOCCO = fifelse(VACC_M.I._PNEUMOCOCCO>0,1,0))

fe_p = fe_wide %>%
  summarise_all(.fun = list(P = function(x){sum(x)/nrow(fe)} )) %>% 
  t()
rownames(fe_p) = str_replace(rownames(fe_p),"_P","")
fe_p = fe_p %>%
  data.frame() %>%
  mutate(VAR = rownames(fe_p))
rownames(fe_p) = NULL
colnames(fe_p) = c("P","VAR")

fwrite(fe_p,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_socioeco_prevalence_selection_all.csv")


#### phis dist ####
library(psych)
#phis = data.frame(NAME1 = c(1), NAME2 = c(1), PHI = c(1))
phis_dist = matrix(data = 0,ncol = ncol(fe_wide), nrow = ncol(fe_wide))

colnames(phis_dist) = colnames(fe_wide)
rownames(phis_dist) = colnames(fe_wide)

for (i in 1:(ncol(fe_wide)-1)) {
  print(i)
  for (j in (i+1):ncol(fe_wide)) {
    print(j)
    mcc = phi(table(as.matrix(fe_wide[,i]), as.matrix(fe_wide[,j])), digits = 5)
    phis_dist[i,j] = mcc
    phis_dist[j,i] = mcc
  }
}
diag(phis_dist) = 1

fwrite(phis_dist,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_socioeco_distance_matrix_all.csv")

#### cluster diag ####

one = function(x){(1 - abs(x))}
phis_as_dist = phis_dist  %>%
  as.data.frame() %>%
  mutate_all(one) %>% 
  as.dist()

library(cluster)
C1 = hclust(phis_as_dist, method = "complete")
C2 = hclust(phis_as_dist, method = "single")
C3 = hclust(phis_as_dist, method = "average")

par(mfrow=c(1,3))
plot(C1,main="Complete Linkage", xlab="", sub="", cex=.9)
plot(C2, main="Single Linkage", xlab="", sub="", cex=.9)
plot((C3),main="Average Linkage", xlab="", sub="", cex=.9)

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/socioeco_trees_all.svg",
    width = 12, height = 6)
par(mfrow=c(1,3))
plot(C1,main="Complete Linkage", xlab="", sub="", cex=.9)
plot(C2, main="Single Linkage", xlab="", sub="", cex=.9)
plot((C3),main="Average Linkage", xlab="", sub="", cex=.9)
dev.off()




P3 = cutree(C1,h =0.4)
table(P3)
library(ggplot2)
library(ggdendro)

#your own labels (now rownames) are supplied in geom_text() and label=label
dis_plot = ggplot() + 
  geom_segment(data=segment(dendro_data(C1)), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendro_data(C1)), aes(x=x, y=y, label=label, hjust=1, angle = 90), size=3) +
  #  coord_flip() + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), expand=c(0.7, 0)) + 
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), expand=c(0.01, 0)) +
  ylab("Distance") +
  geom_hline(yintercept = 0.4, col = "red") +
  theme(axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(angle = 90),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/socioeco_selection_all.svg", width = 12, height = 6)
dis_plot
dev.off()

#### remove duplicates or very similar endpoints ####

fe_ps = fe_p %>%
  mutate(G = P3) %>%
  group_by(G) %>%
  filter(P == max(P))

#gp data ha too much missingness,  cannot be used
prop.table(table((med$SESSO)))





fe_wide_s = ID  %>% cbind(fe_wide[,names(fe_wide) %in% fe_ps$VAR]) 
colnames(fe_wide_s)[1] = "COD_SOGGETTO"
fwrite(fe_ps,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_socioeco_final_selection_all.csv")

#### add demo ####
demo = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/demo_data_all.csv.gz")



demo[, ATS_NAME := fifelse(ATS_NAME == "", NA, ATS_NAME)]
demo[, URBAN := fifelse(URBAN == 1, 1, 0)]


demo = demo[,c("COD_SOGGETTO","SESSO","ANNO_NASCITA","AGE","CONT","ATS_NAME","URBAN")] %>%
  rename(ATS = ATS_NAME,
         CONTINENT = CONT,
         SEX = SESSO,
         YOB = ANNO_NASCITA)


fwrite(list(colnames(demo)[-1]),"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_demo_final_selection_all.csv")




fe_wide_s = fe_wide_s %>%
  merge(demo, by = "COD_SOGGETTO")

fe_wide_s$VACC_N = nvacc

fwrite(fe_wide_s,"N:/output/data/acorbetta/STATINS/ADHERENCE/socioeco_covariates_all.csv")


rm(list=ls())
gc(full = T)
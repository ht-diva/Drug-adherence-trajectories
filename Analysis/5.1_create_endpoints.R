#### create endpoints file ####
rm(list=ls())
gc(full = T)

setwd("N:/output/scripts/acorbetta/")


library(feather)
library(data.table)
setDTthreads(0)
library(dplyr)
library(stringr)
library(tidyr)

fe = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/diseases_data_all.csv")
ch = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv", 
           col.names= c("COD_SOGGETTO","x"))[,1]

fe = merge(ch,fe, all.x = T)
fe[is.na(fe)] <- 0

ID = fe$COD_SOGGETTO

fwrite(list(colnames(fe)[-1]),"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_health_all.csv")

# remove low frequancy diseases
fe_wide = fe %>%
  select(-COD_SOGGETTO) %>%
  select_if(function(x){sum(x)/length(x) >= 0.005}) %>%
  as.data.frame()


fe_p = fe_wide %>%
  summarise_all(.fun = list(P = function(x){sum(x)/nrow(fe)} )) %>% 
  t()

rownames(fe_p) = str_replace(rownames(fe_p),"_P","")
fe_p = fe_p %>%
  data.frame() %>%
  mutate(VAR = rownames(fe_p))
rownames(fe_p) = NULL
colnames(fe_p) = c("P","VAR")

fwrite(fe_p,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_helath_prevalence_selection_all.csv")


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

fwrite(phis_dist,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_health_distance_matrix_all.csv")

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

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/health_trees_all.svg", width = 12, height = 6)
par(mfrow=c(1,3))
plot(C1,main="Complete Linkage", xlab="", sub="", cex=.9)
plot(C2, main="Single Linkage", xlab="", sub="", cex=.9)
plot((C3),main="Average Linkage", xlab="", sub="", cex=.9)
dev.off()




P3 = cutree(C3,h =0.4)
table(P3)
library(ggplot2)
library(ggdendro)

#your own labels (now rownames) are supplied in geom_text() and label=label
dis_plot = ggplot() + 
  geom_segment(data=segment(dendro_data(C3)), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendro_data(C3)), aes(x=x, y=y, label=label, hjust=1, angle = 90), size=3) +
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

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/diseases_selection_all.svg",
    width = 12, height = 6)
dis_plot
dev.off()


#### remove duplicates or very similar endpoints ####

fe_ps = fe_p %>%
  mutate(G = P3) %>%
  group_by(G) %>%
  filter(P == max(P))
MCS = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/diseases_data_all.csv.gz", select = c("COD_SOGGETTO","MCS"))

fe_wide_s = ID %>% cbind(fe_wide[,names(fe_wide) %in% fe_ps$VAR]) 
colnames(fe_wide_s)[1] = "COD_SOGGETTO"
fe_wide_s = fe_wide_s %>% merge(MCS, by = "COD_SOGGETTO")


fwrite(fe_ps,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/covariates_helath_final_selection_all.csv")
fwrite(fe_wide_s,"N:/output/data/acorbetta/STATINS/ADHERENCE/diseases_covariates_all.csv")


rm(list=ls())
gc(full = T)
#### create endpoints file ####
rm(list=ls())
gc(full = T)

setwd("/media/volume/mferro/scripts/")



.libPaths("/shared-directory/sd-tools/apps/R/lib/")
library(feather)
library(data.table)
setDTthreads(0)
library(dplyr)
library(stringr)
library(tidyr)

fe = fread("/media/volume/mferro/data/ENDPOINTS/ADHERENCE_ANALYSIS_ENDPOINT_HISTORY.EndPt")

fe = fe %>%
  select(-date_of_birth, -end_of_followup, -start_of_followup)
ID = fe$FINREGISTRYID

fwrite(list(colnames(fe)[-1]),"/media/volume/mferro/summary_stats/covariates_health.csv")

fe_wide = fe %>%
  select(-FINREGISTRYID) %>%
  select_if(function(x){sum(x)/nrow(fe) >= 0.01}) %>%
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

fwrite(fe_p,"/media/volume/mferro/summary_stats/covariates_helath_prevalence_selection.csv")


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

fwrite(phis_dist,"/media/volume/mferro/summary_stats/covariates_health_distance_matrix.csv")

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

svg("/media/volume/mferro/plots/health_trees.svg", width = 12, height = 6)
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

svg("/media/volume/mferro/plots/diseases_selection.svg", width = 12, height = 6)
dis_plot
dev.off()


#### remove duplicates or very similar endpoints ####

fe_ps = fe_p %>%
  mutate(G = P3) %>%
  group_by(G) %>%
  filter(P == max(P))

fe_wide_s = ID %>% cbind(fe_wide[,names(fe_wide) %in% fe_ps$VAR])
colnames(fe_wide_s)[1] = "FINREGISTRYID"

fwrite(fe_ps,"/media/volume/mferro/summary_stats/covariates_helath_final_selection.csv")
fwrite(fe_wide_s,"/media/volume/mferro/data/ADHERENCE/diseases_covariates.csv")


rm(list=ls())
gc(full = T)
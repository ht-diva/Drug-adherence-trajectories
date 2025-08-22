

setwd("/media/volume/mferro/scripts")

# clean up
rm(list=ls())
gc()

# utilities
library(data.table)
library(tidyr)
library(dplyr)
library(parallel)
library(ggplot2)
library(fda)
library(pbmcapply)
library(furrr)
library(foreach)
library(lubridate)
.libPaths("/shared-directory/sd-tools/apps/R/lib/")


#### read and set ####
smoothed = fread('/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_WIDE_SMOOTHED.csv.gz')

ids =  smoothed[,1]


# prepare for common fda basis/knots
knots = seq(4,(52*5)-4,4)
n_knots = length(knots)
n_order = 4
n_basis = n_knots + n_order - 2
basis = create.bspline.basis(rangeval = c(0,52*5), breaks  = knots)

argvals = seq(0,52*5,4)
ymat = smoothed[,2:ncol(smoothed)] %>%
  mutate_if(is.character, as.numeric) %>%
  t()
str(ymat)
rm(smoothed)
invisible(gc(full = T))


#### create global fda object ####
Wobj = Data2fd(argvals = argvals, y= ymat, basisobj = basis, nderiv = 2)
#tmp = Wobj
#tmp$coefs = tmp$coefs[,100:200]
#plot(tmp)
#takes a long time

save(Wobj, file='/media/volume/mferro/data/ADHERENCE/STATINS/fda_object_5Y_MPR.RData')


load('/media/volume/mferro/data/ADHERENCE/STATINS/fda_object_5Y_MPR.RData')


####  FPCA  ####
PCA = pca.fd(Wobj, nharm = 6)

saveRDS(PCA$meanfd, "/media/volume/mferro/summary_stats/FPCA_mean_function.rds" )
saveRDS(PCA$harmonics, "/media/volume/mferro/summary_stats/FPCA_harmonics.rds" )
fwrite(data.frame(VALUES = PCA$values), "/media/volume/mferro/summary_stats/FPCA_values.csv" )
fwrite(data.frame(VARPROP = PCA$varprop), "/media/volume/mferro/summary_stats/FPCA_varprop.csv" )


svg("/media/volume/mferro/plots/FPCA_components_1to3.svg", width = 12, height = 6)
par(mfrow=c(1,3))
plot(PCA, harm = c(1,2,3), xlab = "Weeks")
dev.off()

svg("/media/volume/mferro/plots/FPCA_components_4to6.svg", width = 12, height = 6)
par(mfrow=c(1,3))
plot(PCA, harm = c(4,5,6), xlab = "Weeks")
dev.off()

sum(PCA$varprop[1:3])

#VMAX <- varmx.pca.fd(PCA)
#par(mfrow=c(1,3))
#plot(VMAX, harm = c(1,2,3))


library(ggplot2)
data_plot = as.data.frame(PCA$harmonics$coefs)
data_plot$WEEKS = seq(0,260,4)


VARPROP = fread("/media/volume/mferro/summary_stats/FPCA_varprop.csv" )
fwrite(data_plot, "/media/volume/mferro/summary_stats/FPCA_components_values.csv" )
data_plot = fread("/media/volume/mferro/summary_stats/FPCA_components_values.csv" )

# plot components
fpca_plot = ggplot(data_plot, aes(x = WEEKS)) +
  geom_line(aes(y= PC1, color =  "FPC1"), linewidth = 1.8) +
  geom_line(aes(y= PC2, color =  "FPC2"), linewidth =1.8) +
  geom_line(aes(y= PC3, color=   "FPC3"), linewidth = 1.8)  +
  scale_color_manual(name = "", values = c("FPC1" = "#E69F00", "FPC2" = "#0072B2", "FPC3" = "#009E73"),
                     labels =  paste0(c("FPC1","FPC2","FPC3"), " : ",round(VARPROP$VARPROP[c(1,2,3)]*100,2),"%")) + 
  ylim(c(-0.1,0.1)) +
  xlab("Weeks")+  
  ylab("Harmonics coefficients")+
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = 20))


svg("/media/volume/mferro/plots/FPCA_components_MAIN.svg", width = 12, height = 6)
fpca_plot
dev.off()


scores = cbind(ids,PCA$scores)
colnames(scores) = c("PATIENT_ID","PC1","PC2","PC3","PC4","PC5","PC6")
dim(scores)


a = ggplot(scores) +
  geom_density(aes(x= PC1), fill =   "#E69F00") +
  ylab("")+
  xlab("FPC1 scores")+
  theme_minimal() +
  theme(legend.position = "none")
b = ggplot(scores) +
  geom_density(aes(x= PC2), fill =  "#0072B2")  +  
  ylab("")+
  xlab("FPC2 scores")+
  theme_minimal() +
  theme(legend.position = "none")

c = ggplot(scores) +
  geom_density(aes(x= PC3), fill=   "#009E73")  +
  ylab("")+
  xlab("FPC3 scores")+
  theme_minimal()+
  theme(legend.position = "none")

library(ggpubr)

scores_plot = ggarrange(a,b,c, ncol =1)


svg("/media/volume/mferro/plots/FPCA_scores_MAIN.svg", width = 10, height = 12)
scores_plot
dev.off()



#### save and clear ####
fwrite(scores,'/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_FPCSCORES.csv',
       col.names = T,row.names = F, quote = F)

rm(list=ls())
gc()

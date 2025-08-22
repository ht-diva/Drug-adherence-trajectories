#clean up
rm(list=ls ())
gc(full=T)
setwd("N:/output/scripts/acorbetta/")

# utilities
library (data.table)
library(tidyr)
library (dplyr)
#library (npreg) #not there library (parallel)
library (ggplot2)
library (fda)
library (pbmcapply)
library (furrr)
library (foreach)
library (lubridate)
#library (kmafda) #not there
#.libPaths ("/shared-directory/sd-tools/apps/R/lib/")
library (feather)
#### read and set ####
smoothed = fread('N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_WIDE_SMOOTHED_ND_5Y_MPR.csv.gz')
ids = smoothed [,1]

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


#takes a long time
save(Wobj, file='N:/output/data/acorbetta/STATINS/ADHERENCE/fda_object_5Y_MPR.RData') 
#rm (Wobi)
invisible(gc())
load('N:/output/data/acorbetta/STATINS/ADHERENCE/fda_object_5Y_MPR.RData')



####  FPCA  ####

PCA = pca.fd(Wobj, nharm = 6)

saveRDS(PCA$meanfd, "N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_mean_function.rds" )
saveRDS(PCA$harmonics, "N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_harmonics.rds" )
fwrite(data.frame(VALUES = PCA$values), "N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_values.csv" )
fwrite(data.frame(VARPROP = PCA$varprop), "N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_varprop.csv" )


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/FPCA_components_1to3.svg", width = 12, height = 6)
par(mfrow=c(1,3))
plot(PCA, harm = c(1,2,3), xlab = "Weeks")
dev.off()

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/FPCA_components_4to6.svg", width = 12, height = 6)
par(mfrow=c(1,3))
plot(PCA, harm = c(4,5,6), xlab = "Weeks")
dev.off()

sum(PCA$varprop[1:3])

library(ggplot2)
data_plot = as.data.frame(PCA$harmonics$coefs)
data_plot$WEEKS = seq(0,260,4)


VARPROP = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_varprop.csv" )
fwrite(data_plot, "N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_components_values.csv" )
data_plot = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/FPCA_components_values.csv" )

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


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/FPCA_components_MAIN.svg", width = 12, height = 6)
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


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/FPCA_scores_MAIN.svg", width = 10, height = 12)
scores_plot
dev.off()



#### save and clear ####
fwrite(scores,'N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_FPCSCORES.csv',
       col.names = T,row.names = F, quote = F)
scores = fread('N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_FPCSCORES.csv')
sum = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_SUMMARY.csv.gz",
            select = c("PATIENT_ID","first"))[PATIENT_ID %in% scores$PATIENT_ID]

fwrite(sum, "N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file_all.csv")

rm(list=ls())
gc(full=T)

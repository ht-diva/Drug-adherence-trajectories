setwd("N:/output/scripts/acorbetta/")


# clean up
rm(list=ls())
gc(full=T)

# utilities
library(data.table)
#library(tidyr)
library(dplyr)
#library(npreg) #not there
library(parallel)
library(ggplot2)
library(fda)
library(foreach)
library(lubridate)

setDTthreads(0)
getDTthreads()


df = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_LONG_FILTERED.csv.gz", nrows = 10000)


#transform to weeks

df = df %>%
  mutate(WEEK = round(DAYS /7)) %>%
  group_by(PATIENT_ID, WEEK) %>%
  summarise(VALUE = mean(POINT_MPR_cap))
gc(full=T)

setDT(df)
df = df[order(PATIENT_ID,WEEK)]


####  widen traj ####
source("N:/output/scripts/acorbetta/adherence/fda_funs.R")
start.time = Sys.time()
df_wide <- widen_traj2(df)
time.taken <- Sys.time() - start.time 
time.taken

rm(df)
invisible(gc(full = T))

# SMOOOOOOOOTH
numCores <- 7
cl <- makeCluster(numCores)
cnames = colnames(df_wide)

gc(full=T)

start.time = Sys.time()
fobjs <- parRapply(cl,df_wide, times = cnames, q=NULL, norder=4,
                   deriv=2, lb = 1e3,
                   FUN =  spline_fda) 
time.taken <- Sys.time() - start.time
time.taken


stopCluster(cl)
showConnections()
rm(cl)



invisible(gc(full=T))


#### plot the worst smoothing examples ####

gcvs = NULL
for (i in c(1:length(fobjs))) {
  gcvs[i] = fobjs[[i]]$gcv
}
data_gcv = data.frame(GCV = gcvs)
gcvs_plot = ggplot(data_gcv, aes(x = (GCV))) +
  geom_density(fill = "#D55E00") +
  scale_x_continuous() +
  xlab("Generalized cross-validatin error") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))


svg("/media/volume/mferro/plots/GCV_smoothing_densities.svg", width = 12, height = 6)
gcvs_plot
dev.off()




####  plot ####

maxo = order(gcvs, decreasing = T)[1:10]
print(maxo)

mino = order(gcvs, decreasing = F)[1:10]
id_ex = df_wide[maxo,]$PATIENT_ID

# plot examples of worst adapting curves

par(mfrow=c(2,5))
for (i in maxo) {
  range = range(fobjs[[i]]$argvals)
  abx = seq(0,365*5/7,4)
  plot(fobjs[[i]]$argvals,fobjs[[i]]$y,xlab="days",ylab="PDC", 
       main=paste0("GCV=",signif(fobjs[[i]]$gcv,2)),
       ylim=c(0,1))
  points(abx,eval.fd(abx, fobjs[[i]]$fd, Lfd=0) ,type="l",col="blue",lwd=2)
}

par(mfrow=c(2,5))
for (i in mino) {
  range = range(fobjs[[i]]$argvals)
  abx = seq(0,365*5/7,4)
  plot(fobjs[[i]]$argvals,fobjs[[i]]$y,xlab="days",ylab="PDC", 
       main=paste0("GCV=",signif(fobjs[[i]]$gcv,2)),
       ylim=c(0,1))
  points(abx,eval.fd(abx, fobjs[[i]]$fd, Lfd=0) ,type="l",col="blue",lwd=2)
}

par(mfrow=c(1,1))
dp = data.frame(x = fobjs[[i]]$argvals, y = fobjs[[i]]$y)
ggplot(dp, aes(x = x, y =y)) +
  geom_point(aes(color = "#D55E00", size = 10)) + 
  geom_smooth(method = NULL,se = F, aes(fill = "#0072B2")) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  ylab("") + xlab("") +
  theme(legend.position = "none")

### plot example ####

i = 18
range = range(fobjs[[i]]$argvals)
abx = seq(0,365*5/7,4)
plot(fobjs[[i]]$argvals,fobjs[[i]]$y,xlab="days",ylab="PDC", 
     main=paste0("GCV=",signif(fobjs[[i]]$gcv,2)),
     ylim=c(0,1))
points(abx,eval.fd(abx, fobjs[[i]]$fd, Lfd=0) ,type="l",col="blue",lwd=2)

#### evaluate data on same axis ####
abx = seq(0,365*5/7,4)
smoothed = array(0,dim=c(length(fobjs),length(abx)))
n = dim(smoothed)[1]
for (i in c(1:length(fobjs))) {
  tmp = eval.fd(abx, fobjs[[i]]$fd, Lfd=0)
  print(i/n)
  smoothed[i,] = ifelse(tmp > 1,
                        1,
                        ifelse(tmp < 0,
                               0,
                               tmp))
}


invisible(gc(full = T))


dim(smoothed)
smoothed = cbind(df_wide$PATIENT_ID,smoothed)
fwrite(as.data.frame(smoothed),'N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_WIDE_SMOOTHED_ND_5Y_MPR.csv.gz')


# clean up
rm(list=ls())
gc(full=T)




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
.libPaths("/shared-directory/sd-tools/apps/R/lib/")
library(feather)

setDTthreads(0)
getDTthreads()

df = fread('/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_LONG_ALL_MPR_FILTERED.csv.gz')

uniqueN(df$PATIENT_ID)

#transform to weeks

df = df %>%
  mutate(WEEK = round(DAYS /7)) %>%
  group_by(PATIENT_ID, WEEK) %>%
  summarise(VALUE = mean(POINT_MPR_cap))
gc(full=T)

setDT(df)
df = df[order(PATIENT_ID,WEEK)]

###  widen traj ####
source("/media/volume/mferro/scripts/fda_funs.R")
df_wide <- widen_traj2(df)
rm(df)
gc(full=T)

# SMOOOOOOOOTH
numCores <- 30
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

#### little comparison on test traj ####
test = df_wide[6839,]
res = pos_smoothing(test,times = cnames, q=seq(0,1,0.25), norder=3,
                    deriv=1)
res4 = pos_smoothing(test,times = cnames, q=seq(0,1,0.25), norder=4,
                    deriv=2)
res$Wfdobj$coefs = exp(res$Wfdobj$coefs)
res4$Wfdobj$coefs = exp(res4$Wfdobj$coefs)
res2 = spline_fda(test,times = cnames, q=seq(0,1,0.25), norder=3,
                  deriv=1)
res3 = spline_fda(test,times = cnames, q=seq(0,1,0.25), norder=4,
                  deriv=2, lb = 500000)

test_id = test$PATIENT_ID
test_set = df[df$PATIENT_ID == test_id,]
par(mfrow=c(1,1))
plot(test_set$DAYS,test_set$POINT_MPR_CAP)
lines(res2, col = "red")
lines(res3, col = "blue")
lines(res$Wfdobj, col = "green")
lines(res4$Wfdobj, col = "orange")


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

#other plot
ggplot(data = df[df$PATIENT_ID %in% id_ex,], aes(x = DAYS, y = POINT_MPR_CAP, group = PATIENT_ID)) +
  geom_step( aes(color = factor(PATIENT_ID)), linetype = "dashed") +
  ylim(c(0,1)) + xlim(c(0,365*4)) +
  theme(legend.position = "none") +
  facet_wrap(~PATIENT_ID, ncol = 5) +
  theme(strip.text.x = element_blank())


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
fwrite(as.data.frame(smoothed), '/media/volume/mferro/data/ADHERENCE/STATINS/ADH_STATINS_WIDE_SMOOTHED.csv.gz')


# clean up
rm(list=ls())
gc(full=T)


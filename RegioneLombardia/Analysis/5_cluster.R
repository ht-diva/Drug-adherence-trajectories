setwd("N:/output/scripts/acorbetta/")

# clean up
rm(list=ls())
gc(full = T)

library(cluster)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(fda)
library(doParallel)
library(parallel)

#### kmeans on pca scores ####
scores = fread('N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_FPCSCORES.csv')

load("N:/output/data/acorbetta/STATINS/ADHERENCE/fda_object_5Y_MPR.RData")

ids = scores[,1]

scores= scores[,-1]
scores_scaled = scale(scores)

#### check optimal k #### 

# Function to perform k-means clustering and calculate total within-cluster sum of squares (explained variance)
kmeans_clustering <- function(k, data) {
  set.seed(123) # For reproducibility
  km_result <- kmeans(data, centers = k, nstart = 25)
  wss <- sum(km_result$withinss) # within-cluster sum of squares
  list(km_result = km_result, wss = wss)
}

# Define the range of clusters to analyze
cluster_range <- 1:10

# Set up parallelization
cl <- makeCluster(10) # Use all but one core
clusterExport(cl, varlist = c("scores", "kmeans_clustering")) # Export variables and functions to the cluster

# Perform k-means clustering in parallel
set.seed(1)
results <- parLapply(cl, cluster_range, kmeans_clustering, data = scores[,1:3])
stopCluster(cl) # Stop the cluster

# Extract WSS (Within-Cluster Sum of Squares) for each k
wss_values <- sapply(results, function(res) res$wss)
wss_values <- wss_values / wss_values[1] # Normalize

# Plot the Explained Variance (Elbow Method)
library(ggplot2)
elbow_plot <- ggplot(data.frame(Clusters = cluster_range, WSS = wss_values), aes(x = Clusters, y = WSS)) +
  geom_line() + 
  geom_point() +
  ylim(0, max(wss_values)) +
  labs(title = "Percentage of SS", x = "Number of Clusters (K)", y = "Within-cluster Sum of Squares (WSS)")

# Display the Elbow Plot
wss_data  = data.frame(Clusters = cluster_range, WSS = wss_values)
fwrite(wss_data, "N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/ADH_STATINS_CLUSTERS_wss.csv",
       col.names = T, row.names = F, quote = F)

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/TRAJ_groups_elbowplot.svg", width = 12, height = 6)
elbow_plot
dev.off()


#### KMEANS ####
set.seed(0)
kmeans_pca = kmeans(scores, centers = 5, nstart=100)

table(kmeans_pca$cluster)

groups = cbind(ids, kmeans_pca$cluster)
colnames(groups) = c("PATIENT_ID","GROUP")
fwrite(groups, "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_CLUSTERS.csv",
       col.names = T, row.names = F, quote = F)

groups = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_CLUSTERS.csv")




tmp = t(aggregate(t(Wobj$coefs),list(groups$GROUP), FUN = median))[-1,]
rownames(tmp) = NULL
fdmeans = Wobj
fdmeans$coefs = tmp


library(ggplot2)
data_plot = as.data.frame(fdmeans$coefs)
data_plot$WEEKS = seq(0,260,4)
colnames(data_plot) = c("G1",  "G2"  ,  "G3"   , "G4"    ,"G5"   , "WEEKS")
fwrite(data_plot, "N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/medioids_trajectories_groups.csv" )
data_plot = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/medioids_trajectories_groups.csv" )
groups[,DESC := fcase(GROUP == 1, "Low",
                      GROUP == 2, "High",
                      GROUP == 3, "Dec",
                      GROUP == 4, "Mid",
                      GROUP == 5, "Inc")]
fwrite(groups, "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_CLUSTERS.csv",
       col.names = T, row.names = F, quote = F)

fwrite(as.data.frame(table(groups$DESC)),
       "N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/groups_numerosity.csv",
       col.names = T, row.names = F, quote = F)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot components
groups_plot = ggplot(data_plot, aes(x = WEEKS)) +
  geom_smooth(aes(y= G1, color =  "G1"), linewidth = 1.5, se = F) +
  geom_smooth(aes(y= G2, color =  "G2"), linewidth = 1.5, se = F) +
  geom_smooth(aes(y= G3, color =  "G3"), linewidth = 1.5, se = F) +
  geom_smooth(aes(y= G4, color =  "G4"), linewidth = 1.5, se = F) +
  geom_smooth(aes(y= G5, color =  "G5"), linewidth = 1.5, se = F) +
  scale_color_manual(name = "", values = cbbPalette,
                     labels =  paste0(c("G1","G2","G3", "G4","G5"), " : ", table(groups$GROUP) )) + 
  ylim(c(0,1)) +
  ylab("Adherence")+
  xlab("Weeks") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = 20))


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/TRAJ_groups_MAIN.svg", width = 12, height = 6)
groups_plot
dev.off()

groups = groups[DESC == "High" | DESC == "Dec"]
sum = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_SUMMARY.csv.gz",
            select = c("PATIENT_ID","first"))[PATIENT_ID %in% groups$PATIENT_ID]

fwrite(sum, "N:/output/data/acorbetta/STATINS/ADHERENCE/cohort_file.csv")

# create grups char
dft = fread( "N:/output/data/acorbetta/STATINS/ADHERENCE/ADH_STATINS_SUMMARY.csv.gz" )

summary_stats <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    Q3 = quantile(x, 0.75, na.rm = TRUE)
  )
}

char = merge(dft,groups) %>%
  mutate(SEX = fifelse(SESSO == "M",1,0)) %>%
  dplyr::select(DESC, n_obs, time, change,lag,sd_lag, seqlen, SEX,AGE) %>%
  group_by(DESC) %>%
  summarise_all(summary_stats)

char$STAT = rep(c("MEAN","SD","1Q","MEDIAN","3Q"), 5)


fwrite(char,
       "N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/groups_characteristics.csv",
       col.names = T, row.names = F, quote = F)


# clean up
rm(list=ls())
gc(full = T)

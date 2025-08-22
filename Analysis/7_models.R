#### SET UP ####
rm(list=ls())
gc(full = T)

setwd("N:/output/scripts/acorbetta/")

library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidymodels)
library(lightgbm)
library(bonsai)
library(finetune)
library(themis)
library(future)
library(future.apply)


#load("N:/output/data/acorbetta/STATINS/ADHERENCE/working_prediction.rds")

#### READ FILE ####
cov = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ALL_COVARIATES_MPR.tsv")
cov = cov[T_GROUPS != ""]
cov = cov[,-c("PC1","PC2","PC3","PC4","PC5","PC6")]
cov[,T_GROUPS := relevel(as.factor(fifelse(T_GROUPS == "High", 0 ,1)), ref = 1)]
cov[,CONTINENT := relevel(as.factor(CONTINENT), ref = "European")]
cov[,ATS := relevel(as.factor(ATS),ref = "ATS DELLA CITTA' METROPOLITANA DI MILANO")]
cov[,COD_SOGGETTO := NULL]
cov[,YOB := NULL]
cov[,ES_FIN := NULL]
cov[,ES_INV := NULL]
cov[,ES_DIS := NULL]


#### SPLIT  ####
set.seed (123)
data_split <- initial_split(cov, prop = 0.7)
train_data <- training(data_split) 
test_data <- testing(data_split)
rm(data_split)
rm(cov)
gc(full=T)




#### MODELS ####

# set up parallelisation
library(doParallel)
cores = detectCores() - 1
registerDoParallel(cores = cores)

#recipe 
recipe <- recipe(T_GROUPS ~ ., data = train_data) %>%
  step_dummy(CONTINENT,ATS) %>%
  step_normalize(all_numeric()) %>%
  step_downsample(T_GROUPS)




#### LOGISTIC ####
logistic_model = logistic_reg()  

logistic_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(logistic_model)

# Fit the logistic model on the training data
final_logistic_model <- fit(logistic_workflow, data = train_data)

# Extract the model fitted with glm
logistic_fit <- extract_fit_parsnip(final_logistic_model)
library(car)
vif(logistic_fit$fit)


library(profvis)
p <- profvis({
  tidy_output <- broom::tidy(logistic_fit$fit, exponentiate = TRUE, conf.int = TRUE)
})



coef_df <- tidy_output %>%
  filter(term != "(Intercept)")  # Remove intercept if not needed

# Add a column for -log10(p-value)
coef_df <- coef_df %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr"),
         log_p_adj = -log10(p_adj_fdr)
  ) %>%
  mutate(
    significant = p_adj_fdr < 0.05,  # Mark as significant if p_adj_fdr < 0.05
    label = ifelse(significant, term, NA)  # Label only significant points
  )


fwrite(coef_df,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coefficients_logistic.csv")
coef_df = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coefficients_logistic.csv")
coef_df[,label := gsub("ATS_|CONTINENT_","",label)]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Volcano plot with labels for significant points
library(ggrepel)
volc_plot = ggplot(coef_df, aes(x = estimate, y = log_p_adj)) +
  geom_point(aes(color = significant), size = 1.5) +  # Color points based on significance
  #coord_trans(x = "log", y = "identity") +
  scale_x_log10(limits = c(0.8,1.25), breaks=seq(0.5, 1.25,0.05), labels = seq(0.5, 1.25,0.05)) +
  #scale_x_continuous(limits = c(0.8,1.25), breaks=seq(0.5, 1.25,0.05), labels = seq(0.5, 1.25,0.05))+ 
  scale_color_manual(values=cbPalette)+
  geom_text_repel(
    aes(label = label),
    size = 5,
    max.overlaps = 10,
    segment.curvature = 0,        # Slightly curved arrows
    segment.ncp = 5,                 # Number of control points for curve smoothness
    segment.angle = 20,              # Angle of segments
    force = 3,                       # Increase force to repel labels further
    force_pull = 0.5                 # Adjust pull force for better separation
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x = "Odds Ratio (Effect Size)",
    y = "-log10(FDR-Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/logistic_volcano_MAIN.svg", width = 10, height = 10)
volc_plot
dev.off()


#### LASSO ####


lasso_model <- logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")
# LASSO is when mixture = 1

# Workflow
lasso_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(lasso_model)

# CV set up cross-validation
set.seed(123)
cv_folds <- vfold_cv(train_data, v = 10)

gc(full = T)

# Define a grid of penalty values to search over
lambda_grid <- grid_regular(penalty(c(-5,0)), levels = 50)

# Tune the model using cross-validation
lasso_tune <- tune_grid(
  lasso_workflow,
  resamples = cv_folds,
  grid = lambda_grid,
  metrics = metric_set(roc_auc)  # Use AUC as the performance metric
)

metr = lasso_tune %>%
  collect_metrics()

fwrite(metr,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/metrics_lasso.csv")

metr %>%
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_errorbar(aes(
    ymin = mean - std_err,
    ymax = mean + std_err
  ),
  alpha = 0.5
  ) +
  geom_line(linewidth = 1.5) +
  scale_x_log10() +
  theme(legend.position = "none", legend.title = element_blank()) +
  theme_minimal()


gc(full=T)

# Extract the best tuning parameter based on ROC AUC
best_lambda <- select_by_one_std_err(lasso_tune,desc(penalty))

# Finalize the workflow with the best lambda
final_lasso_workflow <- finalize_workflow(lasso_workflow, best_lambda)

# Fit the final LASSO model on the training data
final_lasso_model <- fit(final_lasso_workflow, data = train_data)

# Extract the model fitted with glmnet
lasso_fit <- extract_fit_parsnip(final_lasso_model)

library(vip)
coef_lasso = lasso_fit %>%
  vi(lambda = best_lambda$penalty)

fwrite(coef_lasso,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coef_lasso.csv")
coef_lasso = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coef_lasso.csv")
setDT(coef_lasso)
imp_lasso_plot = coef_lasso[abs(Importance) > 0] %>%
  mutate(
    Coefficient = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x = Coefficient, y = Variable, fill = Sign)) +
  geom_col() +
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) +
  theme_minimal()

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/lasso_coefs_MAIN.svg", width = 8, height = 12)
imp_lasso_plot
dev.off()


gc(full = T)


#### LIGHT-GBM ####

# LightGBM model
lightgbm_model <- boost_tree(
  trees = 500,
  mtry = tune(),
  tree_depth = tune(),
  loss_reduction = tune(),
  learn_rate = tune()
) %>%
  set_engine("lightgbm") %>%
  set_mode("classification")

# Cross-validation
set.seed(123)
cv_folds <- vfold_cv(train_data, v = 10)

# Parameter grid
lightgbm_grid <- grid_regular(
  mtry(c(12,60)),
  tree_depth(c(2,50)),
  loss_reduction(),
  learn_rate(),
  levels = 5
)

# Workflow
lightgbm_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(lightgbm_model)

# Tune the model
gc(full=T)
set.seed(0)
t0 = Sys.time()
lightgbm_tune <- tune_grid(
  lightgbm_workflow,
  resamples = cv_folds,
  grid = lightgbm_grid,
  metrics = metric_set(roc_auc)
  )
t =  Sys.time() - t0
t

metr_lgb = lightgbm_tune %>%
  collect_metrics()
fwrite(metr_lgb,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/metrics_lgb.csv")


# Select best parameters
best_params_lgb <- select_best(lightgbm_tune,metric = "roc_auc")


# Finalize workflow
final_lightgbm_workflow <- finalize_workflow(lightgbm_workflow, best_params_lgb)

# Fit final model
final_lightgbm_model <- fit(final_lightgbm_workflow, data = train_data)

# Extract the model fitted with glmnet
lightgbm_fit <- extract_fit_parsnip(final_lightgbm_model)

save.image("N:/output/data/acorbetta/STATINS/ADHERENCE/working_prediction.rds")

# Extract the fitted LightGBM model object
lgb_model <- lightgbm_fit$fit

# Get feature importance
importance_matrix <- lgb.importance(lgb_model, percentage = T)

fwrite(importance_matrix,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/imp_lightgbm.csv")
importance_matrix = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/imp_lightgbm.csv")
importance_matrix[,Feature := gsub("ATS_|CONTINENT_","",Feature)]

imp_lgb_plot = importance_matrix[Gain>=0.01] %>%
  ggplot(aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_col(fill = "#D55E00") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) +
  theme_minimal() 


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/lgb_import_MAIN.svg", width = 8, height = 10)
imp_lgb_plot
dev.off()





#### PREDICTION ####

#### LOGISTIC 

pred_log_p <- predict(final_logistic_model, test_data, type = "prob") %>%
  bind_cols(test_data)

# Predict on test set
pred_log_b <- predict(final_logistic_model, test_data, type = "class") %>%
  bind_cols(test_data)

# Evaluate performance
auc_log = roc_auc(pred_log_p, truth = T_GROUPS, .pred_0)
acc_log = accuracy(pred_log_b, truth = T_GROUPS, .pred_class)
sens_log = sensitivity(pred_log_b, truth = T_GROUPS, .pred_class)
spec_log = specificity(pred_log_b, truth = T_GROUPS, .pred_class)
conf_log = conf_mat(pred_log_b, truth = T_GROUPS, .pred_class)
roc_log <- roc_curve(pred_log_p, truth = T_GROUPS, .pred_0) %>%
  mutate(model = "LOGISTIC")


#### LASSO
pred_lasso_p <- predict(final_lasso_model, test_data, type = "prob") %>%
  bind_cols(test_data)

# Predict on test set
pred_lasso_b <- predict(final_lasso_model, test_data, type = "class") %>%
  bind_cols(test_data)

# Evaluate performance
auc_l = roc_auc(pred_lasso_p, truth = T_GROUPS, .pred_0)
acc_l = accuracy(pred_lasso_b, truth = T_GROUPS, .pred_class)
sens_l = sensitivity(pred_lasso_b, truth = T_GROUPS, .pred_class)
spec_l = specificity(pred_lasso_b, truth = T_GROUPS, .pred_class)
conf_l = conf_mat(pred_lasso_b, truth = T_GROUPS, .pred_class)
roc_l <- roc_curve(pred_lasso_p, truth = T_GROUPS, .pred_0) %>%
  mutate(model = "LASSO")

#### LGB
pred_lgb_p <- predict(final_lightgbm_model, test_data, type = "prob") %>%
  bind_cols(test_data)

# Predict on test set
pred_lgb_b <- predict(final_lightgbm_model, test_data, type = "class") %>%
  bind_cols(test_data)

# Evaluate performance
auc_gb = roc_auc(pred_lgb_p, truth = T_GROUPS, .pred_0)
acc_gb = accuracy(pred_lgb_b, truth = T_GROUPS, .pred_class)
sens_gb = sensitivity(pred_lgb_b, truth = T_GROUPS, .pred_class)
spec_gb = specificity(pred_lgb_b, truth = T_GROUPS, .pred_class)
conf_gb = conf_mat(pred_lgb_b, truth = T_GROUPS, .pred_class)
roc_gb <- roc_curve(pred_lgb_p, truth = T_GROUPS, .pred_0) %>%
  mutate(model = "LIGHT-GB")

#### ROC CURVES ####
# Combine the two ROC datasets into one
combined_roc_data <- bind_rows(roc_log, roc_l, roc_gb)
combined_performances = data.frame(AUC = c(auc_log$.estimate, auc_l$.estimate,auc_gb$.estimate),
                                   ACC = c(acc_log$.estimate, acc_l$.estimate,acc_gb$.estimate),
                                   SENS = c(sens_log$.estimate, sens_l$.estimate,sens_gb$.estimate),
                                   SPEC = c(spec_log$.estimate, spec_l$.estimate,spec_gb$.estimate),
                                   MODEL = c("LOGISTIC","LASSO","LGB"))

fwrite(combined_roc_data,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/performances_roc.csv")
fwrite(combined_performances,"N:/output/data/acorbetta/STATINS/ADHERENCE/performances_summary.csv")


# Plot the ROC curves using ggplot2
# The palette with grey:

# The palette with black:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

roc = ggplot(combined_roc_data, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed", color = "gray" ) +
  scale_color_manual(values=cbPalette) +
  labs(title = "ROC Curves for Multiple Models", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal() +
  theme(legend.position = "bottom")



svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/ROC_CURVES_MAIN.svg", width = 10, height = 10)
roc
dev.off()




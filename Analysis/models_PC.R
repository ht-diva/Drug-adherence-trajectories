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

save.image("N:/output/data/acorbetta/STATINS/ADHERENCE/working_prediction.rds")

#### READ FILE ####
cov = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/ALL_COVARIATES_MPR.tsv")
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
data_split <- initial_split(cov[,-c("T_GROUPS"                                   
  ,"PC3"                                               
  ,"PC4"                                               
  ,"PC5"                                               
  ,"PC6")], prop = 0.7)
train_data <- training(data_split) 
test_data <- testing(data_split)
rm(cov)
gc(full=T)

train_data_1 = train_data[,-c("PC2")]
test_data_1 = test_data[,-c("PC2")]

train_data_2 = train_data[,-c("PC1")]
test_data_2 = test_data[,-c("PC1")]

#### MODEL PC1 ####

# set up parallelisation
library(doParallel)
cores = detectCores() - 1
registerDoParallel(cores = cores)

#recipe 
recipe1 <- recipe(PC1 ~ ., data = train_data_1) %>%
  step_dummy(CONTINENT,ATS) %>%
  step_normalize(all_numeric_predictors()) 




#### REGRESSION ####
regression_model = linear_reg()  %>%
  set_engine("lm")

regression_workflow <- workflow() %>%
  add_recipe(recipe1) %>%
  add_model(regression_model)

# Fit the logistic model on the training data
final_regression_model <- fit(regression_workflow, data = train_data_1)

# Extract the model fitted with glm
regression_fit <- extract_fit_parsnip(final_regression_model)
library(car)
vif(regression_fit$fit)


tidy_output <- broom::tidy(regression_fit$fit, conf.int = TRUE)



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


fwrite(coef_df,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coefficients_regression_PC1.csv")
coef_df = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coefficients_regression_PC1.csv")
coef_df[,label := gsub("ATS_|CONTINENT_","",label)]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Volcano plot with labels for significant points
library(ggrepel)
volc_plot = ggplot(coef_df, aes(x = estimate, y = log_p_adj)) +
  geom_point(aes(color = significant), size = 1.5) +  # Color points based on significance
  #coord_trans(x = "log", y = "identity") +
  #scale_x_log10(limits = c(0.8,1.25), breaks=seq(0.5, 1.25,0.05), labels = seq(0.5, 1.25,0.05)) +
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
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Beta (Effect Size)",
    y = "-log10(FDR-Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/logistic_volcano_MAIN_regression_PC1.svg", width = 10, height = 10)
volc_plot
dev.off()


#### LASSO ####


lasso_model <- linear_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")
# LASSO is when mixture = 1

# Workflow
lasso_workflow <- workflow() %>%
  add_recipe(recipe1) %>%
  add_model(lasso_model)

# CV set up cross-validation
set.seed(123)
cv_folds <- vfold_cv(train_data_1, v = 10)

gc(full = T)

# Define a grid of penalty values to search over
lambda_grid <- grid_regular(penalty(c(-5,0)), levels = 50)

# Tune the model using cross-validation
lasso_tune <- tune_grid(
  lasso_workflow,
  resamples = cv_folds,
  grid = lambda_grid,
  metrics = metric_set(rmse)  # Use AUC as the performance metric
)

metr = lasso_tune %>%
  collect_metrics()

fwrite(metr,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/metrics_lasso_regression_PC1.csv")

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

fwrite(coef_lasso,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coef_lasso_regression_PC1.csv")
coef_lasso = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coef_lasso_regression_PC1.csv")
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

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/lasso_coefs_MAIN_regression_PC1.svg", width = 8, height = 12)
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
  set_mode("regression")

# Cross-validation
set.seed(123)
cv_folds <- vfold_cv(train_data_1, v = 10)

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
  add_recipe(recipe1) %>%
  add_model(lightgbm_model)

# Tune the model
gc(full=T)
set.seed(0)
t0 = Sys.time()
lightgbm_tune <- tune_grid(
  lightgbm_workflow,
  resamples = cv_folds,
  grid = lightgbm_grid,
  metrics = metric_set(rmse)
)
t =  Sys.time() - t0
t

metr_lgb = lightgbm_tune %>%
  collect_metrics()
fwrite(metr_lgb,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/metrics_lgb_regression_PC1.csv")


# Select best parameters
best_params_lgb <- select_best(lightgbm_tune,metric = "rmse")


# Finalize workflow
final_lightgbm_workflow <- finalize_workflow(lightgbm_workflow, best_params_lgb)

# Fit final model
final_lightgbm_model <- fit(final_lightgbm_workflow, data = train_data_1)

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


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/lgb_import_MAIN_regression_PC1.svg", width = 8, height = 10)
imp_lgb_plot
dev.off()





#### PREDICTION ####

#### REGRESSION 

# Make predictions on the test data
pred_reg1 <- predict(final_regression_model, new_data = test_data_1)

# Combine the predictions with actual values
results_reg1 <- bind_cols(PC1 = test_data_1[, PC1], PRED = pred_reg1)

# Calculate performance metrics
metrics_reg1 <- results_reg1 %>%
  metrics(truth = PC1, estimate = .pred) %>%
  mutate(model = "LSR")


#### LASSO
# Make predictions on the test data
pred_lasso1 <- predict(final_lasso_model, new_data = test_data_1)

# Combine the predictions with actual values
results_lasso1 <- bind_cols(PC1 = test_data_1[, PC1], PRED = pred_lasso1)

# Calculate performance metrics
metrics_lasso1 <- results_lasso1 %>%
  metrics(truth = PC1, estimate = .pred) %>%
  mutate(model = "LASSO")

#### LGB
# Make predictions on the test data
pred_lgb1 <- predict(final_lightgbm_model, new_data = test_data_1)

# Combine the predictions with actual values
results_lgb1 <- bind_cols(PC1 = test_data_1[, PC1], PRED = pred_lgb1)

# Calculate performance metrics
metrics_lgb1 <- results_lgb1 %>%
  metrics(truth = PC1, estimate = .pred) %>%
  mutate(model = "LGB")


#### COMPARISON ####

combined_performances1 = bind_rows(metrics_reg1, metrics_lasso1, metrics_lgb1)

fwrite(combined_performances1,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/performances_regression_PC1.csv")


#### MODEL PC2 ####

#recipe 
recipe2 <- recipe(PC2 ~ ., data = train_data_2) %>%
  step_dummy(CONTINENT,ATS) %>%
  step_normalize(all_numeric_predictors()) 




#### REGRESSION ####
regression_model = linear_reg()  %>%
  set_engine("lm")

regression_workflow_2 <- workflow() %>%
  add_recipe(recipe2) %>%
  add_model(regression_model)

final_regression_model_2 <- fit(regression_workflow_2, data = train_data_2)

# Extract the model fitted with glm
regression_fit_2 <- extract_fit_parsnip(final_regression_model_2)
library(car)
vif(regression_fit_2$fit)



tidy_output_2 <- broom::tidy(regression_fit_2$fit, conf.int = TRUE)



coef_df_2 <- tidy_output_2 %>%
  filter(term != "(Intercept)")  # Remove intercept if not needed

# Add a column for -log10(p-value)
coef_df_2 <- coef_df_2 %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr"),
         log_p_adj = -log10(p_adj_fdr)
  ) %>%
  mutate(
    significant = p_adj_fdr < 0.05,  # Mark as significant if p_adj_fdr < 0.05
    label = ifelse(significant, term, NA)  # Label only significant points
  )


fwrite(coef_df_2,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coefficients_regression_PC2.csv")
coef_df_2 = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coefficients_regression_PC2.csv")
coef_df_2[,label := gsub("ATS_|CONTINENT_","",label)]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Volcano plot with labels for significant points
library(ggrepel)
volc_plot = ggplot(coef_df_2, aes(x = estimate, y = log_p_adj)) +
  geom_point(aes(color = significant), size = 1.5) +  # Color points based on significance
  #coord_trans(x = "log", y = "identity") +
  #scale_x_log10(limits = c(0.8,1.25), breaks=seq(0.5, 1.25,0.05), labels = seq(0.5, 1.25,0.05)) +
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
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Beta (Effect Size)",
    y = "-log10(FDR-Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/regression_volcano_MAIN_regression_PC2.svg", width = 10, height = 10)
volc_plot
dev.off()


#### LASSO ####


lasso_model <- linear_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")
# LASSO is when mixture = 1

# Workflow
lasso_workflow_2 <- workflow() %>%
  add_recipe(recipe2) %>%
  add_model(lasso_model)

# CV set up cross-validation
set.seed(123)
cv_folds <- vfold_cv(train_data_2, v = 10)

gc(full = T)

# Define a grid of penalty values to search over
lambda_grid <- grid_regular(penalty(c(-5,0)), levels = 50)

# Tune the model using cross-validation
lasso_tune_2 <- tune_grid(
  lasso_workflow_2,
  resamples = cv_folds,
  grid = lambda_grid,
  metrics = metric_set(rmse)  # Use AUC as the performance metric
)

metr_2 = lasso_tune_2 %>%
  collect_metrics()

fwrite(metr_2,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/metrics_lasso_regression_PC2.csv")

metr_2 %>%
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
best_lambda_2 <- select_by_one_std_err(lasso_tune_2,desc(penalty))

# Finalize the workflow with the best lambda
final_lasso_workflow_2 <- finalize_workflow(lasso_workflow_2, best_lambda_2)

# Fit the final LASSO model on the training data
final_lasso_model_2 <- fit(final_lasso_workflow_2, data = train_data_2)

# Extract the model fitted with glmnet
lasso_fit_2 <- extract_fit_parsnip(final_lasso_model_2)

library(vip)
coef_lasso_2 = lasso_fit_2 %>%
  vi(lambda = best_lambda_2$penalty)

fwrite(coef_lasso_2,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coef_lasso_regression_PC2.csv")
coef_lasso_2 = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/coef_lasso_regression_PC2.csv")
setDT(coef_lasso_2)
imp_lasso_plot_2 = coef_lasso_2[abs(Importance) > 0] %>%
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

svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/lasso_coefs_MAIN_regression_PC2.svg", width = 8, height = 12)
imp_lasso_plot_2
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
  set_mode("regression")

# Cross-validation
set.seed(123)
cv_folds <- vfold_cv(train_data_2, v = 10)

# Parameter grid
lightgbm_grid <- grid_regular(
  mtry(c(12,60)),
  tree_depth(c(2,50)),
  loss_reduction(),
  learn_rate(),
  levels = 5
)

# Workflow
lightgbm_workflow_2 <- workflow() %>%
  add_recipe(recipe2) %>%
  add_model(lightgbm_model)

# Tune the model
gc(full=T)
set.seed(0)
t0 = Sys.time()
lightgbm_tune_2 <- tune_grid(
  lightgbm_workflow_2,
  resamples = cv_folds,
  grid = lightgbm_grid,
  metrics = metric_set(rmse)
)
t =  Sys.time() - t0
t

metr_lgb_2 = lightgbm_tune_2 %>%
  collect_metrics()
fwrite(metr_lgb_2,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/metrics_lgb_regression_PC2.csv")

# Select best parameters
best_params_lgb_2 <- select_best(lightgbm_tune_2,metric = "rmse")


# Finalize workflow
final_lightgbm_workflow_2 <- finalize_workflow(lightgbm_workflow_2, best_params_lgb_2)

# Fit final model
final_lightgbm_model_2 <- fit(final_lightgbm_workflow_2, data = train_data_2)

# Extract the model fitted with glmnet
lightgbm_fit_2 <- extract_fit_parsnip(final_lightgbm_model_2)

save.image("N:/output/data/acorbetta/STATINS/ADHERENCE/working_prediction_regression.rds")

# Extract the fitted LightGBM model object
lgb_model_2 <- lightgbm_fit_2$fit

# Get feature importance
importance_matrix_2 <- lgb.importance(lgb_model_2, percentage = T)

fwrite(importance_matrix_2,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/imp_lightgbm_regression_PC2.csv")
importance_matrix_2 = fread("N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/imp_lightgbm_regression_PC2.csv")
importance_matrix_2[,Feature := gsub("ATS_|CONTINENT_","",Feature)]


imp_lgb_plot_2 = importance_matrix_2[Gain>=0.01] %>%
  ggplot(aes(x = Gain, y = reorder(Feature, Gain))) +
  geom_col(fill = "#D55E00") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) +
  theme_minimal() 


svg("N:/output/data/acorbetta/STATINS/ADHERENCE/plots/lgb_import_MAIN_regression_PC2.svg", width = 8, height = 10)
imp_lgb_plot_2
dev.off()




#### PREDICTION ####

#### REGRESSION 

# Make predictions on the test data
pred_reg2 <- predict(final_regression_model_2, new_data = test_data_2)

# Combine the predictions with actual values
results_reg2 <- bind_cols(PC2 = test_data_2[, PC2], PRED = pred_reg2)

# Calculate performance metrics
metrics_reg2 <- results_reg2 %>%
  metrics(truth = PC2, estimate = .pred) %>%
  mutate(model = "LSR")


#### LASSO
# Make predictions on the test data
pred_lasso2 <- predict(final_lasso_model_2, new_data = test_data_2)

# Combine the predictions with actual values
results_lasso2 <- bind_cols(PC2 = test_data_2[, PC2], PRED = pred_lasso2)

# Calculate performance metrics
metrics_lasso2 <- results_lasso2 %>%
  metrics(truth = PC2, estimate = .pred) %>%
  mutate(model = "LASSO")

#### LGB
# Make predictions on the test data
pred_lgb2 <- predict(final_lightgbm_model_2, new_data = test_data_2)

# Combine the predictions with actual values
results_lgb2 <- bind_cols(PC2 = test_data_2[, PC2], PRED = pred_lgb2)

# Calculate performance metrics
metrics_lgb2 <- results_lgb2 %>%
  metrics(truth = PC2, estimate = .pred) %>%
  mutate(model = "LGB")


#### COMPARISON ####

combined_performances2 = bind_rows(metrics_reg2, metrics_lasso2, metrics_lgb2)

fwrite(combined_performances2,"N:/output/data/acorbetta/STATINS/ADHERENCE/summary_stats/performances_regression_PC2.csv")



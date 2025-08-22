#### functions for splines ####


#f for smoothing the trajectories
# 4 parameters: 1)time series values 2)times 
# 3)quantiles for the knots (equally spaced)
# 4) degree of derivates 
# 5) custum quantiles
spline_fda <- function(time_values,times,q = NULL, deriv=2 , norder = 4, lb = NULL) {
  # calculate observations and time values
  library(fda)
  #print(time_values)
  time_values = as.numeric(time_values[-1])
  #print(time_values)
  times = as.numeric(times[-1])
  y <- as.numeric(time_values[which(!is.na(time_values))])
  #print(y)
  x <- as.numeric(times[which(!is.na(time_values))])
  #print(x)
  x = x[order(x)]
  print(x)
  # select knots with quantiles 
  if (is.null(q)) {breaks = x}
  else {breaks = round(unname(quantile(x,probs=q, type=1))) }
  #create basis
  basis <- create.bspline.basis(breaks = breaks, norder=norder)
  
  #select lambda
  #if (is.null(lb)) {
  #  lambda <- c(1e6,5e5,1e5,5e4,1e4,1e3,1e2,1)
  #  gcv <- numeric(length(lambda))
  #  for (i in 1:length(lambda)){
  #    functionalPar <- fdPar(fdobj=basis, Lfdobj=deriv, lambda=lambda[i])  
  #    gcv[i] <- smooth.basis(x, y, functionalPar)$gcv
  #    lb <- lambda[which.min(gcv)] } }
  # create parameter object to fit spline with 1st and 2nd derivates
  # weight the effect of the minalization w.r.t the least square tearm
  # functional parameter, having arguments: 
  # basis, order of the derivative to be penalized, smoothing parameter.
  functionalPar <- fdPar(fdobj=basis, Lfdobj=int2Lfd(deriv), lambda=lb)
  Xss <- smooth.basis(x, y, functionalPar)
  return(Xss)
  }


wide_to_spline <- function(pat, df, q = NULL, deriv=2 , nodes = 5) {
  # calculate observations and time values
  library(fda)
  library(dplyr)
  
  df_w = df %>%
    dplyr::select(PATIENT_ID, PDC, DAYS) %>% 
    filter(PATIENT_ID == pat) %>%
    pivot_wider(names_from = DAYS, values_from = PDC, id_cols = PATIENT_ID) %>%
    as.data.frame() 
  
  #assign index
  rownames(df_w) = df_w$PATIENT_ID
  df_w = df_w[2:ncol(df_w)]

  #order the columns
  new_order = as.character(sort(as.numeric(colnames(df_w)),decreasing=F))
  df_w <- df_w[, new_order]
  
  y <- as.numeric(df_w[1,])
  x <- as.numeric(colnames(df_w))
  # select knots with quantiles 
  if (is.null(q)) {breaks = x}
  else { breaks = round(unname(quantile(x,probs=q, type=1))) }
  print(length(breaks))
  #create basis
  basis <- create.bspline.basis(breaks, norder=3)
  
  #select lambda
  lambda <- c(1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
  gcv <- numeric(length(lambda))
  for (i in 1:length(lambda)){
    #print(lambda[i])
    functionalPar <- fdPar(fdobj=basis, Lfdobj=deriv, lambda=lambda[i])  
    #print(functionalPar)
    gcv[i] <- smooth.basis(x, y, functionalPar)$gcv
  }
  lb <- lambda[which.min(gcv)]
  #create parameter object to fit spline with 1st and 2nd derivates
  functionalPar <- fdPar(fdobj=basis, Lfdobj=deriv, lambda=lb) #weight the effect of the minalization w.r.t the least square tearm
  # functional parameter, having arguments: 
  # basis, order of the derivative to be penalized, smoothing parameter.
  Xss <- smooth.basis(x, y, functionalPar)
  # estimate obtain with the spline 
  return(Xss)
}


widen_traj <- function(df) {
  df_wide <- df %>%
    dplyr::select(PATIENT_ID, POINT_MPR_CAP, DAYS) %>% 
    tidyr::pivot_wider(names_from = DAYS, values_from = POINT_MPR_CAP, id_cols = PATIENT_ID,
                       values_fn = mean) %>%
    as.data.frame() 
  
  # 1 000 samples is ok
  # 10 000 samples is ok
  # 100 000 samples is ok
  # 200 000 samplpes is not ok w 8gb ram (Error: cannot allocate vector of size 7.8 Gb)
  # maybe we could increase ram or divide in batches, smooth separately, then merge
  
  #order the columns
  new_order = c("PATIENT_ID",(as.character(sort(as.numeric(colnames(df_wide)[-1]),decreasing=F))))
  df_wide <- df_wide[, new_order]
  return(df_wide)
}


widen_traj2 <- function(dt) {
  setDT(dt)
  dt_wide = dcast(dt,
                  PATIENT_ID ~ WEEK,
                  value.var = "VALUE",
                  fun.aggregate = mean)
  setcolorder(dt_wide,c("PATIENT_ID",sort(as.numeric(names(dt_wide)[-1]))))
  gc(full=T)
  return(dt_wide)
}
  


pos_smoothing <- function(time_values,times,q = NULL, deriv=2 , norder = 4) {
  library(fda)
  time_values = time_values[-1]
  y <- as.numeric(time_values[which(!is.na(time_values))])
  #print(y)
  x <- as.numeric(times[which(!is.na(time_values))])
  #print(x)
  # select knots with quantiles 
  if (is.null(q)) {breaks = x}
  else { breaks = round(unname(quantile(x,probs=q, type=1))) }
  print(length(breaks))
  #create basis
  basis <- create.bspline.basis(breaks, norder=norder)
  lambda <- c(1e6,5e5,1e5)
  gcv <- numeric(length(lambda))
  for (i in 1:length(lambda)){
    #print(lambda[i])
    functionalPar <- fdPar(fdobj=basis, Lfdobj=deriv, lambda=lambda[i])  
    #print(functionalPar)
    gcv[i] <- smooth.basis(x, y, functionalPar)$gcv
  }
  lb <- lambda[which.min(gcv)]
  #select lambda
  #create parameter object to fit spline with 1st and 2nd derivates
  functionalPar <- fdPar(fdobj=basis, Lfdobj=deriv, lambda=lb) #weight the effect of the minalization w.r.t the least square tearm
  # functional parameter, having arguments: 
  # basis, order of the derivative to be penalized, smoothing parameter.
  Xss <- smooth.pos(x, y, functionalPar)
  # estimate obtain with the spline
  return(Xss)
}


###Mixed-effects Weighted Quantile Sum regression models

library(future)
library(future.apply)
library(Rsolnp)
library(lme4)
             

gwqs_re <- function(formula, data, weights, mix_name, stratified, valid_var, b = 100,
                 b1_pos = TRUE, b1_constr = FALSE, q = 4, validation = 0.6,
                 family = gaussian, id=NULL, seed = NULL, pred = 0, plan_strategy = "sequential",
                 control = list(rho = 1, outer.iter = 400, inner.iter = 800, delta = 1.0e-7, tol = 1e-8, trace = 0),
                 lp = 0, ln = 0)
{

	if(is.character(family)){
	  if(!family %in% c("gaussian", "binomial") ){
		stop("Only gaussian and binomial family supported")
	  } 
	  family <- get(family, mode = "function", envir = parent.frame())  
	}

	if(is.function(family)) {
	  tmp <- family()
	  if(!tmp$family %in% c("gaussian", "binomial")) {
		stop("Only gaussian and binomial family supported")
	  }
	  family <- tmp  
	}

	if(is.null(family$family)) {
	  stop("Family not recognized")  
	}

  if(!any(grepl("wqs", rownames(attr(terms(formula), "factors"))))) stop("'formula' must contain 'wqs' term: e.g. y ~ wqs + ...")
 
  allvars = all.vars(formula)
  data$wqs = rnorm(nrow(data), 0, 1)  # 初始化wqs
  data <- data[, c(allvars, mix_name, id)]
  dtf <- na.omit(data)
  
  # select rows of the original dataset based on rows selected with na.action

  # defining quantile variables
  if (is.null(q)) {Q = as.matrix(dtf[, mix_name]); q_i = NULL}
  else {
    q_f = quantile_f(dtf, mix_name, q)
    Q = q_f$Q
    q_i = q_f$q_i
  }

  # create stratified variables if stratified is not NULL
  m <- 0
  if (m){
    if(is.null(q)) stop("q must be different from NULL if you want to stratify for a categorical variable")
    strtfd_out = stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name = strtfd_out$mix_name
  }

  if(!is.null(seed) & !is.numeric(seed)) stop("seed must be numeric or NULL")
  if(!is.null(seed)) set.seed(seed)

  N <- nrow(dtf)

  # splitting the dataset
  m <- 0
  rindex = create_rindex(dtf, N, validation, pred, valid_var, m, family)

  # parameters estimation and model fitting
  m <- 0
  if(m[1]) weights <- dtf[, weights, drop = FALSE]
  else  dtf$weights <- weights <- rep(1, N)
  par_model <- par.modl.est(dtf[rindex$it,], Q[rindex$it,], formula, weights[rindex$it], b, b1_pos, b1_constr, family, id, plan_strategy, control, lp, ln)
  
  bres = par_model$bres
  conv = par_model$conv
  bindex = par_model$bindex
  nfuneval = par_model$n_funeval
  n_levels = par_model$n_levels
  strata_names = par_model$strata_names

  mean_weight = mean_weight_f(mix_name, bres, conv, b1_pos, family, n_levels, strata_names)

  # fit the final model with the estimated weights
  wqs_model = model.fit(mean_weight, dtf[rindex$iv,], Q[rindex$iv,], family, id, formula, weights[rindex$iv], b1_pos)  # id

  # prediction
  if(!is.null(rindex$ip)) df_pred = predict_f(Q[rindex$ip,], dtf[rindex$ip,], mean_weight, wqs_model$m_f, formula)
  else df_pred = NULL

  # Plots
  data_plot <- data.frame(mix_name, mean_weight)

  # data_plot = data_plot[order(data_plot$mean_weight, decreasing = TRUE),]
  wqs_index = as.numeric(unlist(wqs_model$wqs))

  # creating the list of elements to return
  results = list(wqs_model$m_f, conv, bres, wqs_index, q_i, bindex, rindex$it, rindex$iv,
                 data_plot, df_pred, rindex$ip)
  names(results) = c("fit", "conv", "bres", "wqs", "q_i", "bindex", "tindex", "vindex",
                     "final_weights", "df_pred", "pindex")  

  return(results)
}


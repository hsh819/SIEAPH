### Mixed-effects WQS Related Functions

# function to create stratified elements in the mixture by levels of factors
stratified_f = function(Q, dtf, stratified, mix_name){
  ls = levels(unlist(dtf[, stratified, drop = FALSE]))
  if(is.null(ls)) stop("'stratified' must be factor")
  llsm = lapply(ls, function(ls){
    mat = diag(as.numeric(dtf[, stratified] == ls))
    sub_dt = mat%*%Q[, mix_name]
    colnames(sub_dt) = paste(mix_name, ls, sep = "_")
    return(sub_dt)
  })
  Q = do.call("cbind", llsm)
  mix_name = colnames(Q)
  strtfd_out = list(Q, mix_name)
  names(strtfd_out) = c("Q", "mix_name")

  return(strtfd_out)
}


# function to create variables with quantile of the components
quantile_f <- function(dtf, mix_name, q){

  if(!is.numeric(q)) stop("'q' must be a number")
  Ql <- lapply(1:length(mix_name), function(i){
    q_i <- unique(quantile(dtf[[mix_name[i]]], probs = seq(0, 1, by = 1/q), na.rm = TRUE))
    if(length(q_i) == 1) q_i = c(-Inf, q_i)
    q <- cut(dtf[[mix_name[i]]], breaks = q_i, labels = FALSE, include.lowest = TRUE) - 1
    return(list(q_i, q))
  })
  q_i <- lapply(Ql, function(x) x[[1]])
  Q <- matrix(unlist(lapply(Ql, function(x) x[[2]])), ncol = length(Ql))
  colnames(Q) <- names(q_i) <- mix_name

  qf_out <- list(Q, q_i)
  names(qf_out) <- c("Q", "q_i")
  return(qf_out)
}


# function to split the dataset
create_rindex <- function(dtf, N, validation, pred, valid_var, m, family){

  if(!m){
    if(!is.numeric(validation) | validation < 0 | validation >= 1) stop("'validation' must be numeric >= 0 and < 1")
    if(!is.numeric(pred) | pred < 0 | pred >= 1) stop("'pred' must be numeric >= 0 and < 1")
    if(pred + validation >= 1) stop("the sum of pred and validation must be between 0 and 1")
    if(pred>0 & family$family == "multinomial") stop("The predictive model is not available for multinomial regression")
    groups <- rep(0, N)
    if(validation > 0) groups[sample(1:N, round(N*validation))] <-1
    if(pred > 0) groups[sample(which(groups!=1), round(N*pred))] <-2
  }
  else{
    groups = unlist(dtf[, valid_var, drop = FALSE])
    if(!any(unique(groups) %in% c(0, 1, 2))) stop("valid_var values must be 0, 1 or 2")
    if(!(0 %in% unique(groups))) stop(("0 must identify test dataset"))
  }
  it = which(groups == 0)
  iv = which(groups == 1)
  ip = which(groups == 2)
  if(length(iv) == 0) iv = it
  if(length(ip) == 0) ip = NULL

  indexl = list(it, iv, ip)
  names(indexl) = c("it", "iv", "ip")
  return(indexl)
}


# parameter names in model matrix
parnames <- function(df, formula, form2){
  if(!is.null(form2)){
    mf <- model.frame(form2, df)
    Y <- model.response(mf, "any")
    df$yz <- ifelse(Y == 0, 0, 1)
  }
  mm <- model.matrix(formula, df)
  colnames(mm)
}


# functtion to define the objective function
objfn <- function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, formula, weights, b1_constr, lp, ln){

  mf <- model.frame(formula, bdtf)
  Y <- model.response(mf, "any")
  if(family$family == "binomial" & class(Y) %in% c("factor", "character")){
    if(class(Y) == "character") Y = factor(Y)
    Y <- as.numeric(Y != levels(Y)[1])
  }
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0

  w <- initp[(kx + 1):(kx + kw)]
  bdtf$wqs <- as.numeric(Q%*%w)
  X <- model.matrix(formula, bdtf)
  b_covs <- initp[1:kx]

  term <- as.numeric(X%*%b_covs) + offset

  f = sum(family$dev.resids(y = Y, mu = family$linkinv(term), wt = weights)) + lp*max(0, initp["wqs"]) - ln*min(0, initp["wqs"])

  return(f)
}


# function to define the equality constraint
linconst = function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, formula, weights, b1_constr, lp, ln){

  w = initp[(kx + 1):(kx + kw)]
  wsum = sum(w)

  return(wsum)
}


# function to determine bounded prameters
bounded_param = function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, formula, weights, b1_constr, lp, ln){

  bp = initp[(kx + 1):(kx + kw*(n_levels-1))]
  if (b1_constr){
    wqs_site = which(Xnames == "wqs")
    bp = c(initp[wqs_site], bp)
  }

  return(bp)
}


# function to define the lower bounds
LBound = function(kw, n_levels, family, b1_pos, b1_constr){

  LB = rep(0, kw)
  if(b1_constr) LB = c(sapply(b1_pos, function(i) ifelse(i, 0, -Inf)), LB)

  return(LB)
}


# function to define the upper bounds
UBound = function(kw, n_levels, family, b1_pos, b1_constr){

  UB = rep(1, kw)
  if(b1_constr) UB = c(sapply(b1_pos, function(i) ifelse(i, Inf, 0)), UB)

  return(UB)
}


# function to define the parameters initial values
values.0 = function(kw, Xnames, kx, n_levels, formula, weights, bdtf, b1_pos, family, id){

  w = rep(1/kw, kw)
  if(is.null(id)) {
	fit = glm(formula, bdtf, family = family, weights = weights)
	bj = coef(fit)
  }
  else {
   
	formula = as.formula(paste(deparse(formula), paste('(1 |', id, ')'), sep="+"))
	if (family$family %in% c("gaussian")) fit = lmer(formula, bdtf, weights = weights)
	else  fit = glmer(formula, bdtf, family = family, weights = weights)
	bj = summary(fit)$coef[, 1]
  
  }   
    
  val.0 = c(bj, w)
  wqs_site <- which(grepl("wqs", Xnames))
  val.0[wqs_site] = sapply(b1_pos, function(i) ifelse(i, 0.0001, -0.0001))

  return(val.0)
}


# optimization function to estimate the weights
optim.f <- function(bdtf, bQ, b1_pos, b1_constr, family, id, formula, weights, control, lp, ln){

  Xnames <- parnames(bdtf, formula, NULL)
  kx <- length(Xnames)

  n_levels <- 2; level_names <- wqsvars <- NULL 
  kw <- dim(bQ)[2]
  initp <- values.0(kw, Xnames, kx, n_levels, formula, weights, bdtf, b1_pos, family, id)
  LowB <- LBound(kw, n_levels, family, b1_pos, b1_constr)
  UpB <- UBound(kw, n_levels, family, b1_pos, b1_constr)
  eq_b <- rep(1, length(b1_pos))

  if(control$trace) cat("start opt")

  opt_res <- tryCatch(solnp(pars = initp, fun = objfn, eqfun = linconst, eqB = eq_b, ineqfun = bounded_param,
                            ineqLB = LowB, ineqUB = UpB, control = control, kw = kw, bdtf = bdtf, Q = bQ, kx = kx,
                            Xnames = Xnames, n_levels = n_levels, level_names = level_names, wqsvars = wqsvars, family = family, 
                            formula = formula, weights = weights, b1_constr = b1_constr, lp = lp, ln = ln), error = function(e) NULL)  # not include id

  if(!is.null(opt_res)) {

	par_opt <- opt_res$pars[(kx + 1):(kx + kw)]
    conv <- opt_res$convergence
    nfuneval <- opt_res$nfuneval
  }
  else {

	par_opt <- initp[(kx + 1):(kx + kw)]
    conv <- 2
    nfuneval <- 0
  }
  out <- list(par_opt, conv, nfuneval)
  names(out) <- c("par_opt", "conv", "nfuneval")

  return(out)
}


# function that fit the wqs model
model.fit <- function(w, bdtf, bQ, family, id, formula, weights, b1_pos){

  if(any(is.na(w)))
  {
    w = rep(1/length(w), length(w))
  }
  
  wqs <- bdtf$wqs <- as.numeric(bQ%*%w)

  if(is.null(id))  m_f = glm(formula, data = bdtf, family = family)
  else {
  
	  formula = as.formula(paste(deparse(formula), paste('(1 |', id, ')'), sep="+"))
	  if (family$family %in% c("gaussian")) m_f = lmer(formula, data = bdtf)
	  else  m_f = glmer(formula, data = bdtf, family = family)

  }

  mf_out = list(wqs, m_f)
  names(mf_out) = c("wqs", "m_f")

  return(mf_out)
}

# function to sample from data rownames
sample_f = function(i, dtf, id){

  if(is.null(id))  bindex = sample(1:nrow(dtf), nrow(dtf), replace=TRUE)
  else { 
	  rownames(dtf) = 1:nrow(dtf)
	  layers <- unique(dtf[[id]])
	  bindex <- c()
	  for(layer in layers){
		layer_data <- dtf[dtf[[id]] == layer, , drop = FALSE]  # 保留原数据框格式不变
		index <- rownames(layer_data)
		index <- sample(index, replace = TRUE)
		bindex <- c(bindex, index)
	  }
	  bindex = as.numeric(bindex)
  }  

  return(bindex)
}

# function that call optim.f to estimate parameters for each bootstrap sample
estimate_param = function(bindex, dtf, Q, b1_pos, b1_constr, family, id, formula, weights, control, lp, ln){

  param = optim.f(dtf[bindex,], Q[bindex,], b1_pos, b1_constr, family, id, formula, weights[bindex], control, lp, ln)

  return(param)
}


# function to be passed to future_lapply to fit the models
model.fit_f = function(i, param, bindex, dtf, Q, b1_pos, family, id, formula, weights){

  b_fit <- model.fit(param[[i]]$par_opt, dtf[bindex[[i]],], Q[bindex[[i]],], family, id, formula, weights, b1_pos)

  return(b_fit)
}

# function that calls the optimization function and the function to fit the model for each bootstrap sample
par.modl.est <- function(dtf, Q, formula, weights, b, b1_pos, b1_constr, family, id, plan_strategy, control, lp, ln){

  if (family$family %in% c("gaussian")) ts = "t"
  else  ts = "z"

  if(!is.numeric(b)) stop("'b' must be a number")

  bindex = lapply(X=1:b, FUN = sample_f, dtf = dtf, id = id)

  plan(plan_strategy)
  param <- future_lapply(X = bindex, FUN = estimate_param, dtf = dtf, Q = Q, b1_pos = b1_pos, b1_constr = b1_constr,
                         family = family, id = id, formula = formula, weights = weights, control = control, lp = lp, ln = ln, future.seed = TRUE)

  plan(plan_strategy)
  b_fit <- future_lapply(X = 1:b, FUN = model.fit_f, param = param, bindex = bindex, dtf = dtf, Q = Q, b1_pos = b1_pos,
                         family = family, id = id, formula = formula, weights = weights, future.seed = TRUE)

  conv <- c(sapply(param, function(i) i$conv))
  nfuneval <- c(sapply(param, function(i) i$nfuneval))

  wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt))

  b1 <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", "Estimate"])
  se <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", "Std. Error"])
  stat <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", paste0(ts, " value")])
  p_val <- sapply(b_fit, function(i) {  # summary(i$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")] 
		stats = summary(i$m_f)$coefficients["wqs", paste0(ts, " value")]
		p = 2*pnorm(-abs(stats))                              	
  })

  n_levels <- 1

  n_non_conv = sum(conv == 2)
  if(n_non_conv == 0 & control$trace) cat(paste0("The optimization function always converged\n"))
  else if(n_non_conv == b) stop("The optimization function never converged\n")
  else if(control$trace) cat(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))

  # estimate mean weight for each component (exclude weights from iterations with failed convergence)

  bres <- as.data.frame(cbind(wght_matrix, b1, se, stat, p_val))
  names(bres) <- c(colnames(Q), "b1", "Std_Error", "stat", "p_val")

  strata_names <- NULL

  par_model_out <- list(bres, conv, bindex, nfuneval, n_levels, strata_names)
  names(par_model_out) <- c("bres", "conv", "bindex", "nfuneval", "n_levels", "strata_names")

  return(par_model_out)
}


# function to estimate mean weights for each component
mean_weight_f = function(mix_name, bres, conv, b1_pos, family, n_levels, strata_names){

  if(b1_pos) mean_weight = apply(bres[bres$b1 > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1 > 0 & conv!=2, "stat"]))
  else mean_weight = apply(bres[bres$b1 < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1 < 0 & conv!=2, "stat"]))

  # if(all(is.nan(mean_weight)))
  #   stop("There are no ", ifelse(b1_pos, "positive", "negative"), " b1 in the bootstrapped models")

  return(mean_weight)
}


# function to build predictive model in case of binomial dist
predict_f = function(Q, dtf, w, m_f, formula){

  dtf$wqs = as.numeric(Q%*%w)
  pred = predict(m_f, newdata = dtf, type = "response")
  y = model.response(model.frame(formula, dtf), "any")
  df_pred = data.frame(y = y, p_y = pred)

  return(df_pred)
}


#' @title Simulate quantized exposures for RFQGC methods
#'
#' @description Simulate quantized exposures for RFQGC methods
#'  
#' @details Simulate right censored survival outcome from Weibull model
#' 
#' @param n Sample size
#' @param corr NULL, or vector of correlations between the first exposure and subsequent exposures (if length(corr) < (length(coef)-1), then this will be backfilled with zeros)
#' @param coef Vector of coefficients for the outcome (i.e. model coefficients for exposures). The length of this determines the number of exposures.
#' @param q Number of levels or "quanta" of each exposure
#' @param shape0 (survival outcomes) baseline shape of weibull distribution \link[stats]{rweibull}
#' @param scale0 (survival outcomes) baseline scale of weibull distribution \link[stats]{rweibull}
#' @param censtime (survival outcomes) administrative censoring time
#' @param mat NULL, or interaction matrix, which must be a symmetric array with elements of 1 or 0, where 1 means there is an interaction between x_i and x_j and 0 means there is no interaction between x_i and x_j.
#' @param ncheck (logical, default=TRUE) adjust sample size if needed so that exposures are exactly evenly distributed (so that qgcomp::quantize(exposure) = exposure)
#'
#' @return a data frame
#' @export
#' @examples
#' set.seed(123)
#' qdat = simdata_surv(
#'   N=10000, corr=c(0.6, 0.3), coef=c(1,1,0,0), 
#'   q = 10, shape0 = 1, scale0 = 10, censtime = 0.1,
#'   mat = matrix(c(0.5, rep(0, 15)), nrow=4, ncol=4, byrow=T)
#' )
#' cor(qdat)
#' 

# qdat = simdata_surv(
#    N=10000, corr=c(.5, .3), coef=c(1,1,0,0), 
#    q = 10, shape0 = 2, scale0 = 10, censtime = 1,
#    mat = matrix(c(0.5, rep(0, 15)), ncol=4, byrow=T) 
# )



simdata_surv = function(
  N = 100,  # sample size
  coef = c(1,0,0,0),  # beta coefficients for X in the outcome model 
  corr = NULL,  # Pearson/spearman (the same here) correlation
  q = 4,
  shape0 = 1,
  scale0 = 10,
  censtime = 1,
  mat = NULL,  # interaction coefficients
  ncheck = TRUE
) {
  rem = N %% q
  if(rem != 0 & ncheck) {
    N = N + rem
    message("Adjusting sample size to have evenly distributed exposure values")
  }
  
  if(!is.null(mat)){
	 if(all(mat == 0)) stop("mat is all zeros")
	 if(any(is.na(mat))) stop("mat contains NA")
	 if(!is_lower_triangular(mat)) stop("mat is not a lower triangular matrix")
	 if(length(coef) != ncol(mat)) stop("The dimension of mat is inconsistent with the length of the coefficients")

  }  

  lst = .quantized_design(N=N, coef=coef, corr=corr, q=q, mat=mat)

  t0 = rweibull(N, shape = shape0, scale = scale0)
  tmg = pmin(censtime, exp(log(t0) - lst$mu/(shape0)))
  
  status = 1.0*(tmg < censtime)
  if(mean(status) < .10) warning("model implies < 10% of observations experience the event of interest")
  dat = data.frame(lst$X, time=tmg, status=status)
  attr(dat, "truecoefs") = list(coef=coef, inter=mat)  
  return(dat)

}

# adjust lower triangular matrix
is_lower_triangular = function(mat) {
  if(nrow(mat) != ncol(mat))
    stop("Matrix rows and columns not equal")
  
  n = nrow(mat)
  
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      if(mat[i, j] != 0) {
        return(FALSE) 
      }
    }
  }
  
  return(TRUE)
}

# design matrix maker for a quantized version of exposures and a single modifier
.quantized_design = function(
  N = 100,                    # sample size
  coef = c(1,0,0,0),          # beta coefficients for X in the outcome model
  corr = 0.75,                # Pearson/spearman (the same here) correlation
  q = 4,
  mat = NULL
) {
  p = length(coef)
  corv = rep(0, length(coef))
  corv[1] = 1.0
  if(is.null(corr)) corr = corv[-1]
  if(length(corr) >= 1) corv[1+(1:length(corr))] = corr
  X = matrix(nrow=N, ncol=p)
  xtemplate = sample(rep(0:(q-1), length.out=N), N, replace=FALSE)
  for(k in 1:p) {
    newx = numeric(N)
    if(corv[k] >= 0) {
      c1 = as.logical(rbinom(N, 1, corv[k]))
      newx[which(c1)] = xtemplate[which(c1)]
      newx[which(!c1)] = sample(xtemplate[which(!c1)])
    }
    if(corv[k] < 0) {
      c1 = as.logical(rbinom(N, 1, -corv[k]))
      newx[which(c1)] = q-1-xtemplate[which(c1)]
      newx[which(!c1)] = sample(q-1-xtemplate[which(!c1)])
    }
    X[,k] = newx
  }
  colnames(X) = paste0("x", 1:p)
  
  if(is.null(mat)) {
    mu = X %*% coef 
  } else {
	  idx = which(mat != 0, arr.ind = TRUE)
	  coef2 = mat[idx]
	  m = nrow(idx)
	  X2 = matrix(nrow=N, ncol=m)
	  for(i in 1:m)
	  {
	    X2[, i] = X[, idx[i,1]] * X[, idx[i, 2]]
	  }  
	  
	  mu = X %*% coef + X2 %*% coef2  
  }
  
  list(mu = mu, X = X)
}



# Correctly identify interaction effects
imat = matrix(c(0.20, 0, 0, 0,
                0, 0, 0, 0,
                0, 0.15, 0, 0,
                0, 0, 0, 0), ncol=4, byrow=TRUE)
# scenario 1
set.seed(123)
dat1 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.0, q=10, shape0=3, scale0=10, censtime=0.01, mat=imat)

# scenario 2
set.seed(123)
dat2 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.75, q=10, shape0=3, scale0=10, censtime=0.01, mat=imat)

# scenario 3
set.seed(123)
dat3 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.0, q=10, shape0=3, scale0=10, censtime=0.02, mat=imat)

# scenario 4
set.seed(123)
dat4 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.75, q=10, shape0=3, scale0=10, censtime=0.02, mat=imat)

# scenario 5
set.seed(123)
dat5 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.0, q=10, shape0=4, scale0=12, censtime=0.06, mat=imat)

# scenario 6
set.seed(123)
dat6 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.75, q=10, shape0=4, scale0=12, censtime=0.06, mat=imat)

# scenario 7
set.seed(123)
dat7 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.0, q=10, shape0=4, scale0=12, censtime=0.1, mat=imat)

# scenario 8
set.seed(123)
dat8 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.75, q=10, shape0=4, scale0=12, censtime=0.1, mat=imat)


# Incorrectly identify interaction effects
imat1 = matrix(c(0, 0, 0, 0,
                0.20, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0), ncol=4, byrow=TRUE)

set.seed(123)
dat9 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.4, mat=imat1)

set.seed(123)
dat10 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.7, mat=imat1)

imat2 = matrix(c(0, 0, 0, 0,
                0.20, 0, 0, 0,
                0, 0.15, 0, 0,
                0, 0, 0, 0), ncol=4, byrow=TRUE)

set.seed(123)
dat11 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.2, mat=imat2)

set.seed(123)
dat12 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.3, mat=imat2)

imat3 = matrix(c(0, 0, 0, 0,
                0.20, 0, 0, 0,
                0, 0.15, 0, 0,
                0, -0.10, 0, 0), ncol=4, byrow=TRUE)
				
set.seed(123)
dat13 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.3, mat=imat3)

set.seed(123)
dat14 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.5, mat=imat3)

imat4 = matrix(c(0, 0, 0, 0,
                0.20, 0, 0, 0,
                0, 0.15, 0, 0,
                0, -0.10, -0.05, 0), ncol=4, byrow=TRUE)

set.seed(123)
dat15 = simdata_surv(N=500, coef=c(0.50, 0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.4, mat=imat4)

set.seed(123)
dat16 = simdata_surv(N=500, coef=c(0.50, -0.25, 0, 0), corr=0.75, q=10, shape0=5, scale0=11, censtime=0.7, mat=imat4)

				 				

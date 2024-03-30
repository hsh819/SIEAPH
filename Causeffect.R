###Causeffect.R

# If the target of inference is the ATE, optimal or generalized full matching, subclassification, or profile matching can be used. 
# If the target of inference is the ATT or ATC, any matching method may be used. When retaining the target estimand is not so important, 
# additional options become available that involve discarding units in such a way that the original estimand is distorted. 
# These include matching with a caliper, matching within a region of common support, cardinality matching, or exact or coarsened exact matching, 
# perhaps on a subset of the covariates.

# Report SMDs before and after matching for each covariate, any prognostically important interactions between covariates, and the prognostic score; this can be reported in a table or in a Love plot.
# Report summaries of balance for other statistics, e.g., the largest mean and maximum eCDF difference among the covariates and the largest SMD among squares, cubes, and interactions of the covariates.

# A marginal effect is a comparison between the expected potential outcome under treatment and the expected potential outcome under control. 
# A conditional effect is the comparison between the expected potential outcomes in the treatment groups within strata.

# The RR, OR, and HR are noncollapsible effect measures, which means the marginal effect on that scale is not a (possibly) weighted average of the conditional effects within strata, even if the stratum-specific effects are of the same magnitude.
# The mean difference and risk difference (RD) are collapsible effect measures, so the same methods can be used to estimate marginal and conditional effects.

# Although there are many possible ways to include covariates (e.g., not just main effects but interactions, smoothing terms like splines, or other nonlinear transformations), 
# it is important not to engage in specification search.
# For this reason, we recommend only including the same terms included in the propensity score model unless there is a strong a priori and justifiable reason to model the outcome differently.


########################################UMAP降维######################################################
# Uniform Manifold Approximation and Projection (UMAP) is an algorithm for dimension reduction based on manifold learning techniques and ideas from topological data analysis. 
# It provides a very general framework for approaching manifold learning and dimension reduction, but can also provide specific concrete realizations.

# Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, 
# but also for general non-linear dimension reduction. The algorithm is founded on three assumptions about the data
# 1.The data is uniformly distributed on Riemannian manifold;
# 2.The Riemannian metric is locally constant (or can be approximated as such);
# 3.The manifold is locally connected.
# From these assumptions it is possible to model the manifold with a fuzzy topological structure. 
# The embedding is found by searching for a low dimensional projection of the data that has the closest possible equivalent fuzzy topological structure.	


compute_min_MMD = function(x, is.synthetic=TRUE, interval=1) {

	if(is.numeric(x)) {  # 模拟数据
		if((x + 99) %% interval == 0) {
			if(is.synthetic) {
				set.seed(x + 99) 
				dat = simdata_synthetic()				
			} else {
				set.seed(x + 99) 
				dat = simdata_causal()			
			}		
			XC = dat[dat$T == 0, 1:(ncol(dat) - 2)]
			XT = dat[dat$T == 1, 1:(ncol(dat) - 2)]
			n = floor(nrow(XC)/ncol(XC))
			if(n < 1) {
			  set.seed(x + 99)
			  M = min_MMD(XC, XT, m=1)
			  warning(paste0('The ', x, 'th iteration: ', 'Control group rows are less than columns'))		
			} else if(n >= 1 & n <= 10) {
				set.seed(x + 99)
				M = min_MMD(XC, XT, m=n)		
			} else if(n > 10) {
				set.seed(x + 99)
				M = min_MMD(XC, XT, m=10)			
			}
		
			return(M)	
		}			
	} else if(is.data.frame(x) | is.matrix(x)) {  # 真实数据
		dat = x
		N = ncol(dat) - 2
		XC = dat[dat$T == 0, 1:N]
		XT = dat[dat$T == 1, 1:N]
		n = floor(nrow(XC)/ncol(XC))
		if(n < 1) {
		  set.seed(N)
		  M = min_MMD(XC, XT, m=1)
		  warning('Control group rows are less than columns')		
		} else if(n >= 1 & n <= 10) {
			set.seed(N)
			M = min_MMD(XC, XT, m=n)		
		} else if(n > 10) {
			set.seed(N)
			M = min_MMD(XC, XT, m=10)			
		}
	
		return(M)	
	
	} else {
		stop("x input error")
	}
		
}


mul_kdrm = function(dat, d, iter=1, replace=TRUE, MMD=list(), interval=1, num=5) {

	N = ncol(dat) - 2
	XC = dat[dat$T == 0, 1:N]
	XT = dat[dat$T == 1, 1:N]
	
	if(is.logical(iter)) {  # 真实数据
		M = MMD
	} else {  # 模拟数据
		k = ceiling(iter / interval)
		M = MMD[[k]]	
	}
	
	# M = min_MMD(XC, XT, m=1)	
	XC_Q = as.matrix(XC) %*% M$min_Q
	colnames(XC_Q) = colnames(XC)

	kernel = gaussian_kernel(XC_Q, XT, kernel_num = num, fix_sigma = M$min_bandwidth)
	Xdist = 2 * (1 - kernel$kernel_val/num)
#	Xdist = sqrt(Xdist)
	
	if(is.logical(iter) == FALSE & replace == TRUE) {  # 模拟数据，重复匹配
		n_neighbors = c(18, 20, 22) # c(16, 18, 20, 22, 24)
		min_dist = c(0.10, 0.30, 0.50) # c(0.10, 0.30, 0.50)		
	} else if(is.logical(iter) == FALSE & replace == FALSE) {  # 模拟数据，不重复匹配
		n_neighbors = c(5, 15, 25)  # c(12, 15, 35)
		min_dist = c(0.10) 		
	} else if(is.logical(iter) == TRUE & replace == TRUE) {  # 真实数据，重复匹配
		n_neighbors = c(10, 18, 22)  # Lalonde: c(10, 18, 22)/c(5, 18, 22)  Air pollution: 
		min_dist = c(0.10)	# LaLonde: c(0.10)  Air pollution:	
	} else {  # 真实数据，不重复匹配
		n_neighbors = c(10, 18, 25)  # LaLonde: c(10, 18, 25)/c(5, 18, 25)  Air pollution: c(10, 18, 22), c(10, 16, 22)
		min_dist = c(0.10)	# LaLonde: c(0.10)/c(0.01, 0.10, 0.20)  Air pollution: c(0.10)
	}

	future::plan(future::multisession)
	# options(future.globals.maxSize = 2 * 1024 ^ 3)
	res = future_lapply(n_neighbors, future.seed=TRUE, FUN=function(i) {
		future_lapply(min_dist, future.seed=TRUE, FUN=function(j) {
		  custom.config = umap.defaults
		  custom.config$random_state = ifelse(is.logical(iter), N, iter + 99)
		  custom.config$n_components = d
		  custom.config$n_neighbors = i
		  custom.config$min_dist = j		
		  umap.dist = umap(Xdist, config=custom.config, input="dist")
		  layout = umap.dist$layout
		  colnames(layout) = paste0("KD", 1:d)
		  return(layout)		 
	  })
	})
	future::plan(future::sequential)  # 关闭进程
	
	return(res)
}


########################################ATT计算函数######################################################

euc_ATT = function(dat, replace=TRUE, formula.select=1, is.m=FALSE) {
	if(formula.select == 1) {
		formula1 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), sep="+"))	
	} else if(formula.select == 2) {
		if( all(c('age', 'education', 're74', 're75', 'un74', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re74^2)', 'I(re75^2)', 'education:re74', 'un74:hispanic', sep="+")		
		} else if( all(c('age', 'education', 're75', 'un75', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re75^2)', 'education:re75', 'un75:hispanic', sep="+")				
		} else {
			newvar = NULL
		}				
		formula1 = as.formula(paste("T ~", paste(paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+"))		
	} else if(formula.select == 3) {
		if( all(c('PM10', 'NO2', 'SO2', 'O3', 'CO', 'age', 'BMI', 'sex1', 'education1', 'education2', 'history1', 'history2') %in% names(dat)) ) {
			newvar = paste('I(PM10^2)', 'I(NO2^2)', 'I(SO2^2)', 'I(O3^2)', 'I(CO^2)', 'I(age^2)', 'I(BMI^2)',   # 'PM10^2', 'NO2^2', 'SO2^2', 'O3^2', 'CO^2', 'age^2', 'BMI^2'
							'age:BMI', 'age:sex1', 'age:education1', 'age:education2', 'BMI:sex1', 'SO2:education1', 'SO2:education2', 'T:age', 'T:BMI', 'T:sex1',
							'O3:CO', 'O3:PM10', 'O3:SO2', sep="+")						
		} else {
			newvar = NULL
		}				
		formula1 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+"))	
	}

	# Matching
	m = matchit(formula1, data = dat, method = "nearest", distance = "euclidean", replace=replace)  # euclidean matching
	# Estimating Treatment Effects
	if(replace) {
		m.data = get_matches(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula2, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula2, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR					
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR					
		}	
	} else {
		m.data = match.data(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula2, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula2, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR	
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR			
		}
	}
	
	# Reporting results
	if(length(unique(dat$Y)) != 2) {		
		res = c(ATT$estimate, ATT$std.error, ATT$conf.low, ATT$conf.high)
	} else {
		lnRR = c(ATT_lnRR$estimate, ATT_lnRR$std.error, ATT_lnRR$conf.low, ATT_lnRR$conf.high)
		lnOR = c(ATT_lnOR$estimate, ATT_lnOR$std.error, ATT_lnOR$conf.low, ATT_lnOR$conf.high)	
		res =  rbind(lnRR, lnOR)
		rownames(res) = paste('euc', c('lnRR', 'lnOR'), sep='_')
	}
	if(is.m) {
		return(list(res=res, m=m))
	} else {
		return(res)		
	}
}

mah_ATT = function(dat, replace=TRUE, formula.select=1, is.m=FALSE) {
	if(formula.select == 1) {
		formula1 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), sep="+"))	
	} else if(formula.select == 2) {
		if( all(c('age', 'education', 're74', 're75', 'un74', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re74^2)', 'I(re75^2)', 'education:re74', 'un74:hispanic', sep="+")		
		} else if( all(c('age', 'education', 're75', 'un75', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re75^2)', 'education:re75', 'un75:hispanic', sep="+")			
		} else {
			newvar = NULL
		}				
		formula1 = as.formula(paste("T ~", paste(paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+"))		
	} else if(formula.select == 3) {
		if( all(c('PM10', 'NO2', 'SO2', 'O3', 'CO', 'age', 'BMI', 'sex1', 'education1', 'education2', 'history1', 'history2') %in% names(dat)) ) {
			newvar = paste('I(PM10^2)', 'I(NO2^2)', 'I(SO2^2)', 'I(O3^2)', 'I(CO^2)', 'I(age^2)', 'I(BMI^2)',   # 'PM10^2', 'NO2^2', 'SO2^2', 'O3^2', 'CO^2', 'age^2', 'BMI^2'
							'age:BMI', 'age:sex1', 'age:education1', 'age:education2', 'BMI:sex1', 'SO2:education1', 'SO2:education2', 'T:age', 'T:BMI', 'T:sex1',
							'O3:CO', 'O3:PM10', 'O3:SO2', sep="+")							
		} else {
			newvar = NULL
		}				
		formula1 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+"))	
	}
	
	# Matching
	m = matchit(formula1, data = dat, method = "nearest", distance = "mahalanobis", replace=replace)  # mahalanobis matching
	# Estimating Treatment Effects
	if(replace) {
		m.data = get_matches(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula2, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula2, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR					
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR					
		}	
	} else {
		m.data = match.data(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula2, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula2, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR	
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR			
		}
	}
	
	# Reporting results
	if(length(unique(dat$Y)) != 2) {
		res = c(ATT$estimate, ATT$std.error, ATT$conf.low, ATT$conf.high)
	} else {
		lnRR = c(ATT_lnRR$estimate, ATT_lnRR$std.error, ATT_lnRR$conf.low, ATT_lnRR$conf.high)
		lnOR = c(ATT_lnOR$estimate, ATT_lnOR$std.error, ATT_lnOR$conf.low, ATT_lnOR$conf.high)	
		res =  rbind(lnRR, lnOR)
		rownames(res) = paste('mah', c('lnRR', 'lnOR'), sep='_')
	}
	if(is.m) {
		return(list(res=res, m=m))
	} else {
		return(res)		
	}	
}

psm_ATT = function(dat, replace=TRUE, formula.select=1, is.m=FALSE) {
	if(formula.select == 1) {
		formula1 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), sep="+"))	
	} else if(formula.select == 2) {
		if( all(c('age', 'education', 're74', 're75', 'un74', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re74^2)', 'I(re75^2)', 'education:re74', 'un74:hispanic', sep="+")		
		} else if( all(c('age', 'education', 're75', 'un75', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re75^2)', 'education:re75', 'un75:hispanic', sep="+")			
		} else {
			newvar = NULL
		}				
		formula1 = as.formula(paste("T ~", paste(paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+"))		
	} else if(formula.select == 3) {
		if( all(c('PM10', 'NO2', 'SO2', 'O3', 'CO', 'age', 'BMI', 'sex1', 'education1', 'education2', 'history1', 'history2') %in% names(dat)) ) {
			newvar = paste('I(PM10^2)', 'I(NO2^2)', 'I(SO2^2)', 'I(O3^2)', 'I(CO^2)', 'I(age^2)', 'I(BMI^2)',   # 'PM10^2', 'NO2^2', 'SO2^2', 'O3^2', 'CO^2', 'age^2', 'BMI^2'
							'age:BMI', 'age:sex1', 'age:education1', 'age:education2', 'BMI:sex1', 'SO2:education1', 'SO2:education2', 'T:age', 'T:BMI', 'T:sex1',
							'O3:CO', 'O3:PM10', 'O3:SO2', sep="+")							
		} else {
			newvar = NULL
		}				
		formula1 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula2 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), newvar, sep="+"))	
	}
	
	# Matching	
	m = matchit(formula1, data = dat, method = "nearest", distance = "glm", replace=replace)  # propensity score matching
	# Estimating Treatment Effects
	if(replace) {
		m.data = get_matches(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula2, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula2, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR					
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR					
		}	
	} else {
		m.data = match.data(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula2, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula2, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR	
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR			
		}
	}
	
	# Reporting results
	if(length(unique(dat$Y)) != 2) {		
		res = c(ATT$estimate, ATT$std.error, ATT$conf.low, ATT$conf.high)
	} else {
		lnRR = c(ATT_lnRR$estimate, ATT_lnRR$std.error, ATT_lnRR$conf.low, ATT_lnRR$conf.high)
		lnOR = c(ATT_lnOR$estimate, ATT_lnOR$std.error, ATT_lnOR$conf.low, ATT_lnOR$conf.high)	
		res =  rbind(lnRR, lnOR)
		rownames(res) = paste('psm', c('lnRR', 'lnOR'), sep='_')
	}	
	if(is.m) {
		return(list(res=res, m=m))
	} else {
		return(res)		
	}		
}

pca_ATT = function(dat, d, replace=TRUE, formula.select=1, distance="euclidean", is.m=FALSE) {
	# PCA
	newnames = paste0("PC", 1:d)
	formula31 = as.formula(paste("T ~", paste0(newnames, collapse="+")))		
	if(formula.select == 1) {
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(newnames, collapse="+"), sep="+"))		
	} else if(formula.select == 2) {
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), sep="+"))	
	} else if(formula.select == 3) {
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), paste0(newnames, collapse="+"), sep="+"))			
	} else if(formula.select == 4) {
		if( all(c('age', 'education', 're74', 're75', 'un74', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re74^2)', 'I(re75^2)', 'education:re74', 'un74:hispanic', sep="+")		
		} else if( all(c('age', 'education', 're75', 'un75', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re75^2)', 'education:re75', 'un75:hispanic', sep="+")			
		} else {
			newvar = NULL
		}
		formula32 = as.formula(paste("T ~", paste(paste0(newnames, collapse="+"), newvar, sep="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), paste0(newnames, collapse="+"), newvar, sep="+"))					
	} else if(formula.select == 5) {
		if( all(c('PM10', 'NO2', 'SO2', 'O3', 'CO', 'age', 'BMI', 'sex1', 'education1', 'education2', 'history1', 'history2') %in% names(dat)) ) {
			newvar = paste('I(PM10^2)', 'I(NO2^2)', 'I(SO2^2)', 'I(O3^2)', 'I(CO^2)', 'I(age^2)', 'I(BMI^2)',  # 'PM10^2', 'NO2^2', 'SO2^2', 'O3^2', 'CO^2', 'age^2', 'BMI^2'
							'age:BMI', 'age:sex1', 'age:education1', 'age:education2', 'BMI:sex1', 'SO2:education1', 'SO2:education2', 'T:age', 'T:BMI', 'T:sex1',
							'O3:CO', 'O3:PM10', 'O3:SO2', sep="+")							
		} else {
			newvar = NULL
		}				
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), paste0(newnames, collapse="+"), newvar, sep="+"))
	}	

	# pca = princomp(dat[, 1:(ncol(dat) - 2)], cor = T)  # cor = T
	# scores = pca$scores[, 1:d]
	pca = prcomp(dat[, 1:(ncol(dat) - 2)], scale. = T)
	scores = pca$x[, 1:d]
	colnames(scores) = newnames
	dat = cbind(dat, scores)		
	
	# Matching: euclidean, scaled_euclidean, mahalanobis, robust_mahalanobis, cosine
	if(formula.select != 5) {
		m = matchit(formula32, data = dat, method = "nearest", distance = distance, replace=replace)			
	} else if(formula.select == 5) {
		if(distance == 'euclidean') {
			Dist = euclidean_dist(formula31, data = dat)
			m = matchit(formula32, data = dat, method = "nearest", distance = Dist, replace=replace)
		} else if(distance == 'scaled_euclidean') {
			Dist = scaled_euclidean_dist(formula31, data = dat)
			m = matchit(formula32, data = dat, method = "nearest", distance = Dist, replace=replace)  	
		} else if(distance == 'mahalanobis') {
			Dist = mahalanobis_dist(formula31, data = dat)
			m = matchit(formula32, data = dat, method = "nearest", distance = Dist, replace=replace)  
		} else if(distance == 'robust_mahalanobis') {
			Dist = robust_mahalanobis_dist(formula31, data = dat)
			m = matchit(formula32, data = dat, method = "nearest", distance = Dist, replace=replace) 
		} else if(distance == 'cosine') {
			Dist = cosine_dist(formula31, data = dat)
			m = matchit(formula32, data = dat, method = "nearest", distance = Dist, replace=replace) 	
		} else {
			stop('distance input error')
		}	
	}
	
	if(replace) {
		m.data = get_matches(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula4, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula4, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR					
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR					
		}	
	} else {
		m.data = match.data(m)
		if(length(unique(dat$Y)) != 2)  # 连续响应变量
		{
			fit = lm(formula4, data = m.data, weights = weights)
			ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights")				
		} else {  # 二分类响应变量
			fit = glm(formula4, data = m.data, weights = weights, family = quasibinomial())
			ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR	
			ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR			
		}
	}
	
	# Reporting results
	if(length(unique(dat$Y)) != 2) {		
		res = c(ATT$estimate, ATT$std.error, ATT$conf.low, ATT$conf.high)
	} else {
		lnRR = c(ATT_lnRR$estimate, ATT_lnRR$std.error, ATT_lnRR$conf.low, ATT_lnRR$conf.high)
		lnOR = c(ATT_lnOR$estimate, ATT_lnOR$std.error, ATT_lnOR$conf.low, ATT_lnOR$conf.high)	
		res =  rbind(lnRR, lnOR)
		rownames(res) = paste('pca', c('lnRR', 'lnOR'), sep='_')
	}
	if(is.m) {
		return(list(res=res, m=m))
	} else {
		return(res)		
	}	
}


# cosine_dist = function(df) {
  
  # mat = as.matrix(df)
  # norms = sqrt(rowSums(mat^2))    # 计算每个向量的模
  # mat = mat / norms    # 标准化每个向量
  # cosine_similarity_matrix = mat %*% t(mat)    # 计算夹角余弦相似度矩阵
  # dist_matrix = 1 - cosine_similarity_matrix    # 将余弦相似度转换为余弦距离

  # return(dist_matrix)
# }

cosine_dist = function(formula, data) {

  if (!inherits(formula, "formula") | !is.data.frame(data)) {    # 检查输入是否正确
    stop("The first argument must be a formula and the second must be a data frame")
  }

  response = all.vars(formula)[1]    # 提取分类变量和协变量
  covariates = all.vars(formula)[-1]
  group_1 = data[data[[response]] == 1, covariates]     # 根据分类变量分割数据
  group_0 = data[data[[response]] == 0, covariates]
  
  mat_1 = as.matrix(group_1)    # 将数据框转换为矩阵
  mat_0 = as.matrix(group_0)
  norms_1 = sqrt(rowSums(mat_1^2))    # 计算每个向量的模
  norms_0 = sqrt(rowSums(mat_0^2))
  mat_1 = mat_1 / norms_1    # 标准化每个向量
  mat_0 = mat_0 / norms_0

  cosine_similarity_matrix = mat_1 %*% t(mat_0)    # 计算两组之间的余弦相似度矩阵
  dist_matrix = 1 - cosine_similarity_matrix    # 将余弦相似度转换为余弦距离

  return(dist_matrix)
}

kdrm_ATT = function(dat, d, iter=1, replace=TRUE, formula.select=1, MMD=list(), interval=1, distance='cosine', is.m=FALSE) {
	# KDRM
	kdrm = mul_kdrm(dat=dat, d=d, iter=iter, replace=replace, MMD=MMD, interval=interval)
	n1 = length(kdrm)
	n2 = length(kdrm[[1]])
	newnames = colnames(kdrm[[1]][[1]])		
	formula31 = as.formula(paste("T ~", paste0(newnames, collapse="+")))		
	if(formula.select == 1) {
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(newnames, collapse="+"), sep="+"))		
	} else if(formula.select == 2) {
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), sep="+"))	
	} else if(formula.select == 3) {
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), paste0(newnames, collapse="+"), sep="+"))			
	} else if(formula.select == 4) {
		if( all(c('age', 'education', 're74', 're75', 'un74', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re74^2)', 'I(re75^2)', 'education:re74', 'un74:hispanic', sep="+")		
		} else if( all(c('age', 'education', 're75', 'un75', 'hispanic') %in% names(dat)) ) {
			newvar = paste('I(age^2)', 'I(age^3)', 'I(education^2)', 'I(re75^2)', 'education:re75', 'un75:hispanic', sep="+")			
		} else {
			newvar = NULL
		}
		formula32 = as.formula(paste("T ~", paste(paste0(newnames, collapse="+"), newvar, sep="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), paste0(newnames, collapse="+"), newvar, sep="+"))					
	} else if(formula.select == 5) {
		if( all(c('PM10', 'NO2', 'SO2', 'O3', 'CO', 'age', 'BMI', 'sex1', 'education1', 'education2', 'history1', 'history2') %in% names(dat)) ) {
			newvar = paste('I(PM10^2)', 'I(NO2^2)', 'I(SO2^2)', 'I(O3^2)', 'I(CO^2)', 'I(age^2)', 'I(BMI^2)',  # 'PM10^2', 'NO2^2', 'SO2^2', 'O3^2', 'CO^2', 'age^2', 'BMI^2'
							'age:BMI', 'age:sex1', 'age:education1', 'age:education2', 'BMI:sex1', 'SO2:education1', 'SO2:education2', 'T:age', 'T:BMI', 'T:sex1',
							'O3:CO', 'O3:PM10', 'O3:SO2', sep="+")							
		} else {
			newvar = NULL
		}				
		formula32 = as.formula(paste("T ~", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+")))	
		formula4 = as.formula(paste("Y ~ T", paste0(setdiff(names(dat), c('T', 'Y')), collapse="+"), paste0(newnames, collapse="+"), newvar, sep="+"))
	}
	
	EstOR = EstRR = Est = array(data=NA, dim=c(n1, n2, 4))
	lstm = replicate(n1, replicate(n2, NULL), simplify = FALSE)
	for(i in 1:n1) {
		for(j in 1:n2) {
			temp = cbind(dat, kdrm[[i]][[j]])
						
			# Matching: euclidean, scaled_euclidean, mahalanobis, robust_mahalanobis, cosine
			if(formula.select != 5) {
				m = matchit(formula32, data = temp, method = "nearest", distance = distance, replace=replace)						
			} else if(formula.select == 5) {
				if(distance == 'euclidean') {
					Dist = euclidean_dist(formula31, data = temp)
					m = matchit(formula32, data = temp, method = "nearest", distance = Dist, replace=replace)  # "data", "random", or "closest"
				} else if(distance == 'scaled_euclidean') {
					Dist = scaled_euclidean_dist(formula31, data = temp)
					m = matchit(formula32, data = temp, method = "nearest", distance = Dist, replace=replace)  	
				} else if(distance == 'mahalanobis') {
					Dist = mahalanobis_dist(formula31, data = temp)
					m = matchit(formula32, data = temp, method = "nearest", distance = Dist, replace=replace)  
				} else if(distance == 'robust_mahalanobis') {
					Dist = robust_mahalanobis_dist(formula31, data = temp)
					m = matchit(formula32, data = temp, method = "nearest", distance = Dist, replace=replace) 
				} else if(distance == 'cosine') {
					Dist = cosine_dist(formula31, data = temp)
					m = matchit(formula32, data = temp, method = "nearest", distance = Dist, replace=replace) 	
				} else {
					stop('distance input error')
				}			
			}
			
			if(replace) {
				m.data = get_matches(m)
				if(length(unique(dat$Y)) != 2)  # 连续响应变量
				{
					fit = lm(formula4, data = m.data, weights = weights)
					ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights")				
				} else {  # 二分类响应变量
					fit = glm(formula4, data = m.data, weights = weights, family = quasibinomial())
					ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR					
					ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass + id, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR					
				}	
			} else {
				m.data = match.data(m)
				if(length(unique(dat$Y)) != 2)  # 连续响应变量
				{
					fit = lm(formula4, data = m.data, weights = weights)
					ATT = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights")				
				} else {  # 二分类响应变量
					fit = glm(formula4, data = m.data, weights = weights, family = quasibinomial())
					ATT_lnRR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnratioavg")  # transform = "exp", logRR	
					ATT_lnOR = avg_comparisons(fit, variables = "T", vcov = ~ subclass, newdata = subset(m.data, T == 1), wts = "weights", comparison = "lnoravg")  # logOR			
				}
			}

			if(length(unique(dat$Y)) != 2) {
				Est[i, j, 1] = ATT$estimate
				Est[i, j, 2] = ATT$std.error
				Est[i, j, 3] = ATT$conf.low
				Est[i, j, 4] = ATT$conf.high			
			} else {
				EstRR[i, j, 1] = ATT_lnRR$estimate
				EstRR[i, j, 2] = ATT_lnRR$std.error
				EstRR[i, j, 3] = ATT_lnRR$conf.low
				EstRR[i, j, 4] = ATT_lnRR$conf.high	
				
				EstOR[i, j, 1] = ATT_lnOR$estimate
				EstOR[i, j, 2] = ATT_lnOR$std.error
				EstOR[i, j, 3] = ATT_lnOR$conf.low
				EstOR[i, j, 4] = ATT_lnOR$conf.high							
			}
			if(is.m) {
				lstm[[i]][[j]] = m
			} 			
		}	
	}
	
	if( (n1*n2) %% 2 == 1) {  # 奇数
		if(length(unique(dat$Y)) != 2) {
			med = median(Est[ , , 1])
			idx = which(med == Est, arr.ind=TRUE)  # 返回行列索引
			se = Est[idx[1], idx[2], 2]
			low = Est[idx[1], idx[2], 3]
			high = Est[idx[1], idx[2], 4]
			m = lstm[[idx[1]]][[idx[2]]]			
		} else {
			med = c(median(EstRR[ , , 1]), median(EstOR[ , , 1]))
			idxRR = which(med[1] == EstRR, arr.ind=TRUE)  # 返回行列索引
			idxOR = which(med[2] == EstOR, arr.ind=TRUE)  # 返回行列索引
			se = c(EstRR[idxRR[1], idxRR[2], 2], EstOR[idxOR[1], idxOR[2], 2])
			low = c(EstRR[idxRR[1], idxRR[2], 3], EstOR[idxOR[1], idxOR[2], 3])
			high = c(EstRR[idxRR[1], idxRR[2], 4], EstOR[idxOR[1], idxOR[2], 4])
			m = lstm[[idxRR[1]]][[idxRR[2]]]				
		}
	} else {  # 偶数
		if(length(unique(dat$Y)) != 2) {
			med = median(Est[ , , 1])
			se = median(Est[ , , 2])
			low = median(Est[ , , 3])
			high = median(Est[ , , 4])
		} else {
			med = c(median(EstRR[ , , 1]), median(EstOR[ , , 1]))
			se = c(median(EstRR[ , , 2]), median(EstOR[ , , 2]))
			low = c(median(EstRR[ , , 3]), median(EstOR[ , , 3]))
			high = c(median(EstRR[ , , 4]), median(EstOR[ , , 4]))				
		}
		m = NULL		
	}
		
	if(length(unique(dat$Y)) != 2) {
		res = c(med, se, low, high)		
	} else {
		res = cbind(med, se, low, high)
		rownames(res) = paste('kdrm', c('lnRR', 'lnOR'), sep='_')
	}
	if(is.m) {
		return(list(res=res, m=m))	
	} else {
		return(res)
	}	
}

# system.time({ kdrm = kdrm_ATT(dat, d=2, iter=1, replace=TRUE, formula.select=1) }); kdrm


#################################################################################################

estimate_ATT = function(method, d=2, replace=TRUE, B=500, N=1000, MMD=list(), interval=1) {

  future::plan(future::multisession)
  ATT = future_lapply(1:B, future.seed=TRUE, FUN=function(x) {	
	set.seed(x + 99)
	dat = simdata_synthetic(N = N)  #
	if(method == 'euclidean') {
		res = euc_ATT(dat, replace=replace, formula.select=1)
		
	} else if(method == 'mahalanobis') {
		res = mah_ATT(dat, replace=replace, formula.select=1)

	} else if(method == 'psm') {
		res = psm_ATT(dat, replace=replace, formula.select=1)  

	} else if(method == 'pca') {
		res = pca_ATT(dat, d=d, replace=replace, formula.select=1)  
		
	} else if(method == 'kdrm') {
		res = kdrm_ATT(dat, d=d, iter=x, replace=replace, formula.select=1, MMD=MMD, interval=interval)  

	} else {
		stop('Method input error')
	}
   
	return(res)   
  
  })
  future::plan(future::sequential)
  
  ATT = do.call(rbind, ATT) 
  MSE = mean((ATT[ , 1] - 1)^2)  # 均方误差
  SE = mean(ATT[ , 2], na.rm=T)
  return(data.frame(method=method, d=d, MSE=MSE, SE=SE))

}

# system.time({ ATT = estimate_ATT(method='kdrm', d=2, replace=TRUE, B=20, MMD=MMD) }); ATT


error_ATT = function(method, d=2, replace=TRUE, B=500, N=2000, MMD=list(), interval=1) {

  future::plan(future::multisession)
  ATT = future_lapply(1:B, future.seed=TRUE, FUN=function(x) {	
	set.seed(x + 99)
	dat = simdata_causal(N = N)  # 
	if(method == 'euclidean') {
		res = euc_ATT(dat, replace=replace, formula.select=1)
		
	} else if(method == 'mahalanobis') {
		res = mah_ATT(dat, replace=replace, formula.select=1)

	} else if(method == 'psm') {
		res = psm_ATT(dat, replace=replace, formula.select=1)  

	} else if(method == 'pca') {
		res = pca_ATT(dat, d=d, replace=replace, formula.select=3)  
		
	} else if(method == 'kdrm') {
		res = kdrm_ATT(dat, d=d, iter=x, replace=replace, formula.select=3, MMD=MMD, interval=interval)  

	} else {
		stop('Method input error')
	}
   
	return(res)   
  
  })
  future::plan(future::sequential)
  
  ATT = do.call(rbind, ATT) 
  MSE = mean((ATT[ , 1] - 2)^2)  # 均方误差
  RMSE = sqrt(MSE)  # 均方根误差
  MAE = mean(abs(ATT[ , 1] - 2))  # 平均绝对误差  
  MAPE = mean(abs(ATT[ , 1] - 2)/2)  # 平均绝对百分比误差
  SMAPE = mean( abs(ATT[ , 1] - 2)/(abs(ATT[ , 1])/2 + 1) )  # 对称平均绝对百分比误差
  SE = mean(ATT[ , 2], na.rm=T)
  return(data.frame(method=method, d=d, MSE=MSE, RMSE=RMSE, MAE=MAE, MAPE=MAPE, SMAPE=SMAPE, SE=SE))

}

# system.time({ ATT = error_ATT(method='kdrm', d=2, replace=FALSE, B=20, MMD=MMD) }); ATT

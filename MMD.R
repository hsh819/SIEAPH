###MMD计算与优化
###MMD.R


###############################################################################################
gaussian_kernel = function(control, treat, kernel_mul = 2.0, kernel_num = 5, fix_sigma = NULL) {
  stopifnot(kernel_num %% 1 == 0) # 检查是否为整数
  stopifnot(kernel_num > 0) # 检查是否为正数
  
  control = as.matrix(control)
  treat = as.matrix(treat)
  N = nrow(control) + nrow(treat)
  X = rbind(control, treat)
  sigma = sd(X)
  dist_sq = as.matrix(dist(X))^2
  
  if (!is.null(fix_sigma)) {
    bandwidth = fix_sigma
  } else {
    bandwidth = sum(dist_sq)/(N*N - N)
  }
  
  bandwidth = bandwidth / (kernel_mul ^ (kernel_num %/% 2))
  bandwidth_list = sapply(0:(kernel_num - 1), function(i) bandwidth * (kernel_mul ^ i))
  kernel_val = lapply(bandwidth_list, function(h) exp(- dist_sq / h))  # exp(- dist_sq / (2*h^2))
  
  return(list(kernel_val=Reduce("+", kernel_val), bandwidth=median(bandwidth_list)))
}

constant_M = function(control, treat) {
  NC = nrow(control)
  NT = nrow(treat)
  
  M11 = matrix(rep(1/NC^2, NC^2), nrow=NC)
  M12 = matrix(rep(-1/(NC*NT), NC*NT), nrow=NC)
  M21 = matrix(rep(-1/(NT*NC), NT*NC), nrow=NT)
  M22 = matrix(rep(1/NT^2, NT^2), nrow=NT)
  M = rbind(cbind(M11, M12), cbind(M21, M22))
  
  return(M)
}

MK_MMD = function(control, treat, kernel_mul = 2.0, kernel_num = 5, fix_sigma = NULL) {

  kernels = gaussian_kernel(control=control, treat=treat, kernel_mul=kernel_mul, kernel_num=kernel_num, fix_sigma=fix_sigma)
  M = constant_M(control=control, treat=treat) 
  KM = as.matrix(kernels$kernel_val) %*% as.matrix(M)
  MMD = sum(diag(KM)) / kernel_num 
  
  return(cbind(MMD=MMD, bandwidth=kernels$bandwidth))
}

TQ = function(control, treat, Q) {
  controlQ = as.matrix(control) %*% Q
  treatQ = as.matrix(treat)
  
  return(list(control=controlQ, treat=treatQ))
}

QR = function(K, M, D) {
  KM = as.matrix(K) %*% as.matrix(M)
  eig = eigen(KM)
  A = eig$vectors[, 1:D]
  if(nrow(A) == ncol(A)) {
     Q = A
  } else {
	  Q = svd(A)$v
	  if( is.complex(Q) ) {
		M = matrix(rnorm(D^2), D, D)
		Q = svd(M)$v
	  }  
  }

  return(Q)
}

data_cut = function(control, treat, m) {

	# 参数
	NC = nrow(control)
	NT = nrow(treat)
	D = ncol(control)
	nC = floor(nrow(control) / m) # 每份大小
	nT = floor(nrow(treat) / m)

	# 检查
	stopifnot( nC >= D)  
	stopifnot( nT >= D)
	stopifnot( D == ncol(treat))  	

	# 随机排列原矩阵行号
	rC = sample(seq_len(NC))
	rT = sample(seq_len(NT))
	
	# 将行号分成m组
	chunks_C = split(rC, cut(rC, breaks = m, labels = FALSE))
	chunks_T = split(rT, cut(rT, breaks = m, labels = FALSE))	

	# 根据行号分割矩阵 
	XC = lapply(chunks_C, function(i) control[i, ])
	XT = lapply(chunks_T, function(i) treat[i, ])
	
	return(list(XC=XC, XT=XT))
}

iter_MMD = function(DX, Q, D, m, niter=20, eps=10-3) {
	
	res = replicate(m, list())
	for(i in 1:m)
	{
	  sub_control = DX$XC[[i]]
	  sub_treat = DX$XT[[i]]
	  CTQ = TQ(sub_control, sub_treat, Q)

	  kernels = gaussian_kernel(control=CTQ$control, treat=CTQ$treat, kernel_num=1)
	  M = constant_M(control=CTQ$control, treat=CTQ$treat)
	  Q = QR(K=kernels$kernel_val, M=M, D=D)				  
	  CTQ = TQ(sub_control, sub_treat, Q)
	  res[[i]]$Q = Q
	  res[[i]]$MMD = MK_MMD(control=CTQ$control, treat=CTQ$treat)			  		  
	}	
	mat = sapply(res, function(x) x$MMD[, 1])
	min1 = min(mat)
	Q = res[[which.min(mat)]]$Q
	flag = TRUE
	iter = 1
	while(flag) 
	{
		for(i in 1:m)
		{
		  sub_control = DX$XC[[i]]
		  sub_treat = DX$XT[[i]]
		  CTQ = TQ(sub_control, sub_treat, Q)

		  kernels = gaussian_kernel(control=CTQ$control, treat=CTQ$treat, kernel_num=1)
		  M = constant_M(control=CTQ$control, treat=CTQ$treat)
		  Q = QR(K=kernels$kernel_val, M=M, D=D)				  
		  CTQ = TQ(sub_control, sub_treat, Q)
		  res[[i]]$Q = Q
		  res[[i]]$MMD = MK_MMD(control=CTQ$control, treat=CTQ$treat)			  		  
		}	
		mat = sapply(res, function(x) x$MMD[, 1])
		min2 = min(mat)		
		if( min2 < min1 ) {
		  if(abs(min2 - min1) < eps) {
		     flag = FALSE			   
		   } else {
		     min1 = min2
		   }		
		} else if( iter > niter) {		
		     flag = FALSE			
		}
		iter = iter + 1	
	}
	
	return(res)
}

min_MMD = function(control, treat, m = 5) {
# control: 控制组协变量矩阵
# treat: 处理组协变量矩阵
# m: 分割数量
     	
	# 数据分割
	if(m > 1) {
	  DX = data_cut(control=control, treat=treat, m=m)    
	} else {
	  DX = list(XC=list(control), XT=list(treat))
	}
	
	# 初始化Q 
	D = ncol(control)	
	M0 = matrix(rnorm(D^2), nrow=D)
	M1 = qr(M0)	
	Q = qr.Q(M1)	
	res = iter_MMD(DX=DX, Q=Q, D=D, m=m)

	# 返回MMD最小值位置
	mat = sapply(res, function(x) x$MMD[, 1])
	min_idx = which(mat == min(mat), arr.ind = TRUE)

	# 获取最小值
	min_Q = res[[min_idx]]$Q
	min_MMD = res[[min_idx]]$MMD

	list( min_Q = min_Q, min_MMD = min_MMD[1], min_bandwidth = min_MMD[2] )
}

# Example
# set.seed(123) 
# NC = 500
# NT = 400
# XC = replicate(20, rlnorm(NC, meanlog = -0.6, sdlog = 0.15))
# XT = replicate(20, rbeta(NT, shape1 = 1, shape2 = 10))


# MMD = sapply(1:10, function(i) {
  # set.seed(123)
  # M = min_MMD(XC, XT, m=i) 
  # c(MMD=M$min_MMD, bandwidth=M$min_bandwidth)

# })
# idx = which.min(MMD[1, ])

# set.seed(123)
# M = min_MMD(XC, XT, m=idx)


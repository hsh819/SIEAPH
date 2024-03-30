###混合物生存分析
###SurvMixture.R




####################################################数据预处理#############################################
setwd("E:/BaiduNetdiskDownload/TAP数据")
setwd("E:/TAPdata")
cohort = read.csv("cohort_TAP_cokrige.csv", header=TRUE)
cohort_pressure = cohort[, c("survpressureday","survpressure","TAP.cokrige.PM2.5.3","agg.cokrige.PM10.3","agg.cokrige.NO2.3","agg.cokrige.SO2.3","TAP.cokrige.O3.3","agg.cokrige.CO.3",                                             						                       
                        "age","BMI","sex","marriage","education","history.pressure","cook.group","duration.group","smoking","exercise","mask","cleaner",  
                        "diastolic","systolic",'survpressurebasetime',"pressure","UnitID","region","retire")]

dat = subset(cohort_pressure, pressure == 0 & (is.na(diastolic) | diastolic < 90) & (is.na(systolic) | systolic < 140)  & retire == 0 & age >= 18 & age <= 65 & !is.na(survpressurebasetime)) # 筛选数据集
# dat = subset(dat, region != 1)						

update_data = function(flag = TRUE, num = 20)
## flag = TRUE，表示删除数据集重复行；flag = FALSE，表示删除暴露变量重复行
## num表示变量的个数，介于11~20之间
{
	if(num > 20 | num < 11)
	{
	  stop("Error: num must be between 11 and 20")
	}
	
	# 去掉缺失值
	dat = na.omit(dat[, 1:num])   	
	if(flag == TRUE)
	{
		dupidx = which(duplicated(dat))  # 数据集重复行下标
		if( length(dupidx) > 0)
		{
			dat = dat[-dupidx, ]  # 去掉重复行
		}
	} else {
		dupidx = which(duplicated(dat[, 3:8]))  # 暴露变量重复行下标
		if( length(dupidx) > 0)
		{
			dat = dat[-dupidx, ]  # 去掉重复行
		}
	}

	if(num >= 19)
	{
		##修改mask类别
		dat$mask[dat$mask == 1] = 0  #不戴或偶尔戴口罩
		dat$mask[dat$mask == 2] = 1  #经常戴口罩
	}

	##单位转换
	colnames(dat)[3:8] = c('PM2.5', 'PM10', 'NO2', 'SO2', 'O3', 'CO')
	dat[, 3:8] = dat[, 3:8] * 10  #单位转换为μg/m³
	dat$NO2 = dat$NO2 / 1.914  #单位转换为ppb
	dat$SO2 = dat$SO2 / 2.660  #单位转换为ppb
	dat$O3 = dat$O3 / 1.9957   #单位转换为ppb, part per billion	
	dat$CO = dat$CO / 1.165  #单位转换为ppb	
	
	# 排除异常样本
	dat[, 3:10] = scale(dat[, 3:10])  # 暴露变量和连续协变量标准化
	outidx = lapply(dat[, 3:10], function(x) { idx = which(x < -3 | x > 3) })
	outidx = do.call(c, outidx)
	outidx = unique(outidx)  # 排除3倍标准差之外样本
	dat = dat[-outidx, ]

	# # 排除出现次数少的类别
	# ID = as.character(dat$UnitID)  # 随机效应
	# tab = table(ID)
	# IDname = names(tab)[tab < 5]
	# IDidx = which(ID %in% IDname)  # 排除出现次数小于5的单位
	# dat = dat[-IDidx, ]
	
	# 离散协变量转虚拟变量
	dat[, 11:num] = lapply(dat[, 11:num], as.factor)  # 离散协变量
	dummy = model.matrix(~., data=dat[, 11:num])  #将分类变量转化为虚拟变量, 去掉了缺失值
	dummy = dummy[, -1]

	colnames(dat)[1:2] = c('time', 'status')
		
	return(list(dat=dat, dummy=dummy))
}

lst = update_data(flag=FALSE, num=14)  # 选择离散变量sex, marriage, education, history.pressure


## 分析数据集
pollutants = c("PM2.5", "PM10", "NO2", "SO2", "O3", "CO")
# covarnames = names(lst$dat)[9:ncol(lst$dat)]
covarnames = c('age', 'BMI', colnames(lst$dummy))

X = as.matrix(lst$dat[, pollutants])  # 暴露变量
time = lst$dat$time  # 生存时间
status = lst$dat$status  # 生存状态
Z = as.matrix(cbind(lst$dat[, c('age', 'BMI')], lst$dummy))  # 协变量

comp_dat = data.frame(cbind(time=time, status=status, X=X, Z=Z))
	

#################################################数据分析#################################################
f_linear = paste(paste(pollutants, collapse = "+"), paste(covarnames, collapse = "+"), sep="+")
formula = as.formula(paste("Surv(time, status) ~", f_linear))

inter = inter_irf(formula, pollutants, comp_dat, q=10)  # "NO2+_O3+"   "O3+_PM2.5+"  "O3+_PM10+"  "PM10+_PM2.5+" "O3+_SO2+"   "CO+_O3+" 

f_quad = paste(c('I(age^2)', 'I(BMI^2)'), collapse = "+")
f_inter = gsub("\\+_", ":", inter$inter) # 替换中间的"+_"为":"
f_inter = sub("\\+$", "", f_inter)  # 替换最后的"+"为""
f_inter = paste(f_inter, collapse = "+")

# fit1 = qgcomp.cox.noboot(formula2, expnms=pollutants, data=comp_dat, q=10)
# fit2 = qgcomp.cox.boot(formula2, expnms=pollutants, data=comp_dat, q=10, B=10, MCsize=1000, parallel=TRUE, parplan=TRUE, degree=2)


finalvar = select_vars(f_linear, f_quad, f_inter)  # 变量筛选
finalvar = paste(names(finalvar), collapse = "+")
formula2 = as.formula(paste("Surv(time, status) ~", paste(f_linear, finalvar, sep="+")))

finalres = rfqgc_surv(formula2, pollutants, comp_dat, q=10, degree=2)


weights = cbind(pos=finalres$qgc$pos.weights, neg=finalres$qgc$neg.weights, finalres$rfqgc$vmp, finalres$rfqgc$rfpos.weights, finalres$rfqgc$rfneg.weights)
coefs = cbind(finalres$qgc$coefs, rbind(finalres$qgc$cindex, c(NA, NA)) )

write.csv(weights, "final_weights.csv")
write.csv(coefs, "final_coefs.csv")

select_vars = function(f_linear, f_quad, f_inter)
{
   formula = as.formula(paste("Surv(time, status) ~", paste(f_linear, f_quad, f_inter, sep="+")))  
   fit = qgcomp.cox.noboot(formula, expnms=pollutants, data=comp_dat, q=10)
   pvalue = summary(fit$fit)$coef[-c(1:15), 5]
   namevar = names(pvalue)
   flag = FALSE
   while(!flag)
   {
	   if(all(pvalue <= 0.1))
	   {
	     flag = TRUE
	     return(pvalue)
	   } else {
	     
		  mid = which.max(pvalue)
		  namevar = namevar[-mid]
		  temp = paste(namevar, collapse = "+")
		  formula = as.formula(paste("Surv(time, status) ~", paste(f_linear, temp, sep="+") ) )
		  fit = qgcomp.cox.noboot(formula, expnms=pollutants, data=comp_dat, q=10)
		  pvalue = summary(fit$fit)$coef[-c(1:15), 5]		  
	   }     
   }
}




#################################################qgcomp#################################################
library(qgcomp)
library(survival)

qdat = simdata_quantized(
     outcometype="survival", 
     n=1000, corr=c(.9,.3), coef=c(1,1,0,0), 
     q = 10, shape0=0.5, scale0=10, censtime=0.01
)  # 左偏


qdat = simdata_quantized(
     outcometype="survival", 
     n=1000, corr=c(.9,.3), coef=c(1,1,0,0), 
     q = 10, shape0=1, scale0=10, censtime=0.1
)  # 指数分布

qdat = simdata_quantized(
     outcometype="survival", 
     n=1000, corr=c(.9,.3), coef=c(1,1,0,0), 
     q = 10, shape0=2, scale0=10, censtime=0.5
)  # Rayleigh 分布

qdat = simdata_quantized(
     outcometype="survival", 
     n=1000, corr=c(.9,.3), coef=c(1,1,0,0), 
     q = 10, shape0=20, scale0=10, censtime=15
)  # 极值分布


data(metals)
head(metals)
expos = c('arsenic', 'barium', 'cadmium', 'calcium', 'chromium', 'copper', 'iron', 
          'lead', 'magnesium', 'manganese', 'mercury', 'selenium', 'silver', 'sodium','zinc')
covars = c('nitrate', 'nitrite', 'sulfate', 'ph', 'total_alkalinity', 'total_hardness', 'mage35')
time = "disease_time"
status = "disease_state"

formulal = paste("Surv(", time,",", status,")")
formular = paste(paste(expos, collapse = "+"), paste(covars, collapse = "+"), sep="+")
formula = as.formula(paste(formulal, formular, sep="~"))

fit1 = qgcomp.cox.noboot(formula, expnms=expos, data=metals, q=10)
fit1
fit2 = qgcomp.cox.boot(formula, expnms=expos, data=metals, q=10, B=5, MCsize=1000, parallel=TRUE, parplan=TRUE)
fit2

# testing (global) proportional hazards
phtest = survival::cox.zph(fit2$fit)
phtest$table[dim(phtest$table)[1], , drop=FALSE]

plot(fit2, suppressprint = TRUE)

# examining the overall hazard ratio as a function of overall exposure
hrs_q = exp(matrix(0:9, ncol=1, byrow=TRUE) %*% fit2$msmfit$coefficients)
colnames(hrs_q) = "Hazard ratio"
print("Hazard ratios by quartiles (0-10%, 10-20%,..., 90%-100%)")

library(randomForestSRC)
obj = rfsrc(formula, data = metals, importance = "permute", seed=123)
vmp1 = vimp(obj, importance = "permute", seed=123)$importance
vmp2 = vimp(obj, importance="permute", block.size=1, seed=123)$importance
find.interaction(obj, method = "maxsubtree", nvar = 8, seed=123)

smp = subsample(obj, B = 100, importance = "permute")
smp = extract.subsample(smp, raw=TRUE, standardize=FALSE)



#################################################模拟数据分析#################################################
library(qgcomp)
library(randomForestSRC)
library(iRF)
library(future)
library(future.apply)


#####################################相关函数#########################################

rfqgc_surv = function(formula, expos, data, q=NULL, seed=123, degree=1)
{
	fit = qgcomp.cox.noboot(f=formula, expnms=expos, data=data, q=q)
	pos.weights = fit$pos.weights
	neg.weights = fit$neg.weights
	pos.weights = norm_weights(pos.weights, expos)
	neg.weights = norm_weights(neg.weights, expos)
	fit2 = qgcomp.cox.boot(f=formula, expnms=expos, data=data, q=q,  B=100, MCsize=10000, parallel=TRUE, parplan=TRUE, degree=degree)
	coefs = summary(fit2)$coef
	cindex = fit2$fit$concordance[6:7]
	cindex = matrix(cindex, ncol=2, dimnames=list('cindex', names(cindex)))
	
	obj = rfsrc(formula=formula, data=data, importance="permute", seed=seed) 
	smp = subsample(obj, importance="permute")
	smp = extract.subsample(smp, raw=TRUE, standardize=FALSE)
	vmp = rbind(vmp=smp$vmp[1:length(expos)], smp$ci[c(1,5), 1:length(expos)])
	vmp = t(swap(vmp))
	rfpos.weights = rf_weights(smp, expos)$pos.weights
	rfneg.weights = rf_weights(smp, expos)$neg.weights

	list(qgc=list(pos.weights=pos.weights, neg.weights=neg.weights, coefs=coefs, cindex=cindex), 
	     rfqgc=list(vmp=vmp, rfpos.weights=rfpos.weights, rfneg.weights=rfneg.weights))
}


norm_weights = function(weights, expos)
{
  	w = rep(NA, length(expos))
  	names(w) = expos
  	w[names(weights)] = weights
  	return(w)
}

swap = function(dat)
{
	for(i in 1:ncol(dat))
	{
	  if(!is.na(dat[2, i]))
	  {
		if(dat[2, i] > dat[3, i])  # 互换第2和第3行
		{
		   temp = dat[2, i]
		   dat[2, i] = dat[3, i]
		   dat[3, i] = temp
		}
	    if(dat[1, i] < dat[2, i])     
	    {
          dat[2, i] = dat[1, i] - 1e-5  # 修改下界   
      } else if(dat[1, i] > dat[3, i])
      {
          dat[3, i] = dat[1, i] + 1e-5  # 修改上界            
      }		  
	  }
	}
	return(dat)
}

rf_weights = function(smp, expos)
{
	vmp = smp$vmp[1:length(expos)]
	ci = smp$ci[c(1,5), 1:length(expos)]

	if(all(vmp >= 0))
	{
  	pos = vmp[vmp >= 0]/sum(vmp[vmp >= 0])
  	pos = norm_weights(pos, expos)
  
  	pos.lower = ci[1, ]/sum(ci[1, ])
  	pos.upper = ci[2, ]/sum(ci[2, ])
  	pos.lower = norm_weights(pos.lower, expos)
  	pos.upper = norm_weights(pos.upper, expos)
  
  	pos.weights = rbind(pos, pos.lower, pos.upper)
  	pos.weights = swap(pos.weights)
  	pos.weights = t(pos.weights)
  	
  	return(list(pos.weights=pos.weights, neg.weights=matrix(NA, nrow=length(expos), ncol=3)))
	
	} else if(all(vmp < 0))
	{
  	neg = vmp[vmp < 0]/sum(vmp[vmp < 0])
  	neg = norm_weights(neg, expos)
  
  	neg.lower = ci[1, ]/sum(ci[1, ])
  	neg.upper = ci[2, ]/sum(ci[2, ])
  	neg.lower = norm_weights(neg.lower, expos)
  	neg.upper = norm_weights(neg.upper, expos)
  	
  	neg.weights = rbind(neg, neg.lower, neg.upper)
  	neg.weights = swap(neg.weights)
  	neg.weights = t(neg.weights)
  
  	list(pos.weights=matrix(NA, nrow=length(expos), ncol=3), neg.weights=neg.weights)	
	
	} else {
  	pos = vmp[vmp >= 0]/sum(vmp[vmp >= 0])
  	neg = vmp[vmp < 0]/sum(vmp[vmp < 0])
  	pos = norm_weights(pos, expos)
  	neg = norm_weights(neg, expos)
  
  	idx1 = which(vmp >= 0)
  	idx2 = which(vmp < 0)
  
  	pos.lower = ci[1, idx1]/sum(ci[1, idx1])
  	pos.upper = ci[2, idx1]/sum(ci[2, idx1])
  	names(pos.lower) = names(pos.upper) = colnames(ci)[idx1]	
  	pos.lower = norm_weights(pos.lower, expos)
  	pos.upper = norm_weights(pos.upper, expos)
  
  	neg.lower = ci[1, idx2]/sum(ci[1, idx2])
  	neg.upper = ci[2, idx2]/sum(ci[2, idx2])
  	names(neg.lower) = names(neg.upper) = colnames(ci)[idx2]
  	neg.lower = norm_weights(neg.lower, expos)
  	neg.upper = norm_weights(neg.upper, expos)
  
  	pos.weights = rbind(pos, pos.lower, pos.upper)
  	pos.weights = swap(pos.weights)
  	pos.weights = t(pos.weights)
  	
  	neg.weights = rbind(neg, neg.lower, neg.upper)
  	neg.weights = swap(neg.weights)
  	neg.weights = t(neg.weights)
  
  	list(pos.weights=pos.weights, neg.weights=neg.weights)	
	
	}
}

inter_irf = function(formula, expos, data, q=NULL, seed=123)
{
	res = qgcomp.cox.noboot(f=formula, expnms=expos, data=data, q=q)
	data$resid =  res$fit$residuals  # Extracting Residuals
	vars = all.vars(formula)
	vars = setdiff(vars, c("time", "status"))
	n = nrow(data)
	plan(multisession, workers=12)
	irf = future_lapply(X=c(0.60, 0.65, 0.70), FUN=fit_irf, vars=vars, data=data, n=n, seed=seed, future.seed=TRUE)	
	# irf = lapply(X=c(0.70, 0.75, 0.80), FUN=fit_irf, expos=expos, data=data, n=n, seed=seed)	
	inter = intersect(intersect(irf[[1]]$int, irf[[2]]$int), irf[[3]]$int)
	return(list(irf=irf, inter=inter))

}

fit_irf = function(ntrain, vars, data, n, seed)
{
	set.seed(seed)
	train.id = sample(seq(1,n), ceiling(n*ntrain))
	if(ntrain == 1)
	{
		test.id = train.id
	} else {	  
		test.id = setdiff(1:n, train.id)
	}

	set.seed(seed)
	fit = iRF(x = data[train.id, vars], 
			  y = data[train.id, "resid"],  # SiRF algorithm with Residuals as the outcomes
			  xtest = data[test.id, vars],
			  ytest = data[test.id, "resid"],
			  n.iter = 10, 
			  n.core = 4,
			  select.iter = TRUE,
			  n.bootstrap = 500
	)

	inter = as.data.frame(fit$interaction[fit$interaction$stability >= 0.5,])
	inter$synergy = rep(NA_real_, nrow(inter))
	for(i in 1:nrow(inter)) { inter$synergy[i] = sum(strsplit(inter$int,"+")[[i]] %in% "-") }
	final = inter[inter$synergy == 0, ]
	return(final)
}


####################################正确识别交互作用####################################
expos = paste0('x', 1:4)
formula = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4)
formula2 = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + I(x1^2) + x2:x3)  # 真实模型


fit = qgcomp.cox.noboot(Surv(time, status) ~ x1 + x2 + x3 + x4 + I(x1^2) + x2:x3 , expnms=expos, data=dat1, q=NULL)
fit2 = qgcomp.cox.boot(Surv(time, status) ~ x1 + x2 + x3 + x4 + I(x1^2) + x2:x3 , expnms=expos, data=dat1, q=NULL,  B=10, MCsize=100, parallel=TRUE, parplan=TRUE, degree = 2)


inter1 = inter_irf(formula, expos, dat1)
inter2 = inter_irf(formula, expos, dat2)
inter3 = inter_irf(formula, expos, dat3)
inter4 = inter_irf(formula, expos, dat4)
inter5 = inter_irf(formula, expos, dat5)
inter6 = inter_irf(formula, expos, dat6)
inter7 = inter_irf(formula, expos, dat7)
inter8 = inter_irf(formula, expos, dat8)


res1 = rfqgc_surv(formula2, expos, dat1)
res2 = rfqgc_surv(formula2, expos, dat2)
res3 = rfqgc_surv(formula2, expos, dat3)
res4 = rfqgc_surv(formula2, expos, dat4)
res5 = rfqgc_surv(formula2, expos, dat5)
res6 = rfqgc_surv(formula2, expos, dat6)
res7 = rfqgc_surv(formula2, expos, dat7)
res8 = rfqgc_surv(formula2, expos, dat8)

res = list(res1, res2, res3, res4, res5, res6, res7, res8)

weights = NULL
coefs = NULL

for(i in 1:8)
{
	temp_w = cbind(pos=res[[i]]$qgc$pos.weights, neg=res[[i]]$qgc$neg.weights, res[[i]]$rfqgc$vmp, res[[i]]$rfqgc$rfpos.weights, res[[i]]$rfqgc$rfneg.weights)
	temp_c = cbind(res[[i]]$qgc$coefs, res[[i]]$qgc$cindex)
	weights = rbind(weights, temp_w)
	coefs = rbind(coefs, temp_c)
}

write.csv(weights, "weights1_8.csv")
write.csv(coefs, "coefs1_8.csv")



####################################错误识别交互作用####################################
expos = paste0('x', 1:4)
formula = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4)
formula3 = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x2)  # 真实模型
formula3_f = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x3) 

formula4 = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3)  # 真实模型
formula4_f = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x3 + x3:x4)  

formula5 = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3 + x2:x4)  # 真实模型
formula5_f = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x2 + x1:x3 + x3:x4)

formula6 = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3 + x2:x4 + x3:x4)  # 真实模型
formula6_f = as.formula(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x3 + x2:x3 + x3:x4)


fit3 = qgcomp.cox.noboot(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x3 + x2:x3 + x3:x4 , expnms=expos, data=dat15, q=NULL)
fit4 = qgcomp.cox.boot(Surv(time, status) ~ x1 + x2 + x3 + x4 + x1:x3 , expnms=expos, data=dat1, q=NULL,  B=10, MCsize=100, parallel=TRUE, parplan=TRUE, degree = 2)



res9 = rfqgc_surv(formula3_f, expos, dat9)
res10 = rfqgc_surv(formula3_f, expos, dat10)
res11 = rfqgc_surv(formula4_f, expos, dat11)
res12 = rfqgc_surv(formula4_f, expos, dat12)
res13 = rfqgc_surv(formula5_f, expos, dat13)
res14 = rfqgc_surv(formula5_f, expos, dat14)
res15 = rfqgc_surv(formula6_f, expos, dat15)
res16 = rfqgc_surv(formula6_f, expos, dat16)

res = list(res9, res10, res11, res12, res13, res14, res15, res16)

weights = NULL
coefs = NULL

for(i in 1:8)
{
	temp_w = cbind(pos=res[[i]]$qgc$pos.weights, neg=res[[i]]$qgc$neg.weights, res[[i]]$rfqgc$vmp, res[[i]]$rfqgc$rfpos.weights, res[[i]]$rfqgc$rfneg.weights)
	temp_c = cbind(res[[i]]$qgc$coefs, res[[i]]$qgc$cindex)
	weights = rbind(weights, temp_w)
	coefs = rbind(coefs, temp_c)
}

write.csv(weights, "weights9_16.csv")
write.csv(coefs, "coefs9_16.csv")

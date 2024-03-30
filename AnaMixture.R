###AnaMixture.R


library(gWQS)
library(qgcomp)
library(gcdnet)
library(bkmr)
library(bkmrhat)
library(bsmim2)


####################################################数据预处理#############################################
setwd("E:/TAPdata")
cohort = read.csv("cohort_TAP_cokrige.csv", header=TRUE)
cohort_glucose = cohort[, c("glucose","TAP.cokrige.PM2.5.3","agg.cokrige.PM10.3","agg.cokrige.NO2.3","agg.cokrige.SO2.3","TAP.cokrige.O3.3","agg.cokrige.CO.3",                                             						                       
                        "age","BMI","sex","marriage","education","history.diabetes","cook.group","duration.group","smoking","exercise","mask","cleaner",  
                        "IFG","UnitID","region","diabetes","retire")]

dat = subset(cohort_glucose, diabetes == 0 & retire == 0 & age <= 65 & age >= 18 & !is.na(IFG))  # 筛选数据集
dat = subset(dat, region != 1)
						

update_data = function(flag = TRUE, num = 21)
## flag = TRUE，表示删除数据集重复行；flag = FALSE，表示删除暴露变量重复行
## num表示变量的个数，介于13~21之间
{
	if(num > 21 | num < 13)
	{
	  stop("Error: num must be between 13 and 21")
	}
	m = num - 2
	
	# 去掉region, diabetes, retire, 缺失值
	dat = na.omit(dat[, c(1:m, 20:21)])   	
	if(flag == TRUE)
	{
		dupidx = which(duplicated(dat))  # 数据集重复行下标
		if( length(dupidx) > 0)
		{
			dat = dat[-dupidx, ]  # 去掉重复行
		}
	} else {
		dupidx = which(duplicated(dat[, 2:7]))  # 暴露变量重复行下标
		if( length(dupidx) > 0)
		{
			dat = dat[-dupidx, ]  # 去掉重复行
		}
	}

	if(m >= 18)
	{
		##修改mask类别
		dat$mask[dat$mask == 1] = 0  #不戴或偶尔戴口罩
		dat$mask[dat$mask == 2] = 1  #经常戴口罩
	}

	##单位转换
	colnames(dat)[2:7] = c('PM2.5', 'PM10', 'NO2', 'SO2', 'O3', 'CO')
	dat[, 2:7] = dat[, 2:7] * 10  #单位转换为μg/m³
	dat$NO2 = dat$NO2 / 1.914  #单位转换为ppb
	dat$SO2 = dat$SO2 / 2.660  #单位转换为ppb
	dat$O3 = dat$O3 / 1.9957   #单位转换为ppb, part per billion	
	dat$CO = dat$CO / 1.165  #单位转换为ppb	
	
	# 排除异常样本
	dat[, 1:9] = scale(dat[, 1:9])  # 结局变量、暴露变量和连续协变量标准化
	outidx = lapply(dat[, 1:9], function(x) { idx = which(x < -3 | x > 3) })
	outidx = do.call(c, outidx)
	outidx = unique(outidx)  # 排除3倍标准差之外样本
	dat = dat[-outidx, ]

	# 排除出现次数少的类别
	ID = as.character(dat$UnitID)  # 随机效应
	tab = table(ID)
	IDname = names(tab)[tab < 5]
	IDidx = which(ID %in% IDname)  # 排除出现次数小于5的单位
	dat = dat[-IDidx, ]
	
	# 离散协变量转虚拟变量
	dat[, c(10:m, ncol(dat))] = lapply(dat[, c(10:m, ncol(dat))], as.factor)  # 离散协变量，包括UnitID
	dummy = model.matrix(~., data=dat[, 10:m])  #将分类变量转化为虚拟变量, 去掉了缺失值
	dummy = dummy[, -1]
	colnames(dat)[ncol(dat)] = 'ID'
		
	return(list(dat=dat, dummy=dummy))
}

#lst = update_data(flag=FALSE, num=21)  
lst = update_data(flag=FALSE, num=15)  # 选择离散变量sex, marriage, education, history.diabetes


## 分析数据集
pollutants = c("PM2.5", "PM10", "NO2", "SO2", "O3", "CO")
covarnames = names(lst$dat)[8:(ncol(lst$dat) - 2)]

X = as.matrix(lst$dat[, pollutants])  # 暴露变量
y = lst$dat$glucose  # 连续响应
ybin = as.numeric(as.character(lst$dat$IFG))  # 离散响应
Z = as.matrix(cbind(lst$dat[, c('age', 'BMI')], lst$dummy))  # 协变量
ID = factor(lst$dat$ID)

comp_dat = list(y=y, ybin=ybin, X=X, Z=Z, ID=ID)  


#################################################gWQS#################################################
library(gWQS)
library(lmerTest)
formulas1 = as.formula("glucose ~ wqs")  # 初始模型，连续
formulas2 = as.formula(paste("glucose ~ wqs", paste(covarnames, collapse = "+"), sep="+"))  # 调整模型，连续
formulas3 = as.formula("IFG ~ wqs")  # 初始模型，离散
formulas4 = as.formula(paste("IFG ~ wqs", paste(covarnames, collapse = "+"), sep="+"))  # 调整模型，离散					   

GWQS = function(formulas, mix_name, data, b1_pos=TRUE, is.gauss=TRUE, seed=123)
{
	if(is.gauss == TRUE)
	{
		results = gwqs(formulas, mix_name=mix_name, data=data, q = 10, validation = 0.6, b = 200, b1_pos=b1_pos, 
						b1_constr = TRUE, family=gaussian, seed=seed, plan_strategy="multisession")  # 并行计算
		vindex = results$vindex
		valdata = cbind(results$fit$data, ID=data[vindex, 'ID'])
		valformulas = as.formula(paste(deparse(formulas), "(1|ID)", sep="+"))
		valresults = lmer(valformulas, data=valdata)
		coefs = summary(valresults)$coef[, -3]  # 去掉df	
	
	} else {
		results = gwqs(formulas, mix_name=mix_name, data=data, q = 10, validation = 0.6, b = 200, b1_pos=b1_pos, 
						b1_constr = TRUE, family=binomial, seed=seed, plan_strategy="multisession")  # 并行计算
		vindex = results$vindex
		valdata = cbind(results$fit$data, ID=data[vindex, 'ID'])		
		valformulas = as.formula(paste(deparse(formulas), "(1|ID)", sep="+"))
		valresults = glmer(valformulas, data=valdata)
		coefs = summary(valresults)$coef
		pvalue = 2 * pnorm(-abs(coefs[, 't value']))  # 近似计算p值
		coefs = cbind(coefs, pvalue)
	
	}
	
	final_weights = results$final_weights
	myorder = c("PM2.5", "PM10", "NO2", "SO2", "O3", "CO")
	final_weights$mix_name = factor(final_weights$mix_name, levels=myorder)
	final_weights = final_weights[order(final_weights$mix_name), ]
	return(list(coefs=coefs, final_weights=final_weights))

}

results1.1 = GWQS(formulas1, mix_name=pollutants, data=lst$dat, b1_pos=TRUE, is.gauss=TRUE)  # 初始模型，连续，正约束
results1.2 = GWQS(formulas1, mix_name=pollutants, data=lst$dat, b1_pos=FALSE, is.gauss=TRUE)  # 初始模型，连续，负约束
results1.3 = GWQS(formulas2, mix_name=pollutants, data=lst$dat, b1_pos=TRUE, is.gauss=TRUE)  # 调整模型，连续，正约束
results1.4 = GWQS(formulas2, mix_name=pollutants, data=lst$dat, b1_pos=FALSE, is.gauss=TRUE)  # 调整模型，连续，负约束

results1.5 = GWQS(formulas3, mix_name=pollutants, data=lst$dat, b1_pos=TRUE, is.gauss=FALSE)  # 初始模型，离散，正约束
results1.6 = GWQS(formulas3, mix_name=pollutants, data=lst$dat, b1_pos=FALSE, is.gauss=FALSE)  # 初始模型，离散，负约束
results1.7 = GWQS(formulas4, mix_name=pollutants, data=lst$dat, b1_pos=TRUE, is.gauss=FALSE)  # 调整模型，离散，正约束
results1.8 = GWQS(formulas4, mix_name=pollutants, data=lst$dat, b1_pos=FALSE, is.gauss=FALSE)  # 调整模型，离散，负约束

lst1 = list(results1.1, results1.2, results1.3, results1.4, results1.5, results1.6, results1.7, results1.8)
lst11 = lapply(lst1, function(x) { x$coefs[1:2, ] })
lst12 = lapply(lst1, function(x) { rownames(x$final_weights) = NULL; x$final_weights })
res11 = do.call(rbind, lst11)
res12 = do.call(cbind, lst12)

wqsname = rep(paste(rep(rep(c('crude', 'adjust'), each=2), 2), rep(c('gauss', 'binom'), each=4), rep(c('pos', 'neg'), 4), sep='_'), each=2)
rownames(res11) = paste(wqsname, rownames(res11), sep='_')
colnames(res12) = paste(wqsname, colnames(res12), sep='_')

write.csv(res11, file="gwqs_coefs.csv")
write.csv(res12, file="gwqs_weights.csv", row.names=FALSE)


##############没有随机效应gwqs##############
reorder_weights = function(results)
{
	final_weights = results$final_weights
	myorder = c("PM2.5", "PM10", "NO2", "SO2", "O3", "CO")
	final_weights$mix_name = factor(final_weights$mix_name, levels=myorder)
	final_weights = final_weights[order(final_weights$mix_name), ]
	return(final_weights)
}

results1.1 = gwqs(formulas2, mix_name=pollutants, data=lst$dat, q = 10, validation = 0.6, b = 100, b1_pos=TRUE, 
						b1_constr = TRUE, family=gaussian, seed=123, plan_strategy="multisession")  #连续，调整模型，正约束 
wqs_weight1 = cbind(summary(results1.1$fit)$coef[2, , drop=FALSE], t(reorder_weights(results1.1)[, 2, drop=FALSE]))


results1.2 = gwqsrh(formulas2, mix_name=pollutants, data=lst$dat, q = 10, validation = 0.6, b = 100, b1_pos=TRUE, rh=20,
						b1_constr = TRUE, family=gaussian, seed=123, plan_strategy="multisession")  #连续，调整模型，正约束 
wqs_weight2 = cbind(results1.2$fit$coef[2, 1:4, drop=FALSE], t(reorder_weights(results1.2)[, 2, drop=FALSE]))


results1.3 = gwqs(formulas2, mix_name=pollutants, data=lst$dat, q = 10, validation = 0, b = 100, b1_pos=TRUE, 
						b1_constr = TRUE, family=gaussian, seed=123, plan_strategy="multisession")  #连续，调整模型，正约束 
wqs_weight3 = cbind(summary(results1.3$fit)$coef[2, , drop=FALSE], t(reorder_weights(results1.3)[, 2, drop=FALSE]))

wqs_weights = rbind(wqs_weight1, wqs_weight2, wqs_weight3)
rownames(wqs_weights) = c("rh1", "rh100", "valid0")

write.csv(t(wqs_weights), "wqs_coef_weights.csv")	


##############随机效应gwqs##############
source("gwqs_re.R")
source("gwqs_related.R")
results1.4 = gwqs_re(formulas2, mix_name = pollutants, data = lst$dat, q = 10, validation = 0.6, b = 100, 
                b1_pos = TRUE, b1_constr = TRUE, family = gaussian, id = 'ID', seed = 123, plan_strategy = "multisession")  #连续，调整模型，正约束
pval1 = 2*pnorm(-abs(summary(results1.4$fit)$coef[, 3]))
wqsre_weight1 = cbind(summary(results1.4$fit)$coef[2, , drop=FALSE], pval=pval1[2], t(results1.4$final_weights[, 2, drop=FALSE]))


gwqsre_rh = function(i)
{
	res = gwqs_re(formulas2, mix_name = pollutants, data = lst$dat, q = 10, validation = 0.6, b = 100, 
					b1_pos = TRUE, b1_constr = TRUE, family = gaussian, id = 'ID', seed = 100+i, plan_strategy = "multisession")  #连续，调整模型，正约束
	final_weights = res$final_weights
	coefs =  summary(res$fit)$coef[1:2, ]
	return(list(coefs=coefs, final_weights=final_weights))
}

future::plan("multisession", workers=10)
lst1 = future_lapply(X = 1:20, FUN = gwqsre_rh, future.seed = TRUE)

lst11 = lapply(lst1, function(x) { x$coefs[2, , drop=FALSE]  })
lst12 = lapply(lst1, function(x) { x$final_weights[, 2, drop=FALSE]  })
lst11 = do.call(rbind, lst11)
lst12 = do.call(cbind, lst12)
pval2 = 2*pnorm(-abs(lst11[, 3]))
wqsre_weight2 = c(apply(lst11, 2, mean), pval=mean(pval2), apply(lst12, 1, mean))


results1.6 = gwqs_re(formulas2, mix_name = pollutants, data = lst$dat, q = 10, validation = 0, b = 100, 
                   b1_pos = TRUE, b1_constr = TRUE, family = gaussian, id = 'ID', seed = 123, plan_strategy = "multisession")  #连续，调整模型，正约束
pval3 = 2*pnorm(-abs(summary(results1.6$fit)$coef[, 3]))
wqsre_weight3 = cbind(summary(results1.6$fit)$coef[2, , drop=FALSE], pval=pval3[2], t(results1.6$final_weights[, 2, drop=FALSE]))				   

wqsre_weights = rbind(wqsre_weight1, wqsre_weight2, wqsre_weight3)
rownames(wqsre_weights) = c("rh1", "rh100", "valid0")

write.csv(t(wqsre_weights), "wqsre_coef_weights.csv")			   


# results1.7 = gwqs_re(formulas2, mix_name = pollutants, data = lst$dat, q = 10, validation = 0.6, b = 100, 
                # b1_pos = FALSE, b1_constr = TRUE, family = gaussian, id = 'ID', seed = 123, plan_strategy = "multisession")  #连续，调整模型，负约束
# wqsre_weight4 = cbind(summary(results1.7$fit)$coef[2, , drop=FALSE], t(results1.7$final_weights[, 2, drop=FALSE]))				

				

#################################################qgcomp#################################################
library(qgcomp)
formulas5 = as.formula(paste("glucose ~", paste(pollutants, collapse = "+")))  # 初始模型，连续
formulas6 = as.formula(paste("glucose ~", paste(paste(pollutants, collapse = "+"), paste(covarnames, collapse = "+"), sep="+")))  # 调整模型，连续
formulas7 = as.formula(paste("IFG ~", paste(pollutants, collapse = "+")))  # 初始模型，离散 
formulas8 = as.formula(paste("IFG ~", paste(paste(pollutants, collapse = "+"), paste(covarnames, collapse = "+"), sep="+")))  # 调整模型，离散


QGCOMP = function(formulas, expnms, data, family=gaussian, is.gauss=TRUE, is.boot=FALSE)
{
   
   if(is.boot)   
   {   
     results = qgcomp.glm.boot(formulas, expnms = expnms, data=data, family=family, q=10, id='ID', B=200, rr=FALSE, seed=123, parallel=TRUE, parplan=TRUE)  # 并行计算
  	 if(is.gauss)
  	 {
  	   coefs = summary(results)$coef
  	 } else {
  	   coefs = summary(results)$coef[, -5]	 # 去掉Z value列
  	 }
       return(list(coefs=coefs, pos.weights=NULL, neg.weights=NULL))	    
   } else {
     results = qgcomp.glm.noboot(formulas, expnms = expnms, data=data, family=family, q=10)
  	 newdata = results$fit$data	
  	 newformulas = as.formula(paste(paste(deparse(formulas), collapse = ""), "+ (1|ID)"))
         
     if(is.gauss)
  	 {
       newresults = lmer(newformulas, data=newdata)
       coeff = summary(newresults)$coef[expnms, 'Estimate']
       pos.weights = coeff[coeff >= 0]/sum(coeff[coeff >= 0]) 
       neg.weights = coeff[coeff < 0]/sum(coeff[coeff < 0])
       coefs = summary(results)$coef	   

  	 } else {
  	   newresults = glmer(newformulas, data=newdata)
  	   coeff = summary(newresults)$coef[expnms, 'Estimate']
  	   pos.weights = coeff[coeff >= 0]/sum(coeff[coeff >= 0]) 
  	   neg.weights = coeff[coeff < 0]/sum(coeff[coeff < 0])
  	   coefs = summary(results)$coef[, -5]	 # 去掉Z value列

  	 }
  	 pos = c(PM2.5=NA, PM10=NA, NO2=NA, SO2=NA, O3=NA, CO=NA)
  	 neg = c(PM2.5=NA, PM10=NA, NO2=NA, SO2=NA, O3=NA, CO=NA)
  	 pos[names(pos.weights)] = pos.weights
  	 neg[names(neg.weights)] = neg.weights
  	 return(list(coefs=coefs, pos.weights=pos, neg.weights=neg))   
   }   
   
}

results2.1 = QGCOMP(formulas5, expnms = pollutants, data=lst$dat, family=gaussian, is.gauss=TRUE, is.boot=FALSE)  # 初始模型，连续，不进行重抽样
results2.2 = QGCOMP(formulas6, expnms = pollutants, data=lst$dat, family=gaussian, is.gauss=TRUE, is.boot=FALSE)  # 调整模型，连续，不进行重抽样
results2.3 = QGCOMP(formulas7, expnms = pollutants, data=lst$dat, family=binomial, is.gauss=FALSE, is.boot=FALSE)  # 初始模型，离散，不进行重抽样
results2.4 = QGCOMP(formulas8, expnms = pollutants, data=lst$dat, family=binomial, is.gauss=FALSE, is.boot=FALSE)  # 调整模型，离散，不进行重抽样

results2.5 = QGCOMP(formulas5, expnms = pollutants, data=lst$dat, family=gaussian, is.gauss=TRUE, is.boot=TRUE)  # 初始模型，连续，进行重抽样
results2.6 = QGCOMP(formulas6, expnms = pollutants, data=lst$dat, family=gaussian, is.gauss=TRUE, is.boot=TRUE)  # 调整模型，连续，进行重抽样  
results2.7 = QGCOMP(formulas7, expnms = pollutants, data=lst$dat, family=binomial, is.gauss=FALSE, is.boot=TRUE)  # 初始模型，离散，进行重抽样
results2.8 = QGCOMP(formulas8, expnms = pollutants, data=lst$dat, family=binomial, is.gauss=FALSE, is.boot=TRUE)  # 调整模型，离散，进行重抽样

lst2 = list(results2.1, results2.2, results2.3, results2.4, results2.5, results2.6, results2.7, results2.8)
lst21 = lapply(lst2, function(x) { x$coefs })
lst22 = lapply(lst2, function(x) { x$pos.weights })
lst23 = lapply(lst2, function(x) { x$neg.weights })

res21 = do.call(rbind, lst21)
qgcname1 = rep(paste( rep(c('crude', 'adjust'), 4), rep(rep(c('gauss', 'binom'), each=2), 2), rep(c('noboot', 'boot'), each=4), sep='_'), each=2)
rownames(res21) = paste(qgcname1, rownames(res21), sep='_')

qgcname2 = paste(rep(c('crude', 'adjust'), 2), rep(c('gauss', 'binom'), each=2), rep('noboot', 4), sep='_')
pos_name = paste(rep(qgcname2, sapply(lst22, length)[1:4]), names(unlist(lst22)), sep='_pos_')
neg_name = paste(rep(qgcname2, sapply(lst23, length)[1:4]), names(unlist(lst23)), sep='_neg_')
res22 = data.frame(pos_name = pos_name, pos_weights = unlist(lst22), neg_name=neg_name, neg_weights = unlist(lst23)) 
                   
write.csv(res21, file="qgcomp_coefs.csv")
write.csv(res22, file="qgcomp_weights.csv", row.names=FALSE)



#################################################bkmr#################################################
# remotes::install_github("jenfb/bkmr")
# remotes::install_github("alexpkeil1/bkmrhat", build_vignettes = TRUE)
library(bkmr)
library(bkmrhat)
library(ggplot2)
library(ggsci)
library(cowplot)


############多链并行估计############################
## 参数设置
R = 10000             ## no. of iterations
burn = 0.5           ## percent burn-in
thin = 40            ## thinning number
sel = seq(burn * R + 1, R, by=thin) 


# set.seed(R)
# samp = sort(sample(length(comp_dat$y), seed))  # 抽取子样本
# kmdat = with(comp_dat, list(y=y[samp], ybin=ybin[samp], Z=X[samp, ], X=Z[samp, ]))  # 暴露变量用z表示，调整协变量用x表示
kmdat = with(comp_dat, list(y=y, ybin=ybin, Z=X, X=Z, ID=ID))  # 暴露变量用z表示，调整协变量用x表示

ncores = 10
future::plan(strategy = future::multisession, workers=ncores)  # 设置并行策略，strategy='sequential', 'multisession', 'multicore', 'cluster' 
start = proc.time()
set.seed(R)
fitkm = kmbayes_parallel(nchains=ncores, y=kmdat$y, Z=kmdat$Z, X=kmdat$X, id=kmdat$ID, iter = ceiling(R/ncores), verbose = TRUE, varsel = TRUE, control.params = list(r.jump2 = 0.5))  # 连续响应，包括协变量，提高M-H算法接收率 
# fitpr = kmbayes_parallel(nchains=ncores, y=kmdat$ybin, Z=kmdat$Z, X=kmdat$X, id=kmdat$ID, iter = ceiling(R/ncores), verbose = TRUE, varsel = TRUE, family="binomial", control.params = list(r.jump2 = 0.5))  # 离散响应，包括协变量，提高M-H算法接收率 


diftime = proc.time() - start
print(paste("Execution time:", round(diftime[3]/3600,2), "hours"))

#########################################################
# load(file="BKMRHAT_fitkmcomb_10000.RData")
load(file="BKMRHAT_fitkmcomb_20000.RData")

## 参数设置
R = 20000             ## no. of iterations
burn = 0.5           ## percent burn-in
thin = 40            ## thinning number
sel = seq(burn * R + 1, R, by=thin) 

## 诊断
multidiag = kmbayes_diagnose(fitkm, warmup=0, digits_summary=2)
# multidiag = kmbayes_diagnose(fitkmmore, warmup=0, digits_summary=2)

## 后验汇总与推断
fitkmcomb = kmbayes_combine(fitkm)
# fitkmcomb = kmbayes_combine(fitkmmore)
summary(fitkmcomb)

## 检查模型收敛性
TracePlot(fit = fitkmcomb, par = "beta", comp=1)  # beta1,...,beta19
TracePlot(fit = fitkmcomb, par = "lambda")
TracePlot(fit = fitkmcomb, par = "sigsq.eps")
TracePlot(fit = fitkmcomb, par = "r", comp = 1)  # r1,...,r6 

## 后验包含概率
# multipips = lapply(fitkm, function(x) t(ExtractPIPs(x)))
pips = ExtractPIPs(fitkmcomb)  # variable selection

## 绘制暴露-响应函数
## i.单变量横截面
pred.resp.univar = PredictorResponseUnivar(fit = fitkmcomb, method = "exact", sel = sel)

tiff(file=paste0("BKMR_univar_gauss_", R, ".tif"), width=12, height=8, units="in", compression="lzw", res=144, family="sans")
# label = c(expression("PM"[2.5]), expression("PM"[10]), expression("NO"[2]), expression("SO"[2]), expression("O"[3]), "CO")
# names(label) = levels(pred.resp.univar$variable)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se))  + theme_bw() + 
    geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
    facet_wrap(~variable, ncol = 3) + 
    # labs(title = "Univariate exposure-response function") + 
    # theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2)) +	
    labs(x = "Exposure", y = "h(Exposure)")
dev.off()


## ii.两变量横截面(等高线)
expos.pairs = subset(data.frame(expand.grid(expos1 = 1:6, expos2 = 1:6)), expos1 < expos2)
expos.pairs
pred.resp.bivar = PredictorResponseBivar(fit = fitkmcomb, min.plot.dist = 1, z.pairs = expos.pairs, sel = sel)

tiff(file=paste0("BKMR_bivar_gauss_", R, ".tif"), width=10, height=8, units="in", compression="lzw", res=144, family="sans")
ggplot(pred.resp.bivar, aes(z1, z2, fill = est)) + geom_raster() + theme_bw() + facet_grid(variable2 ~ variable1) +
    scale_fill_gradientn(colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF")) +
    xlab("Exposure 1") + ylab("Exposure 2") + labs(title="h(Exposure 1, Exposure 2)", fill="Estimate") +
    theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))
dev.off()

# Because we specified min.plot.dist = 0.5 as an argument in the PredictorResponseBivar function, the exposure-response surface is only estimated for points that are within 0.5 units from an observed data point. 
# Points farther than this distance are grayed out in the plot.
# Since it can be hard to see what’s going on in these types of plots, an alternative approach is to investigate the exposure-response function of a single exposure where the second exposure is fixed at various quantiles. 
# This can be done using the PredictorResponseBivarLevels function, which takes as input the bivariate exposure-response function outputted from the previous command, where the argument qs specifies a sequence of quantiles at which to fix the second exposure

## iii.两变量横截面(折线图)
pred.resp.bivar.levels = PredictorResponseBivarLevels(pred.resp.bivar, Z=fitkmcomb$Z, qs = c(0.10, 0.50, 0.90))
myorder = c("PM2.5", "PM10","NO2", "SO2", "O3", "CO")
pred.resp.bivar.levels[, 1:2] = lapply(pred.resp.bivar.levels[, 1:2], function(x) factor(x, levels = myorder))
pred.resp.bivar.levels = pred.resp.bivar.levels[order(pred.resp.bivar.levels$variable1, pred.resp.bivar.levels$variable2),]

tiff(file=paste0("BKMR_bivarlevels_gauss_", R, ".tif"), width=10, height=8, units="in", compression="lzw", res=144, family="sans")
ggplot(pred.resp.bivar.levels, aes(z1, est)) + theme_bw() + geom_smooth(aes(col = quantile), stat = "identity") + 
    facet_grid(variable2 ~ variable1) + xlab("Exposure 1")  + labs(title="h(Exposure 1 | Quantiles of Exposure 2)") + labs(col = "Quantile", y = "Estimate") + 
    theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))
dev.off()
	
# Calculate summary statistics
# In addition to visually inspecting the estimated exposure-response function, one may also wish to calculate a range of summary statistics that highlight specific features of the high-dimensional surface.	
# Cumulative effects
# One potential summary measure of interest is to compute the overall effect of the mixture, by comparing the value of the exposure-response function when all of exposures are at a particular quantile as compared to when all of them are at their median value. 
# The function OverallRiskSummaries allows one to specify a sequence of quantiles using the argument qs and the fixed quantile using the argument q.fixed (the default is the median value).

## iv.混合物总体效应
risks.overall = OverallRiskSummaries(fit = fitkmcomb, qs = seq(0.05, 0.95, by = 0.05), q.fixed = 0.5, method = "exact", sel = sel)  # method='approx' or 'exact'
risks.overall = risks.overall[4:16, ]
tiff(file=paste0("BKMR_overall_gauss_", R, ".tif"), width=10, height=8, units="in", compression="lzw", res=144, family="sans")
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + theme_bw() + 
       geom_pointrange(size = 0.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "brown") + 
       labs(x = "Quantile", y = "Estimates")
dev.off()      
	   
# We can see that as cumulative levels across all exposures increases, the health risks increase. There is also a suggestion of a ceiling effect, which indicates a nonlinear exposure-response function.

# Single-exposure effects
# Another summary that may be of interest is to estimate the contribution of individual exposures to the cumulative effect. For example, we may wish to compare risk when a single exposure is at the 75th percentile as compared to when that exposure is at its 25th percentile, 
# where all of the remaining exposures are fixed to a particular quantile. We refer to this as the single-exposure health risks, and these can be computed using the function SingVarRiskSummaries. The two different quantiles at which to compare the risk are specified using the 
# qs.diff argument, and a sequence of values at which to fix the remaining exposures can be specified using the q.fixed argument.

## v.单个暴露效应
risks.singvar = SingVarRiskSummaries(fit = fitkmcomb, qs.diff = c(0.05, 0.95), q.fixed = c(0.10, 0.50, 0.90), method = "exact", sel = sel)  # method='approx' or 'exact'
risks.singvar
tiff(file=paste0("BKMR_singvar_gauss_", R, ".tif"), width=10, height=8, units="in", compression="lzw", res=144, family="sans")
ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, col = q.fixed)) + theme_bw() + 
	 geom_hline(aes(yintercept = 0), linetype = "dashed", color = "brown") +  
	 geom_pointrange(position = position_dodge(width = 0.75)) + coord_flip() + 
	 labs(x = "", y = "Estimates", col = "Fixed Quantile") + theme(legend.position = "right")
# theme(legend.text = element_text(size = 10), legend.title = element_text(size = 15)) 
dev.off()

# Interactions
# We may wish to compute specific ‘interaction’ parameters. For example, we may which to compare the single-exposure health risks when all of the other exposures are fixed to their 75th percentile to when all of the other exposures are fixed to their 25th percentile. 
# In the previous plot, this corresponds to substracting the estimate represented by the red circle from the estimate represented by the blue circle. This can be done using the function SingVarIntSummaries.

## vi.交互作用
risks.int = SingVarIntSummaries(fit = fitkmcomb, qs.diff = c(0.05, 0.95), qs.fixed = c(0.10, 0.90), method = "exact", sel = sel)  # method='approx' or 'exact'
risks.int
tiff(file=paste0("BKMR_singvarint_gauss_", R, ".tif"), width=10, height=8, units="in", compression="lzw", res=144, family="sans")
ggplot(risks.int, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + theme_bw() + 
    geom_pointrange(position = position_dodge(width = 0.75)) + geom_hline(yintercept = 0, lty = 2, col = "brown") + coord_flip() +
    labs(x = "", y = "Estimates") 
dev.off()

# tiff(file="BMKR_univar_overall_gauss.tif", width=18, height=4, units="in", compression="lzw", res=144, family="sans")
# plot_grid(g1, g2, g3, ncol=3, labels="AUTO", byrow=T, label_colour="black", align="hv")
# dev.off()


############继续采样############################
library(bkmr)
library(bkmrhat)

load(file="BKMRHAT_fitkmcomb_10000.RData")

## 参数设置
R = 20000             ## no. of iterations
burn = 0.5           ## percent burn-in
thin = 40            ## thinning number
sel = seq(burn * R + 1, R, by=thin) 


ncores = 10
future::plan(strategy = future::multisession, workers=ncores)  # 设置并行策略，strategy='sequential', 'multisession', 'multicore', 'cluster' 
start = proc.time()
set.seed(R)
fitkmmore = kmbayes_parallel_continue(fitkm, iter=ceiling(R/ncores))
# fitprmore = kmbayes_parallel_continue(fitpr, iter=ceiling(R/ncores))   
diftime = proc.time() - start
print(paste("Execution time:", round(diftime[3]/3600, 2), "hours")) #输出循环时间

save.image(file=paste0("BKMRHAT_fitkmmore_", R+10000, ".RData"))  #保存工作目录特定文件

################################################


########################################BMIM###################################
# library(devtools)
# devtools::install_github("glenmcgee/bsmim2") 
library(bsmim2)
library(ggplot2)
library(ggsci)
library(cowplot)
# library(knitr)  # RMarkdown文档中嵌入表格
# library(xtable)  # 输出LaTeX或HTML表格

## 参数设置
R = 10000              ## no. of iterations
burn = 0.5           ## percent burn-in
thin = 40            ## thinning number
sel = seq(burn * R + 1, R, by=thin) 


## 模型预测数据
## Construct new grid points, quantiles and weights
getGrid = function(qtl=0.5,pts=20,qtl_lims=c(0.05,0.95))
{
  Xq = c()
  for(ll in 1:ncol(X)){
    tempXq = matrix(apply(X,2,function(x) quantile(x,qtl)),nrow=pts,ncol=ncol(X),byrow=T)
    tempXq[,ll] = seq(quantile(X[,ll],qtl_lims[1]),quantile(X[,ll],qtl_lims[2]),length=pts)
    Xq = rbind(Xq,tempXq)
  }
  return(Xq)
}
X25 = getGrid(qtl=0.25)
X50 = getGrid(qtl=0.50)
X75 = getGrid(qtl=0.75)


## 模型拟合数据
subset_data = function(ids)
{ 
  ## exposure 
  X = as.matrix(X[ids,])
  ##  exposures for model fits
  X_bkmr = list(as.matrix(X[,1]),as.matrix(X[,2]),as.matrix(X[,3]),as.matrix(X[,4]),as.matrix(X[,5]),as.matrix(X[,6]))
  X_SIM = list(X[,1:6])
  X_bsmim = list(X[,1:2],X[,3:6]) 
  ## covariate matrix
  Z = Z[ids,] 
  ## outcomes
  y = y[ids]
  ybin = ybin[ids] 
  ## random effects
  ID = factor(ID[ids])
  ## return data
  return(list(X=X,X_bkmr=X_bkmr,X_SIM=X_SIM,X_bsmim=X_bsmim,Z=Z,y=y,ybin=ybin,ID=ID))
}

prep_data_full = function()
{
  sdat = subset_data(1:nrow(lst$dat))

  ### exposures for predictions
  X25_bkmr = list(as.matrix(X25[,1]), as.matrix(X25[,2]), as.matrix(X25[,3]), as.matrix(X25[,4]), as.matrix(X25[,5]), as.matrix(X25[,6]))
  X50_bkmr = list(as.matrix(X50[,1]), as.matrix(X50[,2]), as.matrix(X50[,3]), as.matrix(X50[,4]), as.matrix(X50[,5]), as.matrix(X50[,6]))
  X75_bkmr = list(as.matrix(X75[,1]), as.matrix(X75[,2]), as.matrix(X75[,3]), as.matrix(X75[,4]), as.matrix(X75[,5]), as.matrix(X75[,6]))
  X25_SIM = list(X25[,1:6])
  X50_SIM = list(X50[,1:6])
  X75_SIM = list(X75[,1:6])
  X25_bsmim = list(X25[,1:2], X25[,3:6])
  X50_bsmim = list(X50[,1:2], X50[,3:6])
  X75_bsmim = list(X75[,1:2], X75[,3:6])
  
  ## return data
  SIM = list(X=sdat$X_SIM,      X25=X25_SIM,   X50=X50_SIM,   X75=X75_SIM)
  bsmim = list(X=sdat$X_bsmim,  X25=X25_bsmim, X50=X50_bsmim, X75=X75_bsmim)
  bkmr = list(X=sdat$X_bkmr,    X25=X25_bkmr,  X50=X50_bkmr,  X75=X75_bkmr)
  return(list(SIM=SIM, bsmim=bsmim, bkmr=bkmr, Z=sdat$Z, y=sdat$y, ybin=sdat$ybin, ID=sdat$ID))
}

pred_twoway = function(obj,qtls=c(0.1,0.9),qtl_lims=c(0.01,0.99),pts=20,include_median=TRUE)
{ 
  ## get predictions at each level (skip 0.5 since it is implicitly computed)
  res = list()
  for(mm in 1:ncol(obj$rho)){
    res_mm = list()
    for(qq in 1:length(qtls)){
      res_mm[[qq]] = predict_hnew_indexwise2(obj,crossM=mm,qtl=qtls[[qq]],qtl_lims=qtl_lims,points=pts)
    }
    names(res_mm) = qtls  ## label
    res[[mm]] = res_mm
  }
  
  ## combine predictions into dataframe for plotting
  df_var1 = df_var2 = df_grid = df_quantile = df_est = c()
  for(xx in 1:ncol(obj$rho)){
    for(yy in 1:ncol(obj$rho)){
      if(xx==yy){
        next
      }
      for(qq in 1:length(qtls)){
        df_var1 = c(df_var1,rep(paste0("Index ",xx),pts))
        df_var2 = c(df_var2,rep(paste0("Index ",yy),pts))
        df_grid = c(df_grid,res[[yy]][[qq]]$grid[(xx-1)*pts+1:pts])
        df_est = c(df_est,res[[yy]][[qq]]$mean[(xx-1)*pts+1:pts])
        df_quantile = c(df_quantile,rep(qtls[qq],pts))
      }
      if(include_median==TRUE){ ## implicitly predicted above
        df_var1 = c(df_var1,rep(paste0("Index ",xx),pts))
        df_var2 = c(df_var2,rep(paste0("Index ",yy),pts))
        df_grid = c(df_grid,res[[xx]][[qq]]$grid[(xx-1)*pts+1:pts]) ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_est = c(df_est,res[[xx]][[qq]]$mean[(xx-1)*pts+1:pts])   ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_quantile = c(df_quantile,rep(0.5,pts))
      }
      
    }
  }
  
  pred_df = data.frame(var1=df_var1,var2=df_var2,grid=df_grid,quantile=df_quantile,est=df_est)
  
  return(pred_df)
}

## 指数模型拟合
bsdat = prep_data_full()
mod_version = 1

if(mod_version==1){  ## fit 1-dim single index model (SIM)

  set.seed(R)
  # fit = bsmim2(y=y,x=bsdat$SIM$X,z=bsdat$Z,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  fit = bsmim2(y=bsdat$y, x=bsdat$SIM$X, z=bsdat$Z, group_id=bsdat$ID, niter=R, nburn=R*burn, nthin=thin, spike_slab=TRUE, gauss_prior=TRUE, stepsize_theta=0.2, nchains = 1)
  pred_assoc = predict_hnew_assoc2(fit)
  pred_overall = predict_hnew_assoc2(fit,overall = TRUE)
  pred_ind = predict_hnew_indexwise2(fit)
  pred_inter = NULL
  SIM_list = list(fit=fit,y=y,pred_assoc=pred_assoc,pred_overall=pred_overall,pred_ind=pred_ind,pred_inter=pred_inter)
  save(SIM_list, file=paste0("BMIM_SIM_list_", R, ".RData") )
} else if(mod_version==2){  ## fit 2-dim multi index model (bsmim)

  set.seed(R)
  # fit = bsmim2(y=y,x=bsdat$bsmim$X,z=bsdat$Z,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  fit = bsmim2(y=y, x=bsdat$bsmim$X, z=bsdat$Z, group_id=bsdat$ID, niter=R, nburn=R*burn, nthin=thin, spike_slab=TRUE, gauss_prior=TRUE, stepsize_theta=0.2, nchains = 1)
  pred_assoc = predict_hnew_assoc2(fit)
  pred_overall = predict_hnew_assoc2(fit,overall = TRUE)
  pred_ind = predict_hnew_indexwise2(fit)
  pred_inter = pred_twoway(fit)
  bsmim_list = list(fit=fit,y=y,pred_assoc=pred_assoc,pred_overall=pred_overall,pred_ind=pred_ind,pred_inter=pred_inter)
  save(bsmim_list, file=paste0("BMIM_bsmim_list_", R, ".RData") )
} else if(mod_version==3){  ## fit 6-dim multi index model (bkmr)

  set.seed(R)
  # fit = bsmim2(y=y,x=bsdat$bkmr$X,z=bsdat$Z,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=0.2,basis.opts=NULL,draw_h=FALSE)
  fit = bsmim2(y=y, x=bsdat$bkmr$X, z=bsdat$Z, group_id=bsdat$ID, niter=R, nburn=R*burn, nthin=thin, spike_slab=TRUE, gauss_prior=TRUE, stepsize_theta=0.2, nchains = 1)
  pred_assoc = predict_hnew_assoc2(fit)
  pred_overall = predict_hnew_assoc2(fit,overall = TRUE)
  pred_ind = predict_hnew_indexwise2(fit)
  pred_inter = NULL
  bkmr_list = list(fit=fit,y=y,pred_assoc=pred_assoc,pred_overall=pred_overall,pred_ind=pred_ind,pred_inter=pred_inter)
  save(bkmr_list, file=paste0("BMIM_bkmr_list_", R, ".RData") )
}


## 结果汇总
#loads an RData file, and returns it with new name
loadRData = function(fileName)
{
    load(fileName)
    get(ls()[ls() != "fileName"])
}


printCI = function(df, col1, col2, dig=4)
{ 
  CI = paste0("(",round(df[,col1],dig),", ",round(df[,col2],dig),")")
  res = CI
  if(col1 > 1){
    res = cbind(df[1:(col1-1)], CI)
  }
  if(col2 < ncol(df)){
    res = cbind(res, df[(col2+1):ncol(df)])
  }
  return(res)
}


mod_names = c("SIM", "bsmim", "bkmr")
mods = list()
for(mm in 1:length(mod_names))
{
  mods[[mm]] = loadRData(paste0("BMIM_", mod_names[mm],"_list_", R, ".RData"))
}
names(mods) = mod_names

## weights
weights_SIM = summarize_thetas(mods$SIM$fit)[[1]]
weights_SIM = printCI(round(weights_SIM, 4),6,7)
weights_bsmim = rbind(summarize_thetas(mods$bsmim$fit)[[1]],summarize_thetas(mods$bsmim$fit)[[2]])
weights_bsmim = printCI(round(weights_bsmim, 4),6,7)
weights_bkmr = round(1 - apply(mods$bkmr$fit$rho, 2, function(x) mean(x==0)),4)

weights = cbind(weights_SIM,weights_bsmim,weights_bkmr)
rownames(weights) = colnames(mods$SIM$fit$x[[1]])
weights = cbind(weights, BKMR_PIP=pips$PIP)
write.csv(weights, file=paste0("BMIM_weights_", R, ".csv"))
# kable(weights, caption="Exposure Weights", booktabs = T)
# print(xtable(weights, digits=4))

## overall
quantile = seq(0.2, 0.8, by=0.05)
overall_SIM = mods$SIM$pred_overall$contrasts[4:16, ]
overall_bsmim = mods$bsmim$pred_overall$contrasts[4:16, ]
overall_bkmr = mods$bkmr$pred_overall$contrasts[4:16, ]


# overall = rbind(overall_SIM, overall_bsmim, overall_bkmr)
# overall = round(overall, 4)
# rownames(overall) = c("SIM", "MIM", "BKMR")
# write.csv(overall, file=paste0("BMIM_overall_", R, ".csv"))
# kable(overall, caption="Overall Effect", booktabs = T)
# print(xtable(overall,digits = 4))

tiff(file=paste0("BMIM_overall_", R, ".tif"), width=18, height=4, units="in", compression="lzw", res=144, family="sans")
SIM = ggplot(overall_SIM, aes(quantile, mean, ymin = lower, ymax = upper)) + theme_bw() + 
       geom_pointrange(size = 0.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "brown") + 
       labs(title = 'Index 1', x = "Quantile", y = "Estimates") + 
       theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))	

bsmim = ggplot(overall_bsmim, aes(quantile, mean, ymin = lower, ymax = upper)) + theme_bw() + 
       geom_pointrange(size = 0.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "brown") + 
       labs(title = 'Index 2', x = "Quantile", y = "Estimates") + 
       theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))	

bkmr = ggplot(overall_bkmr, aes(quantile, mean, ymin = lower, ymax = upper)) + theme_bw() + 
       geom_pointrange(size = 0.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "brown") + 
       labs(title = 'Index 6', x = "Quantile", y = "Estimates") + 
       theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))	   

plot_grid(SIM, bsmim, bkmr, ncol=3, labels="AUTO", byrow=T, label_colour="black", align="hv")

dev.off()  


## i. 单变量(univariate)
pp_SIM = plot_univar_hnew2(mods$SIM$pred_assoc, assoc=F, ylims=NULL)
pp_bsmim = plot_univar_hnew2(mods$bsmim$pred_assoc, assoc=F, ylims=NULL)
pp_bkmr = plot_univar_hnew2(mods$bkmr$pred_assoc, assoc=F, ylims=NULL)


pp_SIM11 = pp_SIM[[1]][[1]]$data  # PM2.5
pp_bsmim11 = pp_bsmim[[1]][[1]]$data  # PM2.5
pp_bkmr11 = pp_bkmr[[1]][[1]]$data  # PM2.5

pp_SIM12 = pp_SIM[[1]][[5]]$data  # O3
pp_bsmim12 = pp_bsmim[[2]][[3]]$data  # O3
pp_bkmr12 = pp_bkmr[[5]][[1]]$data  # O3

ylimit = function(x, y, z)
{
	df = rbind(x, y, z)
	lower = min(df[, 'lower']) - 0.1
	upper = max(df[, 'upper']) + 0.1
	return(c(lower, upper))
}
ylim11 = ylimit(pp_SIM11, pp_bsmim11, pp_bkmr11)
ylim12 = ylimit(pp_SIM12, pp_bsmim12, pp_bkmr12)


tiff(file=paste0("BMIM_univar_comp1_", R, ".tif"), width=18, height=4, units="in", compression="lzw", res=144, family="sans")
# SIM11 = pp_SIM[[1]][[1]]
SIM11 = ggplot(pp_SIM11, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim11) + 
  labs(title="Index 1, Component 1", x="Exposure Component", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2)) 

# bsmim11 = pp_bsmim[[1]][[1]] 
bsmim11 = ggplot(pp_bsmim11, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim11) +   
  labs(title="Index 2, Component 1", x="Exposure Component", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))

# bkmr11 = pp_bkmr[[1]][[1]]
bkmr11 = ggplot(pp_bkmr11, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() +
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim11) +   
  labs(title="Index 6, Component 1", x="Exposure Component", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))

plot_grid(SIM11, bsmim11, bkmr11, ncol=3, labels="AUTO", byrow=T, label_colour="black", align="hv")

dev.off()

tiff(file=paste0("BMIM_univar_comp5_", R, ".tif"), width=18, height=4, units="in", compression="lzw", res=144, family="sans")
# SIM12 = pp_SIM[[1]][[5]]
SIM12 = ggplot(pp_SIM12, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim12) + 
  labs(title="Index 1, Component 5", x="Exposure Component", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2)) 

# bsmim12 = pp_bsmim[[2]][[3]] 
bsmim12 = ggplot(pp_bsmim12, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim12) +   
  labs(title="Index 2, Component 5", x="Exposure Component", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))

# bkmr12 = pp_bkmr[[5]][[1]]
bkmr12 = ggplot(pp_bkmr12, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() +
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim12) +   
  labs(title="Index 6, Component 5", x="Exposure Component", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))

plot_grid(SIM12, bsmim12, bkmr12, ncol=3, labels="AUTO", byrow=T, label_colour="black", align="hv")

dev.off()

tiff(file=paste0("BMIM_univar_comp1_5_", R, ".tif"), width=18, height=8, units="in", compression="lzw", res=144, family="sans")

plot_grid(SIM11, bsmim11, bkmr11, SIM12, bsmim12, bkmr12, ncol=3, labels="AUTO", byrow=T, label_colour="black", align="hv")

dev.off()

## ii.指数(indexwise)
pp_SIM_ind = plot_univar_hnew_indexwise2(mods$SIM$pred_ind)
pp_bsmim_ind = plot_univar_hnew_indexwise2(mods$bsmim$pred_ind)

pp_SIM_ind1 =  pp_SIM_ind[[1]]$data
pp_bsmim_ind1 = pp_bsmim_ind[[1]]$data
pp_bsmim_ind2 = pp_bsmim_ind[[2]]$data

ylim_ind = ylimit(pp_SIM_ind1, pp_bsmim_ind1, pp_bsmim_ind2)

tiff(file=paste0("BMIM_indexwise_", R, ".tif"), width=18, height=4, units="in", compression="lzw", res=144, family="sans")
# SIM_ind1 = pp_SIM_ind[[1]]
SIM_ind1 = ggplot(pp_SIM_ind1, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim_ind) + 
  labs(title="Index 1", x="Exposure", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2)) 

# bsmim_ind1 = pp_bsmim_ind[[1]]
bsmim_ind1 = ggplot(pp_bsmim_ind1, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim_ind) + 
  labs(title="1 of Index 2", x="Exposure", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2)) 

# bsmim_ind2 = pp_bsmim_ind[[2]]
bsmim_ind2 = ggplot(pp_bsmim_ind2, aes(grid, mean, ymin = lower, ymax = upper)) + theme_bw() + 
  geom_smooth(stat = "identity", linetype="solid", size=0.8, color='blue4', alpha=0.7, fill='lightblue') + 
  scale_y_continuous(limits = ylim_ind) + 
  labs(title="2 of Index 2", x="Exposure", y="Estimated exposure-response (h)") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2)) 

plot_grid(SIM_ind1, bsmim_ind1, bsmim_ind2, ncol=3, labels="AUTO", byrow=T, label_colour="black", align="hv")

dev.off()


## iii.指数交互(Indexwise interactions)
pp_bsmim_inter = mods$bsmim$pred_inter

tiff(file=paste0("BMIM_indexwise_inter_", R, ".tif"), width=12, height=8, units="in", compression="lzw", res=144, family="sans")
ggplot(pp_bsmim_inter, aes(grid, est)) + theme_bw() + 
  geom_smooth(aes(col = as.factor(quantile)), stat = "identity", fill="white") + 
  facet_grid(var2 ~ var1, scales = "free_x") + # free_x allows different x axis limits
  labs(title="Indexwise interactions", x="", y="", col="Quantile") +
  theme(plot.title=element_text(size=12, face="plain", color="black", hjust=0.5, lineheight=1.2))
  
dev.off()


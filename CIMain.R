###CIMain.R


####################################Synthetic Dataset###############################################
setwd("~/Rscript")
source("simdata_causal.R")
source("MMD.R")
source("Causeffect.R")

library(umap)
library(future.apply)
library(MatchIt)
library(marginaleffects)


start = proc.time()

future::plan(future::multisession, workers = 20)
span = 1 # Interval
MMD = future_lapply(1:100, future.seed=TRUE, FUN=function(x) {
	compute_min_MMD(x, is.synthetic=TRUE, interval=span)
})
future::plan(future::sequential)

MMD = MMD[!sapply(MMD, is.null)]  # Remove empty elements


# future_lapply parallel
dimension = c(2, 6, seq(10, 100, by=10))

future::plan(future::multisession)
result1 = future_lapply(c('euclidean', 'mahalanobis', 'psm'), future.seed=TRUE, FUN=function(m) {
    estimate_ATT(method = m, replace=TRUE, B=100)  
})

result2 = future_lapply(dimension, future.seed=TRUE, FUN=function(d) {
    estimate_ATT(method = 'pca', d = d, replace=TRUE, B=100)  
})

result3 = future_lapply(dimension, future.seed=TRUE, FUN=function(d) {
    estimate_ATT(method = 'umap', d = d, replace=TRUE, B=100)  
})

result1 = do.call(rbind, result1)
result2 = do.call(rbind, result2)
result3 = do.call(rbind, result3)
result = rbind(result1, result2, result3)
write.csv(result, "~/Rdata/Causeffect_Synthetic123.csv", row.names=F)

result4 = future_lapply(dimension, future.seed=TRUE, FUN=function(d) {
    estimate_ATT(method = 'kdrm', d = d, replace=TRUE, B=100, MMD=MMD)  
})

future::plan(future::sequential)  # Close process

result4 = do.call(rbind, result4)
result = rbind(result, result4)
write.csv(result, "~/Rdata/Causeffect_Synthetic.csv", row.names=F)

diftime = proc.time() - start
print(paste("Execution time:", round(diftime[3]/60, 2), "minutes"))


####################################Austin2009 Dataset###############################################
setwd("~/Rscript")
source("simdata_causal.R")
source("MMD.R")
source("Causeffect.R")

library(umap)
library(future.apply)
library(MatchIt)
library(marginaleffects)

start = proc.time()

future::plan(future::multisession, workers = 20)
span = 1 # Interval
MMD = future_lapply(1:100, future.seed=TRUE, FUN=function(x) {
  compute_min_MMD(x, is.synthetic=FALSE, interval=span)
})
future::plan(future::sequential)

MMD = MMD[!sapply(MMD, is.null)]  # Remove empty elements


# future_lapply parallel
dimension = 2:7

future::plan(future::multisession)
result1 = future_lapply(c('euclidean', 'mahalanobis', 'psm'), future.seed=TRUE, FUN=function(m) {
    error_ATT(method = m, replace=FALSE, B=100)  
})

result2 = future_lapply(dimension, future.seed=TRUE, FUN=function(d) {
    error_ATT(method = 'pca', d = d, replace=FALSE, B=100)  
})

result3 = future_lapply(dimension, future.seed=TRUE, FUN=function(d) {
    error_ATT(method = 'umap', d = d, replace=FALSE, B=100)  
})

result1 = do.call(rbind, result1)
result2 = do.call(rbind, result2)
result3 = do.call(rbind, result3)
result = rbind(result1, result2, result3)
write.csv(result, "~/Rdata/Causeffect_Austin123.csv", row.names=F)

result4 = future_lapply(dimension, future.seed=TRUE, FUN=function(d) {
    error_ATT(method = 'kdrm', d = d, replace=FALSE, B=100, MMD=MMD)  
})

future::plan(future::sequential)  # Close process

result4 = do.call(rbind, result4)
result = rbind(result, result4)
write.csv(result, "~/Rdata/Causeffect_Austin.csv", row.names=F)

diftime = proc.time() - start
print(paste("Execution time:", round(diftime[3]/60, 2), "minutes"))




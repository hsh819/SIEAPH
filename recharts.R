#recharts包在线安装
install.packages("remotes")
library(remotes)
remotes::install_github("cosname/recharts")#从Github上获取recharts包

#recharts包本地安装
install.packages("remotes")
library(remotes)
remotes::install_local("recharts-master.zip")#从本地安装recharts-master.zip

install.packages('recharts',repos = c('http://yihui.name/xran', 'http://cran.rstudio.com'))
library(recharts)

mapData<-head(mapTestData_chs,5)
map=eMap(mapData,namevar=~stdName,datavar=~val1+val2,color=c("red","yellow"),region="china")
map1=map+eTitle(title="Part of province on China")#标题
map2=map1+eLegend(show=TRUE,orient="vertical",x="left",y="top")#图例
map3=map2+eToolbox(show=TRUE,orient="vertical",x="right",y="top")#工具箱
map4=map3+eDataRange(show=TRUE,orient="vertical",text=c("高","低"))#数据区间


provinceMapData<-mapTestData_chs[6:15,]
pmap=eMap(provinceMapData,namevar=~stdName,datavar=~val1+val2,region=unique(provinceMapData$motherName)[1],color=c("red","yellow"))
pmap+eTitle(title="Hebei Province")
pmap+eLegend(orient="vertical")
pmap+eToolbox(orient="vertical")
pmap+eDataRange(orient="vertical")




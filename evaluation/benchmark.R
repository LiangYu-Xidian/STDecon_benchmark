    ### calculation of JSD and RMSE for seqFISH+ and MERFISH

library(stringr)
library(philentropy)

source('D:/seqFISH+/evaluation/our_metrics.R')

### MERFISH
data = readRDS('D:/MERFISH/simMERFISH_100.RDS')
ground_truth = data$gtSpotTopics
ground_truth = as.matrix(ground_truth)
setwd('D:/MERFISH/result/100/')
ls_method = list.files()
Results = data.frame(matrix(ncol = ncol(ground_truth) +2, nrow = length(ls_method)))
for(i in 1:length(ls_method)){

  
  method = read.csv(ls_method[i],row.names = 1)
  colnames(method)[str_detect(colnames(method),'_')] = gsub('_','.',colnames(method)[str_detect(colnames(method),'_')])
  method = method[,colnames(ground_truth)]
  method = as.matrix(method)
  res = benchmark_performance(method,ground_truth)
    Results[i,1] = res$JSD
  Results[i,2] = res$Sum_RMSE
  Results[i,3:8] = res$RMSE
  Results[i,9] = res$corr$estimate
  Results[i,10] = res$kendall$estimate
  Results[i,11] = res$spearman$estimate
  # For MERFISH_50 & MERFISH_50
  # Results[i,1] = res$JSD
  # Results[i,2] = res$Sum_RMSE
  # Results[i,3:ncol(Results)] = res$RMSE

}
colnames(Results) = c('JSD','total_RMSE',colnames(ground_truth),'PCC','kendall','Spearman')
# colnames(Results) = c('JSD','total_RMSE',colnames(ground_truth)) # For MERFISH_50 & MERFISH_50
rownames(Results) = gsub('.csv','',ls_method)
write.csv(Results,'D:/MERFISH/MERFISH_100.csv')



### seqFISH+

ground_truth = read.csv('D:/seqFISH+/datasets/seqFISH+/Out_cell_ratio_1x.csv',row.names=1)
delete_cells_index = which(is.na(ground_truth[,1]))
ground_truth = ground_truth[-delete_cells_index,]
ground_truth = as.matrix(ground_truth)
setwd('D:/seqFISH+/Results/seqFISH/3000/')
ls_method = list.files()
Results = data.frame(matrix(ncol = ncol(ground_truth) +5, nrow = length(ls_method)))

for(i in 1: length(ls_method)){
  method = read.csv(ls_method[i],row.names = 1)
  colnames(method)[str_detect(colnames(method),'_')] = gsub('_','.',colnames(method)[str_detect(colnames(method),'_')])
  method = method[,colnames(ground_truth)]
  method = as.matrix(method)
  res = benchmark_performance(method,ground_truth)
  Results[i,1] = res$JSD
  Results[i,2] = res$Sum_RMSE
  Results[i,3:8] = res$RMSE
  Results[i,9] = res$corr$estimate
  Results[i,10] = res$kendall$estimate
  Results[i,11] = res$spearman$estimate
}
colnames(Results) = c('JSD','total_RMSE',colnames(ground_truth),'PCC','kendall','Spearman')
rownames(Results) = gsub('.csv','',ls_method)

write.csv(Results,'D:/seqFISH+/evaluation/seqFISH_10000.csv')

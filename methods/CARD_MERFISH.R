consumption_time<-system.time({
library(Matrix)
library(data.table)
library(Seurat)
library(CARD)

library("reticulate")
use_python("D:/pycharm/DestVI/venv/Scripts")
sc <- import("scanpy")
pd <- import("pandas")

sc_exp = read.table("D:/MERFISH/raw_somatosensory_sc_exp.txt",header = T,row.names = 1)
sc_anno = read.table("D:/MERFISH/somatosensory_sc_labels.txt",header = F)

simMERFISH_100<-readRDS("D:\\MERFISH\\simMERFISH_100.rds")



data_count<-simMERFISH_100$sim
sparse_mat <- sparseMatrix(i = data_count$i, j = data_count$j, x = data_count$v, dims = c(data_count$nrow, data_count$ncol))
count<-as.matrix(sparse_mat)
rownames(count)<-simMERFISH_100$sim$dimnames[[1]]
colnames(count)<-simMERFISH_100$sim$dimnames[[2]]
rm(data_count)
rm(sparse_mat)

st_location<-simMERFISH_100$st_location
colnames(st_location) = c('x','y')
#dim(st_location)

sc_count = Matrix(as.matrix(sc_exp),sparse = TRUE)
sc_meta = data.frame(matrix(ncol = 3,nrow = ncol(sc_exp)))
colnames(sc_meta) = c('cellID','cellType','sampleInfo')
sc_meta$sampleInfo = "sample1"
sc_meta$cellType = sc_anno$V1
sc_meta$cellID = colnames(sc_count)
rownames(sc_meta) = sc_meta$cellID
count <-t(as.matrix(count))
#dim(count)
CARD_obj = createCARDObject(
  sc_count = as.matrix(sc_exp),
  sc_meta = sc_meta,
  spatial_count = count,
  spatial_location = st_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 0,
  minCountSpot = 0) 
#dim(CARD_obj@spatial_countMat)
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
#dim(CARD_obj@Proportion_CARD)
write.csv(CARD_obj@Proportion_CARD, 'D:/MERFISH/CARD/CARD_MERFISH_100.csv')
})
print(consumption_time)

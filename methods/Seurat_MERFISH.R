consumption_time<-system.time({
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

#preparing ST data
simMERFISH_100<-readRDS("D:\\MERFISH\\simMERFISH_100.rds")
data_count<-simMERFISH_100$sim
sparse_mat <- sparseMatrix(i = data_count$i, j = data_count$j, x = data_count$v, dims = c(data_count$nrow, data_count$ncol))
count<-as.matrix(sparse_mat)
rownames(count)<-simMERFISH_100$sim$dimnames[[1]]
colnames(count)<-simMERFISH_100$sim$dimnames[[2]]
rm(data_count)
rm(sparse_mat)
#st_location<-simMERFISH_100$st_location
count <-t(as.matrix(count))

#preparing scRNA data
sc_exp = read.table('D:/MERFISH/raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('D:/MERFISH/somatosensory_sc_labels.txt',header = F)
cell_type = sc_anno[,1]

sc_seu <- CreateSeuratObject(counts = sc_exp)
cell_type = sc_anno[,1]
sc_seu <- AddMetaData(
  object = sc_seu,
  metadata = cell_type,
  col.name = 'cell_type'
)

sp_seu <- CreateSeuratObject(counts = count)

sc_seu <- SCTransform(sc_seu, verbose = FALSE) %>% RunPCA(verbose = FALSE)
sp_seu <- SCTransform(sp_seu, verbose = FALSE) %>% RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = sc_seu, query = sp_seu, normalization.method = "SCT",reduction='cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = factor(sc_seu$cell_type),
                                  weight.reduction = 'cca')

#decon_mtrx <- as.data.frame(predictions.assay[ , ])
decon_mtrx <- predictions.assay[, 2:7]
rownames(decon_mtrx) = colnames(sp_seu)

write.csv(decon_mtrx, "D:/MERFISH/Seurat/Seurat_MERFISH_100.csv")
})
print(consumption_time)


consumption_time<-system.time({
library(spacexr)
library(Matrix)

simMERFISH_100<-readRDS("D:\\MERFISH\\simMERFISH_100.rds")
data_count<-simMERFISH_100$sim
sparse_mat <- sparseMatrix(i = data_count$i, j = data_count$j, x = data_count$v, dims = c(data_count$nrow, data_count$ncol))
count<-as.matrix(sparse_mat)
rownames(count)<-simMERFISH_100$sim$dimnames[[1]]
colnames(count)<-simMERFISH_100$sim$dimnames[[2]]
rm(data_count)
rm(sparse_mat)
st_location<-simMERFISH_100$st_location
count <-t(as.matrix(count))

sc_counts = read.table("D:/MERFISH/raw_somatosensory_sc_exp.txt",header = T,row.names = 1)
sc_labels = read.csv("D:/MERFISH/somatosensory_sc_labels.txt",header = FALSE)$V1  
names(sc_labels)=colnames(sc_counts)
sc_labels <- as.factor(sc_labels)

sc_reference=Reference(
  counts=sc_counts,
  cell_types=sc_labels
)


st_data=SpatialRNA(
  counts=count,
  coords=st_location,
  require_int=FALSE
)


myRCTD <- create.RCTD(st_data, sc_reference, max_cores = 1, CELL_MIN_INSTANCE = 0)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

weights=myRCTD@results$weights
norm_weights=normalize_weights(weights)

write.csv(as.matrix(norm_weights),'D:/MERFISH/RCTD/RCTD_MERFISH_100.csv')
})
print(consumption_time)
consumption_time<-system.time({
  library(DropletUtils)
  library(Seurat)
  source('D:/DSTG-main/DSTG-main/DSTG/ours_utils.r')

org_st_count = read.csv('D:/MERFISH/SpaDecon/count_matrix.csv',header = T, row.names = 1)

sc_exp = read.table('D:/MERFISH/raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('D:/MERFISH/somatosensory_sc_labels.txt',header = F)
st_count = t(as.matrix(org_st_count))
cell_type = sc_anno[,1]
sc_exp<-as.matrix(sc_exp)
## DSTG
process_data(sc_exp,st_count,cell_type,spot_num = 300, 
             dropout_extract = F,combine_feature = F)
})
print(consumption_time)
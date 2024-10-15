consumption_time<-system.time({
library(STdeconvolve)
library(Matrix)
st.data.file = 'D:/MERFISH/simMERFISH_100.RDS'
simMERFISH_100<-readRDS(st.data.file)
data_count<-simMERFISH_100$sim
sparse_mat <- sparseMatrix(i = data_count$i, j = data_count$j, x = data_count$v, dims = c(data_count$nrow, data_count$ncol))
count<-as.matrix(sparse_mat)
rownames(count)<-simMERFISH_100$sim$dimnames[[1]]
colnames(count)<-simMERFISH_100$sim$dimnames[[2]]
rm(data_count)
rm(sparse_mat)

count <-t(as.matrix(count))

## Sanitize spots
st.counts <- cleanCounts(count, min.lib.size = 150)

## Sanitize genes
st.corpus <- restrictCorpus(st.counts, removeAbove=1.0, removeBelow = 0.001)

## Delete spots who are all zero in gene expressions
col.sum = apply(st.corpus, 2, FUN=sum)
st.corpus.rm = st.corpus[,col.sum != 0]

## Record those spots that are deleted to restore in the future
st.corpus.dl = colnames(st.corpus)[col.sum == 0]
dl.index = (1:length(colnames(st.corpus)))[col.sum == 0]
st.corpus
dl.index


data_values <- st.corpus.rm[, -1]


data_values <- as.data.frame(apply(data_values, 1:2, as.integer))

st.corpus.rm <- data_values

## Fitting model and get the optimum
st.ldas <- fitLDA(t(as.matrix(st.corpus.rm)), Ks = 6)
st.optLDA <- optimalModel(models = st.ldas, opt = "min")

## Extract deconvolved cell-type proportions (theta) 
## and transcriptional profiles (beta)
st.results <- getBetaTheta(st.optLDA, perc.filt = 0.05, betaScale = 1000)

## Results
st.deconProp <- st.results$theta
st.deconGexp <- st.results$beta

# write down proportions
write.csv(st.deconProp, "D:/MERFISH/STdeconvolve/STdeconvolve_MERFISH_100.csv")
})
print(consumption_time)

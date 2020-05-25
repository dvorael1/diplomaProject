# library(SingleCellExperiment)
# library(scater)
# options(stringsAsFactors = FALSE)
# umi <- readRDS("tung/umi.rds")
# umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
# endog_genes <- !rowData(umi.qc)$is_feature_control

M <- matrix(1:9, nrow=3, byrow=TRUE)

cpm_normalize = function(col){
  s = sum(col)
  col = 1000000*col
  col = col/s
  return(col)
}

rpk_normalize = function(col){
  s = sum(col)
  col = 1000*col
  col = col/s
  return(col)
}

t = apply(M,1,rpk_normalize)
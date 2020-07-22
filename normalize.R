# library(SingleCellExperiment)
# library(scater)
# options(stringsAsFactors = FALSE)
# umi <- readRDS("tung/umi.rds")
# umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
# endog_genes <- !rowData(umi.qc)$is_feature_control



calc_cpm <- function (expr_mat, spikes = NULL) # spikes must be integer, not logical vector
{ 
  norm_factor <- colSums( expr_mat ) 
  return( t( t( expr_mat ) / norm_factor ) * 10^6 ) 
}

calc_sf <- function (expr_mat) { 
  gm <- function(cnts) { 
    exp( mean( log( cnts[!(cnts==0)] ) ) ) 
  }
  geomeans <- apply( expr_mat, 1, gm) 
  SF <- function(cnts) { 
    tmp = ( ( cnts / geomeans )[ (is.finite(geomeans) & geomeans > 0) ] )
    median( tmp[tmp>0] ) 
  }
  norm_factor <- apply( expr_mat, 2, SF) 
  keep_cols = norm_factor > 0
  expr_mat[,keep_cols] = t( t( expr_mat[,keep_cols] ) / norm_factor ) 
  return(expr_mat)
}

calc_uq <- function (expr_mat, quantile=0.99) 
{ 
  UQ <- function(x) { quantile( x[ x > 0 ], quantile ) }
  keep_cols = colSums(expr_mat) > 0
  non_zero = expr_mat[,keep_cols]
  uq <- unlist( apply( non_zero, 2, UQ ) ) 
  
  norm_factor <- uq / median(uq) 
  result = ( t( t( non_zero ) / norm_factor) ) 
  expr_mat[,keep_cols] = result
  return(expr_mat)
} 

Down_Sample_Matrix <- function(expr_mat) {
  keep_cols = colSums(expr_mat) > 0
  non_zero = expr_mat[,keep_cols]
  min_lib_size <- min(colSums(non_zero))
  
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(
      sapply(x, function(y) { rbinom(1, y, prob) } )
    ) }
  
  down_sampled_mat <- apply(non_zero, 2, down_sample)
  expr_mat[,keep_cols] = down_sampled_mat
  return(expr_mat)
}

Up_Down_Sample_Matrix <- function(expr_mat,scale=1000) {

  keep_cols = colSums(expr_mat) > 0
  non_zero = expr_mat[,keep_cols]
  
  non_zero@x = non_zero@x + rnorm(length(non_zero@x),mean = 0, sd = 0.5)
  
  
  min_lib_size <- min(colSums(non_zero))
  
  up_down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(
      sapply(x, function(y) { 
        if(y<0.5){
          return(0)
        }
        if(prob >= 1){
          prob = 1
        }
        return(rbinom(1, round(y*scale ), prob))
        
        # if( is.na(t)){
        #   print(paste("prob:",prob,"y:",y, "min",min_lib_size,"sum",sum(x)))
        # }
        #   return(t)
          } )
    ) }
  
  up_down_sampled_mat <- apply(non_zero, 2, up_down_sample)
  expr_mat[,keep_cols] = up_down_sampled_mat
  return(expr_mat)
}

scran_normalize <- function(expr_mat){
  require(scran, quietly = TRUE)
  
  qclust <- quickCluster(expr_mat, min.size = 30)
  expr_mat <- calculateSumFactors(expr_mat, sizes = 15, clusters = qclust)
  
  return(expr_mat)
}
# 
# Up_Down_Sample_Matrix(raw)

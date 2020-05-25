
files <- list.files(path = "C:/Users/Ela/Desktop/owncloud/documents/documents/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)
# files <- list.files(path = "/media/eliska/FC089B5C089B152C/Users/Ela/Desktop/owncloud/documents/documents/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)
r <- lapply(files, readRDS)
# # v = 1:(length(r)/2)
v = 1:1
i=1
dataset = list()
for (i in v) {
  name = substr(files[2*i], 90, 98)
  t1 = r[[2*i-1]]
  t2 = r[[2*i]]

  t1[is.na(t1)] <- 0
  t2[is.na(t2)] <- 0
  
  t3 <- SingleCellExperiment(
    assays = list(counts = t1), 
  )
  t4 <- SingleCellExperiment(
    assays = list(counts = t2), 
  )
  
  keep_feature <- rowSums(counts(t3) > 0) > 0
  t3 <- t3[keep_feature, ]
  
  keep_feature <- rowSums(counts(t4) > 0) > 0
  t4 <- t4[keep_feature, ]
  
  t3 = scater::runPCA(t1,use_coldata = TRUE, detect_outliers = TRUE)
  
  t4 = runPCA(t2, use_coldata = TRUE, detect_outliers = TRUE)
  
  
  reducedDimNames(t3)
  reducedDimNames(t4)

#   for (row1 in row.names(t1)){
#     for (row2 in row.names(t2)[101:500]){
#       tbl_data = list()
#       tbl_data[[1]] = t1[row1,]
#       tbl_data[[2]] = t2[row2,]
#       dataset[[ paste( row1, row2, name, sep=" x ")]] = tbl_data
#       if(i>=100){
#         break
#       }
#     }
#   }

  }

for( i in names(dataset)){
  t = strsplit(i,"x")
  g1 = t[[1]][1] 
  g2 = t[[1]][2]
  if(length(dataset[[i]][[1]][dataset[[i]][[1]]!=0])>=1 && length(dataset[[i]][[2]][dataset[[i]][[2]]!=0])>=1){
    
    pdf(paste( i, ".pdf", sep=""))
    plot(dataset[[i]][[1]],dataset[[1]][[2]], main = i, xlab = g1, ylab = g2 )
    dev.off()
  }
}

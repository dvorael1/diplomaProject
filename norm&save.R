source("normalize.R")
require("Matrix")
set.seed(123)

preprocess = function(data){
  
  keep_feature <- rowSums(data) > 0
  data <- data[keep_feature, ]
  
  keep_cell <- colSums(data) > 0
  data <- data[ ,keep_cell]
  return(data)
}

files <- list.files(path = "C:/Users/elisk/OneDrive/Plocha/Owncloud/documents/documents/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)


v = 1:(length(files)/2)
i=1
size = c(50,100,500,1000)
setwd("rds")

for (s in size) {
  print(paste(i,16,sep = "/"))
  raw = readRDS(files[2*i])
  
  parts = strsplit(files[2*i],"/")
  name = strsplit(parts[[1]][length(parts[[1]])],"\\.")[[1]][1]
  
  raw[is.na(raw)] <- 0
  
  raw = preprocess(raw)
  
  udsm = Up_Down_Sample_Matrix(raw,s)
  if(sum(is.na(udsm))>0){
    stop("NA in udsm")
  }
  saveRDS(udsm,paste(name,"_udsm_",s,".rds",sep=""))
  
}
setwd("..")

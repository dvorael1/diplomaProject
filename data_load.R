source("normalize.R")

plot_scRNA = function(x,y=NULL,xlab = "",ylab = "",main = ""){
  mx = max(x)+1
  my = max(y)+1
  if(!is.null(y)){
    exclude = x == y & x == 0
    if(sum(exclude)>0){
      i = match(TRUE,exclude)
      exclude[i[1]] = FALSE
    }
    plot(x[!exclude],y[!exclude],main=main,xlab = xlab,ylab = ylab, xlim = c(0,mx), ylim=c(0,my))
  }else{
    plot(x,main=main,xlab = xlab,ylab = ylab, xlim = c(0,mx), ylim=c(0,mx))
  }
}

preprocess = function(data, csize = 600, rsize = 10){
  keep_feature <- rowSums(data) > 0
  data <- data[keep_feature, ]
  
  keep_cell <- colSums(data) > 0
  data <- data[ ,keep_cell]
  
  cols = sample.int(dim(data)[2], csize)
  data <- data[ ,cols]
  
  rows = sample.int(dim(data)[1], rsize)
  data <- data[rows, ]

  keep_feature <- rowSums(data) > 0
  data <- data[keep_feature, ]
  
  keep_cell <- colSums(data) > 0
  data <- data[ ,keep_cell]
  return(data)
}

files <- list.files(path = "C:/Users/Ela/Desktop/owncloud/documents/documents/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)
# files <- list.files(path = "/media/eliska/FC089B5C089B152C/Users/Ela/Desktop/owncloud/documents/documents/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)
r <- lapply(files, readRDS)
# # v = 1:(length(r)/2)
v = 1:1

i=1
dataset = list()
setwd("graphs")
for (i in v) {
  name = substr(files[2*i], 90, 98)
  t1 = r[[2*i-1]]
  t2 = r[[2*i]]
  
  t1[is.na(t1)] <- 0
  t2[is.na(t2)] <- 0
  
  raw = preprocess(t2,csize=600,rsize=20)
  
  pdf("normalize_example.pdf",5,4 )
  plot(raw,main="raw",xlab = "cells",ylab = "genes")

  cmp = calc_cpm(raw)
  plot(cmp,main="cmp",xlab = "cells",ylab = "genes")

  sf = calc_sf(raw)
  plot(sf,main="sf",xlab = "cells",ylab = "genes")

  uq = calc_uq(raw)
  plot(uq,main="uq",xlab = "cells",ylab = "genes")

  dsm = Down_Sample_Matrix(raw)
  plot(dsm,main="Down_Sample_Matrix",xlab = "cells",ylab = "genes")

  # scran_norm = scran_normilize(raw)
  # plot(log2(scran_norm),main="scran_normilize",xlab = "cells",ylab = "genes")


  dev.off()

  
  for (row1 in row.names(raw)){
    for (row2 in row.names(raw)){
    
      pdf(paste(paste(row1,row2,sep="_"),"pdf",sep="."),5,4 )
      
      x = raw[row1,]
      y = raw[row2,]
      
      
      plot_scRNA(x,y,main="raw",xlab = row1,ylab = row2)
      
      x = cmp[row2,]
      y = cmp[row2,]
      plot_scRNA(x,y,main="cmp",xlab = row1,ylab = row2)
      
      x = sf[row2,]
      y = sf[row2,]
      plot_scRNA(x,y,main="sf",xlab = row1,ylab = row2)
      
      x = uq[row2,]
      y = uq[row2,]
      plot_scRNA(x,y,main="uq",xlab = row1,ylab = row2)
      
      x = dsm[row2,]
      y = dsm[row2,]
      plot_scRNA(x,y,main="Down_Sample_Matrix",xlab = row1,ylab = row2)
      
      dev.off()
    }
  }
  
}
setwd("..")

# for( i in names(dataset)){
#   t = strsplit(i,"x")
#   g1 = t[[1]][1] 
#   g2 = t[[1]][2]
#   if(length(dataset[[i]][[1]][dataset[[i]][[1]]!=0])>=1 && length(dataset[[i]][[2]][dataset[[i]][[2]]!=0])>=1){
#     
#     pdf(paste( i, ".pdf", sep=""))
#     plot(dataset[[i]][[1]],dataset[[1]][[2]], main = i, xlab = g1, ylab = g2 )
#     dev.off()
#   }
# }

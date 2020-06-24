source("normalize.R")
library("ggpubr")

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
  require("matrixStats")
  
  keep_feature <- rowSums(data) > 0
  data <- data[keep_feature, ]
  
  keep_cell <- colSums(data) > 0
  data <- data[ ,keep_cell]
  
  set.seed(1)
  cols = sample.int(dim(data)[2], csize)
  data <- data[ ,cols]
  
  # rows = sample.int(dim(data)[1], rsize)
  # data <- data[rows,]
  

  keep_feature <- rowMaxs(as.matrix(data)) > 1
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
setwd("filter")
for (i in v) {
  name = substr(files[2*i], 90, 98)
  t1 = r[[2*i-1]]
  t2 = r[[2*i]]
  
  t1[is.na(t1)] <- 0
  
  t2[is.na(t2)] <- 0
  
  
  raw = preprocess(t2,csize=600,rsize=20)
  # 
  # pdf("normalize_example.pdf",5,4 )
  # plot(raw,main="raw",xlab = "cells",ylab = "genes")
  # 
  cmp = calc_cpm(raw)
  # plot(cmp,main="cmp",xlab = "cells",ylab = "genes")
  # 
  sf = calc_sf(raw)
  # plot(sf,main="sf",xlab = "cells",ylab = "genes")
  # 
  uq = calc_uq(raw)
  # plot(uq,main="uq",xlab = "cells",ylab = "genes")
  # 
  dsm = Down_Sample_Matrix(raw)
  # plot(dsm,main="Down_Sample_Matrix",xlab = "cells",ylab = "genes")

  # scran_norm = scran_normilize(raw)
  # plot(scran_norm,main="scran_normilize",xlab = "cells",ylab = "genes")
# 
# 
#   dev.off()

  
  for (row1 in row.names(raw)){
    x = raw[row1,]
    omit = x <= 0
    
    # orderx = order(x)
    if(length(unique(x))<5){
      next
    }
    for (row2 in row.names(raw)){
    
      filename = paste(paste(row1,row2,sep="_"),"pdf",sep=".")
      filenamerev = paste(paste(row2,row1,sep="_"),"pdf",sep=".")
      if(file.exists(filename) || file.exists(filenamerev)){
        next
      }
      
      
      y = raw[row2,]
      omit = y<= 0 | omit | (x<=1 & y<=1)
      
      if(sum(!omit)>=4 && ( length( unique( x)) /length( unique( y)) > 0.5 && length( unique( x))/ length( unique( y)) < 1.5)){
        
        stestraw = cor.test( x, y, method = "spearman", exact=F ) [["p.value"]]
        prgrs = FALSE
        
        cmpx = cmp[row1,]
        cmpy = cmp[row2,]
        stestcmp = cor.test( cmpx, cmpy, method = "spearman", exact=F ) [["p.value"]]
        prgrs = prgrs || stestcmp >= stestraw + 0.1 || stestcmp <= stestraw - 0.1
          
        sfx = sf[row1,]
        sfy = sf[row2,]
        stestsf = cor.test( sfx, sfy, method = "spearman", exact=F ) [["p.value"]]
        prgrs = prgrs || stestsf >= stestraw + 0.1 || stestsf <= stestraw - 0.1
        
        uqx = uq[row1,]
        uqy = uq[row2,]
        stestuq = cor.test( uqx, uqy, method = "spearman", exact=F ) [["p.value"]]
        prgrs = prgrs || stestuq >= stestraw + 0.1 || stestuq <= stestraw - 0.1
        
        dsmx = dsm[row1,]
        dsmy = dsm[row2,]
        stestdsm = cor.test( dsmx, dsmy, method = "spearman", exact=F ) [["p.value"]]
        prgrs = prgrs || stestdsm >= stestraw + 0.1 || stestdsm <= stestraw - 0.1
        
        # orderedy = y[orderx]
        # difsum = orderedy-c(orderedy[-1],0)
        # if(length(unique(y)) > 7){
        # if(min(difsum) == 0 & max(difsum)<3){
    
      
        if(prgrs){
        
          pdf( filename, 5, 4 )
          pname = paste( "raw", "p.value =", stestraw, sep=" ")
          plot_scRNA( x[ !omit ], y[ !omit ], main = pname, xlab = row1, ylab = row2)
          
          pname = paste( "cmp", "p.value =", stestcmp,sep=" ")
          plot_scRNA( cmpx[ !omit ], cmpy[ !omit ], main = pname, xlab = row1, ylab = row2)
          
          pname = paste( "sf", "p.value =", stestdsm, sep=" ")
          plot_scRNA( sfx[ !omit ], sfy[ !omit ], main = pname, xlab = row1, ylab = row2)
          
          pname = paste( "uq", "p.value =", stestuq, sep=" ")
          plot_scRNA( uqx[ !omit ], uqy[ !omit ], main = pname, xlab = row1, ylab = row2)
          
          pname = paste( "Down_Sample_Matrix", "p.value =", stestdsm, sep=" ")
          plot_scRNA( dsmx[! omit ], dsmy[! omit ], main = pname, xlab = row1, ylab = row2)
          
          dev.off()
        }
      }
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

source("normalize.R")
require("Matrix")

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

filter_rows = function(raw){
  lenuniq = function(x){
    return(length(unique(x)))
  }
  
  rsums = rowSums(raw)>50
  uniqrow = apply(raw, 1, lenuniq)
  return(rsums & uniqrow>5 & uniqrow<50)
}


preprocess = function(data, csize = 600, rsize = 10){
  
  keep_feature <- rowSums(data) > 0
  data <- data[keep_feature, ]
  
  keep_cell <- colSums(data) > 0
  data <- data[ ,keep_cell]
  return(data)
}

files <- list.files(path = "C:/Users/elisk/OneDrive/Plocha/Owncloud/documents/documents/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)


v = 1:(length(files)/2)
v = 1:1
i=1

setwd("rds")
  print(paste(i,16,sep="/"))
  scores = seq(0, 0, length.out=25)
  ind1 = seq(0, 0, length.out=25)
  ind2 = seq(0, 0, length.out=25)
  
  pvals.raw = c()
  pvals.50 = c()
  pvals.100 = c()
  pvals.500 = c()
  pvals.1000 = c()
  
  est.raw = c()
  est.50 = c()
  est.100 = c()
  est.500 = c()
  est.1000 = c()
  
  raw = readRDS(files[2*i])
  
  parts = strsplit(files[2*i],"/")
  name = strsplit(parts[[1]][length(parts[[1]])],"\\.")[[1]][1]
  
  
  
  raw[is.na(raw)] <- 0
  
  
  
  raw = preprocess(raw)
  keep_rows = filter_rows(raw)
  raw = raw[keep_rows,]
  
  udsm50 = readRDS(paste(name,"_udsm_50",".rds",sep=""))
  udsm50 = udsm50[keep_rows,]
  
  udsm100 = readRDS(paste(name,"_udsm_100",".rds",sep=""))
  udsm100 = udsm100[keep_rows,]
  
  udsm500 = readRDS(paste(name,"_udsm_500",".rds",sep=""))
  udsm500 = udsm500[keep_rows,]
  
  udsm1000 = readRDS(paste(name,"_udsm_1000",".rds",sep=""))
  udsm1000 = udsm1000[keep_rows,]
  

  
  for (row1 in row.names(raw)){
    x = raw[row1,]
    omit = x <= 0
    
    for (row2 in row.names(raw)){
      
      
      y = raw[row2,]
      omit = y<= 0 | omit | (x<=1 & y<=1)
      
      
      if( which( row.names( raw ) == row1 )[ 1 ] > which( row.names( raw ) == row2 )[1]){
        next  
      }
      
      stestraw = cor.test( x, y, method = "spearman", exact=F )
      prgrs = FALSE
      
      udsm50x = udsm50[row1,]
      udsm50y = udsm50[row2,]
      stest50 = cor.test( udsm50x, udsm50y, method = "spearman", exact=F )
      
      udsm100x = udsm100[row1,]
      udsm100y = udsm100[row2,]
      stest100 = cor.test( udsm100x, udsm100y, method = "spearman", exact=F )
      
      udsm500x = udsm500[row1,]
      udsm500y = udsm500[row2,]
      stest500 = cor.test( udsm500x, udsm500y, method = "spearman", exact=F )
      
      udsm1000x = udsm1000[row1,]
      udsm1000y = udsm1000[row2,]
      stest1000 = cor.test( udsm1000x, udsm1000y, method = "spearman", exact=F )
      
      
      pvals.raw = c( pvals.raw, stestraw [["p.value"]] )
      pvals.50 = c( pvals.50, stest50 [["p.value"]] )
      pvals.100 = c( pvals.100, stest100 [["p.value"]] )
      pvals.500 = c( pvals.500, stest500 [["p.value"]] )
      pvals.1000 = c( pvals.1000, stest1000 [["p.value"]] )
      
      est.raw = c( est.raw, stestraw [["estimate"]] )
      est.50 = c( est.50, stest50 [["estimate"]] )
      est.100 = c( est.100, stest100 [["estimate"]] )
      est.500 = c( est.500, stest500 [["estimate"]] )
      est.1000 = c( est.1000, stest1000 [["estimate"]] )
      
      
      
      pvals = data.frame(pvals.raw, pvals.cpm, pvals.sf, pvals.uq, pvals.dms, pvals.udms)
      est = data.frame(est.raw, est.cpm, est.sf, est.uq, est.dms, est.udms)
      
      
      
      maxsc = max(abs(stestraw[["estimate"]]-stest50[["estimate"]]),
                  abs(stestraw[["estimate"]]-stest100[["estimate"]]),
                  abs(stestraw[["estimate"]]-stest500[["estimate"]]),
                  abs(stestraw[["estimate"]]-stest1000[["estimate"]]))
      
      prgrs = maxsc > min(scores)
      
      if(prgrs){
        ind = which.min(scores)
        scores[ind] = maxsc
        ind1[ind] = row1
        ind2[ind] = row2
        
      }
      
    }
  }
  
  dir.create(paste(name,"sizes",sep = "_"))
  setwd(paste(name,"sizes",sep = "_"))
  
  osc = order(scores, decreasing = TRUE)
  saveRDS(pvals,"pvals.rds")
  saveRDS(est,"est.rds")
  
  for(j in 1:25){
    
    gene1 = ind1[osc[j]]
    gene2 = ind2[osc[j]]
    
    filename = paste( paste(j, gene1, gene2, sep = "_" ), "pdf", sep = "." )
    
    
    x = raw[gene1, ]
    omit = x <= 0
    
    y = raw[gene2, ]
    omit = y<= 0 | omit | (x<=1 & y<=1)
    
    stestraw = cor.test( x, y, method = "spearman", exact=F )[["estimate"]]
    
    
    udsm50x = udsm50[row1,]
    udsm50y = udsm50[row2,]
    stest50 = cor.test( udsm50x, udsm50y, method = "spearman", exact=F ) [["estimate"]]
    
    udsm100x = udsm100[row1,]
    udsm100y = udsm100[row2,]
    stest100 = cor.test( udsm100x, udsm100y, method = "spearman", exact=F ) [["estimate"]]
    
    udsm500x = udsm500[row1,]
    udsm500y = udsm500[row2,]
    stest500 = cor.test( udsm500x, udsm500y, method = "spearman", exact=F ) [["estimate"]]
    
    udsm1000x = udsm1000[row1,]
    udsm1000y = udsm1000[row2,]
    stest1000 = cor.test( udsm1000x, udsm1000y, method = "spearman", exact=F ) [["estimate"]]
    
    pdf( filename, 5, 4 )
    
    pname = paste( "raw", "estimate =", stestraw, sep=" ")
    plot_scRNA( x[ !omit ], y[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    pname = paste( "50", "estimate =", stest50,sep=" ")
    plot_scRNA( udsm50x[ !omit ], udsm50y[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    pname = paste( "100", "estimate =", stest100, sep=" ")
    plot_scRNA( udsm100x[ !omit ], udsm100y[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    
    pname = paste( "500", "estimate =", stest500, sep=" ")
    plot_scRNA( udsm500x[ !omit ], udsm500y[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    pname = paste( "1000", "estimate =", stest1000, sep=" ")
    plot_scRNA( udsm1000x[! omit ], udsm1000y[! omit ], main = pname, xlab = row1, ylab = row2)
    
    
    
    dev.off()
  }
  
  setwd("..")
  
  
}

setwd("..")

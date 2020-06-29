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

files <- list.files(path = "C:/Users/Ela/ownCloud/BIN/project/GSE139369", pattern = "\\.rds$", full.names = TRUE)


# # v = 1:(length(files)/2)
v = 1:1
i=1


for (i in v) {
  scores = seq(0, 0, length.out=25)
  ind1 = seq(0, 0, length.out=25)
  ind2 = seq(0, 0, length.out=25)
  
  parts = strsplit(files[2*i],"/")
  name = strsplit(parts[[1]][length(parts[[1]])],"\\.")[[1]][1]
  
  # t1 = r[[2*i-1]]
  t2 = readRDS(files[2*i])
  
  # t1[is.na(t1)] <- 0
  
  t2[is.na(t2)] <- 0
  
 
  
  raw = preprocess(t2,csize=600,rsize=20)
  
  cmp = calc_cpm(raw)
  
  sf = calc_sf(raw)
  
  uq = calc_uq(raw)
  
  dsm = Down_Sample_Matrix(raw)
  
  
  udsm = Up_Down_Sample_Matrix(raw)

  
  for (row1 in row.names(raw)){
    x = raw[row1,]
    omit = x <= 0
    
    # orderx = order(x)
    if(length(unique(x))<6){
      next
    }
    for (row2 in row.names(raw)){
    
      
      
      filename = paste(paste(row1,row2,sep="_"),"pdf",sep=".")

      y = raw[row2,]
      omit = y<= 0 | omit | (x<=1 & y<=1)
      
      
      if(length(unique(y))<6 || which( row.names( raw ) == row1 )[ 1 ] > which( row.names( raw ) == row2 )[1]){
        next  
      }
      
      # if(sum(!omit)>=4 && ( length( unique( x)) /length( unique( y)) > 0.5 && length( unique( x ))/ length( unique( y )) < 2)){
      if(sum(!omit)>=4){  
        stestraw = cor.test( x, y, method = "spearman", exact=F ) [["estimate"]]
        prgrs = FALSE
        
        cmpx = cmp[row1,]
        cmpy = cmp[row2,]
        stestcmp = cor.test( cmpx, cmpy, method = "spearman", exact=F ) [["estimate"]]
          
        sfx = sf[row1,]
        sfy = sf[row2,]
        stestsf = cor.test( sfx, sfy, method = "spearman", exact=F ) [["estimate"]]
        
        uqx = uq[row1,]
        uqy = uq[row2,]
        stestuq = cor.test( uqx, uqy, method = "spearman", exact=F ) [["estimate"]]
        
        dsmx = dsm[row1,]
        dsmy = dsm[row2,]
        stestdsm = cor.test( dsmx, dsmy, method = "spearman", exact=F ) [["estimate"]]
        
        
        maxsc = max(abs(stestraw-stestcmp),abs(stestraw-stestsf),abs(stestraw-stestuq),abs(stestraw-stestdsm))
      
        prgrs = maxsc > min(scores)
        
        if(prgrs){
          ind = which.min(scores)
          scores[ind] = maxsc
          ind1[ind] = row1
          ind2[ind] = row2
          
        }
      }
    }
  }
  
  setwd("filter")
  
  osc = order(scores, decreasing = TRUE)
  
  for(i in 1:25){
    
    gene1 = ind1[osc[i]]
    gene2 = ind2[osc[i]]
      
    filename = paste( paste(i, gene1, gene2, sep = "_" ), "pdf", sep = "." )
    
    
    x = raw[gene1, ]
    omit = x <= 0
    
    y = raw[gene2, ]
    omit = y<= 0 | omit | (x<=1 & y<=1)
    
    cmpx = cmp[gene1,]
    cmpy = cmp[gene2,]
    stestcmp = cor.test( cmpx, cmpy, method = "spearman", exact=F ) [["estimate"]]
    
    sfx = sf[gene1,]
    sfy = sf[gene2,]
    stestsf = cor.test( sfx, sfy, method = "spearman", exact=F ) [["estimate"]]
    
    uqx = uq[gene1,]
    uqy = uq[gene2,]
    stestuq = cor.test( uqx, uqy, method = "spearman", exact=F ) [["estimate"]]
    
    dsmx = dsm[gene1,]
    dsmy = dsm[gene2,]
    stestdsm = cor.test( dsmx, dsmy, method = "spearman", exact=F ) [["estimate"]]
    
    udsmx = dsm[gene1,]
    udsmy = dsm[gene2,]
    stestudsm = cor.test( udsmx, udsmy, method = "spearman", exact=F ) [["estimate"]]
    
    
    
    
    pdf( filename, 5, 4 )
    
    pname = paste( "raw", "estimate =", stestraw, sep=" ")
    plot_scRNA( x[ !omit ], y[ !omit ], main = pname, xlab = row1, ylab = row2)
  
    pname = paste( "cmp", "estimate =", stestcmp,sep=" ")
    plot_scRNA( cmpx[ !omit ], cmpy[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    pname = paste( "sf", "estimate =", stestdsm, sep=" ")
    plot_scRNA( sfx[ !omit ], sfy[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    
    pname = paste( "uq", "estimate =", stestuq, sep=" ")
    plot_scRNA( uqx[ !omit ], uqy[ !omit ], main = pname, xlab = row1, ylab = row2)
    
    pname = paste( "Down_Sample_Matrix", "estimate =", stestdsm, sep=" ")
    plot_scRNA( dsmx[! omit ], dsmy[! omit ], main = pname, xlab = row1, ylab = row2)
    
    pname = paste( "Up-Down_Sample_Matrix", "estimate =", stestudsm, sep=" ")
    plot_scRNA( udsmx[! omit ], udsmy[! omit ], main = pname, xlab = row1, ylab = row2)
    
    
    dev.off()
  }
  
  
  setwd("..")
  
}



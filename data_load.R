
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
for (i in v) {
  print(paste(i,16,sep="/"))
  # scores = seq(0, 0, length.out=25)
  # ind1 = seq(0, 0, length.out=25)
  # ind2 = seq(0, 0, length.out=25)
  # 
  # pvals.raw = c()
  # pvals.cpm = c()
  # pvals.sf = c()
  # pvals.uq = c()
  # pvals.dms = c()
  # pvals.udms = c()
  # 
  # est.raw = c()
  # est.cpm = c()
  # est.sf = c()
  # est.uq = c()
  # est.dms = c()
  # est.udms = c()
  
  raw = 0
  cpm = 0
  sf = 0
  uq = 0
  dsm = 0
  udsm50 = 0
  udsm100 = 0
  udsm500 = 0
  udsm1000 = 0
  namex = "test"
  namey = "test"
  estsd = data.frame(raw, cpm, sf, uq, dsm, udsm50, udsm100, udsm500,  udsm1000, namex, namey)
  raw = readRDS(files[2*i])
  
  parts = strsplit(files[2*i],"/")
  name = strsplit(parts[[1]][length(parts[[1]])],"\\.")[[1]][1]
  
  
  raw[is.na(raw)] <- 0
  
  raw = preprocess(raw)
  keep_rows = filter_rows(raw)
  raw = raw[keep_rows,]
  
  cpm = readRDS(paste(name,"_cmp",".rds",sep=""))
  cpm = cpm[keep_rows,]
  
  sf = readRDS(paste(name,"_sf",".rds",sep=""))
  sf = sf[keep_rows,]
  
  uq = readRDS(paste(name,"_uq",".rds",sep=""))
  uq = uq[keep_rows,]
  
  dsm = readRDS(paste(name,"_dsm",".rds",sep=""))
  dsm = dsm[keep_rows,]
  
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
    
    print(row1)
    for (row2 in row.names(raw)){
      
      
      y = raw[row2,]
      omit = y<= 0 | omit | (x<=1 & y<=1)
      
      
      if( which( row.names( raw ) == row1 )[ 1 ] > which( row.names( raw ) == row2 )[1]){
        next  
      }
      
      stestraw = cor.test( x, y, method = "spearman", exact=F )[["estimate"]]
      prgrs = FALSE
      
      cpmx = cpm[row1,]
      cpmy = cpm[row2,]
      stestcpm = cor.test( cpmx, cpmy, method = "spearman", exact=F )[["estimate"]]
      
      sfx = sf[row1,]
      sfy = sf[row2,]
      stestsf = cor.test( sfx, sfy, method = "spearman", exact=F )[["estimate"]]
      
      uqx = uq[row1,]
      uqy = uq[row2,]
      stestuq = cor.test( uqx, uqy, method = "spearman", exact=F )[["estimate"]]
      
      dsmx = dsm[row1,]
      dsmy = dsm[row2,]
      stestdsm = cor.test( dsmx, dsmy, method = "spearman", exact=F )[["estimate"]]
      
      
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
      
      # pvals.raw = c( pvals.raw, stestraw [["p.value"]] )
      # pvals.cpm = c( pvals.cpm, stestcpm [["p.value"]] )
      # pvals.sf = c( pvals.sf, stestsf [["p.value"]] )
      # pvals.uq = c( pvals.uq, stestuq [["p.value"]] )
      # pvals.dms = c( pvals.dms, stestdsm [["p.value"]] )
      # pvals.udms = c( pvals.udms, stestudsm [["p.value"]] )
      # 
      # est.raw = c( est.raw, stestraw [["estimate"]] )
      # est.cpm = c( est.cpm, stestcpm [["estimate"]] )
      # est.sf = c( est.sf, stestsf [["estimate"]] )
      # est.uq = c( est.uq, stestuq [["estimate"]] )
      # est.dms = c( est.dms, stestdsm [["estimate"]] )
      # est.udms = c( est.udms, stestudsm [["estimate"]] )
      
      estsd = rbind(estsd,c(stestraw, stestcpm, stestsf, stestuq, stestdsm, 
                            stest50, stest100, stest500, stest1000, row1, row2))
      
      # 
      # pvals = data.frame(pvals.raw, pvals.cpm, pvals.sf, pvals.uq, pvals.dms, pvals.udms)
      # est = data.frame(est.raw, est.cpm, est.sf, est.uq, est.dms, est.udms)
      
      
      # 
      # maxsc = max(abs(stestraw[["estimate"]]-stestcpm[["estimate"]]),
      #             abs(stestraw[["estimate"]]-stestsf[["estimate"]]),
      #             abs(stestraw[["estimate"]]-stestuq[["estimate"]]),
      #             abs(stestraw[["estimate"]]-stestdsm[["estimate"]]),
      #             abs(stestraw[["estimate"]]-stestudsm[["estimate"]]))
      # 
      # prgrs = maxsc > min(scores)
      # 
      # if(prgrs){
      #   ind = which.min(scores)
      #   scores[ind] = maxsc
      #   ind1[ind] = row1
      #   ind2[ind] = row2
      #   
      # }
      
    }
  }
  estsd = estsd[-1,]
  
  dir.create("filter")
  setwd("filter")
  ests.names = c("cpm", "sf", "uq", "dsm", "udsm50", "udsm100", "udsm500", "udsm1000")
  
  for(n in names){
    
    ordered = order(abs(estsd[,n]-estsd[,"raw"]), decreasing = TRUE)
    orderedF = order(abs(estsd[,n]-estsd[,"raw"]))
  
    for(j in 1:5){
      
      gene1 = ind1[orderedF[j], "namex"]
      gene2 = ind2[orderedF[j], "namey"]
      
      filename = paste( paste(n, "best", j, gene1, gene2, sep = "_" ), "pdf", sep = "." )
      
      
      x = raw[gene1, ]
      omit = x <= 0
      
      y = raw[gene2, ]
      omit = y<= 0 | omit | (x<=1 & y<=1)
      
      stestraw = cor.test( x, y, method = "spearman", exact=F )[["estimate"]]
      
      
      cpmx = cpm[gene1,]
      cpmy = cpm[gene2,]
      stestcpm = cor.test( cpmx, cpmy, method = "spearman", exact=F ) [["estimate"]]
      
      sfx = sf[gene1,]
      sfy = sf[gene2,]
      stestsf = cor.test( sfx, sfy, method = "spearman", exact=F ) [["estimate"]]
      
      uqx = uq[gene1,]
      uqy = uq[gene2,]
      stestuq = cor.test( uqx, uqy, method = "spearman", exact=F ) [["estimate"]]
      
      dsmx = dsm[gene1,]
      dsmy = dsm[gene2,]
      stestdsm = cor.test( dsmx, dsmy, method = "spearman", exact=F ) [["estimate"]]
      
      
      udsm50x = udsm50[gene1,]
      udsm50y = udsm50[gene2,]
      stest50 = cor.test( udsm50x, udsm50y, method = "spearman", exact=F ) [["estimate"]]
      
      udsm100x = udsm100[gene1,]
      udsm100y = udsm100[gene2,]
      stest100 = cor.test( udsm100x, udsm100y, method = "spearman", exact=F ) [["estimate"]]
      
      udsm500x = udsm500[gene1,]
      udsm500y = udsm500[gene2,]
      stest500 = cor.test( udsm500x, udsm500y, method = "spearman", exact=F ) [["estimate"]]
      
      udsm1000x = udsm1000[gene1,]
      udsm1000y = udsm1000[gene2,]
      stest1000 = cor.test( udsm1000x, udsm1000y, method = "spearman", exact=F ) [["estimate"]]
      
      pdf( filename, 5, 4 )
      
      pname = paste( "raw", "estimate =", stestraw, sep=" ")
      plot_scRNA( x[ !omit ], y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "cpm", "estimate =", stestcpm,sep=" ")
      plot_scRNA( cpmx[ !omit ], cpmy[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "sf", "estimate =", stestdsm, sep=" ")
      plot_scRNA( sfx[ !omit ], sfy[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      
      pname = paste( "uq", "estimate =", stestuq, sep=" ")
      plot_scRNA( uqx[ !omit ], uqy[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Down_Sample_Matrix", "estimate =", stestdsm, sep=" ")
      plot_scRNA( dsmx[! omit ], dsmy[! omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Up-Down_Sample_Matrix 50", "estimate =", stest50,sep=" ")
      plot_scRNA( udsm50x[ !omit ], udsm50y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Up-Down_Sample_Matrix 100", "estimate =", stest100, sep=" ")
      plot_scRNA( udsm100x[ !omit ], udsm100y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      
      pname = paste( "Up-Down_Sample_Matrix 500", "estimate =", stest500, sep=" ")
      plot_scRNA( udsm500x[ !omit ], udsm500y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Up-Down_Sample_Matrix 1000", "estimate =", stest1000, sep=" ")
      plot_scRNA( udsm1000x[! omit ], udsm1000y[! omit ], main = pname, xlab = gene1, ylab = gene2)
      
      
      dev.off()
      
      gene1 = ind1[orderedF[j], "namex"]
      gene2 = ind2[orderedF[j], "namey"]
      
      filename = paste( paste(n, "worst", j, gene1, gene2, sep = "_" ), "pdf", sep = "." )
      
      
      x = raw[gene1, ]
      omit = x <= 0
      
      y = raw[gene2, ]
      omit = y<= 0 | omit | (x<=1 & y<=1)
      
      stestraw = cor.test( x, y, method = "spearman", exact=F )[["estimate"]]
      
      
      cpmx = cpm[gene1,]
      cpmy = cpm[gene2,]
      stestcpm = cor.test( cpmx, cpmy, method = "spearman", exact=F ) [["estimate"]]
      
      sfx = sf[gene1,]
      sfy = sf[gene2,]
      stestsf = cor.test( sfx, sfy, method = "spearman", exact=F ) [["estimate"]]
      
      uqx = uq[gene1,]
      uqy = uq[gene2,]
      stestuq = cor.test( uqx, uqy, method = "spearman", exact=F ) [["estimate"]]
      
      dsmx = dsm[gene1,]
      dsmy = dsm[gene2,]
      stestdsm = cor.test( dsmx, dsmy, method = "spearman", exact=F ) [["estimate"]]
      
      
      udsm50x = udsm50[gene1,]
      udsm50y = udsm50[gene2,]
      stest50 = cor.test( udsm50x, udsm50y, method = "spearman", exact=F ) [["estimate"]]
      
      udsm100x = udsm100[gene1,]
      udsm100y = udsm100[gene2,]
      stest100 = cor.test( udsm100x, udsm100y, method = "spearman", exact=F ) [["estimate"]]
      
      udsm500x = udsm500[gene1,]
      udsm500y = udsm500[gene2,]
      stest500 = cor.test( udsm500x, udsm500y, method = "spearman", exact=F ) [["estimate"]]
      
      udsm1000x = udsm1000[gene1,]
      udsm1000y = udsm1000[gene2,]
      stest1000 = cor.test( udsm1000x, udsm1000y, method = "spearman", exact=F ) [["estimate"]]
      
      pdf( filename, 5, 4 )
      
      pname = paste( "raw", "estimate =", stestraw, sep=" ")
      plot_scRNA( x[ !omit ], y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "cpm", "estimate =", stestcpm,sep=" ")
      plot_scRNA( cpmx[ !omit ], cpmy[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "sf", "estimate =", stestdsm, sep=" ")
      plot_scRNA( sfx[ !omit ], sfy[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      
      pname = paste( "uq", "estimate =", stestuq, sep=" ")
      plot_scRNA( uqx[ !omit ], uqy[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Down_Sample_Matrix", "estimate =", stestdsm, sep=" ")
      plot_scRNA( dsmx[! omit ], dsmy[! omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Up-Down_Sample_Matrix 50", "estimate =", stest50,sep=" ")
      plot_scRNA( udsm50x[ !omit ], udsm50y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Up-Down_Sample_Matrix 100", "estimate =", stest100, sep=" ")
      plot_scRNA( udsm100x[ !omit ], udsm100y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      
      pname = paste( "Up-Down_Sample_Matrix 500", "estimate =", stest500, sep=" ")
      plot_scRNA( udsm500x[ !omit ], udsm500y[ !omit ], main = pname, xlab = gene1, ylab = gene2)
      
      pname = paste( "Up-Down_Sample_Matrix 1000", "estimate =", stest1000, sep=" ")
      plot_scRNA( udsm1000x[! omit ], udsm1000y[! omit ], main = pname, xlab = gene1, ylab = gene2)
      
      
      dev.off()
      
    }
  }
  setwd("..")
  
  
}

setwd("..")

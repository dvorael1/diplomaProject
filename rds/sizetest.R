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
  
  cpm = readRDS(paste(name,"_cmp",".rds",sep=""))
  
  sf = readRDS(paste(name,"_sf",".rds",sep=""))
  
  uq = readRDS(paste(name,"_uq",".rds",sep=""))
  
  dsm = readRDS(paste(name,"_dsm",".rds",sep=""))
  
  udsm50 = readRDS(paste(name,"_udsm_50",".rds",sep=""))
  
  udsm100 = readRDS(paste(name,"_udsm_100",".rds",sep=""))
  
  udsm500 = readRDS(paste(name,"_udsm_500",".rds",sep=""))
  
  udsm1000 = readRDS(paste(name,"_udsm_1000",".rds",sep=""))

  setwd("GSM4138872_scRNA_BMMC_D1T1_sizes")
    
  gene1 = "COPS3"
  gene2 = "PPIE"
  
  filename = paste( paste("notper", gene1, gene2, sep = "_" ), "pdf", sep = "." )
  
  
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
  

  setwd("..")
  
  l = dim(raw)[2]
  
  raw[gene1,] = raw[gene1,sample.int(l, l)]
  raw[gene2,] = raw[gene2,sample.int(l, l)]
  
  
  print("cpm")
  cpm = calc_cpm(raw)
  
  print("sf")
  sf = calc_sf(raw)
  
  print("uq")
  uq = calc_uq(raw)
  
  print("dsm")
  dsm = Down_Sample_Matrix(raw)
  
  print("50")
  udsm50 = Up_Down_Sample_Matrix(raw,scale = 50)
  
  print("100")
  udsm100 = Up_Down_Sample_Matrix(raw,scale = 100)
  
  print("500")
  udsm500 = Up_Down_Sample_Matrix(raw,scale = 500)
  
  print("1000")
  udsm1000 = Up_Down_Sample_Matrix(raw,scale = 1000)
  
  setwd("GSM4138872_scRNA_BMMC_D1T1_sizes")

  
  filename = paste( paste("per", gene1, gene2, sep = "_" ), "pdf", sep = "." )
  
  
  x = raw[gene1, ]
  omit = x <= 0
  
  y = raw[gene2, ]
  omit = y<= 0 | omit | (x<=1 & y<=1)
  
  stestraw = round(cor.test( x, y, method = "spearman", exact=F )[["estimate"]],digits = 4)
  
  
  cpmx = cpm[gene1,]
  cpmy = cpm[gene2,]
  stestcpm = round(cor.test( cpmx, cpmy, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  sfx = sf[gene1,]
  sfy = sf[gene2,]
  stestsf = round(cor.test( sfx, sfy, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  uqx = uq[gene1,]
  uqy = uq[gene2,]
  stestuq = round(cor.test( uqx, uqy, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  dsmx = dsm[gene1,]
  dsmy = dsm[gene2,]
  stestdsm = round(cor.test( dsmx, dsmy, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  
  udsm50x = udsm50[gene1,]
  udsm50y = udsm50[gene2,]
  stest50 = round(cor.test( udsm50x, udsm50y, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  udsm100x = udsm100[gene1,]
  udsm100y = udsm100[gene2,]
  stest100 = round(cor.test( udsm100x, udsm100y, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  udsm500x = udsm500[gene1,]
  udsm500y = udsm500[gene2,]
  stest500 = round(cor.test( udsm500x, udsm500y, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
  udsm1000x = udsm1000[gene1,]
  udsm1000y = udsm1000[gene2,]
  stest1000 = round(cor.test( udsm1000x, udsm1000y, method = "spearman", exact=F ) [["estimate"]],digits = 4)
  
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
  
  
  setwd("..")
  
}

setwd("..")

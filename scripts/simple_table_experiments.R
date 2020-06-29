source("data_generation.R")
source("methods.R")

#Pearsons Chisqrt
run_Chisq = function(dataset, output_file){
  
  fileConn <- file(output_file, open="w+")
  g1 = "X"
  g2 = "Y"
  for(tbl in dataset){
    pVals = chisq.test(tbl)$p.value
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
    pVals = chisq.test(t(tbl))$p.value
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
  }

  close(fileConn)
}

#correlation test
run_corr_test = function(dataset, output_file){
  require(DescTools)
  fileConn <- file(output_file, open="w+")
  g1 = "X"
  g2 = "Y"
  for(tbl in dataset){
    untbl = Untable(tbl)
    pVals =  cor.test( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )$p.value
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
    
    untbl = Untable(t(tbl))
    pVals =  cor.test( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )$p.value
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
  }

  
  close(fileConn)
}

#mutual information
run_muti = function(dataset, output_file){
  require(infotheo)
  fileConn <- file(output_file, open="w+")
  g1 = "X"
  g2 = "Y"
  for(tbl in dataset){
    untbl = Untable(tbl)
    pVals =  mutinformation( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
    
    untbl = Untable(t(tbl))
    pVals =  mutinformation( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
  }
  
  
  close(fileConn)
}


#conditional entropy
run_contropy = function(dataset, output_file){
  require(infotheo)
  fileConn <- file(output_file, open="w+")
  g1 = "X"
  g2 = "Y"
  for(tbl in dataset){
    untbl = Untable(tbl)
    pVals =  condentropy( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
    
    untbl = Untable(t(tbl))
    pVals =  condentropy( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
  }
  
  
  close(fileConn)
}

#FunChisq
run_FunChisq = function(dataset, output_file){
  
  require(FunChisq)
  fileConn <- file(output_file, open="w+")
  g1 = "X"
  g2 = "Y"
  
  for(tbl in dataset){
    pVals = fun.chisq.test(tbl)$p.value
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn)
    
    pVals = fun.chisq.test(t(tbl))$p.value
    s = paste(g1, g2, pVals, sep = "\t")
    writeLines(s, fileConn) 
  }
  
  close(fileConn)
}

#Loads the results from files and plot the AUROC
run_experiment = function(dataset, edges_gt, input_files, names, d = 0.5){
  source("AUC-multiple.R")
  list.stats = list()
  name = 0
  for( f in input_files){
    name = name + 1
    tp = 0
    fn = 0
    fp = 0
    print(paste(names[name]))
    edges_ex = numeric(NROW(x = edges_gt))
    conn <- file(f,open="r")
    lines <-readLines(conn)
    for (i in 1:length(lines)){
      lcont = unlist(strsplit(lines[i], "[\t]"))
      p = as.numeric(lcont[3])
      edges_ex[i] = -p
    }
    close(conn)
    list.stats = append(list.stats,list(as.vector(edges_ex)))
    
    tp = sum(which(edges_gt == 1) %in% which(edges_gt == edges_ex))
    fn = sum(which(edges_gt == 1) %in% which(edges_gt != edges_ex))
    fp = sum(which(edges_gt == 0) %in% which(edges_gt != edges_ex))
    tn = sum(which(edges_gt == 0) %in% which(edges_gt == edges_ex))
    tpr = tp/(tp+fn)
    fpr = fp/(fp+tp)
    print( paste("TPR: ", tpr, sep=" ") )
    print( paste("FPR: ", fpr, sep=" ") )
  }
  
  names(list.stats) = names
  plot.multiple.ROC.PR.curves(list.stats, list(edges_gt), plot=TRUE, main = d)
}


print("data generation")
dropout_rate = c(0.2, 0.9, 0.99)
data = gen_dataset()
# 
fileFC = "Fuchisq.txt"
fileChi = "Chi.txt"
fileCor = "cor.txt"
fileMuti = "muti.txt"
fileContropy = "contropy.txt"
pdf("performance.pdf",5,5)
for( d in dropout_rate){
  d_data = gen_data_dropout(data = data$data, d = d) 
  print(paste("dropout rate:", d, sep = " "))
  
  print("Running Correlation Test")
  experiment = run_corr_test(d_data, fileCor)
  print("Running ChiSqrt Test")
  experiment = run_Chisq(d_data, fileChi)
  print("Running Funchisq")
  experiment = run_FunChisq(d_data, fileFC)
  print("Running Mutual Information")
  experiment = run_muti(d_data, fileMuti)
  print("Running Conditional Entropy")
  experiment = run_contropy(d_data, fileContropy)
  
  print("Plotting results")
  run_experiment(data, data$gt, c(fileFC, fileChi, fileCor, fileMuti, fileContropy), 
                 c("FunChisq", "Chisq", "Corr","Muti","CondEtrp"),d = d)
}




# mutual information (symetrical) / conditional entropy (directional)
#non uniform distribution 1 2

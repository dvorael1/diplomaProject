libraries <- c("dynverse/dyngen")
#GITHUB_PAT=ca722b2d8907939730c3dc2f3b147ff5394350dc

for(lib in libraries) {
  if(grepl('/', lib)) { # github packages
    pkg <- strsplit(lib, '/')[[1]][2]
    if(! require(pkg, character.only = TRUE, quietly = TRUE)) {
      devtools::install_github(lib)
      require(pkg, character.only = TRUE)
    }
  } else {
    if(! require(lib, character.only = TRUE, quietly = TRUE)) {
      if(! require(BiocManager, quietly = TRUE)) {
        install.packages("BiocManager")
        require(BiocManager)
      }
      if(length(BiocManager::available(lib)) > 0) {
        BiocManager::install(lib)
      } else {
        install.packages(lib)
      }
      require(lib, character.only = TRUE)
    }
  }
}


#generates the dataset
generate_dataset = function(){
  require(dyngen)
  require(tidyverse)
  require(SingleCellExperiment)
  
  params <- simple_params
  options(ncores = 1)
  model <- invoke(generate_model_from_modulenet, params$model)
  
  
  simulation <- invoke(simulate_multiple, params$simulation, model$system)
  
  gs <- invoke(extract_goldstandard, params$gs, simulation, model)
  
  experiment <- invoke(run_experiment, params$experiment, simulation, gs)
  
  normalisation <- invoke(dynnormaliser::normalise_filter_counts, 
                          params$normalisation, experiment$counts)
  
  task <- wrap_dyngen_dataset("readme_dataset", params, model, simulation, 
                              gs, experiment, normalisation)
  
  gt = list(from=model$net$from, to=model$net$to,
            strength= model$net$strength, effect=model$net$effect)
  
  dataset = list(data = task$counts,
                 genes=task$feature_info$gene_id[startsWith(task$feature_info$gene_id, "G")],
                 gt=gt)
  
  return(dataset)
}

#discretization
data_discr = function(dataset){
  require(Ckmeans.1d.dp)
  
  discrete = 0
  genes = dataset$genes
  
  for( g1 in dataset$genes){
    
    x = dataset$data[, g1]
    res = Ckmeans.1d.dp(x, k = c(2, 9))
    xd <- res$cluster
  
    if (discrete == 0){
      discrete = data.frame(xd)
    }else{
      discrete = data.frame(discrete, xd)
    }
    
    
  }
  colnames(discrete) =  genes
  return(discrete)
}


#Pearsons Chisqrt
run_Chisq = function(dataset, output_file, discrete){
  
  fileConn <- file(output_file, open="w+")
  
  for( g1 in dataset$genes){
    
    for(g2 in dataset$genes){
      
      print( paste("Chisq",g1, g2, sep=" ") )
      
      tbl = table(discrete[, g1], discrete[, g2])
      pVals = chisq.test(tbl)$p.value
      
      s = paste(g1, g2, "1",pVals, sep = "\t")
      writeLines(s, fileConn)
    }
  }
  
  
  close(fileConn)
}

#correlation test
run_corr_test = function(dataset, output_file, discrete){
  
  fileConn <- file(output_file, open="w+")
  
  for( g1 in dataset$genes){
    
    for(g2 in dataset$genes){
      
      print( paste("correlation",g1, g2, sep=" ") )
      pVals = cor.test(dataset$data[, g1], dataset$data[, g2])$p.value #pouzit discrete data 
      
      s = paste(g1, g2, "1", (pVals), sep = "\t")
      writeLines(s, fileConn)
    }
  }
  
  
  close(fileConn)
}


#FunChisq
run_FunChisq = function(dataset, output_file, discrete){
  
  require(FunChisq)
  fileConn <- file(output_file, open="w+")
  
  for( g1 in dataset$genes){
    
    for(g2 in dataset$genes){
      
      print( paste("FunChisq", g1, g2, sep=" ") )
      
      tbl = table(discrete[, g1], discrete[, g2])
      
      result = fun.chisq.test(tbl)
      
      s = paste(g1, g2, result$estimate, result$p.value, sep = "\t")
      writeLines(s, fileConn)
    }
  }
  
  
  close(fileConn)
}

#Loads the results from files and plot the AUROC
plot_results = function(dataset,input_files,names){
  source("AUC-multiple.R")
  genes = dataset$genes
  size = NROW(genes)
  edges_gt = matrix(0, nrow=size, ncol=size)
  rownames(edges_gt) =  genes
  colnames(edges_gt) =  genes
  list.stats = c()
  for( i in 1:NROW(dataset$gt$from)){
    g1 = dataset$gt$from[i]
    g2 = dataset$gt$to[i]
    if(g1 %in% genes && g2 %in% genes){
      edges_gt[g1,g2] = 1
    }
  }
  list.stats = list()
  name = 0
  for( f in input_files){
    name = name + 1
    tp = 0
    fn = 0
    fp = 0
    print(paste(names[name]))
    edges_ex = matrix(0, nrow=size, ncol=size)
    rownames(edges_ex) =  genes
    colnames(edges_ex) =  genes
    conn <- file(f,open="r")
    lines <-readLines(conn)
    for (i in 1:length(lines)){
      lcont = unlist(strsplit(lines[i], "[\t]"))
      x = as.numeric(lcont[3])
      p = as.numeric(lcont[4])
      # if(x>0.5 && p<0.05){
      #   edges_ex[lcont[1],lcont[2]] = 1
      # }
      edges_ex[lcont[1],lcont[2]] = -p
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
  plot.multiple.ROC.PR.curves(list.stats, list(as.vector(edges_gt)), plot=TRUE)
}

#example for evaluation funtion run_FunChisq take a lot of time an data ganerated by generate_dataset function contains
#diffetent genes each run so for debugging purposes I suggest to execute each command outside the function.
evaluate_example = function(){
  print("data preparation") 
  data = generate_dataset()
  
  file = "Fuchisq.txt"
  print("Running Funchisq")
  experiment = run_FunChisq(data,file)
  print("Plotting results")
  plot_results(data,file)
  
}

#evaluate_example()

print("data preparation")
# data = generate_dataset()
# 
fileFC = "Fuchisq.txt"
fileChi = "Chi.txt"
fileCor = "cor.txt"
# discrete = data_discr(data)
# print("Running Correlation Test")
# experiment = run_corr_test(data,fileCor,discrete)
# print("Running ChiSqrt Test")
# experiment = run_Chisq(data,fileChi,discrete)
# print("Running Funchisq")
# experiment = run_FunChisq(data,fileFC,discrete)

print("Plotting results")
plot_results(data,c(fileFC,fileChi,fileCor),c("FunChisq","Chisq","Corr"))

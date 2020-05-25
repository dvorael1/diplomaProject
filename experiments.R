source("data_generation.R")
source("methods.R")
fileChi = "Chi.txt"
fileCor = "cor.txt"
fileMuti = "muti.txt" 
fileContropy = "contropy.txt"
fileFC = "Funchisq.txt"
files = c(fileChi,fileCor, fileMuti, fileContropy, fileFC)
# files = fileMuti
mtd_names = c( "Chisq", "Corr","MutI","CondEtrp","FunChisq")
# mtd_names = "MutI"

delete_files = function(files){
  for (f in files) {
    file.remove(f)
  }
}


run_experiment = function(dataset, edges_gt, input_files, names, title = 0.5){
  source("AUC-multiple.R")
  list.stats = list()
  name = 0
  list_tables = list()
  
  for( f in input_files){
    
    name = name + 1
    list_tables[[name]] = 1
    print(paste(names[name]))
    edges_ex = numeric(NROW(x = edges_gt))
    conn <- file(f,open="r")
    lines <-readLines(conn)
    min = 1
    for (i in 1:length(lines))
    {
      
      lcont = unlist(strsplit(lines[i], "[\t]"))
      p = as.numeric(lcont[3])
      
      # if(min>p){
      #   list_tables[[name]] = i
      #   min = p
      # }
        
      edges_ex[i] = -p
    }
    close(conn)
    list.stats = append(list.stats,list(as.vector(edges_ex)))

  }
  # for(i in 1:NROW(names))
  # {
  #   plot_table(table = dataset$data[[list_tables[[i]]]] , main = paste(names[i],dataset$names[[list_tables[[i]]]]))
  # 
  # }
  print("results are:")
  names(list.stats) = names
  ggplot.ROC.PR.curves(list.stats, list(edges_gt), plot=TRUE, title = title)
}

# for 1 config use only differentially expressed genes
process_real_data = function(pdf_name = "real.pdf"){
  delete_files(  files = files )
  
  pdf(pdf_name, 5,4 )


  print("data generation")
  data = gen_real_dataset(expressed = FALSE,sizeG = -1,sizeP = -1)

  for( tbl in data$data)
  {
    run_corr_test(tbl, fileCor)
    run_Chisq(tbl, fileChi)
    run_FunChisq(tbl, fileFC)
    run_muti(tbl, fileMuti)
    run_contropy(tbl, fileContropy)
  }


  print("Plotting results")
  run_experiment(data, data$gt, files, mtd_names, title = "1st configuration" )
  delete_files(  files = files )
  
  data = gen_real_dataset(expressed = TRUE,sizeG = -1, sizeP = -1)


  for( tbl in data$data)
  {
    run_corr_test(tbl, fileCor)
    run_Chisq(tbl, fileChi)
    run_FunChisq(tbl, fileFC)
    run_muti(tbl, fileMuti)
    run_contropy(tbl, fileContropy)
  }


  print("Plotting results")
  run_experiment(data, data$gt, files, mtd_names,title = "2nd configuration")
  delete_files(  files = files )

  
  dev.off()
  
}

process_simulated_data = function( dropout_rate = c(0.2, 0.8, 0.9, 0.99), noise_levels = c(0, 0.1, 0.2, 0.5, 1) ){
  delete_files(  files = files )
  
  pdf("simulated.pdf",7,5)
  
  print("data generation")
  for( noise in noise_levels){
    
    
    data = gen_simulated_dataset(directional = FALSE,noise = noise)
    
    for( d in dropout_rate)
    {
      
      d_data = gen_data_dropout(data = data$data, d = d) 
      print(paste("dropout rate:", d, sep = " "))
      
      for( tbl in d_data)
      {
        run_corr_test(tbl, fileCor)
        run_Chisq(tbl, fileChi)
        run_FunChisq(tbl, fileFC)
        run_muti(tbl, fileMuti)
        run_contropy(tbl, fileContropy)
      }
      
      
      print("Plotting results")
      run_experiment(data, data$gt, files, mtd_names,title = paste("2nd configuration dropout:",d,"noise:",noise,sep = " "))
      delete_files(  files = files )
    }
    
    print("data generation")
    data = gen_simulated_dataset(directional = TRUE,noise = noise)
    
    for( d in dropout_rate)
    {
      
      d_data = gen_data_dropout(data = data$data, d = d) 
      print(paste("dropout rate:", d, sep = " "))
      
      for( tbl in d_data)
      {
        run_corr_test(tbl, fileCor)
        run_Chisq(tbl, fileChi)
        run_FunChisq(tbl, fileFC)
        run_muti(tbl, fileMuti)
        run_contropy(tbl, fileContropy)
      }
      
      
      print("Plotting results")
      run_experiment(data, data$gt, files, mtd_names,title = paste("1st configuration dropout:",d,"noise:",noise,sep = " "))
      delete_files(  files = files )
    }
  }
  dev.off()
}

test_sample_size = function(samples = c(100,1000,10000,100000,1000000)){
  delete_files(  files = files )
  pdf("samples.pdf",7,5)
  
  for( s in samples)
  {
    
    data = gen_simulated_dataset(directional = FALSE,noise = 0.2,n = s)
    print(paste("sample size:", s, sep = " "))
    d_data = gen_data_dropout(data = data$data, d = 0.9) 
    
    for( tbl in d_data)
    {
      run_corr_test(tbl, fileCor)
      run_Chisq(tbl, fileChi)
      run_FunChisq(tbl, fileFC)
      run_muti(tbl, fileMuti)
      run_contropy(tbl, fileContropy)
    }
    
    
    print("Plotting results")
    run_experiment(data, data$gt, files, mtd_names,title = paste("2nd configuration sample size:",s,sep = " "))
    delete_files(  files = files )
  }
  
  
  for( s in samples)
  {
    
    data = gen_simulated_dataset(directional = TRUE,noise = 0.2,n = s)
    print(paste("sample size:", s, sep = " "))
    d_data = gen_data_dropout(data = data$data, d = 0.9) 
    
    for( tbl in d_data)
    {
      run_corr_test(tbl, fileCor)
      run_Chisq(tbl, fileChi)
      run_FunChisq(tbl, fileFC)
      run_muti(tbl, fileMuti)
      run_contropy(tbl, fileContropy)
    }
    
    
    print("Plotting results")
    run_experiment(data, data$gt, files, mtd_names,title = paste("1st configuration sample size:",s,sep = " "))
    delete_files(  files = files )
  }
  
  dev.off()
  
}


# process_simulated_data(dropout_rate = c(0.2,0.8,0.9,0.99),noise_levels = c(0, 0.1, 0.2, 0.5, 1))
process_real_data()
test_sample_size()

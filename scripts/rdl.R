load_GSE136369 = function(sizeP = -1, sizeG = -1, expressed = TRUE){
  DE <- read.table("tung\\TPs.txt")
  
  molecules <- read.table("tung\\molecules_normalized.txt", sep = "\t", header = TRUE)
  anno <- read.table("tung\\annotation.txt", sep = "\t", header = TRUE)
  p_names = unique(anno$individual)
  
  collection = filter_dataset(molecules, p_names, anno, DE, sizeG, sizeP, expressed = expressed)
  
  datadsc = data_discr(collection$dataset, dim = 1)
  
  gt = c()
  data = list()
  names = list()
  
  if(expressed)
  {
    for(g in 1:NROW(datadsc))
    { 
      tbl_data = list()
      l = 0
      for(n in 1:NROW(p_names))
      {
        f = l + 1
        l = l + NROW(collection$p_indexes[[n]])
        tbl_data[[n]] = datadsc[g,c(f:l)]
        
      }
      data[[g]] = create_table(tbl_data)
      names[[g]] = paste("Person to",collection$names[g],sep = " ")
    }
    gt = c(rep.int(1, collection$sizeGDE), rep.int(0, collection$sizeG))
  }
  else
  {
    i = 1
    for(g in 1:NROW(datadsc))
    { 
      tbl_data = list()
      l = 0
      for(n in 1:NROW(p_names))
      {
        f = l + 1
        l = l + NROW(collection$p_indexes[[n]])
        tbl_data[[n]] = datadsc[g,c(f:l)]
        
      }
      tbl = create_table(tbl_data)
      data[[ 2 * g - 1]] = tbl
      data[[ 2 * g ]] = t(tbl)
      names[[ 2 * g - 1]] = paste("Person","to",collection$names[g],sep = " ")
      names[[ 2 * g ]] = paste(collection$names[g],"to","Person",sep = " ")
      gt = c(gt, 1, 0)
    }
  }
  
  dataset = list(data = data ,gt = gt, names = names)
  # print_dataset(dataset)
  return(dataset)
}

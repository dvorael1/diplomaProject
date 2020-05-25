require(FunChisq)

gen_simulated_dataset = function(n = 1000, edges = 200, noise = 0.0, directional = TRUE , nrow = 3)
{
  set.seed(123457)
  gt = c()
  data = list()
  for( i in 1:edges)
  {
    if(directional)
    {
      if(i < 0.5*edges)
      {
        tbl = simulate_tables(n = n, type = "functional", noise = noise)
        gt = c(gt, 1)
        
      }
      else
      {
        tbl = simulate_tables(n = n, type = "independent", noise = noise)
        gt = c(gt, 0)
      }
      
      m = matrix(c(unlist(tbl$sample.list)),nrow=nrow)
      data[[ i ]] = as.table(m)
      
    }
    else
    {
      tbl = simulate_tables(n = n, type = "many.to.one",noise = noise)
      
      m = matrix(c(unlist(tbl$sample.list)), nrow = nrow)
      data[[ 2 * i - 1]] = as.table(m)
      data[[ 2 * i  ]] = as.table( t(m) )
      
      gt = c(gt,1,0)
    }
    
  }
  dataset = list(data = data ,gt = gt)
  return(dataset)
}

dropout = function(x, d = 0.0, n = 1000)
{
  indexes = sample(1:n, as.integer( d * n + 1) , replace=F)
  x[1:n %in% indexes] = 'A'
  return(x)
}

gen_data_dropout = function(data, n = 1000, d = 0.2 )
{
  require(DescTools)
  i = 1
  dataset = list()
  
  for( tbl in data)
  {
    untbl = Untable(tbl)
    x = dropout(untbl[,1],d = d, n = n)
    y = dropout(untbl[,2],d = d, n = n)
    dataset[[i]] = table(x,y)
    # dataset[[i]][1,1] = rbinom(1,dataset[[i]][1,1],1-d)
    dataset[[i]] = dataset[[i]] + 0.1
    i = i + 1
  } 
  return(dataset)
}

filter_dataset = function(dataset, p_names, anno, DE, sizeG = -1, sizePin = -1, expressed = FALSE)
{
  
  names_array = names(dataset)
  nc = c()
  for(n in names_array)
  {
    nc = c(nc,substr(n,1,7))
  }
  
  if(sizePin < 0)
  {
    sizePin = NCOL(dataset)
  }
  if (sizeG < 0)
  {
    sizeG = NROW(dataset)
  }
  
  p_indexes = c()
  person_indexes = list()
  i = 1
  for( n in p_names)
  {
    sizeP = min(sum( (nc %in% n)), sizePin)
    person_indexes[[i]] = which(nc %in% n)[1:sizeP]
    p_indexes = c(p_indexes, person_indexes[[i]])
    i = i + 1
  }
  
  if(expressed){
    sizeGDE = min(NROW(DE),0.5*sizeG)
    sizeG = sizeG - sizeGDE
    g_indexes = c( which(row.names.data.frame(dataset) %in% DE[ ,1])[1:sizeGDE], 
                   which(!(row.names.data.frame(dataset) %in% DE[ ,1]))[1:sizeG])
  }
  else
  {
    sizeGDE = min(NROW(DE),sizeG)
    sizeG = 0
    g_indexes = c(which(row.names.data.frame(dataset) %in% DE[ ,1])[1:(sizeGDE+sizeG)])
  }
  
  collection = list()
  collection$dataset = dataset[g_indexes,p_indexes]
  collection$sizeGDE = sizeGDE
  collection$p_indexes = person_indexes
  collection$sizeG = sizeG
  collection$names = row.names.data.frame(dataset)[g_indexes]
  
  return(collection)
}

gen_real_dataset = function(sizeP = -1, sizeG = -1, expressed = TRUE){
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

print_dataset = function(dataset){
  for( t in dataset$data){
    
    plot_table(table = t)
  } 
}

#discretization
data_discr = function(dataset, size = -1, dim = 2){
  require(Ckmeans.1d.dp)
  
  discrete = dataset
  
  if(dim == 1)
  {
    size = NROW(dataset)
  }
  
  if(dim == 2)
  {
    size = NCOL(dataset)
  }
  
  for( i in 1:size){
    
    if(dim == 1)
    {
      x = dataset[i, ]
    }
    else
    {
      x = dataset[, i]
    }
    
    res = Ckmeans.1d.dp(x, k = c(2,9))
    xd <- res$cluster
    
    if(dim == 1)
    {
      discrete[i, ] = xd
    }
    else
    {
      discrete[, i] = xd
    }
  }
  
  return(discrete)
}

create_table = function(dataset)
{
  size = NROW(dataset)
  m = max(as.numeric(dataset[[1]]))
  for( i in 2:size){
    m2 = max(as.numeric(dataset[[i]]))
    if(m < m2){
      m = m2
    }
  }
  t = matrix(nrow = size, ncol = m)
  
  for( i in 1:size){
    for( j in 1:m){
      t[i,j] = sum(dataset[[i]]==j)
    }
    
  }
  return(as.table(t))
}

# use boxplot
# data = gen_real_dataset(expressed = FALSE,sizeG = -1,sizeP = -1)

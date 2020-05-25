#Pearsons Chisqrt
run_Chisq = function(tbl, output_file , g1 = "X", g2 = "Y"){
  
  pVals = chisq.test(tbl)$p.value
  s = paste(g1, g2, pVals, sep = "\t")
  write(s, file = output_file, append = TRUE)

}

#correlation test
run_corr_test = function(tbl, output_file , g1 = "X", g2 = "Y"){
  require(DescTools)

  untbl = Untable(tbl)
  pVals =  cor.test( as.numeric(untbl[,1]), as.numeric(untbl[,2]), method ="pearson" )$p.value
  s = paste(g1, g2, pVals, sep = "\t")
  write(s, file = output_file, append = TRUE)

}

#mutual information
run_muti = function(tbl, output_file , g1 = "X", g2 = "Y"){
  require(infotheo)

  untbl = Untable(tbl)
  pVals =  mutinformation( as.numeric(untbl[,1]), as.numeric(untbl[,2]) )
  pVals = 1-pVals
  # print(paste("Muti value: ", pVals))
  s = paste(g1, g2, pVals, sep = "\t")
  write(s, file = output_file, append = TRUE)
    
}


#conditional entropy
run_contropy = function(tbl, output_file , g1 = "X", g2 = "Y"){
  require(infotheo)

  untbl = Untable(tbl)
  # rename the pVals variable because it is not p-value of conditional entropy it is result of the test
  pVals =  condentropy( as.numeric(untbl[,2]), as.numeric(untbl[,1]) )
  s = paste(g1, g2, pVals, sep = "\t")
  write(s, file = output_file, append = TRUE)
  
}

#FunChisq
run_FunChisq = function(tbl, output_file , g1 = "X", g2 = "Y"){
  require(FunChisq)

  pVals = fun.chisq.test(tbl)$p.value
  s = paste(g1, g2, pVals, sep = "\t")
  write(s, file = output_file, append = TRUE)
  
}

run_all = function(tbl, chisqr = "Chi.txt", cor = "cor.txt", muti = "muti.txt", 
                   con_ent = "contropy.txt", funchisq = "Funchisq.txt" ){
  
  run_Chisq(tbl = tbl, output_file = chisqr)
  run_corr_test(tbl = tbl, output_file = cor)
  run_muti(tbl = tbl, output_file = muti)
  run_contropy(tbl = tbl, output_file = con_ent)
  run_FunChisq(tbl = tbl, output_file = funchisq)
  
}

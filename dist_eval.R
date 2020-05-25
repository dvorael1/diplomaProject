# libraries <- c("devtools", "hemberg-lab/scRNA.seq.funcs", "edgeR", "monocle", 
#                "MAST", "ROCR", "nghiavtr/BPSC", "scde")
# 
# for(lib in libraries) {
#   
#   if(grepl('/', lib)) { # github packages
#     pkg <- strsplit(lib, '/')[[1]][2]
#     if(! require(pkg, character.only = TRUE, quietly = TRUE)) {
#       devtools::install_github(lib)
#       require(pkg, character.only = TRUE)
#     } 
#   } else {
#     if(! require(lib, character.only = TRUE, quietly = TRUE)) {
#       if(! require(BiocManager, quietly = TRUE)) {
#         install.packages("BiocManager")
#         require(BiocManager)
#       }
#       if(length(BiocManager::available(lib)) > 0) {
#         BiocManager::install(lib)
#       } else {
#         install.packages(lib)      
#       }
#       require(lib, character.only = TRUE)
#     }
#   }
# }
set.seed(1)

DE <- read.table("C:\\Users\\Ela\\Desktop\\owncloud\\BIN\\project\\TPs.txt")
notDE <- read.table("TNs.txt")
GroundTruth <- list(DE = as.character(unlist(DE)), notDE = as.character(unlist(notDE)))

molecules <- read.table("C:\\Users\\Ela\\Desktop\\owncloud\\BIN\\project\\molecules.txt", sep = "\t", header = TRUE)
anno <- read.table("C:\\Users\\Ela\\Desktop\\owncloud\\BIN\\project\\annotation.txt", sep = "\t", header = TRUE)
keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"
data <- molecules[,keep]
group <- anno[keep,1]
batch <- anno[keep,4]
# remove genes that aren't expressed in at least 6 cells
gkeep <- rowSums(data > 0, na.rm = TRUE) > 5
counts <- data[gkeep,]
# Library size normalization
lib_size = colSums(counts)
norm <- t(t(counts)/lib_size * median(lib_size)) 
# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules

pVals <- apply(norm, 1, function(x) {
  ks.test(x[group == "NA19101"], x[group == "NA19239"])$p.value } )
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")

sigDE <- names(pVals)[pVals < 0.05]
length(sigDE) 
# Number of KS-DE genes
sum(GroundTruth$DE %in% sigDE) 
sum(GroundTruth$notDE %in% sigDE)

tp <- sum(GroundTruth$DE %in% sigDE)
fp <- sum(GroundTruth$notDE %in% sigDE)
tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05])
fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05])
tpr <- tp/(tp + fn)
fpr <- fp/(fp + tn)
cat(c(tpr, fpr))

pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                 names(pVals) %in% GroundTruth$notDE] 
truth <- rep(1, times = length(pVals));
truth[names(pVals) %in% GroundTruth$DE] = 0;
pred <- ROCR::prediction(pVals, truth)
perf <- ROCR::performance(pred, "tpr", "fpr")
ROCR::plot(perf)
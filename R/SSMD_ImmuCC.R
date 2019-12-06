setwd('/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/New Server')

source("/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/New Server/SSMD find module.R")
source('/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/New Server/RMThreshold.R')
source("/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/New Server/CIBERSORT_modified.R")


##############################
#ImmuCC
ImmuCC <- function(sig_matrix, mixture_file){
  
  # Load the function of CIBERSORT
  #source("/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/New Server/CIBERSORT_modified.R")
  # Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/
  perm <- 100
  results <- CIBERSORT(sig_matrix, mixture_file, perm)
  results
}
##############################



# mixture_file <- read.table(expression,header=T,sep="\t",check.names=F)
# mixture_file=as.matrix(mixture_file)
# rownames(mixture_file) <- mixture_file[,1]
# mixture_file=mixture_file[,-1]
# mixture_file=as.data.frame(mixture_file)
# training_data='/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/testing data/ImmunCC/immuCC_old/srep40508-s1.txt'
# original_sig_matrix<- read.table(training_data,header=T,sep="\t",row.names=1,check.names=F)
# save(original_sig_matrix,file='ImmuCC_sig_matrix.RData')

expression='/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/New Server/SmallIntestine.txt'
expression_data=read.table(expression)

SSMD=SSMD_find_module(expression_data)

SSMD_module=SSMD

###############
SSMD_modules_plain <- list()
nn <- c()
N <- 0
for (i in 1:length(SSMD_module)) {
  for (j in 1:length(SSMD_module[[i]])) {
    N <- N + 1
    SSMD_modules_plain[[N]] <- SSMD_module[[i]][[j]]
  }
  nn <- c(nn, names(SSMD_module[[i]]))
}
names(SSMD_modules_plain) <- nn

#########
SSMD_gene=c()
for (i in 1:length(SSMD_modules_plain)) {
  SSMD_gene=c(SSMD_gene,SSMD_modules_plain[[i]])
}

# library(VennDiagram)
# 
# venn.diagram(x = list(ImmuCC = row.names(original_sig_matrix),SSMD = SSMD_core_marker), filename = "venn",
#              main = 'Genes', cex = 2, cat.cex = 1.3,fill = c("orange","violet"))

aaa=intersect(rownames(original_sig_matrix),SSMD_gene)
modify_sig_matrix=original_sig_matrix[aaa,]

ImmuCC.proportion_modify <- ImmuCC (modify_sig_matrix,expression_data)

#ImmuCC.proportion_original <- ImmuCC (original_sig_matrix,expression_data)
#save(Immune.proportion,file='/Users/xiaoyulu/Documents/RA/NEW/201811_NMF/Pipeline/testing data/ImmunCC/new_ImmuCC_1003/newImmuCC_SmallIntestine.RData')


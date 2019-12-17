.onLoad <- function(libname, pkgname) {
  utils::data(SSMD_ImmuCC_core_markers, package = pkgname, envir = parent.env(environment()))
  SSMD_ImmuCC_core_markers <- SSMD::SSMD_ImmuCC_core_markers
  assign("SSMD_ImmuCC_core_markers", SSMD_ImmuCC_core_markers, envir = parent.env(environment()))

  utils::data(SSMD_ImmuCC_labeling_genes, package = pkgname, envir = parent.env(environment()))
  SSMD_ImmuCC_labeling_genes <- SSMD::SSMD_ImmuCC_labeling_genes
  assign("SSMD_ImmuCC_labeling_genes", SSMD_ImmuCC_labeling_genes, envir = parent.env(environment()))

  utils::data(original_sig_matrix, package = pkgname, envir = parent.env(environment()))
  original_sig_matrix <- SSMD::original_sig_matrix
  assign("original_sig_matrix", original_sig_matrix, envir = parent.env(environment()))
}

ImmuCC_modify <- function(mixture_file,sig_matrix){
  
  # Load the function of CIBERSORT
  #source("CIBERSORT_modified.R")
  # Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/
  perm <- 100
  results <- CIBERSORT(sig_matrix, mixture_file, perm)
  results
}

SSMD_ImmuCC <- function(data11) {
  
  SSMD_module=SSMD_find_module_ImmuCC(data11)[[1]]
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
  
  SSMD_matrix=as.matrix(data11[which(rownames(data11)%in%SSMD_gene),])
  ImmuCC_original_matrix=data11[which(rownames(data11)%in% rownames(original_sig_matrix) ),]
  
  cor_matrix=cor(t(ImmuCC_original_matrix),t(SSMD_matrix))
  cor_matrix[,1]=apply(cor_matrix, 1, function(x) max(x))
  modify_sig_gene=rownames(cor_matrix)[which(cor_matrix[,1]>0.8)]
  
  #aaa=intersect(rownames(original_sig_matrix),SSMD_gene)
  modify_sig_matrix=original_sig_matrix[modify_sig_gene,]
  
  ImmuCC.proportion_modify <- ImmuCC_modify (data11,modify_sig_matrix)
  predict_p <- ImmuCC.proportion_modify
  predict_sig = modify_sig_matrix
  return(predict_p,predict_sig)
}



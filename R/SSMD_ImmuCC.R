.onLoad <- function(libname, pkgname) {
  utils::data(labeling_matrix, package = pkgname, envir = parent.env(environment()))
  SSMD_ImmuCC_core_markers <- SSMD::SSMD_ImmuCC_core_markers
  assign("SSMD_ImmuCC_core_markers", SSMD_ImmuCC_core_markers, envir = parent.env(environment()))

  utils::data(SSMD_ImmuCC_labeling_genes, package = pkgname, envir = parent.env(environment()))
  SSMD_ImmuCC_labeling_genes <- SSMD::SSMD_ImmuCC_labeling_genes
  assign("SSMD_ImmuCC_labeling_genes", SSMD_ImmuCC_labeling_genes, envir = parent.env(environment()))

  utils::data(original_sig_matrix, package = pkgname, envir = parent.env(environment()))
  original_sig_matrix <- SSMD::original_sig_matrix
  assign("original_sig_matrix", original_sig_matrix, envir = parent.env(environment()))
}

ImmuCC_modify <- function(sig_matrix,mixture_file){
  
  # Load the function of CIBERSORT
  #source("CIBERSORT_modified.R")
  # Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/
  perm <- 100
  results <- CIBERSORT(sig_matrix, mixture_file, perm)
  results
}

SSMD_ImmuCC <- function(data11) {
  
  SSMD_module=SSMD_find_module(data11)[[1]]
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
  
  aaa=intersect(rownames(original_sig_matrix),SSMD_gene)
  modify_sig_matrix=original_sig_matrix[aaa,]
  
  ImmuCC.proportion_modify <- ImmuCC_modify (modify_sig_matrix,data11)
  return(ImmuCC.proportion_modify)
}



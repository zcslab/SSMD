#ImmuCC
.onLoad <- function(libname, pkgname) {
  utils::data(original_sig_matrix, package = pkgname, envir = parent.env(environment()))
  original_sig_matrix <- SSMD::original_sig_matrix
  assign("original_sig_matrix", original_sig_matrix, envir = parent.env(environment()))
}

ImmuCC <- function(data11){

  # Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/
  perm <- 100
  results <- CIBERSORT(sig_matrix=original_sig_matrix, mixture_file=data11, perm)
  
  predict_p <- results
  predict_sig = as.matrix(original_sig_matrix)
  LM <- cal_Zscore_small(predict_sig)
  sig_gene_list <- find_1_genelist(LM)
  
  #E-Score
  e_mat <- cal_escore(predict_sig, t(predict_p), data11)
  #list(predict_p = proportion_matrix,sig_gene_list = module_keep_plain)
  return(list(SigMat=predict_sig, ProMat=t(predict_p), mk_gene=sig_gene_list,Escore_vector=e_mat))
  
}

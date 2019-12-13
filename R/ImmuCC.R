#ImmuCC
.onLoad <- function(libname, pkgname) {
  utils::data(original_sig_matrix, package = pkgname, envir = parent.env(environment()))
  original_sig_matrix <- SSMD::original_sig_matrix
  assign("original_sig_matrix", original_sig_matrix, envir = parent.env(environment()))
}

ImmuCC <- function(mixture_file){
  
  # Load the function of CIBERSORT
  #source("CIBERSORT_modified.R")
  # Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/
  perm <- 100
  results <- CIBERSORT(original_sig_matrix, mixture_file, perm)
  results
}
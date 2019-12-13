#ImmuCC
ImmuCC <- function(sig_matrix, mixture_file){
  
  # Load the function of CIBERSORT
  #source("CIBERSORT_modified.R")
  # Note: the scirpts of CIBERSORT.R is a method developed by Newman et al.and can be accesssed upon an request from https://cibersort.stanford.edu/
  perm <- 100
  results <- CIBERSORT(sig_matrix, mixture_file, perm)
  results
}
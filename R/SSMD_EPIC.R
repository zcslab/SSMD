#hint
#when we use Epic tool, we should divided them into two parts blood and immune cells
#we can write two functions for keep one parameter to distinguish them

#EPIC
.onLoad <- function(libname, pkgname) {
  utils::data(Mouse_human_mapping, package = pkgname, envir = parent.env(environment()))
  Mouse_human_mapping <- SSMD::Mouse_human_mapping
  assign("Mouse_human_mapping", Mouse_human_mapping, envir = parent.env(environment()))
}

SSMD_EPIC <- function(data11,tissue) {
  SSMD_module=SSMD_find_module_EPIC(data11,tissue)[[1]]
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
  
  library(EPIC)
  
  if (tissue=='Blood')
    reference=EPIC::BRef
  if (tissue=='Inflammatory')
    reference=EPIC::TRef
  
  data_c1=data11
  SSMD_matrix=as.matrix(data_c1[which(rownames(data_c1)%in%SSMD_gene),])
  common_gene <- intersect(rownames(data_c1), Mouse_human_mapping[,2])
  data_c2 <- data_c1[common_gene,]
  human_mouse_mapping2 <- Mouse_human_mapping
  rownames(human_mouse_mapping2) <- human_mouse_mapping2[,2]
  ccc <- human_mouse_mapping2[rownames(data_c2),]
  data_c2=as.matrix(data_c2)
  rownames(data_c2) <- ccc[,1]
  data_c3 <- data_c2[unique(rownames(data_c2)),]
  EPIC_original_matrix=data_c3[which(rownames(data_c3)%in%reference$sigGenes),]
  
  cor_matrix=cor(t(EPIC_original_matrix),t(SSMD_matrix))
  cor_matrix[,1]=apply(cor_matrix, 1, function(x) max(x))
  reference$sigGenes=rownames(cor_matrix)[which(cor_matrix[,1]>0.8)]
  
  out <- EPIC(data_c3,reference)
  predict_p <- as.matrix(out[[2]])
  gene_name=intersect(rownames(Mouse_human_mapping),reference$sigGenes)
  predict_sig = as.matrix(reference$refProfiles[gene_name,])
  rownames(predict_sig)=Mouse_human_mapping[rownames(predict_sig),2]
  LM <- cal_Zscore_small(predict_sig)
  sig_gene_list <- find_1_genelist(LM)
  

  #E-Score
  e_mat <- SSMD_cal_escore(predict_sig, t(predict_p), data11)
  #list(predict_p = proportion_matrix,sig_gene_list = module_keep_plain)
  return(list(SigMat=predict_sig, ProMat=t(predict_p), mk_gene=sig_gene_list,Escore_vector=e_mat))
  
}

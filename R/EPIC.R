#EPIC
.onLoad <- function(libname, pkgname) {
  utils::data(Mouse_human_mapping, package = pkgname, envir = parent.env(environment()))
  Mouse_human_mapping <- SSMD::Mouse_human_mapping
  assign("Mouse_human_mapping", Mouse_human_mapping, envir = parent.env(environment()))
}

my_EPIC <- function(data11,tissue){
  library(EPIC)
  
  if (tissue=='Blood')
    reference=EPIC::BRef
  if (tissue=='Inflammatory')
    reference=EPIC::TRef

  data_c1=data11
  common_gene <- intersect(rownames(data_c1), Mouse_human_mapping[,2])
  data_c2 <- data_c1[common_gene,]
  human_mouse_mapping2 <- Mouse_human_mapping
  rownames(human_mouse_mapping2) <- human_mouse_mapping2[,2]
  ccc <- human_mouse_mapping2[rownames(data_c2),]
  data_c2=as.matrix(data_c2)
  rownames(data_c2) <- ccc[,1]
  data_c3 <- data_c2[unique(rownames(data_c2)),]
  
  out <- EPIC(data_c3,reference)
  predict_p <- as.matrix(out[[2]])
  gene_name=intersect(rownames(Mouse_human_mapping),reference$sigGenes)
  predict_sig = as.matrix(reference$refProfiles[gene_name,])
  rownames(predict_sig)=Mouse_human_mapping[rownames(predict_sig),2]
  LM <- cal_Zscore_small(predict_sig)
  sig_gene_list <- find_1_genelist(LM)
  
  list(predict_p,predict_sig,sig_gene_list)
}

.onLoad <- function(libname, pkgname) {
  utils::data(Mouse_Cancer_labeling_matrix, package = pkgname, envir = parent.env(environment()))
  Mouse_Cancer_labeling_matrix <- SSMD::Mouse_Cancer_labeling_matrix
  assign("Mouse_Cancer_labeling_matrix", Mouse_Cancer_labeling_matrix, envir = parent.env(environment()))
  
  utils::data(Mouse_Cancer_core_marker, package = pkgname, envir = parent.env(environment()))
  Mouse_Cancer_core_marker <- SSMD::Mouse_Cancer_core_marker
  assign("Mouse_Cancer_core_marker", Mouse_Cancer_core_marker, envir = parent.env(environment()))
  
  utils::data(Mouse_hematopoietic_labeling_matrix, package = pkgname, envir = parent.env(environment()))
  Mouse_hematopoietic_labeling_matrix <- SSMD::Mouse_hematopoietic_labeling_matrix
  assign("Mouse_hematopoietic_labeling_matrix", Mouse_hematopoietic_labeling_matrix, envir = parent.env(environment()))
  
  utils::data(Mouse_hematopoietic_core_marker, package = pkgname, envir = parent.env(environment()))
  Mouse_hematopoietic_core_marker <- SSMD::Mouse_hematopoietic_core_marker
  assign("Mouse_hematopoietic_core_marker", Mouse_hematopoietic_core_marker, envir = parent.env(environment()))
  
  utils::data(Mouse_Brain_labeling_matrix, package = pkgname, envir = parent.env(environment()))
  Mouse_Brain_labeling_matrix <- SSMD::Mouse_Brain_labeling_matrix
  assign("Mouse_Brain_labeling_matrix", Mouse_Brain_labeling_matrix, envir = parent.env(environment()))
  
  utils::data(Mouse_Brain_core_marker, package = pkgname, envir = parent.env(environment()))
  Mouse_Brain_core_marker <- SSMD::Mouse_Brain_core_marker
  assign("Mouse_Brain_core_marker", Mouse_Brain_core_marker, envir = parent.env(environment()))
}


SSMD <- function(data11,tissue) {
  
  BCV_ttest2 <- function(data0, rounds = 20, slice0 = 2, maxrank0 = 4, msep_cut = 0.01) {
    x <- data0
    fff_cc <- c()
    for (kk in 1:rounds) {
      cv_result <- bcv::cv.svd.gabriel(x, slice0, slice0, maxrank = maxrank0)
      fff_cc <- rbind(fff_cc, cv_result$msep)
    }
    fff_cc[is.na(fff_cc)]=0
    pp <- c()
    ddd <- apply(fff_cc, 2, mean)
    ddd <- ddd/sum(ddd)
    for (kk in 1:(ncol(fff_cc) - 1)) {
      pp_c <- 1
      if (mean(fff_cc[, kk], na.rm = T) > mean(fff_cc[, kk + 1], na.rm = T)) {
        if (ddd[kk] > msep_cut) {
          pp_c <- stats::t.test(fff_cc[, kk], fff_cc[, kk + 1])$p.value
        }
      }
      pp <- c(pp, pp_c)
    }
    return(list(pp, fff_cc))
  }
  
  
  ############################
  
  # caculate the base in selected list
  Compute_Rbase_SVD <- function(bulk_data, tg_R1_lists_selected) {
    tg_R1_lists_st_ccc <- tg_R1_lists_selected
    data_c <- bulk_data
    Base_all <- c()
    for (i in 1:length(tg_R1_lists_st_ccc)) {
      tg_data_c <- data_c[tg_R1_lists_st_ccc[[i]], ]
      cc <- svd(tg_data_c)$v[, 1]
      ccc <- stats::cor(cc, t(tg_data_c))
      if (mean(ccc) < 0) {
        cc <- -cc
      }
      Base_all <- rbind(Base_all, cc)
    }
    rownames(Base_all) <- 1:nrow(Base_all)
    if (length(names(tg_R1_lists_selected)) > 0) {
      rownames(Base_all) <- names(tg_R1_lists_selected)
    }
    return(Base_all)
  }
  
  
  #################
  if (tissue=='Inflammatory'){
    tg_core_marker_set=Mouse_Cancer_core_marker
    marker_stats1_uni=Mouse_Cancer_labeling_matrix
  }
  if (tissue=='Central Nervous System'){
    tg_core_marker_set=Mouse_Brain_core_marker
    marker_stats1_uni=Mouse_Brain_labeling_matrix
  }  
  if (tissue=='Hematopoietic System'){
    tg_core_marker_set=Mouse_hematopoietic_core_marker
    marker_stats1_uni=Mouse_hematopoietic_labeling_matrix
  }
  cell_type = names(tg_core_marker_set)
  i = 1
  intersect_marker1 = vector("list")
  for (cell in cell_type) {
    name = cell
    if (cell == "T") {
      tg_marker <- names(which(marker_stats1_uni[, "CD4T"] >= 1 | marker_stats1_uni[, "CD8T"] >= 1 | marker_stats1_uni[, "T"] >= 1))
    } else {
      tg_marker <- names(which((marker_stats1_uni[, cell] >= 1)))
    }
    intersect_marker1[[i]] = intersect(tg_marker, rownames(data11))
    names(intersect_marker1)[[i]] = name
    i = i + 1
  }
  
  #######################
  intersect_marker1_choose = vector("list", length = length(intersect_marker1))
  for (i in 1:length(intersect_marker1)) {
    intersect_marker1_choose[i] = intersect_marker1[i]
  }
  names(intersect_marker1_choose) = names(intersect_marker1)
  intersect_marker1_choose[sapply(intersect_marker1_choose, is.null)] = NULL
  
  
  ################
  
  infor.list = vector("list", length(intersect_marker1_choose))
  marker_modules = vector("list")
  
  for (i in 1:length(intersect_marker1_choose)) {
    
    name=names(intersect_marker1_choose)[i]
    data_c<-data11[intersect(rownames(data11),intersect_marker1_choose[[i]]),]
    
    corr=cor(t(data_c))
    corr[sapply(corr, is.na)] = 0
    
    
    if(dim(corr)[1]<100){
      thr=0.6
    }else{
      ###gene size must be large enough
      res <- rm.get.threshold(corr,interactive =F,plot.spacing =F,plot.comp =F,save.fit=F,interval=c(0.4,max(abs(corr[which(corr!=1)]))))
      
      #suppressWarnings()
      dis=res$tested.thresholds[which(res$dist.Expon>res$dist.Wigner & res$tested.thresholds>0.6)][1]
      if ( is.na(dis) ){
        dis=0
      }
      p.ks=res$tested.thresholds[which(res$p.ks>0.05)][1]
      if ( is.na(p.ks) ){
        p.ks=0
      }
      thr=max(dis,p.ks,0.6)
    }
    ######
    print('##################')
    print(thr)
    print('##################')
    
    cleaned.matrix <- rm.denoise.mat(corr, threshold = thr, keep.diag = TRUE)
    clust=hclust(dist(cleaned.matrix))
    
    written_list=rep(0, dim(corr)[1])
    names(written_list)=row.names(corr)
    n=1
    cut_value=2
    t=cutree(clust, k = cut_value)
    keep_k=vector("list",cut_value)
    marker_modules_cell_type=vector("list")
    marker_modules_length=1
    
    while(cut_value<length(clust$order))
    {
      t=cutree(clust, k = cut_value)
      for (k in 1:cut_value) {
        d=t[which(t==k)]
        mean=mean(abs(corr[names(d),names(d)]))
        #print(mean)
        
        if ( mean>=thr & length(d) >= 6 ){
          if (sum(written_list[names(d)])==0){
            keep_sample=names(d)
            written_list[keep_sample]=n
            n=n+1
            # print(d)
            # print(mean)
            keep_corr=mean
            keep_sample=names(d)
            marker_modules_cell_type[[marker_modules_length]]=names(d)
            marker_modules_length=marker_modules_length+1
            keep_k[[k]]=keep_sample
            print(mean)
          }
        }
      }  
      
      #t=cutree(clust, k = cut_value)
      cut_value=cut_value+1
    }
    
    keep_k[sapply(keep_k, is.null)] = NULL
    
    infor.list[[i]]=list(name, keep_k)
    marker_modules[[i]]=marker_modules_cell_type
    if(length(marker_modules_cell_type)!=0){
      for (t in 1:length(marker_modules[[i]])) {
        names(marker_modules[[i]])[t]=paste(name,t,sep = '_')
      }
    }else{
      marker_modules[[i]]=NULL
    }
    
  }
  marker_modules[sapply(marker_modules, is.null)] = NULL
  
  ###############
  marker_modules_plain <- list()
  nn <- c()
  N <- 0
  for (i in 1:length(marker_modules)) {
    for (j in 1:length(marker_modules[[i]])) {
      N <- N + 1
      marker_modules_plain[[N]] <- marker_modules[[i]][[j]]
    }
    nn <- c(nn, names(marker_modules[[i]]))
  }
  names(marker_modules_plain) <- nn
  
  ################
  Stat_all <- as.data.frame(nn)
  aa <- c()
  for (i in 1:length(nn)) {
    aa <- rbind(aa, unlist(strsplit(nn[i], "_")))
  }
  Stat_all$CT <- aa[, 1]
  Stat_all$CTN <- as.numeric(aa[, 2])
  colnames(Stat_all)[1] <- "ID"
  
  mean_value = list()
  Core_overlap_number = list()
  Core_overlap_rate = list()
  BCV_rank = list()
  
  for (i in 1:length(marker_modules_plain)) {
    ############
    data0 = data11[marker_modules_plain[[i]], ]
    corr = stats::cor(t(data0))
    mean_value[[i]] = mean(corr)
    #####################
    gene_name = sapply(names(marker_modules_plain)[i], function(y) strsplit(y, split = "_")[[1]][[1]])
    core_marker = tg_core_marker_set[[which(names(tg_core_marker_set) == gene_name)]]
    interaction_marker = intersect(core_marker, marker_modules_plain[[i]])
    Core_overlap_number[[i]] = length(interaction_marker)
    Core_overlap_rate[[i]] = (length(interaction_marker)/length(marker_modules_plain[[i]]))
    BCV = BCV_ttest2(data0, maxrank0 = 10)
    BCV_rank[[i]] = mean(BCV[[2]][, 1]/apply(BCV[[2]], 1, sum))
  }
  
  Stat_all$mean = mean_value
  Stat_all$Core_overlap_number = Core_overlap_number
  Stat_all$Core_overlap_rate = Core_overlap_rate
  Stat_all$BCV_rank = BCV_rank
  
  j = 1
  module_keep = vector("list")
  for (module_cell in unique(Stat_all$CT)) {
    aaa = Stat_all[which(Stat_all$CT == module_cell), ]
    bbb = aaa$ID[which((aaa$Core_overlap_number == max(max(unlist(aaa$Core_overlap_number)),2))|(aaa$Core_overlap_number>=10)
                       |((aaa$Core_overlap_number>=5)&(aaa$Core_overlap_rate>=0.5)))]
    module_keep[[j]] =marker_modules_plain[which(names(marker_modules_plain) %in% bbb)]
    j = j + 1
  }
  #module_keep[sapply(module_keep, is.null)] = NULL
  
  
  # aaa <- Compute_Rbase_SVD(data11, module_keep) get propotion for selected modules
  proportion = vector("list")
  for (i in 1:length(module_keep)) {
    # name=colnames(marker_stats1_uni)[i] print(name)
    if (length(module_keep[[i]]) > 0) {
      my_list <- list()
      my_list <- module_keep[[i]]
      ###### aaa is base: estimated propotion
      aaa <- Compute_Rbase_SVD(data11, my_list)
      rownames(aaa) = sapply(rownames(aaa), function(y) strsplit(y, split = "_")[[1]][[1]])
      proportion[[i]] = aaa
      names(proportion)[i] = sapply(rownames(aaa), function(y) strsplit(y, split = "_")[[1]][[1]])
      # dim(aaa) dim(tProp) correlation between estimated and true propotion ccc <- cor(t(aaa),t(tProp)) modules_cor_tPro[[i]]=ccc print(ccc)
    } else {
      print("NO Marker")
    }
  }
  #list(Stat_all = Stat_all, module_keep = module_keep, proportion = proportion)
  list(predict_p = proportion,sig_gene_list = module_keep)
  
}




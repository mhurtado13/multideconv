
##### Functions for running and processing cell-type deconvolution for pipeline multideconv: https://github.com/mhurtado13/multideconv

# author: Marcelo Hurtado
# email: marcelo.hurtado@inserm.fr
# organization: INSERM CRCT - Pancaldi team 21
# place: Toulouse, France


libraries_set <- function(){
  suppressMessages(library("BiocManager"))
  suppressMessages(library("devtools"))
  suppressMessages(library("pak"))
  suppressMessages(library("remotes"))
  suppressMessages(library("tidyr"))
  suppressMessages(library("dplyr"))
  suppressMessages(library("matrixStats"))
  suppressMessages(library("reshape2"))
  suppressMessages(library("purrr"))
  suppressMessages(library("tidygraph"))
  suppressMessages(library("stringr"))
  suppressMessages(library("tibble"))
  suppressMessages(library("gplots"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("AnnotationDbi"))
  suppressMessages(library("RColorBrewer"))
  suppressMessages(library("pheatmap"))
  suppressMessages(library("ggfortify"))
  suppressMessages(library("Hmisc"))
  suppressMessages(library("ggpubr"))
  suppressMessages(library("ggstatsplot"))
  suppressMessages(library("dendextend"))
  suppressMessages(library("stats"))
  suppressMessages(library("rms"))
  suppressMessages(library("uuid"))
  suppressMessages(library("parallel"))
  suppressMessages(library("factoextra"))
  suppressMessages(library("doParallel"))
  suppressMessages(library("foreach"))
  suppressMessages(library("omnideconv"))
}

libraries_set()

dir.create(file.path(getwd(), "Results"))

#' Compute deconvolution processing
#'
#' @param deconvolution A matrix with unprocessed cell deconvolution results
#' @param corr A numeric value with the minimum correlation allowed to group cell deconvolution features
#' @param zero A numeric value with the maximum proportion of zeros values allowed the deconvolution features to have. Features with higher number of zeros across samples (>zero) will be discarded. Default is 0.9
#' @param high_corr A numeric value with the threshold for pairwise correlation. Pair of features correlated with more than this threshold are classify as 'high correlated' and choose randomly one of them. Default is 0.9
#' @param seed A numeric value to specificy the seed. This ensures reproducibility during the choice step of high correlated features.
#' @param return Boolean value to whether return and saved the plot and csv files of deconvolution generated during the run inside the Results/ directory.
#'
#'
#'
#' @return A list containing
#'
#' - A matrix with the deconvolution after processing
#' - The deconvolution subgroups per cell type
#' - The deconvolution groups discarded caused they are all belonging to the same method
#' - The discarded features because they contain a high number of zeros across samples (> zero)
#' - Discarded features due to low variance across samples
#' - Discarded cell types because they are not supported in the pipeline
#' - High correlated deconvolution pairs (>high_corr)
#'
#' @export


#' @examples
#'
#' dt = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123, return = T)
#'
#'
compute.deconvolution.analysis <- function(deconvolution, corr, zero = 0.9, high_corr = 0.9, seed = NULL, cells_extra = NULL, file_name = NULL, return = T){
  deconvolution.mat = deconvolution
  
  #####Unsupervised filtering 
  
  #Remove high zero number features
  cat(paste0("Removing features with high zero number ", round(zero*100,2), "%...............................................................\n\n"))
  deconvolution.mat = deconvolution.mat[, colSums(deconvolution.mat == 0, na.rm=TRUE) < round(zero*nrow(deconvolution.mat)) , drop=FALSE]
  diff_colnames <- setdiff(colnames(deconvolution), colnames(deconvolution.mat))
  zero_features <- deconvolution[, diff_colnames]
  
  #Remove low_variance features
  variance = remove_low_variance(deconvolution.mat, plot = return)
  deconvolution.mat = variance[[1]]
  low_variance_features = variance[[2]]
  
  # #Scale deconvolution features by columns for making them comparable between cell types (0-1). 
  # cat("Scaling deconvolution features for comparison between cell types...............................................................\n\n")
  # for (i in 1:ncol(deconvolution.mat)) {
  #   deconvolution.mat[,i] = deconvolution.mat[,i]/max(deconvolution.mat[,i])
  # } 
  
  
  #####Cell types split 
  cat("Splitting deconvolution features per cell type...............................................................\n\n")
  cells_types = compute.cell.types(deconvolution.mat, cells_extra)
  cells = cells_types[[1]]
  cells_discarded = cells_types[[2]]
  
  ######Pairwise correlation filtering (Highly correlated variables >0.9) within cell types
  cat("Finding group of features with high correlation between each other...............................................................\n\n")
  features_high_corr = list()
  j = 1
  for (i in 1:length(cells)) {
    data = cells[[i]]
    if(is.null(ncol(data))==T){
      cells[[i]] = data
    }else if(ncol(data)>1){
      data = removeCorrelatedFeatures(data, high_corr, names(cells)[i], seed)
      cells[[i]] = data[[1]]
      if(length(data[[2]])>0 && is.null(data[[3]])==F){
        features_high_corr[[j]] = data[[2]]
        names(features_high_corr)[j] = data[[3]]
        j = j+1
      }
    }
  }
  
  #####Subgrouping of deconvolution features
  res = list()
  groups = list()
  #groups_similarity = list()
  groups_discard = list()
  for (i in 1:length(cells)) {
    x = compute_subgroups(cells[[i]], file_name = names(cells)[i], thres_corr = corr)
    res = c(res, x[1])
    groups = c(groups, x[2])
    #groups_similarity = c(groups_similarity, x[3])
    groups_discard = c(groups_discard, x[3])
  }

  names_cells = names(cells)
  
  names(res) = names_cells
  names(groups) = names_cells
  #names(groups_similarity) = names_cells
  names(groups_discard) = names_cells
  
  #####Preparing output
  dt = c()
  for (i in 1:length(res)) {
    dt = c(dt, res[[i]])
  }
  dt = data.frame(dt)
  rownames(dt) = rownames(deconvolution.mat)
  
  #####Create and export table with subgroups
  
  #Count number of subgroups - Linear-based
  idx = c()
  for (i in 1:length(groups)){
    if(length(groups[[i]])>0){
      for (j in 1:length(groups[[i]])){
        idx = c(idx, names(groups[[i]])[[j]])
      } 
    }
  }  
  data.groups = data.frame(matrix(nrow = length(idx), ncol = 2)) #Create table
  colnames(data.groups) = c("Cell_subgroups", "Methods-signatures")
  data.groups$Cell_subgroups = idx #Assign subgroups 
  
  #Save methods corresponding to each subgroup
  contador = 1
  for (i in 1:length(groups)){
    if(length(groups[[i]])>0){
      for (j in 1:length(groups[[i]])){
        data.groups[contador,2] = paste(groups[[i]][[j]], collapse ="\n")
        contador = contador + 1
      } 
    }
  }
  
  
  #Count number of subgroups - Proportionality-based
  # idy = c()
  # for (i in 1:length(groups_similarity)){
  #   if(length(groups_similarity[[i]])>0){
  #     for (j in 1:length(groups_similarity[[i]])){
  #       idy = c(idy, names(groups_similarity[[i]])[[j]])
  #     } 
  #   }
  # }
  # data.groups.similarity = data.frame(matrix(nrow = length(idy), ncol = 2)) #Create table
  # colnames(data.groups.similarity) = c("Cell_subgroups", "Methods-signatures")
  # data.groups.similarity$Cell_subgroups = idy #Assign subgroups 
  # 
  # #Save methods corresponding to each subgroup
  # contador = 1
  # for (i in 1:length(groups_similarity)){
  #   if(length(groups_similarity[[i]])>0){
  #     for (j in 1:length(groups_similarity[[i]])){
  #       data.groups.similarity[contador,2] = paste(groups_similarity[[i]][[j]], collapse ="\n")
  #       contador = contador + 1
  #     } 
  #   }
  # }
  # 
  #Save data to export
  if(return){
    data.output = data.groups
    write.csv(dt, 'Results/Deconvolution_after_subgrouping.csv')
    write.csv(data.output, 'Results/Cell_subgroups.csv', row.names = F) 
  }
  
  message("Deconvolution features subgroupped")
  
  dend_column = as.dendrogram(hclust(dist(dt), method = "ward.D2"))
  
  ht1 = Heatmap(t(scale(dt)), border = T, cluster_columns = dend_column, 
                column_gap = unit(8, "mm"), name = "Deconvolution scores",
                clustering_method_rows = "ward.D2", 
                column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm"), 
                column_dend_reorder = T, row_dend_reorder = F,
                show_row_names = T, 
                show_heatmap_legend = T, 
                row_names_gp = gpar(fontsize = 10), 
                column_names_gp = gpar(fontsize =10), 
                width = unit(40, "cm"), height = unit(40, "cm"),
                heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), legend_width = unit(12, "cm"), 
                                            legend_heigh = unit(12, "cm"), title_gp = gpar(fontsize = 12)))
  
  pdf(paste0("Results/Heatmap_deconvolution_after_groupping_", file_name), height = 20, width = 25)
  draw(ht1, show_heatmap_legend = T, heatmap_legend_side = "left", annotation_legend_side = 'left')
  dev.off()
  
  results = list(dt, res, groups, groups_discard, zero_features, low_variance_features, cells_discarded, features_high_corr)
  names(results) = c("Deconvolution matrix", "Deconvolution groups scores per cell types", "Deconvolution groups - Linear-based correlation",
                     "Discarded groups with equal method", "Discarded features with high number of zeros", "Discarded features with low variance", "Discarded cell types",
                     "High correlated deconvolution groups (>0.9) per cell type")
  return(results)
  
}

compute.deconvolution.preprocessing = function(deconv){
  cat("Preprocessing deconvolution features...............................................................\n\n")

  #Remove NA (this need to be check -- not possible to have NAs values in deconv)
  deconv <- deconv %>%
    mutate(across(everything(), ~ replace_na(.x, 0)))
  
  #Convert mcp and xcell features to proportions by row-scaling 
  # for (i in 1:nrow(deconv)) {
  #   idx = grep("MCP", colnames(deconv))
  #   if(length(idx)>0){deconv[,idx][i,] = deconv[,grep("MCP", colnames(deconv))][i,]/sum(deconv[,grep("MCP", colnames(deconv))][i,])}
  #   idx = grep("XCell", colnames(deconv))
  #   if(length(idx)>0){deconv[,idx][i,] = deconv[,grep("XCell", colnames(deconv))][i,]/sum(deconv[,grep("XCell", colnames(deconv))][i,])}
  # } 
  # 
  ##### Edit cell names for consistency across features
  
  ##### Macrophages (M0, M1, M2)
  Macrophages = deconv[, grep("acrophage", colnames(deconv)), drop = F]
  M0 = deconv[,grep("M0", colnames(deconv)), drop = F]
  M1 = deconv[,grep("M1", colnames(deconv)), drop = F]
  M2 <- deconv[,grep("M2", colnames(deconv)), drop = F]
  if(length(grep("LM22", colnames(M2)))>0){M2 <- M2[,-grep("LM22", colnames(M2)), drop = F]}else{M2 <- M2} 
  test = deconv[,grep("LM22", colnames(deconv)), drop = F]
  test = test[,grep("Macrophages.M2", colnames(test)), drop = F]
  M2 = cbind(M2, test)
  
  idx = which(colnames(Macrophages)%in%c(colnames(M0), colnames(M1), colnames(M2)))
  if(length(idx)>0){
    Macrophages = Macrophages[,-idx, drop = F]
  }
  
  if(ncol(Macrophages)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Macrophages)), drop = F]
    colnames(Macrophages) = stringr::str_replace(colnames(Macrophages), "Macrophages", "Macrophages.cells") 
    colnames(Macrophages) = stringr::str_replace(colnames(Macrophages), "Macrophage(?!.)", "Macrophages.cells")  
  }

  if(ncol(M0)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(M0)), drop = F]
    colnames(M0) = stringr::str_replace(colnames(M0), "Macrophages_M0", "Macrophages.M0") 
    colnames(M0) = stringr::str_replace(colnames(M0), "_M0", "_Macrophages.M0")
  }
  if(ncol(M1)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(M1)), drop = F]
    colnames(M1) = stringr::str_replace(colnames(M1), "Macrophages_M1", "Macrophages.M1") 
    colnames(M1) = stringr::str_replace(colnames(M1), "Macrophages_M1", "Macrophages.M1")
    colnames(M1) = stringr::str_replace(colnames(M1), "Macrophage_M1", "Macrophages.M1")
    colnames(M1) = stringr::str_replace(colnames(M1), "_M1", "_Macrophages.M1")
  }
  if(ncol(M2)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(M2)), drop = F]
    colnames(M2) = stringr::str_replace(colnames(M2), "Macrophage_M2", "Macrophages.M2") 
    colnames(M2) = stringr::str_replace(colnames(M2), "Macrophages_M2", "Macrophages.M2")
    colnames(M2) = stringr::str_replace(colnames(M2), "_M2", "_Macrophages.M2")
  }
  
  ##### Monocytes
  Monocytes = deconv[,grep("Mono|mono", colnames(deconv)), drop = F]
  if(ncol(Monocytes)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Monocytes)), drop = F]
    colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Monocytic_lineage", "Monocytes") 
    colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Monocyte(?!s)", "Monocytes") 
    colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Mono(?!cytes)", "Monocytes")
    colnames(Monocytes) = stringr::str_replace(colnames(Monocytes), "Mono(?!cytes)", "Monocytes")
  }

  ##### Neutrophils
  Neutrophils <- deconv[,grep("Neu", colnames(deconv)), drop = F]
  if(ncol(Neutrophils)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Neutrophils)), drop = F]
    colnames(Neutrophils) = stringr::str_replace(colnames(Neutrophils), "Neutrophil(?!s)", "Neutrophils") 
    colnames(Neutrophils) = stringr::str_replace(colnames(Neutrophils), "Neu(?!trophils)", "Neutrophils") 
  }

  ### NK cells
  NK = deconv[,grep("NK", colnames(deconv)), drop = F]
  NKT = NK[,grep("NKT", colnames(NK)), drop = F]
  NK.activated <- NK[,grep("activated", colnames(NK), value = TRUE), drop = F]
  NK.resting <- NK[,grep("resting", colnames(NK), value = TRUE), drop = F]
  
  idx = which(colnames(NK)%in%c(colnames(NK.activated), colnames(NK.resting), colnames(NKT)))
  if(length(idx)>0){
    NK = NK[,-idx, drop = F] 
  }
   
  if(ncol(NK)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(NK)), drop = F]
    colnames(NK) = stringr::str_replace(colnames(NK), "NK(?!.)", "NK.cells")
    colnames(NK) = stringr::str_replace(colnames(NK), "NK_cells", "NK.cells")
    colnames(NK) = stringr::str_replace(colnames(NK), "NK_cell", "NK.cells")
  }

  if(ncol(NKT)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(NKT)), drop = F]
    colnames(NKT) = stringr::str_replace(colnames(NKT), "NKT_", "NKT.")
  }
  
  if(ncol(NK.activated)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(NK.activated)), drop = F]
    colnames(NK.activated) = stringr::str_replace(colnames(NK.activated), "NK.cells.activated", "NK.activated")
    colnames(NK.activated) = stringr::str_replace(colnames(NK.activated), "NK.cells_activated", "NK.activated")
  }

  if(ncol(NK.resting)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(NK.resting)), drop = F]
    colnames(NK.resting) = stringr::str_replace(colnames(NK.resting), "NK.cells.resting", "NK.resting")
    colnames(NK.resting) = stringr::str_replace(colnames(NK.resting), "NK.cells_resting", "NK.resting")
  }

  
  ### CD4 cells
  CD4 <- deconv[,grep("CD4", colnames(deconv)), drop = F]
  CD4.memory.activated = CD4[,grep("activated", colnames(CD4)), drop = F]
  CD4.memory.resting = CD4[,grep("resting", colnames(CD4)), drop = F]
  CD4.naive = CD4[,grep("naive", colnames(CD4)), drop = F]
  CD4.non.regulatory = CD4[,grep("regulatory", colnames(CD4)), drop = F]
  
  idx = which(colnames(CD4)%in%c(colnames(CD4.memory.activated), colnames(CD4.memory.resting), colnames(CD4.naive), colnames(CD4.non.regulatory)))
  if(length(idx)>0){CD4 = CD4[,-idx, drop = F]}
   
  if(ncol(CD4)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4)), drop = F]
    colnames(CD4) = stringr::str_replace(colnames(CD4), "CD4", "CD4.cells")
    colnames(CD4) = stringr::str_replace(colnames(CD4), "T.cells.CD4.cells", "CD4.cells")
  }

  if(ncol(CD4.memory.activated)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4.memory.activated)), drop = F]
    colnames(CD4.memory.activated) = stringr::str_replace(colnames(CD4.memory.activated), "CD4_memory_activated", "CD4.memory.activated")
  }
  
  if(ncol(CD4.memory.resting)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4.memory.resting)), drop = F]
    colnames(CD4.memory.resting) = stringr::str_replace(colnames(CD4.memory.resting), "CD4_memory_resting", "CD4.memory.resting")
  }
  
  if(ncol(CD4.naive)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4.naive)), drop = F]
    colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "CD4_naive", "CD4.naive")
    colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "CD4._naive", "CD4.naive")
    colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "T.cells.CD4.naive", "CD4.naive")
    colnames(CD4.naive) = stringr::str_replace(colnames(CD4.naive), "T_cells_CD4.naive", "CD4.naive")
  }

  if(ncol(CD4.non.regulatory)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4.non.regulatory)), drop = F]
    colnames(CD4.non.regulatory) = stringr::str_replace(colnames(CD4.non.regulatory), "T_cell_CD4._.non.regulatory.", "T.cells.non.regulatory")
  }

  #### CD8
  CD8 <- deconv[,grep("CD8", colnames(deconv)), drop = F]
  if(ncol(CD8)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD8)), drop = F]
    colnames(CD8) = stringr::str_replace(colnames(CD8), "T_cells_CD8", "CD8.cells")
    colnames(CD8) = stringr::str_replace(colnames(CD8), "T_cell_CD8", "CD8.cells")
    colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8_T_cells", "CD8.cells")
    colnames(CD8) = stringr::str_replace(colnames(CD8), "T.cells.CD8", "CD8.cells")
    colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8(?!.)", "CD8.cells")
    colnames(CD8) = stringr::str_replace(colnames(CD8), "CD8.cells.", "CD8.cells")
  }

  ##### Regulatory T cells 
  Tregs = deconv[,grep("regs", colnames(deconv)), drop = F]
  if(ncol(Tregs)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Tregs)), drop = F]
    colnames(Tregs) = stringr::str_replace(colnames(Tregs), "T_cell_regulatory_.Tregs.", "T.cells.regulatory")
    colnames(Tregs) = stringr::str_replace(colnames(Tregs), "T.cells.regulatory..Tregs.", "T.cells.regulatory")
    colnames(Tregs) = stringr::str_replace(colnames(Tregs), "Tregs", "T.cells.regulatory")
  }

  ##### Helper T cells 
  Thelper = deconv[,grep("helper", colnames(deconv)), drop = F]
  if(ncol(Thelper)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Thelper)), drop = F]
    colnames(Thelper) = stringr::str_replace(colnames(Thelper), "T.cells.follicular.helper", "T.cells.helper")
    colnames(Thelper) = stringr::str_replace(colnames(Thelper), "T_cells_follicular_helper", "T.cells.helper")
  }

  ##### Gamma delta T cells 
  Tgamma = deconv[,grep("gamma", colnames(deconv)), drop = F]
  if(ncol(Tgamma)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Tgamma)), drop = F]
    colnames(Tgamma) = stringr::str_replace(colnames(Tgamma), "T_cells_gamma_delta", "T.cells.gamma.delta")
    colnames(Tgamma) = stringr::str_replace(colnames(Tgamma), "T_cell_gamma_delta", "T.cells.gamma.delta")
  }

  ##### Dendritic cells (activated, resting)
  Dendritic = deconv[,grep("endritic", colnames(deconv)), drop = F]
  Dendritic.activated = Dendritic[,grep("activated", colnames(Dendritic)), drop = F]
  Dendritic.resting = Dendritic[,grep("resting", colnames(Dendritic)), drop = F]
  
  idx = which(colnames(Dendritic)%in%c(colnames(Dendritic.activated), colnames(Dendritic.resting)))
  if(length(idx)>0){Dendritic = Dendritic[,-idx, drop = F]}

  if(ncol(Dendritic)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Dendritic)), drop = F]
    colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Myeloid_dendritic_cells", "Dendritic.cells")
    colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Myeloid_dendritic_cell", "Dendritic.cells")
    colnames(Dendritic) = stringr::str_replace(colnames(Dendritic), "Dendritic_cells", "Dendritic.cells")
  }

  if(ncol(Dendritic.activated)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Dendritic.activated)), drop = F]
    colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "dendritic_cell_activated", "Dendritic.activated.cells")
    colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "Dendritic.cells.activated", "Dendritic.activated.cells")
    colnames(Dendritic.activated) = stringr::str_replace(colnames(Dendritic.activated), "Dendritic_cells_activated", "Dendritic.activated.cells")
  }
  
  if(ncol(Dendritic.resting)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Dendritic.resting)), drop = F]
    colnames(Dendritic.resting) = stringr::str_replace(colnames(Dendritic.resting), "Dendritic.cells.resting", "Dendritic.resting.cells")
    colnames(Dendritic.resting) = stringr::str_replace(colnames(Dendritic.resting), "Dendritic_cells_resting", "Dendritic.resting.cells")
  }

  ##### CAF cells 
  CAF = deconv[,grep("CAF|Cancer_associated_fibroblast", colnames(deconv)), drop = F]

  if(ncol(CAF)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CAF)), drop = F]
    colnames(CAF) = stringr::str_replace(colnames(CAF), "Cancer_associated_fibroblast", "CAF")
    colnames(CAF) = stringr::str_replace(colnames(CAF), "CAFs", "CAF")
  }

  ##### Cancer cells
  Cancer = deconv[,grep("ancer", colnames(deconv)), drop = F]

  if(ncol(Cancer)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Cancer)), drop = F]
    colnames(Cancer) = stringr::str_replace(colnames(Cancer), "cancer", "Cancer")
    colnames(Cancer) = stringr::str_replace(colnames(Cancer), "Cancer.cells", "Cancer")
  }
  
  malignant = deconv[,grep("alignant", colnames(deconv)), drop = F]

  if(ncol(malignant)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(malignant)), drop = F]
    colnames(malignant) = stringr::str_replace(colnames(malignant), "Malignant", "Cancer")
    colnames(malignant) = stringr::str_replace(colnames(malignant), "Cancer_cells", "Cancer")
    colnames(malignant) = stringr::str_replace(colnames(malignant), "Cancer.cells", "Cancer")
    if(ncol(Cancer)>0){
      Cancer = cbind(Cancer, malignant)
    }else{
      Cancer = malignant
    }
  }

  ##### Endothelial cells
  Endothelial = deconv[,grep("dothelial", colnames(deconv)), drop = F]

  if(ncol(Endothelial)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Endothelial)), drop = F]
    colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial_cells", "Endothelial")
    colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial.cells", "Endothelial")
    colnames(Endothelial) = stringr::str_replace(colnames(Endothelial), "Endothelial_cell", "Endothelial")
  }

  ##### Eosinophils cells
  Eosinophils = deconv[,grep("osinophil", colnames(deconv)), drop = F]
  if(ncol(Eosinophils)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Eosinophils)), drop = F]
    colnames(Eosinophils) = stringr::str_replace(colnames(Eosinophils), "Eosinophil(?!.)", "Eosinophils")
  }
  
  ##### Plasma cells
  Plasma = deconv[,grep("lasma", colnames(deconv)), drop = F]

  if(ncol(Plasma)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Plasma)), drop = F]
    colnames(Plasma) = stringr::str_replace(colnames(Plasma), "plasma(?!.)", "Plasma.cells")
    colnames(Plasma) = stringr::str_replace(colnames(Plasma), "Plasma_cells", "Plasma.cells")
  }

  ##### Myocytes cells
  Myocytes = deconv[,grep("yocytes", colnames(deconv)), drop = F]
  if(ncol(Myocytes)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Myocytes)), drop = F]
  }
  
  ##### Fibroblasts cells 
  Fibroblasts = deconv[,grep("ibroblast", colnames(deconv)), drop = F]
  if(ncol(Fibroblasts)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Fibroblasts)), drop = F]
  }
  
  ##### Mast cells
  Mast = deconv[,grep("Mast", colnames(deconv)), drop = F]
  Mast.activated = Mast[,grep("activated", colnames(Mast)), drop = F]
  Mast.resting = Mast[,grep("resting", colnames(Mast)), drop = F]
  
  idx = which(colnames(Mast)%in%c(colnames(Mast.activated), colnames(Mast.resting)))
  if(length(idx)>0){Mast = Mast[,-idx, drop = F]}  

  if(ncol(Mast)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Mast)), drop = F]
    colnames(Mast) = stringr::str_replace(colnames(Mast), "Mast_cell(?!.)", "Mast.cells")
    colnames(Mast) = stringr::str_replace(colnames(Mast), "Mast_cells", "Mast.cells")
  }

  if(ncol(Mast.activated)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Mast.activated)), drop = F]
    colnames(Mast.activated) = stringr::str_replace(colnames(Mast.activated), "Mast.cells.activated", "Mast.activated.cells")
    colnames(Mast.activated) = stringr::str_replace(colnames(Mast.activated), "Mast_cells_activated", "Mast.activated.cells")
  }

  if(ncol(Mast.resting)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(Mast.resting)), drop = F]
    colnames(Mast.resting) = stringr::str_replace(colnames(Mast.resting), "Mast.cells.resting", "Mast.resting.cells")
    colnames(Mast.resting) = stringr::str_replace(colnames(Mast.resting), "Mast_cells_resting", "Mast.resting.cells")
  }

  ##### B cells
  B.naive = deconv[,grep("naive", colnames(deconv)), drop = F] #Can be a problem in the future if more cells are add and have the pattern of "naive"
  if(ncol(B.naive)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(B.naive)), drop = F]
    colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B.cells.naive", "B.naive.cells")
    colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B_cells_naive", "B.naive.cells")
    colnames(B.naive) = stringr::str_replace(colnames(B.naive), "B_cell_naive", "B.naive.cells")
  }
  
  B.memory = deconv[,grep("memory", colnames(deconv)), drop = F] #Can be a problem in the future if more cells are add and have the pattern of "memory"
  if(ncol(B.memory)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(B.memory)), drop = F]
    colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B.cells.memory", "B.memory.cells")
    colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B_cells_memory", "B.memory.cells")
    colnames(B.memory) = stringr::str_replace(colnames(B.memory), "B_cell_memory", "B.memory.cells")
  }
  
  idx = which(colnames(deconv)%in%c(colnames(B.naive), colnames(B.memory)))
  if(length(idx)>0){
    deconv = deconv[,-idx, drop = F]
  } #"deconv" will include also the cells haven't been re-named, as it is the last cell + B cells
  
  if(ncol(deconv)>0){
    colnames(deconv) = stringr::str_replace(colnames(deconv), "B_cells", "B.cells")
    colnames(deconv) = stringr::str_replace(colnames(deconv), "B_cell", "B.cells")
    colnames(deconv) = stringr::str_replace(colnames(deconv), "B_lineage", "B.cells")
    colnames(deconv) = stringr::str_replace(colnames(deconv), "_B(?!.)", "_B.cells")
    B = deconv[,grep("B.cells", colnames(deconv)), drop = F]
    if(ncol(B)>0){
      deconv = deconv[,-which(colnames(deconv)%in%colnames(B)), drop = F]
    }
  }

  ## All features no renamed are 'extra' cell types because they do not belong to current nomenclature
  extra = deconv

  cell_types = cbind(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated,
                     CD4.memory.resting, CD4.naive, CD4.non.regulatory, CD8, Tregs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, 
                     Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated, Mast.resting, CAF, extra)
  
  cat("Checking consistency in deconvolution cell fractions across patients...............................................................\n\n")
  
  combinations = c("Quantiseq", "Epidish_BPRNACan_",  "Epidish_BPRNACanProMet", "Epidish_BPRNACan3DProMet", "Epidish_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                   "Epidish_CBSX.NSCLC.PBMCs.scRNAseq", "Epidish_CCLE.TIL10", "Epidish_TIL10", "Epidish_LM22", "DeconRNASeq_BPRNACan_", "DeconRNASeq_BPRNACanProMet",
                   "DeconRNASeq_BPRNACan3DProMet", "DeconRNASeq_CBSX.HNSCC.scRNAseq", "DeconRNASeq_CBSX.Melanoma.scRNAseq", "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq", "DeconRNASeq_CCLE.TIL10",
                   "DeconRNASeq_TIL10", "DeconRNASeq_LM22", "CBSX_BPRNACan_", "CBSX_BPRNACanProMet", "CBSX_BPRNACan3DProMet", "CBSX_CBSX.HNSCC.scRNAseq",
                   "CBSX_CBSX.Melanoma.scRNAseq", "CBSX_CBSX.NSCLC.PBMCs.scRNAseq", "CBSX_CCLE.TIL10", "CBSX_TIL10", "CBSX_LM22", "DWLS_BPRNACan_",  "DWLS_BPRNACanProMet", "DWLS_BPRNACan3DProMet", "DWLS_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                   "DWLS_CBSX.NSCLC.PBMCs.scRNAseq", "DWLS_CCLE.TIL10", "DWLS_TIL10", "DWLS_LM22")

  error = F
  for(i in combinations){
    idx = grep(i, colnames(cell_types))
    if(length(idx)>0){
      mat = cell_types[,idx] #A matrix of samples as rows and features only with combination[i] as columns
      sums = round(rowSums(mat), 2)
      if(all(sums == 1) == F){
        cat(paste("\n\nTotal sum across samples of combination", i, "is not 1! Remember these are proportions and the total should be 1\n"))
        cat("Samples which sum with combination", i, "is not 1:\n\n", paste0(names(sums)[sums != 1], collapse = "\n"), "\n")
        error = T
      }else{
        cat(paste("\nTotal sum across samples of combination", i, "is", round(sum(mat[1, ]), 2))) #Print only sum of 1st row 
      }
    }
  }
  
  if(error){
    warning("\nPlease verify your matrix")
  }
  
  return(cell_types)
}

#' Splitting deconvolution features by cell types
#' 
#' \code{compute.cell.types} Split deconvolution features into matrices for each cell type (features X samples)
#' 
#' @param data Deconvolution matrix from GEMDeCan output (samples X features).
#' @return list of cell types matrices (samples X features).
#' 
#' @details Cell types included: B, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.activated, CD4.resting, CD8, Tregs, Dendritic, 
#' Dendritic.activated, Dendritic.resting, Cancer, Endothelial, CAF
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
compute.cell.types = function(data, cells_extra = NULL){
  ##### B cells
  B = grep("B.cells", colnames(data))
  B = data[, B, drop = FALSE]
  ##### B naive
  B.naive = grep("B.naive", colnames(data))
  B.naive = data[, B.naive, drop = FALSE]
  ##### B memory
  B.memory = grep("B.memory", colnames(data))
  B.memory = data[, B.memory, drop = FALSE]  
  ##### Macrophages (M0, M1, M2)
  Macrophages = grep("Macrophages.cells", colnames(data))
  Macrophages = data[, Macrophages, drop = FALSE]  
  M0 = grep("Macrophages.M0", colnames(data))
  M0 = data[, M0, drop = FALSE]  
  M1 = grep("Macrophages.M1", colnames(data))
  M1 = data[, M1, drop = FALSE]  
  M2 = grep("Macrophages.M2", colnames(data))
  M2 = data[, M2, drop = FALSE]  
  ##### Monocytes
  Monocytes = grep("Monocytes", colnames(data))
  Monocytes = data[, Monocytes, drop = FALSE]   
  ##### Neutrophils
  Neutrophils = grep("Neutrophils", colnames(data))
  Neutrophils = data[, Neutrophils, drop = FALSE]     
  ##### NK cells (activated, resting)
  NK = grep("NK.cells", colnames(data))
  NK = data[, NK, drop = FALSE]   
  NK.activated = grep("NK.activated", colnames(data))
  NK.activated = data[, NK.activated, drop = FALSE]   
  NK.resting = grep("NK.resting", colnames(data))
  NK.resting = data[, NK.resting, drop = FALSE]   
  ##### NKT cells
  NKT = grep("NKT.cells", colnames(data))
  NKT = data[, NKT, drop = FALSE]   
  ##### CD4 cells (activated, resting)
  CD4 = grep("CD4.cells", colnames(data))
  CD4 = data[, CD4, drop = FALSE]   
  #memory = CD4[,grep("memory", colnames(CD4))]
  #helper = CD4[,grep("Th", colnames(CD4))]
  #CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(memory), colnames(helper)))]
  CD4.memory.activated = grep("CD4.memory.activated", colnames(data))
  CD4.memory.activated = data[, CD4.memory.activated, drop = FALSE]   
  CD4.memory.resting = grep("CD4.memory.resting", colnames(data))
  CD4.memory.resting = data[, CD4.memory.resting, drop = FALSE]   
  CD4.naive = grep("CD4.naive", colnames(data))
  CD4.naive = data[, CD4.naive, drop = FALSE]
  ##### CD8 cells
  CD8 = grep("CD8.cells", colnames(data))
  CD8 = data[, CD8, drop = FALSE]
  #naive = grep("naive", colnames(CD8))
  #naive = CD8[, naive, drop = FALSE]
  #memory = grep("memory", colnames(CD8))
  #memory = CD8[, memory, drop = FALSE]
  #CD8 = CD8[,-which(colnames(CD8)%in%c(colnames(memory), colnames(naive)))]
  ##### Regulatory T cells 
  Tregs = grep("T.cells.regulatory", colnames(data))
  Tregs = data[, Tregs, drop = FALSE]
  ##### Non regulatory T cells 
  T.non.regs = grep("T.cells.non.regulatory", colnames(data))
  T.non.regs = data[, T.non.regs, drop = FALSE]
  ##### Helper T cells 
  Thelper = grep("T.cells.helper", colnames(data))
  Thelper = data[, Thelper, drop = FALSE]
  ##### Gamma delta T cells 
  Tgamma = grep("T.cells.gamma.delta", colnames(data))
  Tgamma = data[, Tgamma, drop = FALSE]
  ##### Dendritic cells (activated, resting)
  Dendritic = grep("Dendritic.cells", colnames(data))
  Dendritic = data[, Dendritic, drop = FALSE]
  Dendritic.activated = grep("Dendritic.activated", colnames(data))
  Dendritic.activated = data[, Dendritic.activated, drop = FALSE]
  Dendritic.resting = grep("Dendritic.resting", colnames(data))
  Dendritic.resting = data[, Dendritic.resting, drop = FALSE]
  ##### Cancer cells
  Cancer = grep("Cancer", colnames(data))
  Cancer = data[, Cancer, drop = FALSE]
  ##### Endothelial cells
  Endothelial = grep("Endothelial", colnames(data))
  Endothelial = data[, Endothelial, drop = FALSE]
  ##### Eosinophils cells
  Eosinophils = grep("Eosinophils", colnames(data))
  Eosinophils = data[, Eosinophils, drop = FALSE]
  ##### Plasma cells
  Plasma = grep("Plasma.cells", colnames(data))
  Plasma = data[, Plasma, drop = FALSE]
  ##### Myocytes cells
  Myocytes = grep("Myocytes", colnames(data))
  Myocytes = data[, Myocytes, drop = FALSE]
  ##### Fibroblasts cells
  Fibroblasts = grep("Fibroblasts", colnames(data))
  Fibroblasts = data[, Fibroblasts, drop = FALSE]
  ##### Mast cells
  Mast = grep("Mast.cells", colnames(data))
  Mast = data[, Mast, drop = FALSE]
  Mast.activated = grep("Mast.activated", colnames(data))
  Mast.activated = data[, Mast.activated, drop = FALSE]
  Mast.resting = grep("Mast.resting", colnames(data))
  Mast.resting = data[, Mast.resting, drop = FALSE]
  ##### CAF cells 
  CAF = grep("CAF", colnames(data))
  CAF = data[, CAF, drop = FALSE]
  
  #####Output list
  cell_types = list(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated, CD4.memory.resting, CD4.naive,
                    CD8, Tregs, T.non.regs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated,
                    Mast.resting, CAF)
  
  names(cell_types) = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated",
                        "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial",
                        "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated", "Mast.resting", "CAF")
  
  cell_types_matrix = cbind(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.memory.activated, CD4.memory.resting, CD4.naive,
                      CD8, Tregs, T.non.regs, Thelper, Tgamma, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Endothelial, Eosinophils, Plasma, Myocytes, Fibroblasts, Mast, Mast.activated,
                      Mast.resting, CAF)
  
  ####Add extra cells (if exist)
  if(is.null(cells_extra) == F){
    extra = list()
    for (i in 1:length(cells_extra)){
      extra[[i]] = grep(cells_extra[i], colnames(data))
      extra[[i]] = data[, extra[[i]], drop = FALSE] 
      names(extra)[i] = cells_extra[[i]]
    }
    extra_df = do.call(cbind, extra)
    
    cell_types = c(cell_types, extra)
    cell_types_matrix = cbind(cell_types_matrix, extra_df)
  }
  
  ####Discarded cell types
  cell_types_discarded = data[,!(colnames(data)%in%colnames(cell_types_matrix)), drop = F]
  
  return(list(cell_types, cell_types_discarded))
  
}

#' Remove highly correlated features
#' 
#' \code{removeCorrelatedFeatures} Remove one of two highly correlated features above a threshold 
#' 
#' @param data Matrix (samples X features).
#' @param threshold Cutoff to define above which corr number you want to consider highly correlated (default = 0.9).
#' @return Matrix with only one feature from a pair of highly correlated variables.
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
removeCorrelatedFeatures <- function(data, threshold, name, n_seed) {
  
  features_high_corr = c()
  cell_name = c()
  # Compute correlation matrix
  corr_matrix <- cor(data)
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    set.seed(n_seed)
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      mutate_all(~ifelse(. > threshold, ., NA)) %>%
      select_if(~all(!is.na(.)))
    
    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features
    
    if(ncol(feature)>1){
      keep = colnames(feature)[sample(ncol(feature), size = 1)] #From high corr group, keep only one feature
      
      print(paste0("Highly correlated features (r>", threshold,"): ", paste(colnames(feature), collapse = ', ')))
      cat(paste0("Keeping only feature: ", keep, "\n\n")) 
      
      if(length(features_high_corr)>0){
        features_high_corr = c(features_high_corr, colnames(feature))
      }else{
        features_high_corr = colnames(feature)
      }
      
      feature = feature[,-which(colnames(feature)%in%keep), drop = F] 

      if(contador==1){
        new_data <- data[, -which(colnames(data)%in%colnames(feature)), drop = F] #Remove rest of the features from original data
      }else{
        new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature)), drop = F] #Remove rest of the features from original data
      }
      contador = contador + 1
      cell_name = name
    }else{
      if(contador == 1){ 
        new_data = data
      }else{ #If it already started the loop
        new_data = new_data
      }
    }
  }
  
  if(length(cell_name)==0){
      cell_name = NULL
  }

  return(list(new_data, features_high_corr, cell_name))
}

#' Post-filtering after subgroups computation
#' 
#' \code{remove_subgroups} Remove subgroups that have the same method across different signatures 
#' 
#' @param groups List of Groups of features within cell types.
#' @return List of position of groups which have features of same method.
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
remove_subgroups = function(groups){
  lis = c()
  for (pos in 1:length(groups)){
    x = c()
    if(length(groups[[pos]])!=0){
      for (i in 1:length(groups[[pos]])) {
        x =  c(x,str_split(groups[[pos]][[i]], "_")[[1]][[1]])
      }
      if(length(unique(x)) == 1){
        lis = c(lis, pos)
      }
    }
  }
  
  return(lis)
} 

#' Post-filtering after subgroups computation
#' 
#' \code{compute_subgroups} Remove subgroups that have the same method across different signatures 
#' 
#' @param data Cell type deconvolution matrix from GEMDeCan output after splitting (samples X features).
#' @param thres_similarity Threshold for grouping features by proportionality.
#' @param thres_corr Threshold for grouping features by Pearson correlation.
#' @param thres_change Accepted variation between two subgroups to decide if keep subgrouping or finish the iteration.
#' @param file_name Name of the file.
#' 
#' @return List of subgroups of features within cell types. (data frame with subgroups, subgroups of corr, subgroups by similarity)
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
compute_subgroups = function(deconvolution, thres_corr, file_name){
  data = data.frame(deconvolution)
  cell_subgroups = list()
  #cell_groups_similarity = list()
  cell_groups_discard = list()
  if (ncol(data) < 2) {
    warning("Deconvolution features with less than two columns for subgrouping (skipping)\n")
    return(list(data, cell_subgroups, cell_groups_discard))
  }else{
    # #################### Proportionality-based correlation   
    # is_similar <- function(value1, value2, threshold) {return(abs(value1 - value2) <= threshold)}
    # similarity_matrix <- matrix(FALSE, nrow = ncol(data), ncol = ncol(data), dimnames = list(names(data), names(data)))
    # for (col1 in names(data)) {
    #   for (col2 in names(data)) {
    #     similarity <- all(mapply(is_similar, data[[col1]], data[[col2]], MoreArgs = list(0.05))) #similarity threshold = 0.05
    #     similarity_matrix[col1, col2] <- similarity
    #   }
    # }
    # get_upper_tri <- function(cormat){
    #   cormat[lower.tri(cormat, diag = T)]<- NA
    #   return(cormat)
    # }
    # upper_tri <- get_upper_tri(similarity_matrix)
    # x <- melt(upper_tri) %>%
    #   na.omit() %>%
    #   mutate_all(as.character)
    # indice = 1
    # subgroup = list()
    # vec = unique(x$Var1)
    # while(length(vec)>0){
    #   sub = x[which(x$Var1%in%vec[1]),] 
    #   sub = sub[which(sub$value==T),]
    #   if(nrow(sub)!=0){
    #     subgroup[[indice]] = c(vec[1], sub$Var2)
    #     x = x[-which(x$Var1%in%subgroup[[indice]]),] #Variable 1
    #     x = x[-which(x$Var2%in%subgroup[[indice]]),] #Variable 2
    #     vec = vec[-which(vec%in%subgroup[[indice]])]
    #     indice = indice + 1
    #   }else{
    #     indice = indice
    #     vec = vec[-1]
    #   }
    # }
    # 
    # if(length(subgroup)!=0){
    #   for (i in 1:length(subgroup)){ 
    #     names(subgroup)[i] = paste0(file_name, "_Subgroup.Similarity.", i) #Name subgroups
    #   }
    #   lis = remove_subgroups(subgroup) #Map subgroups with same method
    #   if(length(lis)>0){
    #     cell_groups_discard = subgroup[lis]
    #     subgroup = subgroup[-lis] #Remove subgroups if all subgroupped features belong to the same method
    #   }
    #   
    #   if(length(subgroup)!=0){  #check if after removal of subgroups with equal method, you still have subgroups
    #     cell_groups_similarity = subgroup
    #     data_sub = c()
    #     for(i in 1:length(cell_groups_similarity)){ #Create data frame with features subgroupped
    #       sub = data.frame(data[,colnames(data)%in%cell_groups_similarity[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
    #       sub$median = rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroups across patients 
    #       data_sub = data.frame(cbind(data_sub, sub$median)) #Save median in a new data frame
    #       colnames(data_sub)[i] = names(cell_groups_similarity)[i]
    #       name = colnames(data)[which(!(colnames(data)%in%cell_groups_similarity[[i]]))]
    #       data = data[,-which(colnames(data)%in%cell_groups_similarity[[i]])] #Remove from deconvolution features that are subgrouped
    #       if(ncol(data.frame(data))==1){
    #         data = as.data.frame(data)
    #         colnames(data)[1] = name
    #       }
    #     }
    #     
    #     rownames(data_sub) = rownames(data) #List of patients
    #     data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_groups_similarity)])
    #     colnames(data_sub) = names(cell_groups_similarity)
    #     
    #     data = cbind(data, data_sub) #Join subgroups in deconvolution file
    #   }else{
    #     cell_groups_similarity = list()
    #   }
    #   k = 2
    # }else{
    #   k = 3
    # }
    # if(ncol(data) == 1){ #everything is already subgroupped
    #   return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))  
    # }
    
    #################### Linear-based correlation
    #if(k==2 | k==3){
      terminate = FALSE
      iteration = 1
      while (terminate == FALSE) {
        corr_df <- correlation(data.matrix(data))
        vec = colnames(data)
        indice = 1
        subgroup = list()
        data_sub = c() 
        while(length(vec)>0){ #Keep running until no features are left
          if(vec[1] %in% corr_df$measure1){ #Check if feature still no-grouped
            tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature against the others
            tab = tab[tab$r>thres_corr,] #Select features corr above the threshold 
            if(nrow(tab)!=0){ #If algorithm found features above corr
              subgroup[[indice]] = c(vec[1], tab$measure2) #Save features as subgroup
              idx = which(corr_df$measure1 %in% subgroup[[indice]])
              if(length(idx)>0){corr_df = corr_df[-idx,]} #Remove features already subgroupped 
              idy = which(corr_df$measure2 %in% subgroup[[indice]])
              if(length(idy)>0){corr_df = corr_df[-idy,]} #Remove features already subgroupped
              vec = vec[-which(vec%in%subgroup[[indice]])] #Remove feature already subgroupped from vector  
              indice = indice + 1
            }else{ #Condition when there is no correlation above the threshold (features no subgroupped)
              corr_df = corr_df[-which(corr_df$measure1 == vec[1]),] #Remove variable from corr matrix to keep subgrouping the others
              if(length(which(corr_df$measure2==vec[1]))>0){corr_df = corr_df[-which(corr_df$measure2 == vec[1]),]}
              vec = vec[-1] #Remove variable from vector to keep analyzing the others 
              indice = indice #Not increase index cause no subgroup appeared
            }
          }else{ #If feature is not in corr matrix it means that there is no any significant correlation against it and other features 
            vec = vec[-1] #Remove variable from vector to keep analyzing the others
            indice = indice  #Not increase index cause no subgroup appeared
          }
        }
        
        if(length(subgroup)!=0){
          for (i in 1:length(subgroup)){ #Name subgroups
            names(subgroup)[i] = paste0(file_name, "_Subgroup.", i, ".Iteration.", iteration)
          }
          ###Check whenever some subgroups belong to the same method
          if(iteration == 1){
            idx = remove_subgroups(subgroup) #Map subgroups with same method
            if(length(idx)>0){
              if(length(cell_groups_discard)>0){
                cell_groups_discard = c(cell_groups_discard, subgroup[idx])
                duplica = which(duplicated(cell_groups_discard)) #check if there are subgroups duplicated discarded
                if(length(duplica)>0){
                  cell_groups_discard = cell_groups_discard[-duplica]
                }
              }
              else{
                cell_groups_discard = subgroup[idx]
              }
              subgroup = subgroup[-idx] #Remove subgroups if all subgroupped features belong to the same method
            }
          }
          
          if(length(subgroup)!=0){ #check if after removal of subgroups with equal method, you still have subgroups (when iteration == 1)
            #Take median expression of subgroups
            for(i in 1:length(subgroup)){ #Create data frame with features subgroupped
              sub = data.frame(data[,colnames(data)%in%subgroup[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
              sub$median = rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroup across patients 
              data_sub = data.frame(cbind(data_sub, sub$median)) #Save median in a new data frame
              colnames(data_sub)[i] = names(subgroup)[i]
              name = colnames(data)[which(!(colnames(data)%in%subgroup[[i]]))]
              data = data.frame(data[,-which(colnames(data)%in%subgroup[[i]])]) #Remove from deconvolution features that are subgrouped
              if(ncol(data.frame(data))==1){
                data = as.data.frame(data)
                colnames(data)[1] = name
              }
            }
            
            rownames(data_sub) = rownames(data) #List of patients
            
            if(iteration == 1){ #Save what is inside the first subgroups
              cell_subgroups = subgroup 
              data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_subgroups)])
              colnames(data_sub) = names(cell_subgroups)
            }else{
              for (i in 1:length(subgroup)) {
                cell_subgroups[[length(cell_subgroups)+1]] = subgroup[[i]]
                names(cell_subgroups)[length(cell_subgroups)] = names(subgroup)[i]
              }
            }
            
            if(ncol(data)!=0){
              data = cbind(data, data_sub)
            }else{
              data = data_sub #data will be 0 if all deconvolution features have been subgroupped
              terminate = TRUE
            }
            iteration = iteration + 1
          }else{
            terminate = TRUE #when the only subgroup that keep grouping is composed from the same method
          }
          
        }else{
          terminate = TRUE
        }
      }
      
      if(is.null(data_sub)==FALSE){
        data = cbind(data, data_sub)
      }
      
      idx = which(duplicated(t(data)))
      if(length(idx)>0){
        names = colnames(data)[idx]
        data = data[,-idx, drop = F]
        colnames(data) = names 
      }
    #}
    
    return(list(data, cell_subgroups, cell_groups_discard))
  }
  
}

#' Pearson correlation for features
#' 
#' \code{correlation} Perform pairwise correlations between two matrices and return organized table of only significant corr
#' 
#' @param data Matrix (samples X features).
#' @return Matrix with corr, p-value, abs_corr and significant correlations.
#' 
#' -------------------------------------------------------------------------------------------------------------
#' 
correlation <- function(data) {
  
  M <- Hmisc::rcorr(as.matrix(data), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))
  
  corr_df = Mdf %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
  
  corr_df = na.omit(corr_df)  #remove the ones that are the same TFs (pval = NA)
  corr_df <- corr_df[which(corr_df$sig_p==T),]  #remove not significant
  corr_df <- corr_df[order(corr_df$r, decreasing = T),]
  corr_df$AbsR =  abs(corr_df$r)
  
  return(corr_df)
  
}

remove_low_variance <- function(data, plot = T) {
  vars <- apply(data, 2, var)
  threshold = summary(vars)[[2]]
  cat("Checking distribution of variances...............................................................\n\n")
  cat("Chosen threshold is:", threshold, "\n\n")
  cat("Saving variance distribution plot in Results/ folder\n\n")
  low_variance <- which(vars < threshold)
  cat("Removing", length(low_variance), "features with variance across samples below this threshold...............................................................\n\n")
  
  if(plot){
    pdf("Results/Distribution_variances_deconvolution.pdf")
    hist(vars, main = "Distribution of deconvolution variances across samples\nRemoving features below threshold (low variance)", xlab = "Variance", col = "skyblue", border = "white", xlim = range(vars))
    lines(density(vars), col = "red", lwd = 2)
    legend("topright", legend = c("Density", paste("Threshold =", round(threshold, 5))), col = c("red", "orange"), lty = c(1, 2), lwd = c(2, 2))
    
    # Shade region below threshold
    abline(v = threshold, col = "orange", lwd = 2, lty = 2)
    x <- density(vars)$x
    y <- density(vars)$y
    polygon(c(min(x[vars < threshold]), x[vars < threshold], max(x[vars < threshold])), 
            c(0, y[vars < threshold], 0), col = adjustcolor("orange", alpha.f = 0.3), border = NA)
    dev.off() 
  }
  
  data_filt = data[,-low_variance, drop = F]
  low_var_features = data[,low_variance, drop = F]
  
  res = list(data_filt, low_var_features)
  return(res)
}

computeQuantiseq <- function(TPM_matrix) {
  require(immunedeconv)
  TPM_matrix = TPM_matrix[rownames(TPM_matrix)%in%rownames(immunedeconv::dataset_racle$expr_mat),] #To avoid problems regarding gene names (quantiseq error)
  
  quantiseq = immunedeconv::deconvolute(TPM_matrix, "quantiseq", tumor = T) %>% 
    column_to_rownames("cell_type") %>%
    t()
  
  colnames(quantiseq) = paste0("Quantiseq_", colnames(quantiseq))
  colnames(quantiseq) <- colnames(quantiseq) %>%
    str_replace_all(., " ", "_")
  
  return(quantiseq)
}

computeMCP <- function(TPM_matrix, genes_path) {
  require(MCPcounter)
  genes <- read.table(paste0(genes_path, "/MCPcounter/MCPcounter-genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  mcp <- MCPcounter.estimate(TPM_matrix, genes = genes, featuresType = "HUGO_symbols", probesets = NULL) %>%
    t()
  
  colnames(mcp) = paste0("MCP_", colnames(mcp))
  colnames(mcp) <- colnames(mcp) %>%
    str_replace_all(., " ", "_")
  
  return(mcp)
}

computeXCell <- function(TPM_matrix) {
  require(immunedeconv)
  
  xcell = immunedeconv::deconvolute(TPM_matrix, "xcell") %>% 
    column_to_rownames("cell_type") %>%
    t()
  
  colnames(xcell) = paste0("XCell_", colnames(xcell))
  colnames(xcell) <- colnames(xcell) %>%
    str_replace_all(., " ", "_")
  
  return(xcell)
}

computeCBSX_parallel = function(TPM_matrix, signatures, name, password, workers){
  cl = parallel::makeCluster(workers)
  registerDoParallel(cl)  
  
  cbsx = foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    source("src/environment_set.R")
    signature <- read.delim(signatures[[i]], row.names=1)
    signature_name = str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeCBSX(TPM_matrix, signature, name, password, signature_name)
  }
  
  parallel::stopCluster(cl)
  unregister_dopar() 
  
  return(cbsx)
}

computeCBSX = function(TPM_matrix, signature_file, name, password, name_signature){
  set_cibersortx_credentials(name, password)
  cbsx = omnideconv::deconvolute_cibersortx(TPM_matrix, signature_file)
  
  colnames(cbsx) = paste0("CBSX_", name_signature, "_", colnames(cbsx))
  colnames(cbsx) <- colnames(cbsx) %>%
    str_replace_all(., " ", "_")
  
  return(cbsx)
}

computeDWLS_parallel = function(TPM_matrix, signatures, workers){
  cl = parallel::makeCluster(workers)
  registerDoParallel(cl)  
  
  dwls = foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    source("src/environment_set.R")
    signature <- read.delim(signatures[[i]], row.names=1)
    signature_name = str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeDWLS(TPM_matrix, signature, signature_name)
  }
  
  parallel::stopCluster(cl)
  unregister_dopar() 
  
  return(dwls)
}

computeDWLS = function(TPM_matrix, signature_file, name_signature){
  require(DWLS)
  genes = rownames(signature_file)
  
  signature_file <- signature_file %>%
    apply(., 2, as.numeric) %>%
    data.frame() %>%
    mutate("Genes" = genes) %>%
    column_to_rownames("Genes") %>%
    as.matrix()
  
  dwls = omnideconv::deconvolute_dwls(TPM_matrix, signature_file, dwls_submethod = "SVR", verbose = T)
  
  colnames(dwls) = paste0("DWLS_", name_signature, "_", colnames(dwls))
  colnames(dwls) <- colnames(dwls) %>%
    str_replace_all(., " ", "_")
  
  return(dwls)
}

computeEpiDISH = function(TPM_matrix, signature_file, name_signature){
  require(EpiDISH)
  epi <- epidish(TPM_matrix, as.matrix(signature_file), method = "RPC", maxit = 500)
  epidish = epi$estF
  
  colnames(epidish) = paste0("Epidish_", name_signature, "_", colnames(epidish))
  colnames(epidish) <- colnames(epidish) %>%
    str_replace_all(., " ", "_")
  
  return(epidish)
}

computeDeconRNASeq = function(TPM_matrix, signature_file, name_signature){
  require(DeconRNASeq)

  decon <- DeconRNASeq(TPM_matrix, data.frame(signature_file))
  deconRNAseq = decon$out.all
  rownames(deconRNAseq) = colnames(TPM_matrix)
  
  colnames(deconRNAseq) = paste0("DeconRNASeq_", name_signature, "_", colnames(deconRNAseq))
  colnames(deconRNAseq) <- colnames(deconRNAseq) %>%
    str_replace_all(., " ", "_")
  
  return(deconRNAseq)
}


compute_methods_variable_signature = function(TPM_matrix, signatures, methods = c("CBSX", "Epidish", "DeconRNASeq", "DWLS"), exclude = NULL, cbsx.name, cbsx.token, doParallel = F, workers = NULL){
  
  db=list.files(signatures, full.names = T, pattern = "\\.txt$")
  name_exclude = c()
  if(is.null(methods)==F){
    cat("\nThe following method-signature combinations are going to be calculated...............................................................\n")
    
    cat("\nMethods\n")
    for (method in methods) {
      cat("* ", method, "\n", sep = "")
    }
    cat("\nSignatures\n")
    for (i in 1:length(db)) {
      name = stringr::str_split(basename(db[[i]]), "\\.")[[1]][1]
      cat("* ", name, "\n", sep = "")
      if(is.null(exclude)==F && name %in% exclude){
        name_exclude = c(name_exclude, name)
      }
    }
    
    if(length(name_exclude)>0){
      cat("\nExcluding signatures: ", paste0(name_exclude, collapse = ", "), "\n")
    }
    
    deconvolution = list()
    
    if("CBSX" %in% methods){
      if(is.null(cbsx.name)==T || is.null(cbsx.token)==T){
        cat("\nYou select to run CBSX but no credentials were found")
        cat("\nPlease set your credentials in the function for running CibersortX")
        stop()
      }
    }
    
    for (i in 1:length(db)) {
      signature <- read.delim(db[[i]], row.names=1)
      signature_name = str_split(basename(db[[i]]), "\\.")[[1]][1]
      
      ###Check whether common genes between counts and signature have values different than 0 to avoid NAs
      # common.data <- rownames(TPM_matrix) %in% rownames(signature)
      # data.check <- TPM_matrix[common.data,]
      # zero = any(rowSums(data.check != 0) == 0)
      # if(zero){
      #   exclude = c(exclude, signature_name)
      #   warning("Common genes between count matrix and signature ", signature_name, " are all zero values")
      # }
      
      if(!is.null(exclude) && signature_name %in% exclude) {
        next
      }else{
        if("DeconRNASeq"%in%methods){
          cat("\nRunning DeconRNASeq...............................................................\n\n")
          deconrnaseq <- computeDeconRNASeq(TPM_matrix, signature, signature_name)}
        if("Epidish"%in%methods){
          cat("\nRunning Epidish...............................................................\n\n")
          epidish <- computeEpiDISH(TPM_matrix, signature, signature_name)}
        if("DWLS"%in%methods){
          if(doParallel == F){
            cat("\nRunning DWLS...............................................................\n\n")
            dwls <- computeDWLS(TPM_matrix, signature, signature_name)}}
        if("CBSX"%in%methods){
          if(doParallel == F){
            cat("\nRunning CBSX...............................................................\n\n")
            cbsx <- computeCBSX(TPM_matrix, signature, cbsx.name, cbsx.token, signature_name)}}
        combined_data <- NULL
        if (exists("deconrnaseq")) {
          combined_data <- deconrnaseq
        }
        if (exists("epidish")) {
          if (is.null(combined_data)) {
            combined_data <- epidish
          } else {
            combined_data <- cbind(combined_data, epidish)
          }
        }
        if (exists("cbsx")) {
          if (is.null(combined_data)) {
            combined_data <- cbsx
          } else {
            combined_data <- cbind(combined_data, cbsx)
          }
        }
        if (exists("dwls")) {
          if (is.null(combined_data)) {
            combined_data <- dwls
          } else {
            combined_data <- cbind(combined_data, dwls)
          }
        }
        
        deconvolution[[i]] <- combined_data
      }
    }
  
    deconv = do.call(cbind, deconvolution)
    
    if("DWLS"%in%methods & doParallel == T){
      cat("\nRunning DWLS in parallel using", workers,"workers...............................................................\n\n")
      dwls <- computeDWLS_parallel(TPM_matrix, db, workers)
      deconv = cbind(deconv, dwls)
    }
    
    if("CBSX"%in%methods & doParallel == T){
      cat("\nRunning CBSX in parallel using", workers,"workers...............................................................\n\n")
      cbsx <- computeCBSX_parallel(TPM_matrix, db, cbsx.name, cbsx.token, workers)
      deconv = cbind(deconv, cbsx)
    }
    
    return(deconv)
  }else{
    cat("\nNo methods to be calculated using variable signatures.")
    return(NULL)
  }
  
}

#' Compute deconvolution
#'
#'The function calculates cell abundance based on cell type signatures using different methods and signatures. Methods available are Quantiseq, MCP, XCell, CibersortX, EpiDISH, DWLS and DeconRNASeq. Provided signatures included signatures based on bulk and methylation data (7 methods and 10 signature in total). Signatures are present in the src/signatures directory, user can add its own signatures by adding the .txt files in this same folder. Second generation methods to perform deconvolution based on single cell data are also available if scRNAseq object is provided.
#'
#' @param raw.counts A matrix with the raw counts (samples as columns and genes symbols as rows)
#' @param methods A character vector with the deconvolution methods to run. Default are "Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"
#' @param signatures_exclude A character vector with the signatures to exclude from the src/signatures folder.
#' @param normalized If raw.counts are not available, user can input its normalized counts. In that case this argument need to be set to False.
#' @param doParallel Whether to do or not parallelization. Only CBSX and DWLS methods will run in parallel.
#' @param workers Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1
#' @param return Whether to save or not the csv file with the deconvolution features
#' @param credentials.mail (Optional) Credential email for running CibersortX. If not provided, cibersortX method will not be run.
#' @param credentials.token (Optional) Credential token for running CibersortX. If not provided, cibersortX method will not be run.
#' @param sc_deconv Whether to run or not deconvolution methods based on single cell.
#' @param ncores If sc_deconv = T, number of cores to use for running the second-generation methods.
#' @param sc_matrix If sc_deconv = T, the matrix of counts across cells from the scRNAseq object is provided.
#' @param cell_annotations If sc_deconv = T, a character vector indicating the cell labels (same order as the count matrix)
#' @param cell_samples If sc_deconv = T, a character vector indicating the cell samples IDs (same order as the count matrix)
#' @param name_sc_signature If sc_deconv = T, the name you want to give to the signature generated
#' @param file_name File name for the csv files and plots saved in the Results/ directory
#'
#' @return
#'
#' A matrix of cell type deconvolution features across samples
#'
#' @export
#'
#' @examples
#'
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", file_name = "Tutorial")
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", methods = c("Quantiseq", "MCP", "XCell", "DWLS"), file_name = "Test")
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", signatures_exclude = "BPRNACan", file_name = "Tutorial")
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", sc_deconv = T, sc_matrix = sc.object, cell_annotations = cell_labels, cell_samples = bath_ids, name_sc_signature = "Signature_test", file_name = "Test")
#'
#'

compute.deconvolution <- function(raw.counts, methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), signatures_exclude = NULL, normalized = T, doParallel = F, workers = NULL, return = T,
                                  credentials.mail = NULL, credentials.token = NULL, sc_deconv = F, ncores = NULL, sc_matrix = NULL, cell_annotations = NULL, cell_samples = NULL, name_sc_signature = NULL, file_name = NULL){

  path_signatures = 'signatures'
  
  if(normalized == T){
    cat("Performing TPM normalization ................................................................................\n\n")
    TPM_matrix = data.frame(ADImpute::NormalizeTPM(raw.counts)) 
  }else{ #If no raw counts are available
    TPM_matrix = data.frame(raw.counts)
  }

  cat("Running deconvolution using the following methods...............................................................\n\n")
  for (method in methods) {
    cat("* ", method, "\n", sep = "")
  }
  
  if("Quantiseq" %in% methods){
    cat("\nRunning Quantiseq...............................................................\n")
    quantiseq = computeQuantiseq(TPM_matrix)}
  # if("MCP" %in% methods){
  #   cat("\nRunning MCPCounter...............................................................\n")
  #   mcp = computeMCP(TPM_matrix, path_signatures)}
  # if("xCell" %in% methods){
  #   xcell = computeXCell(TPM_matrix)
  #   cat("\nRunning XCell...............................................................\n")}
  # 
  default_sig = "Quantiseq"
  methods = methods[!(methods %in% default_sig)]
  if(length(methods) == 0){
    methods = NULL
  }
  deconv_sig = compute_methods_variable_signature(TPM_matrix, signatures = path_signatures, method = methods, exclude = signatures_exclude, cbsx.name = credentials.mail, cbsx.token = credentials.token, doParallel, workers)
  
  deconv_default <- NULL
  if (exists("quantiseq")) {
    deconv_default <- quantiseq
  }
  if (exists("mcp")) {
    if (is.null(deconv_default)) {
      deconv_default <- mcp
    } else {
      deconv_default <- cbind(deconv_default, mcp)
    }
  }
  if (exists("xcell")) {
    if (is.null(deconv_default)) {
      deconv_default <- xcell
    } else {
      deconv_default <- cbind(deconv_default, xcell)
    }
  }
  
  if(is.null(deconv_sig)){
    all_deconvolution_table = deconv_default
  }else{
    all_deconvolution_table = cbind(deconv_default, deconv_sig)
  }
  
  if(sc_deconv){
    message("Running second generation cell-type deconvolution methods using scRNAseq\n")
    if(is.null(sc_object)==T){
      stop("No single cell object has been provided for deconvolution.")
    }else{
      deconv_sc = compute_sc_deconvolution_methods(raw.counts, sc_matrix, cell_annotations, cell_samples, name_sc_signature, normalized = normalized,
                                                   n_cores = ncores, cbsx_name = credentials.mail, cbsx_token = credentials.token)
      all_deconvolution_table = cbind(data.frame(all_deconvolution_table), deconv_sc)
    }
  }
  
  deconvolution = compute.deconvolution.preprocessing(data.frame(all_deconvolution_table))
  
  if(return){
    write.csv(deconvolution, paste0("Results/Deconvolution_", file_name, ".csv"))
  }

  return(deconvolution)
  
}

compute_sc_deconvolution_methods = function(raw_counts, sc_object, cell_annotations, samples_ids, name_object, normalized = T,
                                            n_cores = NULL, cbsx_name = NULL, cbsx_token = NULL){
  
  if(normalized){
    bulk_counts = ADImpute::NormalizeTPM(raw_counts) 
  }else{ 
    bulk_counts = raw_counts
  }
  
  
  if(is.null(n_cores)){
    n_cores = parallel::detectCores() - 1
    message("\nUsing ",n_cores, " cores available for running.........................................................\n")
  }
  
  message("\nRunning AutogeneS...............................................................\n")
  autogenes = deconvolute_autogenes(bulk_gene_expression = bulk_counts, single_cell_object = sc_object,
                                    cell_type_annotations = as.character(cell_annotations), verbose = T)$proportions
  
  message("\nRunning BayesPrism...............................................................\n")
  bayesprism = deconvolute_bayesprism(bulk_gene_expression = raw_counts, single_cell_object = sc_object,
                                      cell_type_annotations = cell_annotations, n_cores = n_cores)$theta 
  
  message("\nRunning Bisque...............................................................\n")
  bisque = deconvolute_bisque(bulk_gene_expression = as.matrix(raw_counts), single_cell_object = sc_object,
                              cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)$bulk_props 
  
  
  message("\nRunning CibersortX...............................................................\n")
  if(is.null(cbsx_name)&is.null(cbsx_token)){
    message("CibersortX credentials not given, method will not be used.\n")
    cibersortx = NULL
  }else{
    set_cibersortx_credentials(cbsx_name, cbsx_token)
    model_cbsx <- omnideconv::build_model(sc_object, as.character(cell_annotations), 
                                          batch_ids = samples_ids, method = "cibersortx") 
    cibersortx = deconvolute_cibersortx(bulk_gene_expression = as.matrix(bulk_counts), 
                                        signature = model_cbsx, verbose = T)
  }
  
  message("\nRunning CPM...............................................................\n")
  cpm = deconvolute_cpm(bulk_gene_expression = data.frame(raw_counts), single_cell_object = sc_object, no_cores = 4, #raw counts
                        cell_type_annotations = as.character(cell_annotations), verbose = T)$cellTypePredictions
  
  message("\nRunning DWLS...............................................................\n")
  model_dwls <- omnideconv::build_model_dwls(sc_object, as.character(cell_annotations), 
                                             batch_ids = samples_ids, dwls_method = "mast_optimized", ncores = n_cores) 
  deconvolution_dwls = deconvolute_dwls(bulk_gene_expression = bulk_counts, signature = model_dwls, 
                                        dwls_submethod = "DampenedWLS", verbose = T)
  
  message("\nRunning MOMF...............................................................\n")
  model_momf <- omnideconv::build_model(bulk_gene_expression = raw_counts, sc_object, 
                                        as.character(cell_annotations), batch_ids = samples_ids, method = "momf") #raw counts 
  momf = deconvolute_momf(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = sc_object,
                          signature = model_momf, method = "KL", verbose = T)
  
  message("\nRunning MuSiC...............................................................\n")
  music = deconvolute_music(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = sc_object,
                            cell_type_annotations = cell_annotations,  batch_ids = samples_ids, verbose = T)$Est.prop.weighted
  
  message("\nRunning SCDC...............................................................\n")
  scdc = deconvolute_scdc(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = sc_object, 
                          cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)$prop.est.mvw  
  
  
  results = list(AutogeneS = autogenes, BayesPrism = bayesprism, CBSX = cibersortx, CPM = cpm, 
                 DWLS = deconvolution_dwls, MOMF = momf, MuSic = music, SCDC = scdc)
  
  
  results <- lapply(names(results), function(method) {
    deconv_method <- results[[method]]
    
    colnames(deconv_method) <- paste0(method, "_", name_object, "_", colnames(deconv_method))
    colnames(deconv_method) <- str_replace_all(colnames(deconv_method), " ", "_")
    
    return(data_element)
  })
  
  results = do.call(cbind, results)
  
  return(results)
  
  ### METHODS NOT YET IMPLEMENTED 
  
  # 1. BSeq-sc: Need CIBERSORT source code
  # message("\nRunning BSeq-sc...............................................................\n") 
  # bseqsc = deconvolute_bseqsc(bulk_gene_expression = as.matrix(bulk.data), single_cell_object = sc_object,
  #                             cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)
  
  # 2. CDSeq: Crash machine
  #message("\nRunning CDSeq...............................................................\n") 
  # cdseq = deconvolute(bulk_gene_expression = raw_counts, single_cell_object = counts.matrix, no_cores = 6, method = "cdseq",
  #                     cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)
  
  # 3. SCADEN: Takes a lot of time
  # message("\nRunning Scaden...............................................................\n")
  # model_scaden <- omnideconv::build_model(bulk_gene_expression = bulk.data, counts.matrix, as.character(cell_annotations), 
  #                                         batch_ids = samples_ids, method = "scaden") 
  # scaden = deconvolute_scaden(bulk_gene_expression = bulk_counts, 
  #                             signature = model_scaden, verbose = T)
  
}

# Recursive function to get all nested subgroup elements
get_all_cells <- function(subgroup_name, cell_subgroups) {
  if (subgroup_name %in% names(cell_subgroups)) { #Check if subgroup_name is key in cell_subgroups
    # If the subgroup contains further subgroups, retrieve their base elements 
    nested_cells <- unlist(lapply(cell_subgroups[[subgroup_name]], get_all_cells, cell_subgroups = cell_subgroups))
    return(unique(nested_cells)) # Only return base elements
  } else {
    # If subgroup_name is not a key in cell_subgroups, it is considered a base composition
    return(subgroup_name)
  }
}


# Deconvolution dictionary
# deconvolution_dictionary = function(deconvolution, progeny){
#   
#   #Remove column with zero counts
#   idx = which(colSums(deconvolution) == 0)
#   if(length(idx)!=0){
#     deconvolution = deconvolution[,-idx]
#   }
#   
#   #Compute module correlation between deconvolution features and PROGENy pathways
#   x = compute.modules.relationship(deconvolution, progeny, return = T, plot = T, file_name = "Deconvolution_Pathways")
#   
#   #Create distance matrix and hierarchical clustering for the PROGENy pathways
#   d <- dist(t(x[[1]]))
#   dendrogram <- hclust(d)
#   
#   pdf("Results/Pathways_clusters")
#   plot(dendrogram)
#   dev.off()
#   
#   #Identify the two pathway clusters
#   clusters <- cutree(dendrogram, k = 2)
#   clusters <- split(names(clusters), clusters) 
#   names(clusters) = c("Cluster_1", "Cluster_2")
#   
#   #Calculate mean correlation for each cluster
#   corr_matrix = data.frame(x[[1]])
#   corr_matrix$Cluster1_Score <- rowMeans(corr_matrix[, clusters[[1]]], na.rm = TRUE)
#   corr_matrix$Cluster2_Score <- rowMeans(corr_matrix[, clusters[[2]]], na.rm = TRUE)
#   
#   #Classify features based on the higher cluster score
#   corr_matrix$Classification <- ifelse(
#     corr_matrix$Cluster1_Score > corr_matrix$Cluster2_Score,
#     "Cluster1", "Cluster2"
#   )
#   
#   #Rename deconvolution features with cluster information
#   deconv_names = paste0(rownames(corr_matrix), "_", corr_matrix$Classification)
#   colnames(deconvolution) = deconv_names
#   
#   return(list(Deconvolution = deconvolution, Clusters = clusters))
# }



utils::globalVariables(c("mcp", "xcell" ,"i", ".", "samples_ids", "multisession", ".data", "Patient", "var", "id", "P", "sig_p", "r", "y", "p", "average", "Cells", "variable", "value", "pval_value"))

#' Compute deconvolution preprocessing
#'
#' Give consistent names and patterns following the method_signature_cell structure to the deconvolution features
#'
#' @param deconv A dataframe with the unprocessed deconvolution features
#'
#' @return A matrix of the preprocessed deconvolution features with fixed and consistent names across the different methods and signatures following the nomenclature specified in multideconv (see Readme)
#' @export
#'
#' @examples
#'
#' data("deconvolution")
#'
#' deconvolution = compute.deconvolution.preprocessing(deconvolution)
#'
compute.deconvolution.preprocessing = function(deconv){
  cat("Preprocessing deconvolution features...............................................................\n\n")

  #Remove NA (this need to be check -- not possible to have NAs values in deconv)
  deconv <- deconv %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, 0)))

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
    colnames(CD4) = stringr::str_replace(colnames(CD4), "T.cells.CD4(?!\\.cells)", "CD4.cells")
    colnames(deconv) <- stringr::str_replace(colnames(deconv), "_CD4$", "_CD4.cells")
    colnames(CD4) = stringr::str_replace(colnames(CD4), "^CD4(?!\\.cells)", "CD4.cells")
  }

  if(ncol(CD4.memory.activated)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4.memory.activated)), drop = F]
    colnames(CD4.memory.activated) = stringr::str_replace(colnames(CD4.memory.activated), "CD4_memory_activated", "CD4.memory.activated")
    colnames(CD4.memory.activated) = stringr::str_replace(colnames(CD4.memory.activated), "T.cells.CD4.memory.activated", "CD4.memory.activated")
  }

  if(ncol(CD4.memory.resting)>0){
    deconv = deconv[,-which(colnames(deconv)%in%colnames(CD4.memory.resting)), drop = F]
    colnames(CD4.memory.resting) = stringr::str_replace(colnames(CD4.memory.resting), "CD4_memory_resting", "CD4.memory.resting")
    colnames(CD4.memory.resting) = stringr::str_replace(colnames(CD4.memory.resting), "T.cells.CD4.memory.resting", "CD4.memory.resting")

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
    colnames(Plasma) = stringr::str_replace(colnames(Plasma), "plasma(?!.)", "Plasma")
    colnames(Plasma) = stringr::str_replace(colnames(Plasma), "Plasma_cells", "Plasma")
    colnames(Plasma) = stringr::str_replace(colnames(Plasma), "Plasma.cells", "Plasma")
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


#' Cell types split from deconvolution
#'
#' @param data A matrix with the deconvolution results
#' @param cells_extra A string specifying the cells names to consider and that are not including in the nomenclature of multideconv (see Readme)
#'
#' @return A list containing:
#' - A sublist with each cell type features as an element recover from the different signatures
#' - Discarded cell types (this will happen if the cell types are not supported. See the READme for more information about this)
#'
#' @export
#'
#' @examples
#'
#' data("deconvolution")
#'
#' cells_types = compute.cell.types(deconvolution)
#' cells = cells_types[[1]]
#' cells_discarded = cells_types[[2]]
#' cells_types = compute.cell.types(deconvolution, cells_extra = c("mesenchymal", "basophils"))
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
                        "Eosinophils", "Plasma", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated", "Mast.resting", "CAF")

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
    extra_df = do.call(cbind, unname(extra))

    cell_types = c(cell_types, extra)
    cell_types_matrix = cbind(cell_types_matrix, extra_df)
  }

  ####Discarded cell types
  cell_types_discarded = data[,!(colnames(data)%in%colnames(cell_types_matrix)), drop = F]

  return(list(cell_types, cell_types_discarded))

}

#' Remove high correlated cell deconvolution features
#'
#' If two deconvolution features within a specific cell type are found to be highly correlated, one feature is kept randomly for further analysis.
#'
#' @param data Deconvolution matrix
#' @param threshold Threshold for defined high correlated features
#' @param name Cell type name corresponding to the given matrix in 'data'
#' @param n_seed Seed to ensure reproducibility regarding the choice of the feature.
#'
#' @return A list containing
#'
#' - Deconvolution matrix with only one deconvolution feature per high-correlated pair.
#' - Highly correlated features found
#' - Cell type name
#'
removeCorrelatedFeatures <- function(data, threshold, name, n_seed) {

  features_high_corr = c()
  cell_name = c()
  # Compute correlation matrix
  corr_matrix <- stats::cor(data)
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    set.seed(n_seed)
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      dplyr::mutate_all(~ifelse(. > threshold, ., NA)) %>%
      dplyr::select_if(~all(!is.na(.)))

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

#' Remove subgroups that have the same method across different signatures
#'
#' @param groups Cell groups of features within cell types.
#'
#' @return List of position of groups which have features of same method.
#'
remove_subgroups = function(groups){
  lis = c()
  for (pos in 1:length(groups)){
    x = c()
    if(length(groups[[pos]])!=0){
      for (i in 1:length(groups[[pos]])) {
        x =  c(x,stringr::str_split(groups[[pos]][[i]], "_")[[1]][[1]])
      }
      if(length(unique(x)) == 1){
        lis = c(lis, pos)
      }
    }
  }

  return(lis)
}

#' Compute deconvolution subgroups
#'
#' @param deconvolution A matrix with unprocessed cell deconvolution results
#' @param thres_corr A numeric value with the minimum correlation allowed to group cell deconvolution features
#' @param file_name Base name for subgroup
#'
#' @return A list containing
#'
#' - A matrix with the processed deconvolution features
#' - Cell subgroups obtained by linear correlation
#' - Cell subgroups obtained by proportionality correlation
#' - Discard cell features either because of low variance or high zero number
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
    #       sub$median = matrixStats::rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroups across patients
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
            sub$median = matrixStats::rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroup across patients
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

#' Perform pairwise correlation across all features
#'
#' @param data Matrix with features to correlate
#'
#' @return Dataframe containing all significant correlations (pvalue < 0.05)
#'
correlation <- function(data) {

  M <- Hmisc::rcorr(as.matrix(data), type = "pearson")
  Mdf <- purrr::map(M[c("r", "P", "n")], ~data.frame(.x))

  corr_df = Mdf %>%
    purrr::map(~tibble::rownames_to_column(.x, var="measure1")) %>%
    purrr::map(~tidyr::pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    dplyr::bind_rows(.id = "id") %>%
    tidyr::pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    dplyr::mutate(sig_p = ifelse(p < .05, T, F),
                  p_if_sig = ifelse(sig_p, p, NA),
                  r_if_sig = ifelse(sig_p, r, NA))

  corr_df = stats::na.omit(corr_df)  #remove the ones that are the same Tfeatures (pval = NA)
  corr_df <- corr_df[which(corr_df$sig_p==T),]  #remove not significant
  corr_df <- corr_df[order(corr_df$r, decreasing = T),]
  corr_df$AbsR =  abs(corr_df$r)

  return(corr_df)

}


#' Remove low variance deconvolution features
#'
#' @param data Deconvolution features
#' @param plot Whether to save or not the plot of variance distribution in the Results/ directory.
#'
#' @return A list containing
#'
#' - Deconvolution matrix after removal of low variance.
#' - Discarded low variance features.
#'
remove_low_variance <- function(data, plot = FALSE) {
  vars <- apply(data, 2, var)
  threshold = summary(vars)[[2]]
  cat("Checking distribution of variances...............................................................\n\n")
  cat("Chosen threshold is:", threshold, "\n\n")
  cat("Saving variance distribution plot in Results/ folder\n\n")
  low_variance <- which(vars < threshold)
  cat("Removing", length(low_variance), "features with variance across samples below this threshold...............................................................\n\n")

  if(plot){
    grDevices::pdf("Results/Distribution_variances_deconvolution.pdf")
    graphics::hist(vars, main = "Distribution of deconvolution variances across samples\nRemoving features below threshold (low variance)", xlab = "Variance", col = "skyblue", border = "white", xlim = range(vars))
    graphics::lines(stats::density(vars), col = "red", lwd = 2)
    graphics::legend("topright", legend = c("Density", paste("Threshold =", round(threshold, 5))), col = c("red", "orange"), lty = c(1, 2), lwd = c(2, 2))

    # Shade region below threshold
    graphics::abline(v = threshold, col = "orange", lwd = 2, lty = 2)
    x <- stats::density(vars)$x
    y <- stats::density(vars)$y
    graphics::polygon(c(min(x[vars < threshold]), x[vars < threshold], max(x[vars < threshold])),
            c(0, y[vars < threshold], 0), col = grDevices::adjustcolor("orange", alpha.f = 0.3), border = NA)
    grDevices::dev.off()
  }

  data_filt = data[,-low_variance, drop = F]
  low_var_features = data[,low_variance, drop = F]

  res = list(data_filt, low_var_features)
  return(res)
}

#' Compute cell type processing
#'
#' @param deconvolution Deconvolution output of compute.deconvolution() with features as columns and samples as rows
#' @param corr Minimum correlation threshold for subgroupping the deconvolution features
#' @param seed A numeric value to specificy the seed. This ensures reproducibility during the choice step of high correlated features.
#' @param cells_extra A string specifying the cells names to consider and that are not including in the nomenclature of multideconv (see Readme)
#' @param file_name A string specifying the file name of the .csv file with the deconvolution subgroups
#' @param return Boolean value to whether return and saved the plot and csv files of deconvolution generated during the run inside the Results/ directory.
#'
#' @return A list containing
#'
#' - A matrix with the deconvolution after processing
#' - The deconvolution subgroups per cell type
#' - The deconvolution subgroups composition
#' - The deconvolution groups discarded caused they are all belonging to the same method
#' - The discarded features because they contain a high number of zeros across samples (> 90%)
#' - Discarded features due to low variance across samples
#' - Discarded cell types because they are not supported in the pipeline
#' - High correlated deconvolution pairs (>high_corr)
#'
#' @export
#'
#' @examples
#'
#' data("deconvolution")
#'
#' processed_deconvolution = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123)
#'
#' processed_deconvolution = compute.deconvolution.analysis(deconvolution, cells_extra = "mesenchymal")
#'
compute.deconvolution.analysis <- function(deconvolution, corr = 0.7, seed = NULL, cells_extra = NULL, file_name = NULL, return = FALSE){
  deconvolution.mat = deconvolution

  #####Unsupervised filtering

  #Remove high zero number features
  cat(paste0("Removing features with high zero number 90%...............................................................\n\n"))
  deconvolution.mat = deconvolution.mat[, colSums(deconvolution.mat == 0, na.rm=TRUE) < round(0.9*nrow(deconvolution.mat)) , drop=FALSE]
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
      data = removeCorrelatedFeatures(data, 0.9, names(cells)[i], seed)
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
  if(return == TRUE){
    data.output = data.groups
    utils::write.csv(dt, paste0('Results/Deconvolution_after_subgrouping_', file_name,'.csv'))
    utils::write.csv(data.output, paste0('Results/Cell_subgroups_', file_name,'.csv'), row.names = F)

    dend_column = stats::as.dendrogram(stats::hclust(stats::dist(dt), method = "ward.D2"))

    ht1 = ComplexHeatmap::Heatmap(t(scale(dt)), border = T, cluster_columns = dend_column,
                                  column_gap = grid::unit(8, "mm"), name = "Deconvolution scores",
                                  clustering_method_rows = "ward.D2",
                                  column_dend_height = grid::unit(2, "cm"), row_dend_width = grid::unit(2, "cm"),
                                  column_dend_reorder = T, row_dend_reorder = F,
                                  show_row_names = T,
                                  show_heatmap_legend = T,
                                  row_names_gp = gpar(fontsize = 10),
                                  column_names_gp = gpar(fontsize =10),
                                  width = grid::unit(40, "cm"), height = grid::unit(40, "cm"),
                                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), legend_width = grid::unit(12, "cm"),
                                                              legend_heigh = grid::unit(12, "cm"), title_gp = gpar(fontsize = 12)))

    grDevices::pdf(paste0("Results/Heatmap_deconvolution_after_groupping_", file_name), height = 20, width = 25)
    ComplexHeatmap::draw(ht1, show_heatmap_legend = T, heatmap_legend_side = "left", annotation_legend_side = 'left')
    grDevices::dev.off()
  }

  message("Deconvolution features subgroupped")

  results = list(dt, res, groups, groups_discard, zero_features, low_variance_features, cells_discarded, features_high_corr)
  names(results) = c("Deconvolution matrix", "Deconvolution subgroups per cell types", "Deconvolution subgroups composition",
                     "Discarded groups with equal method", "Discarded features with high number of zeros", "Discarded features with low variance", "Discarded cell types",
                     "High correlated deconvolution groups (>0.9) per cell type")
  return(results)

}


#' Computes QuanTIseq
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#'
#' @return A matrix with cell abundance deconvolve with QuanTIseq
#'
computeQuantiseq <- function(TPM_matrix) {
  TPM_matrix = TPM_matrix[rownames(TPM_matrix)%in%rownames(immunedeconv::dataset_racle$expr_mat),] #To avoid problems regarding gene names (quantiseq error)

  quantiseq = immunedeconv::deconvolute(TPM_matrix, "quantiseq", tumor = T) %>%
    tibble::column_to_rownames("cell_type") %>%
    t()

  colnames(quantiseq) = paste0("Quantiseq_", colnames(quantiseq))
  colnames(quantiseq) <- colnames(quantiseq) %>%
    stringr::str_replace_all(., " ", "_")

  return(quantiseq)
}

#' Computes MCPcounter
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param genes_path Path containing the MCP genes
#
#'
#' @return A matrix with cell enrichment scores from MCP
#'
computeMCP <- function(TPM_matrix, genes_path) {
  genes <- utils::read.table(paste0(genes_path, "/MCPcounter/MCPcounter-genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  mcp <- MCPcounter::MCPcounter.estimate(TPM_matrix, genes = genes, featuresType = "HUGO_symbols", probesets = NULL) %>%
    t()

  colnames(mcp) = paste0("MCP_", colnames(mcp))
  colnames(mcp) <- colnames(mcp) %>%
    stringr::str_replace_all(., " ", "_")

  return(mcp)
}

#' Computes XCell
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#'
#' @return A matrix with cell enrichment scores from XCell.
#'
computeXCell <- function(TPM_matrix) {
  #t(xCell::xCellAnalysis(counts))  # try this: more cells
  xcell = immunedeconv::deconvolute(TPM_matrix, "xcell") %>%
    tibble::column_to_rownames("cell_type") %>%
    t()

  colnames(xcell) = paste0("XCell_", colnames(xcell))
  colnames(xcell) <- colnames(xcell) %>%
    stringr::str_replace_all(., " ", "_")

  return(xcell)
}

#' Compute CIBERSORTx (CBSX) in parallel across multiple signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signatures Path where signatures files are located
#' @param name Credential email for running CIBERSORTx.
#' @param password Credential token for running CIBERSORTx.
#' @param workers Number of processes available to run on parallel.
#'
#' @return A matrix with cell abundance deconvolve with CBSX
#'
computeCBSX_parallel = function(TPM_matrix, signatures, name, password, workers){
  cl = parallel::makeCluster(workers)
  doParallel::registerDoParallel(cl)

  cbsx = foreach::foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    library(multideconv)
    signature <- utils::read.delim(signatures[[i]], row.names=1)
    signature_name = stringr::str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeCBSX(TPM_matrix, signature, name, password, signature_name)
  }

  parallel::stopCluster(cl)
  unregister_dopar()

  return(cbsx)
}

#' Computes CIBERSORTx (CBSX) using one signature
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name Credential email for running CIBERSORTx.
#' @param password Credential token for running CIBERSORTx.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with CBSX
#'
computeCBSX = function(TPM_matrix, signature_file, name, password, name_signature){
  omnideconv::set_cibersortx_credentials(name, password)
  cbsx = omnideconv::deconvolute_cibersortx(TPM_matrix, signature_file)

  colnames(cbsx) = paste0("CBSX_", name_signature, "_", colnames(cbsx))
  colnames(cbsx) <- colnames(cbsx) %>%
    stringr::str_replace_all(., " ", "_")

  return(cbsx)
}

#' Compute DWLS in parallel across multiple signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signatures Path where signatures files are located
#' @param workers Number of processes available to run on parallel.
#'
#' @return A matrix with cell abundance deconvolve with DWLS
#'
computeDWLS_parallel = function(TPM_matrix, signatures, workers){
  cl = parallel::makeCluster(workers)
  doParallel::registerDoParallel(cl)

  dwls = foreach::foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    library(multideconv)
    signature <- utils::read.delim(signatures[[i]], row.names=1)
    signature_name = stringr::str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeDWLS(TPM_matrix, signature, signature_name)
  }

  parallel::stopCluster(cl)
  unregister_dopar()

  return(dwls)
}

#' Computes DWLS
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with DWLS
#'
computeDWLS = function(TPM_matrix, signature_file, name_signature){
  genes = rownames(signature_file)

  signature_file <- signature_file %>%
    apply(., 2, as.numeric) %>%
    data.frame() %>%
    dplyr::mutate("Genes" = genes) %>%
    tibble::column_to_rownames("Genes") %>%
    as.matrix()

  dwls = omnideconv::deconvolute_dwls(TPM_matrix, signature_file, dwls_submethod = "SVR", verbose = T)

  colnames(dwls) = paste0("DWLS_", name_signature, "_", colnames(dwls))
  colnames(dwls) <- colnames(dwls) %>%
    stringr::str_replace_all(., " ", "_")

  return(dwls)
}

#' Computes MOMF
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param sc_object A matrix with the counts from scRNAseq object (genes as rows and cells as columns)
#' @param signature_file The signature file to use.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with MOMF
#'
computeMOMF = function(TPM_matrix, sc_object, signature_file, name_signature){

  genes = rownames(signature_file)

  signature_file <- signature_file %>%
    apply(., 2, as.numeric) %>% #rownames are removed here
    data.frame() %>%
    dplyr::mutate("Genes" = genes) %>% #set original rownames
    tibble::column_to_rownames("Genes") %>%
    as.matrix()

  momf = omnideconv::deconvolute_momf(bulk_gene_expression = TPM_matrix, single_cell_object = as.matrix(sc_object),
                                      signature = signature_file, method = "KL", verbose = T)$cell.prop

  colnames(momf) = paste0("MOMF_", name_signature, "_", colnames(momf))
  colnames(momf) <- colnames(momf) %>%
    stringr::str_replace_all(., " ", "_")

  return(momf)
}

#' Compute MOMF in parallel across multiple signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param sc_object A matrix with the counts from scRNAseq object (genes as rows and cells as columns)
#' @param signatures Path where signatures files are located
#' @param workers Number of processes available to run on parallel.
#'
#' @return A matrix with cell abundance deconvolve with MOMF
#'
computeMOMF_parallel = function(TPM_matrix, sc_object, signatures, workers){
  cl = parallel::makeCluster(workers)
  doParallel::registerDoParallel(cl)

  momf = foreach::foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    library(multideconv)
    signature <- utils::read.delim(signatures[[i]], row.names=1)
    signature_name = stringr::str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeMOMF(TPM_matrix, sc_object, signature, signature_name)
  }

  parallel::stopCluster(cl)
  unregister_dopar()

  return(momf)
}

#' Computes EpiDISH
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with EpiDISH
#'
computeEpiDISH = function(TPM_matrix, signature_file, name_signature){
  epi <- EpiDISH::epidish(TPM_matrix, as.matrix(signature_file), method = "RPC", maxit = 500)
  epidish = epi$estF

  colnames(epidish) = paste0("Epidish_", name_signature, "_", colnames(epidish))
  colnames(epidish) <- colnames(epidish) %>%
    stringr::str_replace_all(., " ", "_")

  return(epidish)
}

#' Computes DeconRNASeq
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with DeconRNASeq
#'
#' @import pcaMethods
computeDeconRNASeq = function(TPM_matrix, signature_file, name_signature){
  library(pcaMethods) #Explicitly loading pcaMethods to access functions like prep(), which are not imported via @import. Problem with DeconRNASeq that uses not imported functions from pcaMethods

  decon <- DeconRNASeq::DeconRNASeq(TPM_matrix, data.frame(signature_file))
  deconRNAseq = decon$out.all
  rownames(deconRNAseq) = colnames(TPM_matrix)

  colnames(deconRNAseq) = paste0("DeconRNASeq_", name_signature, "_", colnames(deconRNAseq))
  colnames(deconRNAseq) <- colnames(deconRNAseq) %>%
    stringr::str_replace_all(., " ", "_")

  return(deconRNAseq)
}

#' Compute deconvolution methods with variable signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (samples as columns and genes symbols as rows)
#' @param signatures A path with a directory where signatures are located
#' @param algos A character vector with the methods to compute (Default methods are CBSX, Epidish, DeconRNASeq and DWLS)
#' @param exclude (Optional) A character vector with the signature to exclude
#' @param cbsx.name CIBERSORTx credential mail if CBSX will be run
#' @param cbsx.token CIBERSORTx credential token if CBSX will be run
#' @param doParallel Boolean value to specify if DWLS and CBSX should run in parallel (default is False)
#' @param workers Number of worker process to run during parallelization (default is NULL)
#' @param sc_obj A matrix with the counts from scRNAseq object (genes as rows and cells as columns) to run MOMF method. If NULL, MOMF is ignored.
#'
#' @return A matrix with the deconvolution features corresponding to all combinations of methods-signatures specified
#' @export
#'
#' @references
#'
#' Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.
#' Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363
#'
#' Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Dietrich, Alexander and Merotto, Lorenzo and Pelz, Konstantin and Eder, Bernhard and Zackl, Constantin and Reinisch, Katharina and
#' Edenhofer, Frank and Marini, Federico and Sturm, Gregor and List, Markus and Finotello, Francesca. (2024) https://doi.org/10.1101/2024.06.10.598226
#'
compute_methods_variable_signature = function(TPM_matrix, signatures, algos = c("CBSX", "Epidish", "DeconRNASeq", "DWLS", "MOMF"), exclude = NULL, cbsx.name, cbsx.token, doParallel = FALSE, workers = NULL, sc_obj = NULL){

  signature_dir = "Results/custom_signatures"
  default_sig = list.files(signatures, full.names = T, pattern = "\\.txt$")
  user_files = list.files(signature_dir, full.names = TRUE, pattern = "\\.txt$")

  db <- c(default_sig, user_files)

  name_exclude = c()
  if(is.null(algos)==F){
    cat("\nThe following method-signature combinations are going to be calculated...............................................................\n")

    cat("\nMethods\n")
    for (deconv_method in algos) {
      cat("* ", deconv_method, "\n", sep = "")
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

    if("CBSX" %in% algos){
      if(is.null(cbsx.name)==T || is.null(cbsx.token)==T){
        cat("\nYou select to run CBSX but no credentials were found")
        cat("\nPlease set your credentials in the function for running CIBERSORTx")
        stop()
      }
    }

    for (i in 1:length(db)) {
      signature <- utils::read.delim(db[[i]], row.names=1)
      signature_name = stringr::str_split(basename(db[[i]]), "\\.")[[1]][1]

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
        if("DeconRNASeq"%in%algos){
          cat("\nRunning DeconRNASeq...............................................................\n\n")
          deconrnaseq <- computeDeconRNASeq(TPM_matrix, signature, signature_name)}
        if("Epidish"%in%algos){
          cat("\nRunning Epidish...............................................................\n\n")
          epidish_res <- computeEpiDISH(TPM_matrix, signature, signature_name)}
        if("DWLS"%in%algos){
          if(doParallel == F){
            cat("\nRunning DWLS...............................................................\n\n")
            dwls <- computeDWLS(TPM_matrix, signature, signature_name)}}
        if("CBSX"%in%algos){
          if(doParallel == F){
            cat("\nRunning CBSX...............................................................\n\n")
            cbsx <- computeCBSX(TPM_matrix, signature, cbsx.name, cbsx.token, signature_name)}}
        if("MOMF"%in%algos & is.null(sc_obj) == F){
          if(doParallel == F){
            cat("\nRunning MOMF...............................................................\n\n")
            momf <- computeMOMF(TPM_matrix, sc_obj, signature, signature_name)}}
        combined_data <- NULL
        if (exists("deconrnaseq")) {
          combined_data <- deconrnaseq
        }
        if (exists("epidish_res")) {
          if (is.null(combined_data)) {
            combined_data <- epidish_res
          } else {
            combined_data <- cbind(combined_data, epidish_res)
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
        if (exists("momf")) {
          if (is.null(combined_data)) {
            combined_data <- momf
          } else {
            combined_data <- cbind(combined_data, momf)
          }
        }
        deconvolution[[i]] <- combined_data
      }
    }

    deconv = do.call(cbind, deconvolution)

    if("DWLS"%in%algos && doParallel == T){
      cat("\nRunning DWLS in parallel using", workers,"workers...............................................................\n\n")
      dwls <- computeDWLS_parallel(TPM_matrix, db, workers)
      deconv = cbind(deconv, dwls)
    }

    if("CBSX"%in%algos && doParallel == T){
      cat("\nRunning CBSX in parallel using", workers,"workers...............................................................\n\n")
      cbsx <- computeCBSX_parallel(TPM_matrix, db, cbsx.name, cbsx.token, workers)
      deconv = cbind(deconv, cbsx)
    }

    if("MOMF"%in%algos && doParallel == T && is.null(sc_obj) == F){
      cat("\nRunning MOMF in parallel using", workers,"workers...............................................................\n\n")
      momf <- computeMOMF_parallel(TPM_matrix, sc_obj, db, workers)
      deconv = cbind(deconv, momf)
    }

    return(deconv)
  }else{
    cat("\nNo methods to be calculated using variable signatures.")
    return(NULL)
  }

}

#' Compute deconvolution
#'
#'The function calculates cell abundance based on cell type signatures using different methods and signatures. Methods available are Quantiseq, MCP, XCell, CIBERSORTx, EpiDISH, DWLS and DeconRNASeq. Provided signatures included signatures based on bulk and methylation data (7 methods and 10 signature in total). Signatures are present in the src/signatures directory, user can add its own signatures by adding the .txt files in this same folder. Second generation methods to perform deconvolution based on single cell data are also available if scRNAseq object is provided.
#'
#' @param raw.counts A matrix with the raw counts (samples as columns and genes symbols as rows)
#' @param methods A character vector with the deconvolution methods to run. Default are "Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"
#' @param signatures_exclude A character vector with the signatures to exclude from the src/signatures folder.
#' @param normalized If raw.counts are not available, user can input its normalized counts. In that case this argument need to be set to False.
#' @param doParallel Whether to do or not parallelization. Only CBSX and DWLS methods will run in parallel.
#' @param workers Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1
#' @param return Whether to save or not the csv file with the deconvolution features
#' @param create_signature Whether to create or not the signatures using the methods MOMF, CBSX, DWLS and BSeq-SC. If TRUE, sc_matrix shuld be provide.
#' @param credentials.mail (Optional) Credential email for running CIBERSORTx. If not provided, CIBERSORTx method will not be run.
#' @param credentials.token (Optional) Credential token for running CIBERSORTx. If not provided, CIBERSORTx method will not be run.
#' @param sc_deconv Whether to run or not deconvolution methods based on single cell.
#' @param sc_matrix If sc_deconv = T, the matrix of counts across cells from the scRNAseq object is provided.
#' @param sc_metadata Dataframe with metadata from the single cell object. The matrix should include the columns cell_label and sample_label.
#' @param cell_label If sc_deconv = T, a character vector indicating the cell labels (same order as the count matrix)
#' @param sample_label If sc_deconv = T, a character vector indicating the cell samples IDs (same order as the count matrix)
#' @param cell_markers Named list with the genes markers names as Symbol per cell types to be used to create the signature using the BSeq-SC method. If NULL, the method will be ignored during the signature creation.
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
#' data("raw_counts")
#' data("cell_labels")
#' data("sample_labels")
#' data("metacells_data")
#' data("metacells_metadata")
#' data("pseudobulk")
#'
#' deconv = compute.deconvolution(raw_counts, normalized = TRUE,
#'                                methods = c("Epidish", "DeconRNASeq"), return = FALSE)
#'
#'
#' @references
#'
#' Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.
#' Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363
#'
#' Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Dietrich, Alexander and Merotto, Lorenzo and Pelz, Konstantin and Eder, Bernhard and Zackl, Constantin and Reinisch, Katharina and
#' Edenhofer, Frank and Marini, Federico and Sturm, Gregor and List, Markus and Finotello, Francesca. (2024) https://doi.org/10.1101/2024.06.10.598226
#'
compute.deconvolution <- function(raw.counts, methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS","MOMF"), signatures_exclude = NULL, normalized = TRUE, doParallel = FALSE, workers = NULL, return = TRUE, create_signature = FALSE,
                                  credentials.mail = NULL, credentials.token = NULL, sc_deconv = FALSE, sc_matrix = NULL, sc_metadata = NULL, cell_label = NULL, sample_label = NULL, cell_markers = NULL, name_sc_signature = NULL, file_name = NULL){

  path_signatures = system.file("signatures", package = "multideconv")

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
  default_sig = "Quantiseq" #This was including MCP and XCell before
  methods = methods[!(methods %in% default_sig)]
  if(length(methods) == 0){
    methods = NULL
  }

  if(create_signature == T){
    message("\nCreating static signatures...............................................................\n")
    if(is.null(sc_matrix)==T || is.null(sc_metadata) == T){
      stop("No single cell object or metadata has been provided for creating signature.")
    }else{
      signatures = create_sc_signatures(sc_matrix, sc_metadata, cell_label, sample_label, credentials.mail = credentials.mail, credentials.token = credentials.token,
                                        bulk_rna = raw.counts, cell_markers, name_signature = name_sc_signature)
    }
  }

  deconv_sig = compute_methods_variable_signature(TPM_matrix, signatures = path_signatures, algos = methods, exclude = signatures_exclude, cbsx.name = credentials.mail, cbsx.token = credentials.token, doParallel, workers, sc_matrix)

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
    if(is.null(sc_matrix)==T){
      stop("No single cell object has been provided for deconvolution.")
    }else{
      deconv_sc = compute_sc_deconvolution_methods(raw.counts, normalized = normalized, sc_matrix, sc_metadata, cell_label,
                                                   sample_label, name_sc_signature, n_cores = workers)
      all_deconvolution_table = cbind(data.frame(all_deconvolution_table), deconv_sc)
    }
  }

  deconvolution = compute.deconvolution.preprocessing(data.frame(all_deconvolution_table))

  if(return == TRUE){
    utils::write.csv(deconvolution, paste0("Results/Deconvolution_", file_name, ".csv"))
  }

  return(deconvolution)

}

#' Compute second-generation deconvolution methods
#'
#' @param raw_counts A matrix with raw counts (samples as columns and genes symbols as rows)
#' @param normalized Boolean value to specify if raw_counts need to be normalized (If no raw_counts are available and argument corresponds to already normalized counts this arguments needs to be set to False)
#' @param sc_object A matrix with the counts from scRNAseq object (genes as rows and cells as columns)
#' @param sc_metadata Dataframe with metadata from the single cell object. The matrix should include the columns cell_label and sample_label.
#' @param cell_annotations A character vector with the cell labels (need to be of the same order as in the sc_object)
#' @param samples_ids A character vector with the samples labels (need to be of the same order as in the sc_object)
#' @param name_object Signature name to use in the generated single cell signature for deconvolving the bulk RNAseq data
#' @param n_cores Number of cores to use for paralellization. If no number is set, detectCores() - 1 will be set as the number.
#' @param return Whether to save or not the csv file with the deconvolution features.
#' @param file_name File name for the .csv file to save with the deconvolution results.
#'
#' @return A matrix of deconvolution features across samples from your bulk counts based on the second generation methods.
#' @export
#'
#' @references
#' Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.
#' Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363
#'
#' Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Dietrich, Alexander and Merotto, Lorenzo and Pelz, Konstantin and Eder, Bernhard and Zackl, Constantin and Reinisch, Katharina and
#' Edenhofer, Frank and Marini, Federico and Sturm, Gregor and List, Markus and Finotello, Francesca. (2024) https://doi.org/10.1101/2024.06.10.598226
#'
compute_sc_deconvolution_methods = function(raw_counts, normalized = TRUE, sc_object, sc_metadata, cell_annotations, samples_ids, name_object, n_cores = NULL, return = FALSE, file_name = NULL){

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
  autogenes = omnideconv::deconvolute_autogenes(bulk_gene_expression = bulk_counts, single_cell_object = as.matrix(sc_object),
                                                cell_type_annotations = as.character(sc_metadata[,cell_annotations]), verbose = T)$proportions

  message("\nRunning BayesPrism...............................................................\n")
  bayesprism = omnideconv::deconvolute_bayesprism(bulk_gene_expression = raw_counts, single_cell_object = as.matrix(sc_object),
                                                  cell_type_annotations = as.character(sc_metadata[,cell_annotations]), n_cores = n_cores)$theta

  message("\nRunning Bisque...............................................................\n")
  bisque = omnideconv::deconvolute_bisque(bulk_gene_expression = as.matrix(raw_counts), single_cell_object = as.matrix(sc_object),
                                          cell_type_annotations = as.character(sc_metadata[,cell_annotations]), batch_ids = as.character(sc_metadata[,samples_ids]), verbose = T)$bulk_props

  message("\nRunning CPM...............................................................\n")
  cpm = omnideconv::deconvolute_cpm(bulk_gene_expression = data.frame(raw_counts), single_cell_object = as.matrix(sc_object), no_cores = 4, #raw counts
                                    cell_type_annotations = as.character(sc_metadata[,cell_annotations]), verbose = T)$cellTypePredictions

  message("\nRunning MuSiC...............................................................\n")
  music = omnideconv::deconvolute_music(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = as.matrix(sc_object),
                                        cell_type_annotations = as.character(sc_metadata[,cell_annotations]),  batch_ids = as.character(sc_metadata[,samples_ids]), verbose = T)$Est.prop.weighted

  message("\nRunning SCDC...............................................................\n")
  scdc = omnideconv::deconvolute_scdc(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = as.matrix(sc_object),
                                      cell_type_annotations = as.character(sc_metadata[,cell_annotations]), batch_ids = as.character(sc_metadata[,samples_ids]), verbose = T)$prop.est.mvw


  results = list(AutogeneS = autogenes, BayesPrism = bayesprism, Bisque = bisque,
                 CPM = cpm, MuSic = music, SCDC = scdc)


  results <- lapply(names(results), function(method) {
    deconv_method <- results[[method]]

    colnames(deconv_method) <- paste0(method, "_", name_object, "_", colnames(deconv_method))
    colnames(deconv_method) <- stringr::str_replace_all(colnames(deconv_method), " ", "_")

    return(deconv_method)
  })

  results = do.call(cbind, results)

  if(return == TRUE){
    utils::write.csv(results, paste0("Results/Deconvolution_sc_", file_name, ".csv"))
  }

  return(results)

  ### METHODS NOT YET IMPLEMENTED

  # 1. BSeq-sc: Need CIBERSORT source code
  # message("\nRunning BSeq-sc...............................................................\n")
  # bseqsc = omnideconv::deconvolute_bseqsc(bulk_gene_expression = as.matrix(bulk_counts), signature = signatures[["BSeqsc"]], verbose = T)

  # 2. CDSeq: Crash machine
  # message("\nRunning CDSeq...............................................................\n")
  # cdseq = omnideconv::deconvolute_cdseq(bulk_gene_expression = raw_counts, single_cell_object = as.matrix(sc_object), no_cores = n_cores,
  #                     cell_type_annotations = as.character(sc_metadata[,cell_annotations]), batch_ids = as.character(sc_metadata[,samples_ids]), verbose = T)

  # 3. SCADEN: Takes a lot of time
  # message("\nRunning Scaden...............................................................\n")
  # model_scaden <- omnideconv::build_model(bulk_gene_expression = bulk.data, counts.matrix, as.character(cell_annotations),
  #                                         batch_ids = samples_ids, method = "scaden")
  # scaden = deconvolute_scaden(bulk_gene_expression = bulk_counts,
  #                             signature = model_scaden, verbose = T)

}


#' Create meta-cells from a single cell object using the KNN algorithm. This function is adapted from the R package hdWGCNA (Morabito et al., 2023)
#'
#' @param sc_object A matrix with the counts from scRNAseq object (genes as rows and cells as columns)
#' @param labels_column A character vector with the cell labels (need to be of the same order as in the sc_object)
#' @param samples_column A character vector with the samples labels (need to be of the same order as in the sc_object)
#' @param exclude_cells Cell types to discard from metacell algorithm.
#' @param min_cells The minimum number of cells in a particular grouping to construct metacells.
#' @param k Number of nearest neighbors to aggregate for KNN algorithm.
#' @param max_shared The maximum number of cells to be shared across two metacells.
#' @param n_workers Number of cores to use for paralellization.
#' @param min_meta Minimum number of metacells allowed. Below this number, metacells of this cell type will be discarded.
#'
#' @return A list with two elements:
#' - The metacell count matrix (genes as rownames and cells as columns)
#' - The metadata matrix corresponding to the metacell object
#'
#' @export
#'
#' @references
#'
#' Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
#'
#' Morabito, S., Reese, F., Rahimzadeh, N., Miyoshi, E., & Swarup, V. (2023). hdWGCNA identifies co-expression networks in high-dimensional transcriptomics data. Cell Reports Methods, 3(6), 100498. https://doi.org/10.1016/j.crmeth.2023.100498
#'
#'
create_metacells = function(sc_object, labels_column, samples_column, exclude_cells = NULL, min_cells = 50, k = 15, max_shared = 15, n_workers = 4, min_meta = 10){

  message("\nCreating metacells...............................................................\n")
  ## Setup sc object
  data <- hdWGCNA::SetupForWGCNA(
    sc_object,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "MetaCells" # the name of the hdWGCNA experiment
  )

  rm(sc_object)
  gc()

  ### Parallelize work
  data@meta.data$cells_labels = data@meta.data[,labels_column]
  data@meta.data$samples_ids = data@meta.data[,samples_column]
  Seurat::Idents(data) = data@meta.data$cells_labels
  subset_data = list()
  contador = 1
  cells_ids = unique(data@meta.data$cells_labels)
  cells = cells_ids[!cells_ids %in% exclude_cells]
  for (cell_type in cells) {
    for (patient in unique(data@meta.data$samples_ids)) {
      x = subset(x = data, subset = samples_ids == patient, idents = cell_type, return.null = T)
      if(is.null(x) == F){
        subset_data[[contador]] <- x
        contador = contador + 1
      }
    }
  }

  # Set up parallelization using the future package
  future::plan(future::multisession, workers = n_workers) # Adjust the number of workers based on your system
  # Run the function in parallel
  results <- future.apply::future_lapply(subset_data, FUN = process_group, min_cells, k, max_shared, labels_column, samples_column, future.seed = TRUE)
  results <- results[!sapply(results, is.null)]

  #Stop parallelization from running in the background
  future::plan(future::sequential)

  rm(data, subset_data)
  gc()

  message("\nMetacells done!...............................................................\n")

  # Combine results into a single Seurat object (if needed)
  counts_sc = do.call(cbind, lapply(results, '[[', 1))
  metadata = do.call(rbind, lapply(results, '[[', 2))

  rm(results)
  gc()

  n_cells = table(metadata[,labels_column])
  low_count_cells <- n_cells[n_cells < min_meta]

  cat("\nNumber of metacells per cell type\n")
  print(n_cells)

  cat("\nRemoving metacells with less than", min_meta, "...............................................................\n")
  print(names(low_count_cells))

  metadata = metadata %>%
    dplyr::filter(!.data[[labels_column]] %in% names(low_count_cells))
  counts_sc = counts_sc[,colnames(counts_sc) %in% rownames(metadata)]

  return(list(Counts = counts_sc, Metadata = metadata))

}

#' Compute deconvolution benchmark
#'
#' @param deconvolution The deconvolution matrix output from compute.deconvolution()
#' @param groundtruth A matrix with the cell type proportions (samples as rows and cell types as columns). Cell types names should correspond to the ones on the deconvolution matrix.
#' @param cells_extra A string specifying the cells names to consider and that are not including in the nomenclature of multideconv (see Readme)
#' @param corr_type Secifies the type of correlations to compute ('spearman' or 'pearson').
#' @param scatter Boolean value to specify if scatter plots should be returned.
#' @param pval A numeric value with the pvalue to use for selecting significant features.
#' @param plot Boolean value to whether save or not the plot of the benchmark in the Results/ directory.
#' @param file_name A string specifying the name of the plot saved in Results/
#' @param width A numeric value with the width for the returned plot.
#' @param height A numeric value with the height for the returned plot.
#'
#' @return A correlation matrix between the cell type deconvolution combinations and the real cell proportions.
#' @export
#'
#' @examples
#'
#' data("deconvolution")
#' data("cells_groundtruth")
#'
#' corr_matrix = compute.benchmark(deconvolution, cells_groundtruth, cells_extra = "Myeloid.cells",
#'                                 corr_type = "pearson", scatter = FALSE)
#'
compute.benchmark = function(deconvolution, groundtruth, cells_extra = NULL, corr_type = "spearman", scatter = TRUE, plot = FALSE, pval = 0.05, file_name = NULL, width = 16, height = 8){

  groundtruth = groundtruth[rownames(deconvolution),] #Order samples to match both features

  cell_types = c("B.cells", "B.naive.cells", "B.memory.cells", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated",
                 "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated.cells", "Dendritic.resting.cells", "Cancer", "Endothelial",
                 "Eosinophils", "Plasma", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated.cells", "Mast.resting.cells", "CAF")

  cell_types = c(cell_types, cells_extra)

  pattern <- paste0("(_", gsub("\\.", "\\\\.", cell_types), ")$", collapse = "|")
  deconvolution_combinations <- unique(gsub(pattern, "", colnames(deconvolution)))

  deconvolution_combinations = gsub("(BPRNACan3DProMet|BPRNACanProMet|BPRNACan)", "\\1_", deconvolution_combinations)

  ###Correlation function
  corr_bench <- function(data, corr = "pearson", pval = 0.05) {
    M <- Hmisc::rcorr(as.matrix(data), type = corr)

    # Only keep the three matrix elements: r, P, n
    Mdf <- purrr::map(M[c("r", "P", "n")], ~data.frame(.x))

    corr_df <- Mdf %>%
      purrr::map(~tibble::rownames_to_column(.x, var = "measure1")) %>%
      purrr::map(~tidyr::pivot_longer(.x, -measure1, names_to = "measure2")) %>%
      dplyr::bind_rows(.id = "id") %>%
      tidyr::pivot_wider(names_from = id, values_from = value) %>%
      dplyr::mutate(
        r = as.numeric(r),
        P = as.numeric(P),
        sig_p = ifelse(P < pval, TRUE, FALSE),
        p_if_sig = ifelse(sig_p, P, NA),
        r_if_sig = ifelse(sig_p, r, NA)
      )

    return(corr_df)
  }

  #####Scatter plot function
  scatter_plots = function(deconv, ground, method){
    plots <- list()
    for (i in 1:ncol(deconv)) {
      data = cbind(ground, deconv[,i])
      colnames(data) = c("x", "y")
      cor_test <- stats::cor.test(data$x, data$y)
      cor_value <- cor_test$estimate  # Correlation coefficient
      p_value <- cor_test$p.value    # p-value

      p <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(color = "blue", size = 0.1, alpha = 0.7) +  # Customize the points
        ggplot2::geom_smooth(method = "lm", se = T, color = "red") +  # Add regression line
        ggplot2::theme_minimal() +  # Apply a minimal theme
        ggplot2::labs(
          x = colnames(ground),
          #title = paste0("Linear correlation - ", colnames(ground)),  # Set the title
          y = colnames(deconv)[i],                 # Set the x-axis label
        ) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 5),
          axis.text.x = ggplot2::element_text(size = 5),
          axis.text.y = ggplot2::element_text(size = 5),
          axis.title.y = ggplot2::element_text(size = 5)  # Adjust title font size and position
        ) +
        ggplot2::geom_text(
          ggplot2::aes(x = mean(data$x), y = max(data$y), label = paste("r =", round(cor_value, 2), ", pval = ", round(p_value, 2))),  # Add the correlation text
          size = 2,  # Adjust the font size of the correlation coefficient text
          hjust = 0.5,  # Adjust the horizontal alignment of the text
          vjust = -1     # Adjust the vertical position of the text
        )

      plots[[i]] <- p
    }
    return(plots)
  }

  cell_clusters = colnames(groundtruth)

  ###Correlation matrix
  corr_matrix = data.frame(matrix(ncol = length(deconvolution_combinations), nrow = length(cell_clusters)))
  pval_matrix = data.frame(matrix(ncol = length(deconvolution_combinations), nrow = length(cell_clusters)))
  rownames(corr_matrix) = cell_clusters
  colnames(corr_matrix) = deconvolution_combinations
  rownames(pval_matrix) = cell_clusters
  colnames(pval_matrix) = deconvolution_combinations
  cells_discard = c()
  plots_all = list()

  ###Correlation computation
  for (i in 1:length(cell_clusters)) {
    idx = grep(paste0("_", cell_clusters[i], "$"), colnames(deconvolution))
    if(length(idx)==0){
      cells_discard = c(cells_discard, cell_clusters[i])
    }

    deconv = deconvolution[,idx, drop = F]

    ground = groundtruth[,cell_clusters[i],drop=F]

    ###Scatter plots
    if(scatter == T){
      if(ncol(deconv)!=0){
        # Edit deconv to match with reference and have the same number of scatter plots per cell type
        deconv_sub = deconv
        colnames(deconv_sub) <- gsub(paste0("_", cell_clusters[i], "$"), "", colnames(deconv_sub))
        colnames(deconv_sub) = gsub("(BPRNACan3DProMet|BPRNACanProMet|BPRNACan)", "\\1_", colnames(deconv_sub))
        missing_cols <- setdiff(colnames(corr_matrix), colnames(deconv_sub)) #Find which columns are missing from reference

        # Fill with 0 no-existing columns (missing from reference)
        for (col in missing_cols) {
          deconv_sub[[col]] <- 0
        }
        deconv_sub <- deconv_sub[, names(corr_matrix)] #reorder columns matching the reference

        plots = scatter_plots(deconv_sub, ground, corr_type)

        plots_all[[i]] = plots
      }
    }

    x = corr_bench(cbind(deconv, ground), corr_type, pval)
    x = x[which(x$measure1==colnames(ground)),] #only taking corr against ground truth

    for (j in 1:ncol(corr_matrix)) {
      idx = grep(colnames(corr_matrix)[j], x$measure2)
      if(length(idx) == 0){
        corr_matrix[i,j] = NA
      }else{
        corr_matrix[i,j] = x$r[idx]
      }
    }

    for (j in 1:ncol(pval_matrix)) {
      idx = grep(colnames(pval_matrix)[j], x$measure2)
      if(length(idx) == 0){
        pval_matrix[i,j] = NaN
      }else{
        pval_matrix[i,j] = x$P[idx]
      }
    }
  }


  ###Benchmarking plot
  if(length(cells_discard)>0){
    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%cells_discard),]
    pval_matrix = pval_matrix[-which(rownames(pval_matrix)%in%cells_discard),]
  }

  corr_matrix[nrow(corr_matrix)+1,] = colMeans(corr_matrix, na.rm = T)
  rownames(corr_matrix)[nrow(corr_matrix)] = "average"

  pval_matrix[nrow(pval_matrix)+1,] = 0
  rownames(pval_matrix)[nrow(pval_matrix)] = "average"

  ##Order methods
  corr_matrix = t(corr_matrix) %>%
    data.frame() %>%
    dplyr::arrange(average) %>%
    t() %>%
    data.frame()

  corr_df <- reshape2::melt(corr_matrix)
  pval_df = reshape2::melt(pval_matrix[,colnames(corr_matrix)]) #Take the same order as corr_matrix

  corr_df = corr_df %>%
    dplyr::mutate(Cells = rep(rownames(corr_matrix), ncol(corr_matrix)),
                  pval_value = pval_df$value)

  g <- corr_df %>%
    ggplot2::ggplot(ggplot2::aes(Cells, variable, fill=value, label=round(value,2))) +
    ggplot2::geom_tile() +
    ggplot2::labs(x = NULL, y = NULL, fill = paste0(corr_type, "'s\nCorrelation"), title=file_name, subtitle = paste0("Only showing significant correlations (<", pval, ")")) +
    ggplot2::scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    ggplot2::geom_text(data = subset(corr_df, pval_value <= pval)) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete(expand=c(0,0)) +
    ggplot2::scale_y_discrete(expand=c(0,0)) +
    ggpubr::rotate_x_text(angle = 45) + ggplot2::theme(axis.text.x=ggtext::element_markdown()) + ggplot2::theme(axis.text.y=ggtext::element_markdown())

  if(plot){
    grDevices::pdf(paste0("Results/Benchmark_plot_", file_name,".pdf"), width = width, height = height)
    plot(g)
    grDevices::dev.off()
  }

  return(corr_matrix)

}

#' Create pseudo bulk from single cell object
#'
#' @param sc_obj A Seurat single cell object
#' @param cells_labels A character vector with the cell labels (need to be of the same order as in the sc_obj)
#' @param sample_labels A character vector with the samples labels (need to be of the same order as in the sc_obj)
#' @param normalized Whether pseudobulk should be or not TPM normalized
#' @param file_name A string specifying the name of the .csv pseudobulk saved in Results/
#'
#' @return A gene count matrix (genes as rows and samples as columns)
#' @export
#'
create_sc_pseudobulk = function(sc_obj, cells_labels, sample_labels, normalized = TRUE, file_name){

  #Convert to SingleCell
  sc_obj@meta.data$Patient = as.factor(sc_obj@meta.data[,sample_labels])
  sc_obj@meta.data$new_annotation = as.factor(sc_obj@meta.data[,cells_labels])
  sce = Seurat::as.SingleCellExperiment(sc_obj)

  ##Aggregating counts
  aggr_counts <- glmGamPoi::pseudobulk(sce, group_by = glmGamPoi::vars(Patient), aggregation_functions = list(counts = "rowMeans2", .default = "rowMeans2"))
  pseudo_counts = data.frame(aggr_counts@assays@data$counts)

  if(normalized == TRUE){
    pseudo_counts = ADImpute::NormalizeTPM(pseudo_counts, log=F) %>%
      data.frame()
  }

  #Save pseudobulk matrix
  utils::write.table(pseudo_counts, file = paste0("Results/", file_name, ".csv"), quote = F, sep = "\t", row.names = T)

  return(pseudo_counts)

}


#' Create cell type signatures from scRNAseq
#'
#'
#' @param sc_obj A matrix with the counts from scRNAseq object (genes as rows and cells as columns)
#' @param sc_metadata Dataframe with metadata from the single cell object. The matrix should include the columns cell_label and sample_label.
#' @param cells_labels A character vector with the cell labels (need to be of the same order as in the sc_object)
#' @param sample_labels A character vector with the samples labels (need to be of the same order as in the sc_object)
#' @param credentials.mail (Optional) Credential email for running CIBERSORTx If not provided, CIBERSORTx method will not be run.
#' @param credentials.token (Optional) Credential token for running CIBERSORTx. If not provided, CIBERSORTx method will not be run.
#' @param bulk_rna A matrix of bulk data. Rows are genes, columns are samples. This is needed for MOMF method, if not given the method will not be run.
#' @param cell_markers Named list with the genes markers names as Symbol per cell types to be used to create the signature using the BSeq-SC method. If NULL, the method will be ignored during the signature creation.
#' @param name_signature A string indicating the signature name. This will be added as a suffix in each method (e.g. CBSX_name_signature, DWLS_name_signature)
#'
#' @return A list containing the cell signatures per method. Signatures are directly saved in Results/custom_signatures folder, these will be used to run deconvolution.
#' @export
#'
#' @references
#' Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.
#' Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363
#'
#' Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Dietrich, Alexander and Merotto, Lorenzo and Pelz, Konstantin and Eder, Bernhard and Zackl, Constantin and Reinisch, Katharina and
#' Edenhofer, Frank and Marini, Federico and Sturm, Gregor and List, Markus and Finotello, Francesca. (2024) https://doi.org/10.1101/2024.06.10.598226
#'
create_sc_signatures = function(sc_obj, sc_metadata, cells_labels, sample_labels, credentials.mail = NULL, credentials.token = NULL, bulk_rna = NULL, cell_markers = NULL, name_signature = NULL){

  signature_dir = "Results/custom_signatures/"

  sc_obj = as.matrix(sc_obj)

  cat("\nRunning DWLS...............................................................\n")
  model_dwls <- omnideconv::build_model_dwls(sc_obj, as.character(sc_metadata[,cells_labels]),
                                             dwls_method = "mast_optimized", ncores = 1) %>%
    data.frame() %>%
    tibble::rownames_to_column("NAME")

  utils::write.table(model_dwls, paste0(signature_dir, "DWLS-", name_signature,"-scRNAseq.txt"), row.names = F, quote = F, sep = "\t")

  signatures = list(DWLS = model_dwls)

  cat("\nRunning CIBERSORTx...............................................................\n")

  if(is.null(credentials.mail)==T || is.null(credentials.token)==T){
    warning("\nNo credentials were found\n")
    cat("\nPlease set your credentials in the function for running CIBERSORTx")
    cat("\nSkipping...............")
  }else{
    omnideconv::set_cibersortx_credentials(credentials.mail, credentials.token)
    model_cbsx = omnideconv::build_model(sc_obj, as.character(sc_metadata[,cells_labels]),
                                         batch_ids = as.character(sc_metadata[,sample_labels]), method = "cibersortx") %>%
      data.frame() %>%
      tibble::rownames_to_column("NAME")

    utils::write.table(model_cbsx, paste0(signature_dir, "CBSX-", name_signature,"-scRNAseq.txt"), row.names = F, quote = F, sep = "\t")
    signatures[[length(signatures) + 1]] = model_cbsx
    names(signatures)[length(signatures)] = "CBSX"
  }

  if(is.null(bulk_rna) == F){
    cat("\nRunning MOMF...............................................................\n")
    model_momf <- omnideconv::build_model_momf(sc_obj, as.character(sc_metadata[,cells_labels]),
                                               bulk_gene_expression = bulk_rna) %>%
      data.frame() %>%
      tibble::rownames_to_column("NAME")

    utils::write.table(model_momf, paste0(signature_dir, "MOMF-", name_signature,"-scRNAseq.txt"), row.names = F, quote = F, sep = "\t")

    signatures[[length(signatures) + 1]] = model_momf
    names(signatures)[length(signatures)] = "MOMF"
  }else{
    warning("\nMOMF requires the matrix of bulk data to deconvolve. As no bulk data has been provided, MOMF is skipped.")
  }

  if(is.null(cell_markers) == F){
    cat("\nRunning BSeq-sc...............................................................\n")
    model_bseq <- omnideconv::build_model_bseqsc(sc_obj, as.character(sc_metadata[,cells_labels]),
                                                 markers = cell_markers, batch_ids = as.character(sc_metadata[,sample_labels])) %>%
      data.frame() %>%
      tibble::rownames_to_column("NAME")

    utils::write.table(model_bseq, paste0(signature_dir, "BSeqSC-", name_signature,"-scRNAseq.txt"), row.names = F, quote = F, sep = "\t")

    signatures[[length(signatures) + 1]] = model_bseq
    names(signatures)[length(signatures)] = "BSeqsc"
  }else{
    warning("\nBSeq-sc requires cell type marker genes. As no cell_markers have been provided, BSeq-sc is skipped")
  }

  return(signatures)

}

unregister_dopar <- function() {
  if (!is.null(foreach::getDoParRegistered())) {
    # switch back to sequential backend
    foreach::registerDoSEQ()
    gc()
  }
}

# Function to process each group in metacells
process_group <- function(data, min_cells = 50, k = 15, max_shared = 15, labels_column, samples_column) {

  if (ncol(data) < min_cells) {
    cat("Skipping group: Less than", min_cells, "cells in this subset\n")
    return(NULL) # Return NULL for groups with insufficient cells
  }

  #Create Metacells by Groups celltype_patient
  seurat_obj = hdWGCNA::MetacellsByGroups(seurat_obj = data,
                                          min_cells = min_cells,
                                          group.by = c(labels_column, samples_column),
                                          reduction = 'pca',
                                          k = k,
                                          max_shared = max_shared,
                                          ident.group = labels_column)

  meta = hdWGCNA::GetMetacellObject(seurat_obj)

  result <- list(
    counts = as.matrix(meta@assays[["RNA"]]@counts),
    metadata = meta@meta.data
  )

  # Clean memory inside worker
  rm(data, seurat_obj, meta)
  gc()

  return(result)

}



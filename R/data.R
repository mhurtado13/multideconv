#' Cell labels
#'
#' A character vector or factor representing the cell type labels for each cell.
#'
#' @format A vector of length equal to the number of cells and in the same order as it is on the single cell object.
#'
#' @examples
#' data(cell_labels)
#' head(cell_labels)
"cell_labels"

#' Sample labels
#'
#' A character vector or factor representing the sample labels for each cell.
#'
#' @format A vector of length equal to the number of cells and in the same order as it is on the single cell object.
#'
#' @examples
#' data(sample_labels)
#' head(sample_labels)
"sample_labels"

#' Cell groundtruth
#'
#' A matrix with the cell type proportions of the single cell object.
#'
#' @format Matrix with samples as rows and cell types as columns
#'
#' @examples
#' data(cells_groundtruth)
#' head(cells_groundtruth)
"cells_groundtruth"

#' Deconvolution matrix
#'
#' A matrix with the cell type deconvolution features obtained from compute.deconvolution()
#'
#' @format Matrix with samples as rows and deconvolution features as columns
#'
#' @examples
#' data(deconvolution)
#' head(deconvolution)
"deconvolution"

#' Metacells data
#'
#' A gene expression matrix corresponding to the single cell object, in this case it corresponds to the reduced single cell object (metacells), not the original.
#'
#' @format Matrix with samples as columns and genes names (SYMBOL) as rows
#'
#' @source Senosain et al. (2023), doi: 10.1158/2767-9764.CRC-22-0373. PMID: 37501683; PMCID: PMC10370362
#'
#' @examples
#' data(metacells_data)
#' head(metacells_data)
"metacells_data"

#' Metacells metadata
#'
#' A dataframe corresponding to the metadata of the single cell object, in this case it corresponds to the reduced single cell object (metacells). It needs to include the columns specified as cell_labels and samples_labels.
#'
#' @format Matrix with samples as rows and meta information as columns
#'
#' @source Senosain et al. (2023), doi: 10.1158/2767-9764.CRC-22-0373. PMID: 37501683; PMCID: PMC10370362
#'
#' @examples
#' data(metacells_metadata)
#' head(metacells_metadata)
"metacells_metadata"

#' Pseudobulk
#'
#' A gene expression matrix after aggregation of the single cell object (pseudobulk).
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @examples
#' data(pseudobulk)
#' head(pseudobulk)
"pseudobulk"

#' Raw counts
#'
#' Raw gene expression matrix from bulk RNAseq.
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(raw_counts)
#' head(raw_counts)
"raw_counts"

#' Cell subgroups
#'
#' Cell subgroups composition
#'
#' @format A list with the cell subgroups
#'
#' @examples
#' data(subgroups)
#' subgroups[[1]]
"subgroups"

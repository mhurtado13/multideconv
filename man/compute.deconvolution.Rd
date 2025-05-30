% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_deconvolution.R
\name{compute.deconvolution}
\alias{compute.deconvolution}
\title{Compute deconvolution}
\usage{
compute.deconvolution(
  raw.counts,
  methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"),
  signatures_exclude = NULL,
  normalized = TRUE,
  doParallel = FALSE,
  workers = NULL,
  return = TRUE,
  create_signature = FALSE,
  credentials.mail = NULL,
  credentials.token = NULL,
  sc_deconv = FALSE,
  sc_matrix = NULL,
  sc_metadata = NULL,
  methods_sc = c("Autogenes", "BayesPrism", "Bisque", "CPM", "MuSic", "SCDC"),
  cell_label = NULL,
  sample_label = NULL,
  cell_markers = NULL,
  methods_sig = c("DWLS", "CIBERSORTx", "MOMF", "BSeqsc"),
  name_sc_signature = NULL,
  file_name = NULL
)
}
\arguments{
\item{raw.counts}{A matrix with the raw counts (samples as columns and genes symbols as rows)}

\item{methods}{A character vector with the deconvolution methods to run. Default are "Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"}

\item{signatures_exclude}{A character vector with the signatures to exclude from the src/signatures folder.}

\item{normalized}{If raw.counts are not available, user can input its normalized counts. In that case this argument need to be set to False.}

\item{doParallel}{Whether to do or not parallelization. Only CBSX and DWLS methods will run in parallel.}

\item{workers}{Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1}

\item{return}{Whether to save or not the csv file with the deconvolution features}

\item{create_signature}{Whether to create or not the signatures using the methods MOMF, CBSX, DWLS and BSeq-SC. If TRUE, sc_matrix shuld be provide.}

\item{credentials.mail}{(Optional) Credential email for running CIBERSORTx. If not provided, CIBERSORTx method will not be run.}

\item{credentials.token}{(Optional) Credential token for running CIBERSORTx. If not provided, CIBERSORTx method will not be run.}

\item{sc_deconv}{Whether to run or not deconvolution methods based on single cell.}

\item{sc_matrix}{If sc_deconv = T, the matrix of counts across cells from the scRNAseq object is provided.}

\item{sc_metadata}{Dataframe with metadata from the single cell object. The matrix should include the columns cell_label and sample_label.}

\item{methods_sc}{A character vector with the sc-deconvolution methods to run. Default are "Autogenes", "BayesPrism", "Bisque", "CPM", "MuSic", "SCDC"}

\item{cell_label}{If sc_deconv = T, a character vector indicating the cell labels (same order as the count matrix)}

\item{sample_label}{If sc_deconv = T, a character vector indicating the cell samples IDs (same order as the count matrix)}

\item{cell_markers}{Named list with the genes markers names as Symbol per cell types to be used to create the signature using the BSeq-SC method. If NULL, the method will be ignored during the signature creation.}

\item{methods_sig}{A character vector specifying which methods to run. Options are "DWLS", "CIBERSORTx", "MOMF", and "BSeqsc". Default runs all available methods.}

\item{name_sc_signature}{If sc_deconv = T, the name you want to give to the signature generated}

\item{file_name}{File name for the csv files and plots saved in the Results/ directory}
}
\value{
A matrix of cell type deconvolution features across samples
}
\description{
The function calculates cell abundance based on cell type signatures using different methods and signatures. Methods available are Quantiseq, MCP, XCell, CIBERSORTx, EpiDISH, DWLS and DeconRNASeq. Provided signatures included signatures based on bulk and methylation data (7 methods and 10 signature in total). Signatures are present in the src/signatures directory, user can add its own signatures by adding the .txt files in this same folder. Second generation methods to perform deconvolution based on single cell data are also available if scRNAseq object is provided.
}
\examples{

data("raw_counts")
data("cell_labels")
data("sample_labels")
data("metacells_data")
data("metacells_metadata")
data("pseudobulk")

deconv = compute.deconvolution(raw_counts, normalized = TRUE,
                               methods = c("Epidish", "DeconRNASeq"), return = FALSE)


}
\references{
Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.
Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363

Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Dietrich, Alexander and Merotto, Lorenzo and Pelz, Konstantin and Eder, Bernhard and Zackl, Constantin and Reinisch, Katharina and
Edenhofer, Frank and Marini, Federico and Sturm, Gregor and List, Markus and Finotello, Francesca. (2024) https://doi.org/10.1101/2024.06.10.598226
}

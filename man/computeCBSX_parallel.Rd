% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_deconvolution.R
\name{computeCBSX_parallel}
\alias{computeCBSX_parallel}
\title{Compute CibersortX (CBSX) in parallel across multiple signatures}
\usage{
computeCBSX_parallel(TPM_matrix, signatures, name, password, workers)
}
\arguments{
\item{TPM_matrix}{A matrix with TPM normalized counts (genes symbols as rows and samples as columns).}

\item{signatures}{Path where signatures files are located}

\item{name}{Credential email for running CibersortX.}

\item{password}{Credential token for running CibersortX.}

\item{workers}{Number of processes available to run on parallel.}
}
\value{
A matrix with cell abundance deconvolve with CBSX
}
\description{
Compute CibersortX (CBSX) in parallel across multiple signatures
}

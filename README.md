# multideconv

Integrative pipeline for cell type deconvolution from bulk RNAseq using first and second generation methods

### Install from source
Clone the repository:
```
git clone https://github.com/mhurtado13/multideconv
```

## Environment

multideconv uses [renv](https://rstudio.github.io/renv/index.html) for creating a reproducible r-environment that can be easily share between users, this, to allow reproducibility and avoid libraries conflicts. Setting it up will install the neccessary packages, along with their specific versions in an isolated environment. 

For this, open the project `multideconv.Rproj` and in the R console run:

```r
# To activate the R environment (if you are using it for the first time)
renv::activate()
# To download and install all the require libraries and packages (if you are using it for the first time)
renv::restore() 
```

Once all packages have been installed, you can start creating your own scripts but be sure to still be inside the .Rproj!

Note that this is a once-step only when running multideconv for the first time. For the following times, you will only need to open the `multideconv.Rproj` and you are ready to go!

If you are doing other analyses that require the installation of extra libraries not present in the environment, you can install them as usual but after that, make sure to run `renv::snapshot()` to update your environment.

Make sure to run `renv::deactivate()` when finishing running multideconv, to avoid conflicts whenever you start a different project.

For more information about how `renv` works, visit https://rstudio.github.io/renv/articles/renv.html.

## Cell type deconvolution

multideconv performs and integrates results of cell type deconvolution using XX methods including first and second generation algorithms combined with different signatures based on gene expression, methylation data and single cell (Figure 1). 

<p align="center">
 <img src="man/Cell_deconvolution.png?raw=true" />
</p>

<p align="center"><i>
   Figure 1. Cell type deconvolution.
</i></p>

## Deconvolution processing

To process our deconvolution features, multideconv applied a combination of unsupervised filtering techniques and iterative linear based correlations to form subgroups across all our cell types. The algorithm is presented in Figure # and as output it returns a simplified deconvolution matrix composed with subgroups of methods-signatures whose expression is similar across samples. Additionally, it returns two matrices of cell-types with high proportion of zeros across samples and a matrix composed of discarded cell types not present in the categories mentioned above. Finally, it also gives to the user the list of subgroups for each type of correlation and a list of high-correlated features after pairwise analysis. 

<p align="center">
 <img src="man/Deconvolution_pipeline.png?raw=true" />
</p>

<p align="center"><i>
  Figure 2. Deconvolution processing algorithm.
</i></p>

## Nomenclature handling

In order to be able to process all deconvolution features as a result of using different methods and signatures from single cell and bulk RNAseq data, multideconv uses an universal nomenclature for 35 cell types (Table #). This to ensure concordance in downstream analysis when comparing and subgrouping results of different methods. We handle the nomenclature automatically when the user it is using the output of compute.deconvolution() (so don't worry!). However if you are using your own deconvolution results or using scRNAseq data please make sure to use the nomenclature provided in Figure 3. 

<p align="center">
 <img src="man/cell_types.png?raw=true" />
</p>

<p align="center"><i>
   Figure 3. Cell types nomenclature.
</i></p>

## How to add cell types other than the ones present in the nomenclature?

If you want multideconv to consider other cells, it is pretty simple! Just use the argument cells_extra in the function compute.deconvolution.analysis(). 

Let's say you want to add mesenchymal and basophils cells:

```r
processed_deconvolution = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123, cells_extra = c("mesenchymal", "basophils")) 
```

And that's it, just make sure the name you are putting in cells_extra is exactly the name of your cells in your deconvolution matrix!

## How do I use my single cell data for deconvolution?

multideconv includes second generation methods that allows users to input single cell data and use it for constructing signatures 'on-the-fly' and deconvolve the bulk RNAseq data.

Because of scRNAseq can have big and sparse matrix, we applied 'preprocessing' steps to avoid crashing the computer. For this, multideconv construct metacells from the single cell data in a parallelize manner and the compute deconvolution in the reduced object (Figure 4).


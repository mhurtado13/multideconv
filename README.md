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

## Deconvolution processing

To process our deconvolution features, multideconv applied a combination of unsupervised filtering techniques and iterative linear based correlations to form subgroups across all our cell types. The algorithm is presented in Figure # and as output it returns a simplified deconvolution matrix composed with subgroups of methods-signatures whose expression is similar across samples. Additionally, it returns two matrices of cell-types with high proportion of zeros across samples and a matrix composed of discarded cell types not present in the categories mentioned above. Finally, it also gives to the user the list of subgroups for each type of correlation and a list of high-correlated features after pairwise analysis. 

## Nomenclature handling



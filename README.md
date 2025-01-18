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

## General usage

The function calculates cell abundance based on cell type signatures using different methods and signatures through the function `compute.deconvolution` which takes as input the bulk RNAseq gene expression matrix either as raw or normalized counts. See `examples/deconvolve.Rmd` to see a detailed explanation of the parameters inside the function.

```r
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", file_name = "Tutorial") 
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", methods = c("Quantiseq", "MCP", "XCell", "DWLS"), file_name = "Test") 
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", signatures_exclude = "BPRNACan", file_name = "Tutorial")
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", sc_deconv = T, sc_matrix = sc.object, cell_annotations = cell_labels, cell_samples = bath_ids, name_sc_signature = "Signature_test", file_name = "Test")
```

**Note**: CIBERSORTx is included in the deconvolution methods, but it's not an open-source program. To run it, please ask for a token in [CibersortX](https://cibersortx.stanford.edu/register.php) and once obtained, provided your username and password on the parameters `credentials.mail` and `credentials.token`. 

Also CIBERSORTx ask to have [docker](https://kinsta.com/blog/install-docker-ubuntu/#installing-docker-desktop-on-ubuntu) install in your computer. Once installed, be sure to concede permission to manage docker as a non-root user. For verify this, go to the terminal and run:

```
docker ps
```

If you don't have any error, congrats you are go to go!

If you receive an error of permission, just run: 
```
sudo groupadd docker
sudo usermod -aG docker ${USER}
```

Then restart your computer and try again
```
docker ps
```

Now, you should be able to run deconvolution without any problems :)

For processing the deconvolution features obtained from `compute.deconvolution`, you can use the `compute.deconvolution.analysis` function.

```r
processed_deconvolution = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123, return = T)
```

Users can also compute second-generation deconvolution methods using their single cell data. For this use the function `compute_sc_deconvolution_methods`. Remember that this function is already included in `compute.deconvolution` when setting `sc_deconv = T`.

```r
deconv_sc = compute_sc_deconvolution_methods(raw.counts, sc_matrix, cell_annotations, cell_samples, name_sc_signature, normalized = T, n_cores = 4, cbsx_name = "XXX", cbsx_token = "XXX")
```

## How to add cell types other than the ones present in the nomenclature?

If you want multideconv to consider other cells, it is pretty simple! Just use the argument cells_extra in the function compute.deconvolution.analysis(). 

Let's say you want to add mesenchymal and basophils cells:

```r
processed_deconvolution = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123, cells_extra = c("mesenchymal", "basophils")) 
```

And that's it, just make sure the name you are putting in cells_extra is exactly the name of your cells in your deconvolution matrix!

## How to add other signatures?

You can include other signatures into the analysis by adding them as .txt into the folder `signatures`.

## How does my single cell data is used for deconvolution?

multideconv includes second generation methods that allows users to input single cell data and use it for constructing signatures 'on-the-fly' and deconvolve the bulk RNAseq data.

Because most of the times scRNAseq data can be big and sparse matrices, we applied 'preprocessing' steps to avoid crashing the computer. For this, multideconv construct metacells per cell type and patient from the single cell data using the KNN algorithm. This is done essentially by executing different tasks per cluster in parallel for N workers and once finished, it computes deconvolution using the reduced single cell object (Figure 4).

<p align="center">
 <img src="man/scRNAseq_deconvolution.png?raw=true" />
</p>

<p align="center"><i>
   Figure 4. Single cell RNAseq metacells construction
</i></p>

## Additional notes
multideconv is fully implemented and adapted as part of [CellTFusion](https://github.com/VeraPancaldiLab/CellTFusion) R package, an integration tool for cell deconvolution couple with transcription factor activities to deconvolve cell states of the tumor microenvironment. If you would like to use it in this context we invite you to visit our github repository.

## Citing multideconv

If you use multideconv in a scientific publication, we would appreciate citation to the :

```
XXXXX
```

## Acknowledgements

This repository was created by [Marcelo Hurtado](https://github.com/mhurtado13) in the [Network Biology for Immuno-oncology (NetB(IO)²)](https://www.crct-inserm.fr/en/netbio2_en/) group at the Cancer Research Center of Toulouse in supervision of [Vera Pancaldi](https://github.com/VeraPancaldi).

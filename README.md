multideconv
================

An integrative pipeline for combining first and second generation cell
type deconvolution results

<!-- badges: start -->
<!-- badges: end -->
<p align="center">
<img src="man/figures/overview.png?raw=true"/>
</p>
<p align="center">
<em>Figure 1. A schematic overview of the `multideconv` pipeline</em>
</p>

## Installation

You can install the development version of `multideconv` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pkg_install("VeraPancaldiLab/multideconv")
```

## General usage

These are basic examples which shows you how to use `multideconv` for
different tasks. For a detailed tutorial, see [Get
started](https://VeraPancaldiLab.github.io/multideconv/articles/multideconv.html)

``` r
library(multideconv)
```

The function calculates cell abundance based on cell type signatures
using different methods and signatures through the function
`compute.deconvolution` which takes as input the bulk RNAseq gene
expression matrix either as raw or normalized counts.

``` r
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", file_name = "Tutorial") 
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", methods = c("Quantiseq", "MCP", "XCell", "DWLS"), file_name = "Test") 
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", signatures_exclude = "BPRNACan", file_name = "Tutorial")
deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", sc_deconv = T, sc_matrix = sc.object, cell_label = cell_labels, sample_label = bath_ids, name_sc_signature = "Signature_test", file_name = "Test")
```

**NOTE**: `CIBERSORTx` is included in the deconvolution methods, but
it’s not an open-source program. To run it, please ask for a token in
[CIBERSORTx](https://cibersortx.stanford.edu/register.php) and once
obtained, provided your username and password on the parameters
`credentials.mail` and `credentials.token`.

Also CIBERSORTx ask to have
[docker](https://kinsta.com/blog/install-docker-ubuntu/#installing-docker-desktop-on-ubuntu)
install in your computer. Once installed, be sure to concede permission
to manage docker as a non-root user. For verify this, go to the terminal
and run:

    docker ps

If you don’t have any error, congrats you are go to go!

If you receive an error of permission, just run:

    sudo groupadd docker
    sudo usermod -aG docker ${USER}

Then restart your computer and try again

    docker ps

Now, you should be able to run deconvolution using CIBERSORTx without
any problems :)

For processing the deconvolution features obtained from
`compute.deconvolution`, you can use the
`compute.deconvolution.analysis` function.

``` r
processed_deconvolution = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123, return = T)
```

Users can also compute second-generation deconvolution methods using
their single cell data. For this use the function
`compute_sc_deconvolution_methods`. Remember that this function is
already included in `compute.deconvolution` when setting
`sc_deconv = T`.

``` r
deconv_sc = compute_sc_deconvolution_methods(raw_counts, sc_object, sc_metadata, cell_annotations, samples_ids, name_object, normalized = T, n_cores = 4, cbsx_name = "XXX", cbsx_token = "XXX")
```

## Cell types nomenclature

`multideconv` works based on established cell naming conventions (Figure
2) to simplify analysis and processing. Thus, if you would like to use
your own deconvolution results or signatures, please make sure to follow
these formats.

<p align="center">
<img src="man/figures/cell_types.png?raw=true"/>
</p>
<p align="center">
<em>Figure 2. Cell types nomenclature for `multideconv`</em>
</p>

## How to add cell types other than the ones present in the nomenclature?

If you want multideconv to consider other cells, it is pretty simple!
Just use the argument cells_extra in the function
`compute.deconvolution.analysis()`.

Let’s say you want to add mesenchymal and basophils cells:

``` r
processed_deconvolution = compute.deconvolution.analysis(deconvolution, corr = 0.7, seed = 123, cells_extra = c("mesenchymal", "basophils")) 
```

And that’s it, just make sure the name you are putting in cells_extra is
exactly the name of your cells in your deconvolution matrix!

## How to add other signatures?

You can include other signatures into the analysis by adding them as
.txt into the folder `Results/custom_signatures`.

## How does my single cell data is used for deconvolution?

`multideconv` includes second generation methods that allows users to
input single cell data and use it for constructing signatures
‘on-the-fly’ and deconvolve the bulk RNAseq data.

Because most of the times scRNAseq data can be big and sparse matrices,
we applied ‘preprocessing’ steps to avoid crashing the computer. For
this, `multideconv` construct metacells per cell type and patient from
the single cell data using the KNN algorithm. This is done essentially
by executing different tasks per cluster in parallel for N workers and
once finished, it computes deconvolution using the reduced single cell
object.

## Authors

`multideconv` was developed by [Marcelo
Hurtado](https://github.com/mhurtado13) in supervision of [Vera
Pancaldi](https://github.com/VeraPancaldi) and is part of the
[Pancaldi](https://github.com/VeraPancaldiLab) team. Currently, Marcelo
is the primary maintainer of this package.

## Citing multideconv

If you use `multideconv` in a scientific publication, we would
appreciate citation to the :

*Hurtado, M., Essabbar, A., Khajavi, L., & Pancaldi, V. (2025).
multideconv – Integrative pipeline for cell type deconvolution from bulk
RNAseq using first and second generation methods. bioRxiv.
<https://doi.org/10.1101/2025.04.29.651220>*

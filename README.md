Introduction to Sincast\!
================
Yidi Deng
JAN, 11th, 2021





























## Query Bian etal.(2020) with Sincast

This is the very first Sincast version\! A R package will be release
soon\! Here I will give a very brife introduction to Sincast using the

1.  Query data from [Bian et
    al.(2020)](https://doi.org/10.1038/s41586-020-2316-7): **Deciphering
    human macrophage development at single-cell resolution.**

2.  Referncce data from [Rajab et
    al.(2021)](https://doi.org/10.1016/j.stemcr.2021.04.010): \*\*An
    integrated analysis of human myeloid cells identifies gaps in in
    vitro models of in vivo biology

The query data can be downloaded at NCBI
[GSE133345](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133345),
and the reference data can be downloaded at the [Stemformatic data
protal](https://www.stemformatics.org/atlas/myeloid). They are also
available at our github page.

Now download [sincast
package.R](https://github.com/meiosis97/Sincast/blob/main/sincast%20package.R)
to get started\!

## Load functions

Put
[matMult.cpp](https://github.com/meiosis97/Sincast/blob/main/matMult.cpp)
into your R working directory and run

``` r
source('YOUR_DIR/sincast package.R')
```

This will automatically dowload the reuired R packages (a propmt will
show up), and load Sincast functions into your global environment.

## Read your data
!(this is a ploy)[./SincastDemo_atlas.png]

``` r
#load the query data
#query.annotation <- read.table('GSE133345_Annotations_of_all_1231_embryonic_cells_updated_0620.txt')
#query.data <- read.table('GSE133345_Quality_controled_UMI_data_of_all_1231_embryonic_cells.txt')
#load the reference data
#reference.data <- read.table('stemformatics_atlas_myeloid.1.0.expression_filtered.txt')
#reference.annotation <- read.delim('stemformatics_atlas_myeloid.1.0.samples.tsv', row.names = 1)
```

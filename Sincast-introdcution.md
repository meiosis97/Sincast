Query Bian etal.(2020) with Sincast
================
Yidi Deng
JAN, 11th, 2021





























## Download data

We first load the query data from [Bian et al.(2017)](): **Deciphering
human macrophage development at single-cell resolution.** The data can
be downloaded at NCBI
[GSE133345](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133345),
as well as our [github page]()

``` r
#load the query data
query.annotation <- read.table('GSE133345_Annotations_of_all_1231_embryonic_cells_updated_0620.txt')
query.data <- read.table('GSE133345_Quality_controled_UMI_data_of_all_1231_embryonic_cells.txt')
#load the reference data
reference.data <- read.table('stemformatics_atlas_myeloid.1.0.expression_filtered.txt')
reference.annotation <- read.delim('stemformatics_atlas_myeloid.1.0.samples.tsv', row.names = 1)
```

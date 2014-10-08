---
layout: page
title: Sample input data files
---

The internal data structure for [R/qtl2](http://kbroman.org/qtl2) is
different from that of [R/qtl](http://www.rqtl.org), and the input
data file format has also changed. Details on the new internal data
structure are in the
[R/qtl2 developer guide](../assets/vignettes/developer.html).  Details
on the new data file format are described in
[a separate vignette](../assets/vignettes/input_files.html).

For simple crosses (such as a backcross or intercross), one can
continue to use the [old R/qtl formats](http://rqtl.org/sampledata/),
load them with `qtl::read.cross()`, and then convert the data to the
new format with `qtl2::convert2cross2()`.

The following are sample input data files in the new R/qtl2 format.

---

### RIL by selfing

Data from
[Moore et al. (2013) Genetics 195:1077-1086](http://www.genetics.org/content/195/3/1077.abstract)
(the second replicate of RILs).

- [`grav2.yaml`](../assets/sampledata/grav2/grav2.yaml), the control file ([YAML format](http://www.yaml.org/))
- [`grav2_geno.csv`](../assets/sampledata/grav2/grav2_geno.csv), genotype data
- [`grav2_gmap.csv`](../assets/sampledata/grav2/grav2_gmap.csv), genetic map
- [`grav2_pheno.csv`](../assets/sampledata/grav2/grav2_pheno.csv), phenotype data
- [`grav2_phenocovar.csv`](../assets/sampledata/grav2/grav2_phenocovar.csv), phenotype covariates
  (times, in hours, for each of the phenotype columns in [`grav2_pheno.csv`](../assets/sampledata/grav2/grav2_pheno.csv))

You can load these data into R as follows:

    library(qtl2)
    grav2 <- read_cross2("http://kbroman.org/qtl2/assets/sampledata/grav2/grav2.yaml")

You can also [peruse the data at GitHub](https://github.com/kbroman/qtl2/tree/gh-pages/assets/sampledata/grav2).

---

### An intercross

Data from [Grant et al. (2006) Hepatology 44:174-185](http://www.ncbi.nlm.nih.gov/pubmed/16799992)

- [`iron.yaml`](../assets/sampledata/iron/iron.yaml), the control file ([YAML format](http://www.yaml.org/))
- [`iron_geno.csv`](../assets/sampledata/iron/iron_geno.csv), genotype data
- [`iron_gmap.csv`](../assets/sampledata/iron/iron_gmap.csv), genetic map
- [`iron_covar.csv`](../assets/sampledata/iron/iron_covar.csv), covariate data (sex and cross direction)
- [`iron_pheno.csv`](../assets/sampledata/iron/iron_pheno.csv), phenotype data (strictly numeric)
- [`iron_phenocovar.csv`](../assets/sampledata/iron/iron_phenocovar.csv), phenotype covariates
  (a bit silly, really; just indicates that the phenotype columns name are
  the tissues that were measured).

You can load these data into R as follows:
  
    library(qtl2)
    iron <- read_cross2("http://kbroman.org/qtl2/assets/sampledata/iron/iron.yaml")

You can also [peruse the data at GitHub](https://github.com/kbroman/qtl2/tree/gh-pages/assets/sampledata/iron).

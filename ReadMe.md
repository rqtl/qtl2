### R/qtl2db

[![Build Status](https://travis-ci.org/rqtl/qtl2db.svg?branch=master)](https://travis-ci.org/rqtl/qtl2db)

[Karl Broman](http://kbroman.org)

[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](https://rqtl.org), to better handle
high-dimensional data and complex cross designs. It is split into
multiple packages:

- [qtl2geno](https://github.com/rqtl/qtl2geno), for calculating genotype
  probabilities, imputations, and genetic maps
- [qtl2scan](https://github.com/rqtl/qtl2scan), for QTL genome scans and
  related calculations
- [qtl2plot](https://github.com/rqtl/qtl2plot), for data visualization
- [qtl2convert](https://github.com/rqtl/qtl2convert),
  for converting data among the R/qtl2,
  [DOQTL](https://www.bioconductor.org/packages/release/bioc/html/DOQTL.html),
  and [R/qtl](https://rqtl.org) formats
- [qtl2db](https://github.com/qtl2db), for connecting to genome databases

---

### Installation

R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), but
it can be installed from a mini-CRAN at [rqtl.org](https://rqtl.org).
Make sure you have the latest version of [R (3.4.2)](https://cran.r-project.org).

    install.packages("qtl2", repos="https://rqtl.org/qtl2cran")

The [qtl2](https://github.com/rqtl/qtl2) package is
inspired by the
[tidyverse package](https://cran.r-project.org/package=tidyverse);
it is basically empty, but when you install it, the
[qtl2geno](https://github.com/rqtl/qtl2geno),
[qtl2scan](https://github.com/rqtl/qtl2scan),
[qtl2plot](https://github.com/rqtl/qtl2plot),
[qtl2convert](https://github.com/rqtl/qtl2convert) and
[qtl2db](https://github.com/rqtl/qtl2db) packages, plus a
bunch of dependencies, will be installed.

Alternatively, you can install R/qtl2 from its source on
[GitHub](https://github.com/rqtl). (But note that compiling the C++
code can be rather slow.)

On _Windows_, you'll need [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

On _Mac OS X_, you'll need the
[command-line developer tools](https://mac-how-to.gadgethacks.com/how-to/install-command-line-developer-tools-without-xcode-0168115/),
as well as [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS).

You then need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
and [RcppEigen](https://github.com/RcppCore/RcppEigen).
(Additional, secondary dependencies will also be installed.)

    install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))

Finally, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github("rqtl/qtl2")

---

### Usage

The [R/qtl2db](https://github.com/rqtl/qtl2db) package contains two
key functions, `create_genes_query_func()` and
`create_variants_query_func()`, for creating functions that access
databases of genes and variants.

For example, we've prepared SQLite databases with mouse genes, and with
variants (SNPs, indels, and structural variants) in the eight mouse
founder lines for the Collaborative Cross.

- [`cc_variants.sqlite` doi:10.6084/m9.figshare.5280229.v1](https://doi.org/10.6084/m9.figshare.5280229.v1)
- [`mouse_genes.sqlite` doi:10.6084/m9.figshare.5280238.v2](https://doi.org/10.6084/m9.figshare.5280238.v2)

If you download those two files locally, you could use the following
code to create "accessor" functions:

```r
library(qtl2db)
query_genes <- create_gene_query_func("mouse_genes.sqlite")
query_variants <- create_variant_query_func("cc_variants.sqlite")
```

You can then use these two functions to grab the records for mouse
genes on the one hand, and CC variants on the other. For example, to
grab the genes and variants in a 2 Mbp region centered at 97.5 Mbp on
mouse chromosome 2, we'd do the following:

```r
genes <- query_genes(2, 96.5, 98.5)
variants <- query_variants(2, 96.5, 98.5)
```

Why the complexity? We want [R/qtl2](http://kbroman.org) functions to
be able to query SNPs or genes, but we don't want to prescribe how the
databases might be set up; they could be sitting in a SQLite database
(as above), or they might be in some other kind of database or accessed via
some Web API. Rather than requiring the SNP and gene databases to be
in a certain form, we'll instead ask users to provide query functions
like `query_genes()` and `query_variants()`.

Finally, note that we've also created a smaller version of the mouse genes
database, containing just the records with `source=="MGI"`. If you
just want full genes, you can use this:

- [`mouse_genes_mgi.sqlite` doi:10.6084/m9.figshare.5286019.v1](https://doi.org/10.6084/m9.figshare.5286019.v1)

---

### Vignettes

- [user guide](http://kbroman.org/qtl2/assets/vignettes/user_guide.html)
- [input file formats](http://kbroman.org/qtl2/assets/vignettes/input_files.html)
  (see also the [sample data files](http://kbroman.org/qtl2/pages/sampledata.html))
- [developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html)
- [HMM benchmarks](http://kbroman.org/qtl2/assets/vignettes/hmm_benchmarks.html)
- [linear regression benchmarks](http://kbroman.org/qtl2/assets/vignettes/linreg_benchmarks.html)

---

#### License

[Licensed](License.md) under [GPL-3](https://www.r-project.org/Licenses/GPL-3).

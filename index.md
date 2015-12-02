---
layout: page
title: R/qtl2
tagline: QTL analysis for high-dimensional data and complex crosses
description: R/qtl2, a reimplementation of R/qtl to better handle high-dimensional data and complex cross designs
---

R/qtl2 (aka qtl2) is a reimplementation of the QTL analysis software
[R/qtl](http://www.rqtl.org), to better handle high-dimensional data
and complex cross designs.  It is split
into [qtl2geno](https://github.com/rqtl/qtl2geno) (for calculating
genotype probabilities, imputations, and genetic maps) and
[qtl2scan](https://github.com/rqtl/qtl2scan) (for QTL genome scans and
related calculations).

---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You first need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
and [RcppEigen](https://github.com/RcppCore/RcppEigen).
(Additional, secondary dependencies will also be installed)

    install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))

Then, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github(c("rqtl/qtl2geno", "rqtl/qtl2scan"))

---


### Vignettes

- [user guide](assets/vignettes/user_guide.html)
- [input file formats](assets/vignettes/input_files.html)
  (see also the [sample data files](pages/sampledata.html))
- [developer guide](assets/vignettes/developer_guide.html)
- [HMM benchmarks](assets/vignettes/hmm_benchmarks.html)
- [linear regression benchmarks](assets/vignettes/linreg_benchmarks.html)

---

### License

The [qtl2geno](https://github.com/rqtl/qtl2geno) and
[qtl2scan](https://github.com/rqtl/qtl2scan) packages are free
software; you can redistribute them and/or modify them under the terms
of the GNU General Public License, version 3, as published by the Free
Software Foundation.

These programs are distributed in the hope that they will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>

---

Sources on [github](http://github.com):

- The [source for the qtl2geno](https://github.com/rqtl/qtl2geno)
- The [source for the qtl2scan](https://github.com/rqtl/qtl2scan)
- The [source for the website](https://github.com/kbroman/qtl2/tree/gh-pages)

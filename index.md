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
genotype probabilities, imputations, and genetic maps),
[qtl2scan](https://github.com/rqtl/qtl2scan) (for QTL genome scans and
related calculations), and
[qtl2plot](https://github.com/rqtl/qtl2plot) (for data visualization).
A further package, [qtl2convert](https://github.com/rqtl/qtl2convert),
contains functions for converting data among the R/qtl2,
[DOQTL](https://www.bioconductor.org/packages/release/bioc/html/DOQTL.html),
and [R/qtl](http://rqtl.org) formats, for example to convert genotype
probabilities produced by DOQTL to the format needed by qtl2scan, or
to convert qtl2scan results to the format produced by `scanone` in
R/qtl, so that they may be graphed with the R/qtl functions.

---

![Warning](assets/pics/warning_icon.png)

In R/qtl2 version 0.5, we made major revisions to some of the
central data structures, and a number of steps in QTL analyses have
changed. See the revised
[user guide](assets/vignettes/user_guide.html), or
[this description of the changes in version 0.5](assets/vignettes/version05_new.html).
A couple of functions for converting objects from the format for
Rqtl2 version 0.4 and the new format are in
[`convert_04_to_05.R`](assets/convert_04_to_05.R).

---

### Installation

R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), but
it can be installed from a mini-CRAN at [rqtl.org](http://rqtl.org).

    install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

The [qtl2](https://github.com/rqtl/qtl2) package is
inspired by the
[tidyverse package](https://cran.r-project.org/package=tidyverse);
it is basically empty, but when you install it, the
[qtl2geno](https://github.com/rqtl/qtl2geno),
[qtl2scan](https://github.com/rqtl/qtl2scan),
[qtl2plot](https://github.com/rqtl/qtl2plot), and
[qtl2convert](https://github.com/rqtl/qtl2convert) packages, plus a
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


### Documentation

- [user guide](assets/vignettes/user_guide.html)
- [input file formats](assets/vignettes/input_files.html)
  (also see the [sample data files](pages/sampledata.html) and the
  [qtl2data repository](https://github.com/rqtl/qtl2data))
- [differences between R/qtl and R/qtl2](assets/vignettes/rqtl_diff.html)
- [developer guide](assets/vignettes/developer_guide.html)
- [HMM benchmarks](assets/vignettes/hmm_benchmarks.html)
- [linear regression benchmarks](assets/vignettes/linreg_benchmarks.html)

---

### License

The [qtl2geno](https://github.com/rqtl/qtl2geno),
[qtl2scan](https://github.com/rqtl/qtl2scan),
[qtl2plot](https://github.com/rqtl/qtl2plot), and
[qtl2convert](https://github.com/rqtl/qtl2convert)
packages are free
software; you can redistribute them and/or modify them under the terms
of the GNU General Public License, version 3, as published by the Free
Software Foundation.

These programs are distributed in the hope that they will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>

---

Sources on [github](https://github.com):

- The [source for the qtl2geno](https://github.com/rqtl/qtl2geno)
- The [source for the qtl2scan](https://github.com/rqtl/qtl2scan)
- The [source for the qtl2plot](https://github.com/rqtl/qtl2plot)
- The [source for the qtl2convert](https://github.com/rqtl/qtl2convert)
- The [source for the website](https://github.com/kbroman/qtl2/tree/gh-pages)

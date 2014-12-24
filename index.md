---
layout: page
title: R/qtl2
tagline: QTL analysis for high-dimensional data and complex crosses
description: R/qtl2, a reimplementation of R/qtl to better handle high-dimensional data and complex cross designs
---

R/qtl2 (aka qtl2) is a reimplementation of the QTL analysis software
[R/qtl](http://www.rqtl.org), to better handle high-dimensional data
and complex cross designs.


---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from its
[GitHub repository](http://github.com/kbroman/qtl2). You first need to
install the [devtools](https://github.com/hadley/devtools) package.

    install.packages("devtools")

Then install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github("kbroman/qtl2")

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

The R/qtl2 package is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>

---

Sources on [github](http://github.com):

- The [source for the package](https://github.com/kbroman/qtl2/tree/master)
- The [source for the website](https://github.com/kbroman/qtl2/tree/gh-pages)

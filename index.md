---
layout: page
title: R/qtl2
tagline: QTL analysis for high-dimensional data and complex crosses
description: R/qtl2, a reimplementation of R/qtl to better handle high-dimensional data and complex cross designs
---

[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](https://rqtl.org), to better handle
high-dimensional data and complex cross designs.

---

![Warning](assets/pics/warning_icon.png)

In R/qtl2 version 0.7, we merged the multiple packages qtl2geno, qtl2scan,
qtl2plot, and qtl2db, into a single package
[qtl2](https://github.com/rqtl/qtl2). The multiple packages proved awkward and confusing.
The [qtl2convert](https://github.com/rqtl/qtl2convert) package (for
converting data among the R/qtl2,
[DOQTL](https://www.bioconductor.org/packages/release/bioc/html/DOQTL.html),
and [R/qtl](https://rqtl.org) formats) will remain a separate package.

---

### Discussion group

For discussion/questions about R/qtl2, join the
[rqtl2-disc google group](https://groups.google.com/forum/#!forum/rqtl2-disc).
Or join [rqtl-announce](https://groups.google.com/forum/#!forum/rqtl-announce)
for announcements about R/qtl and R/qtl2.
(We'll try to keep the original
[rqtl-disc group](https://groups.google.com/forum/#!forum/rqtl-disc)
for the discussion/questions about the original R/qtl only.)

---

### Installation

R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), but
it can be installed from a mini-CRAN at [rqtl.org](https://rqtl.org).
Make sure you have the latest version of [R (3.4.2)](https://cran.r-project.org).

    install.packages("qtl2", repos="https://rqtl.org/qtl2cran")

_Alternatively_, you can install R/qtl2 from its source on
[GitHub](https://github.com/rqtl). (But note that compiling the C++
code can be rather slow.)

On _Windows_, you'll need [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

On _Mac OS X_, you'll need the
[command-line developer tools](https://mac-how-to.gadgethacks.com/how-to/install-command-line-developer-tools-without-xcode-0168115/).

You then need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
[RcppEigen](https://github.com/RcppCore/RcppEigen),
[RSQLite](https://github.com/rstats-db/RSQLite), and
[qtl](https://rqtl.org).
(Additional, secondary dependencies will also be installed.)

    install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl"))

Finally, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github("rqtl/qtl2")

---


### Documentation

- [user guide](assets/vignettes/user_guide.html)
- [input file formats](assets/vignettes/input_files.html)
  (also see the [sample data files](pages/sampledata.html) and the
  [qtl2data repository](https://github.com/rqtl/qtl2data))
- [preparing DO mouse data for R/qtl2](pages/prep_do_data.html)
- [differences between R/qtl and R/qtl2](assets/vignettes/rqtl_diff.html)
- [developer guide](assets/vignettes/developer_guide.html)
- [HMM benchmarks](assets/vignettes/hmm_benchmarks.html)

---

### License

[Licensed](LICENSE) under [GPL-3](https://www.r-project.org/Licenses/GPL-3).

---

Sources on [github](https://github.com):

- The [source for the qtl2](https://github.com/rqtl/qtl2)
- The [source for the qtl2convert](https://github.com/rqtl/qtl2convert)
- The [source for the website](https://github.com/kbroman/qtl2/tree/gh-pages)

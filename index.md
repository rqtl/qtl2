---
layout: page
title: R/qtl2
tagline: QTL analysis for high-dimensional data and complex crosses
description: R/qtl2, a reimplementation of R/qtl to better handle high-dimensional data and complex cross designs
---

[R/qtl2](https://kbroman.org/qtl2) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](https://rqtl.org), to better handle
high-dimensional data and complex cross designs.

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
Make sure you have the latest version of [R](https://cran.r-project.org).

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
- [categorized list of functions in R/qtl2](pages/rqtl2_functions.html)
- [input file formats](assets/vignettes/input_files.html)
  (also see the [sample data files](pages/sampledata.html) and the
  [qtl2data repository](https://github.com/rqtl/qtl2data))
- [preparing DO mouse data for R/qtl2](pages/prep_do_data.html)
- [genotype diagnostics for diversity outbred mice](assets/vignettes/do_diagnostics.html)
- [differences between R/qtl and R/qtl2](assets/vignettes/rqtl_diff.html)
- [developer guide](assets/vignettes/developer_guide.html)
- [HMM benchmarks](assets/vignettes/hmm_benchmarks.html)
- [Tutorial on R/qtl2](https://smcclatchy.github.io/mapping/) by [Susan McClatchy](https://github.com/smcclatchy) and
  [Dan Gatti](https://github.com/dmgatti)

---

### Citation

To cite R/qtl2 in publications, use:

> Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen &#346;,
> Yandell BS, Churchill GA (2018)
> R/qtl2: software for mapping quantitative trait loci with
> high-dimensional data and multi-parent populations.
> [Genetics](http://genetics.org) 211:495-502
> [doi:10.1534/genetics.118.301595](https://doi.org/10.1534/genetics.118.301595)
> [![pdf](https://kbroman.org/pages/icons16/pdf-icon.png)](http://www.genetics.org/content/genetics/211/2/495.full.pdf)

---

### License

[Licensed](LICENSE) under [GPL-3](https://www.r-project.org/Licenses/GPL-3).

---

Sources on [github](https://github.com):

- The [source for the qtl2](https://github.com/rqtl/qtl2)
- The [source for the qtl2convert](https://github.com/rqtl/qtl2convert)
- The [source for the website](https://github.com/kbroman/qtl2/tree/gh-pages)

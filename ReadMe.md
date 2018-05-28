### [R/qtl2](https://kbroman.org/qtl2)

[![Build Status](https://travis-ci.org/rqtl/qtl2.svg?branch=master)](https://travis-ci.org/rqtl/qtl2)

[Karl Broman](https://kbroman.org)

[R/qtl2](https://kbroman.org/qtl2) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](https://rqtl.org), to better handle
high-dimensional data and complex cross designs.

It was previously split into the multiple packages qtl2geno, qtl2scan,
qtl2plot, and qtl2db, but that proved awkward and confusing and so the
packages have been combined. But also see the related package,
[qtl2convert](https://github.com/rqtl/qtl2convert), for converting
data among the R/qtl2,
[DOQTL](https://www.bioconductor.org/packages/release/bioc/html/DOQTL.html),
and [R/qtl](https://rqtl.org) formats.
[qtl2convert](https://github.com/rqtl/qtl2convert) will remain a
separate package.

---

### Installation

R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), but
it can be installed from a mini-CRAN at [rqtl.org](https://rqtl.org).
Make sure you have the latest version of [R (3.4.4)](https://cran.r-project.org).

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

### Vignettes

- [user guide](https://kbroman.org/qtl2/assets/vignettes/user_guide.html)
- [input file formats](https://kbroman.org/qtl2/assets/vignettes/input_files.html)
  (see also the
  [sample data files](https://kbroman.org/qtl2/pages/sampledata.html)
  and the [qtl2data repository](https://github.com/rqtl/qtl2data))
- [preparing DO mouse data for R/qtl2](https://kbroman.org/qtl2/pages/prep_do_data.html)
- [differences between R/qtl and R/qtl2](https://kbroman.org/qtl2/assets/vignettes/rqtl_diff.html)
- [developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html)
- [HMM benchmarks](https://kbroman.org/qtl2/assets/vignettes/hmm_benchmarks.html)

---

### License

[Licensed](LICENSE) under [GPL-3](https://www.r-project.org/Licenses/GPL-3).

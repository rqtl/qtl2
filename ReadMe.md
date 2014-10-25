### R/qtl2

[![Build Status](https://travis-ci.org/kbroman/qtl2.png?branch=master)](https://travis-ci.org/kbroman/qtl2)

[Karl W Broman](http://kbroman.org)

[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](http://www.rqtl.org), to better
handle high-dimensional data and complex cross designs.

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

- [user guide](http://kbroman.org/qtl2/assets/vignettes/user_guide.html)
- [input file formats](http://kbroman.org/qtl2/assets/vignettes/input_files.html)
  (see also the [sample data files](http://kbroman.org/qtl2/pages/sampledata.html))
- [developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html)
- [HMM benchmarks](http://kbroman.org/qtl2/assets/vignettes/hmm_benchmarks.html)
- [linear regression benchmarks](http://kbroman.org/qtl2/assets/vignettes/linreg_benchmarks.html)

---

#### License

[Licensed](LICENSE) under [GPL-3](http://www.r-project.org/Licenses/GPL-3).

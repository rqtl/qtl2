### [R/qtl2](https://kbroman.org/qtl2/)

[![R-CMD-check](https://github.com/rqtl/qtl2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rqtl/qtl2/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/qtl2)](https://cran.r-project.org/package=qtl2)

[Karl Broman](https://kbroman.org)

[R/qtl2](https://kbroman.org/qtl2/) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](https://rqtl.org), to better handle
high-dimensional data and complex cross designs.

Also see the related packages,
[qtl2convert](https://github.com/rqtl/qtl2convert) (for converting
data among the R/qtl2, DOQTL, and R/qtl formats),
[qtl2fst](https://github.com/rqtl/qtl2fst) (for storing genotype
probabilities on disk), and [qtl2ggplot](https://github.com/byandell/qtl2ggplot)
(for [ggplot2](https://ggplot2.tidyverse.org/)-based data visualizations).

---

### Installation

Install R/qtl2 from [CRAN](https://cran.r-project.org):

    install.packages("qtl2")

---

### Documentation

- [user guide](https://kbroman.org/qtl2/assets/vignettes/user_guide.html)
- [categorized list of functions in R/qtl2](https://kbroman.org/qtl2/pages/rqtl2_functions.html)
- [input file formats](https://kbroman.org/qtl2/assets/vignettes/input_files.html)
  (see also the
  [sample data files](https://kbroman.org/qtl2/pages/sampledata.html)
  and the [qtl2data repository](https://github.com/rqtl/qtl2data))
- [using qtl2fst for on-disk genotype probabilities](https://kbroman.org/qtl2/assets/vignettes/qtl2fst.html)
- [preparing diversity outbred mouse data for R/qtl2](https://kbroman.org/qtl2/pages/prep_do_data.html)
- [genotype diagnostics for diversity outbred mice](https://kbroman.org/qtl2/assets/vignettes/do_diagnostics.html)
- [identifying sample mix-ups in diversity outbred mice](https://kbroman.org/qtl2/assets/vignettes/do_mixups.html)
- [differences between R/qtl and R/qtl2](https://kbroman.org/qtl2/assets/vignettes/rqtl_diff.html)
- [developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html)
- [HMM benchmarks](https://kbroman.org/qtl2/assets/vignettes/hmm_benchmarks.html)
- [Tutorial on R/qtl2](https://smcclatchy.github.io/mapping/) by [Susan McClatchy](https://github.com/smcclatchy) and
  [Dan Gatti](https://github.com/dmgatti)

---

### Citation

To cite R/qtl2 in publications, use:

> Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen &#346;,
> Yandell BS, Churchill GA (2019)
> R/qtl2: software for mapping quantitative trait loci with
> high-dimensional data and multi-parent populations.
> [Genetics](https://academic.oup.com/genetics) 211:495-502
> [doi:10.1534/genetics.118.301595](https://doi.org/10.1534/genetics.118.301595)

---

### License

Licensed under [GPL-3](https://www.r-project.org/Licenses/GPL-3).

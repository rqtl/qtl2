
all: vignettes
.PHONY: vignettes

VIGNETTES = assets/vignettes/linreg_benchmarks.html
vignettes: ${VIGNETTES}

assets/vignettes/linreg_benchmarks.html: assets/vignettes/linreg_benchmarks.Rmd
	cd $(<D);R -e 'library(knitr);knit2html("$(<F)")'


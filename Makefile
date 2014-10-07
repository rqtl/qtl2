
all: vignettes data
.PHONY: vignettes data

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

VIGNETTES = assets/vignettes/linreg_benchmarks.html
vignettes: ${VIGNETTES}

assets/vignettes/linreg_benchmarks.html: assets/vignettes/linreg_benchmarks.Rmd
	cd $(<D);R -e 'library(knitr);knit2html("$(<F)")'

data: assets/sampledata/grav2/grav2.yaml assets/sampledata/iron/iron.yaml

assets/sampledata/grav2/grav2.yaml: assets/sampledata/scripts/grav2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} grav2cross2.R

assets/sampledata/iron/iron.yaml: assets/sampledata/scripts/iron2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} iron2cross2.R

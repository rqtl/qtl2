
all: vignettes data external_vignettes
.PHONY: vignettes data external_vignettes

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

VIGNETTES = assets/vignettes/linreg_benchmarks.html
vignettes: ${VIGNETTES}

EXTERNAL_VIGNETTES = assets/vignettes/developer.html assets/vignettes/input_files.html
external_vignettes: ${EXTERNAL_VIGNETTES}

assets/vignettes/linreg_benchmarks.html: assets/vignettes/linreg_benchmarks.Rmd
	cd $(<D);R -e 'library(knitr);knit2html("$(<F)")'
	rm $(@D)/linreg_benchmarks.md

data: assets/sampledata/grav2/grav2.yaml assets/sampledata/iron/iron.yaml

assets/sampledata/grav2/grav2.yaml: assets/sampledata/scripts/grav2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} grav2cross2.R

assets/sampledata/iron/iron.yaml: assets/sampledata/scripts/iron2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} iron2cross2.R

assets/vignettes/developer.html: ../qtl2/vignettes/developer.Rmd
	cd $(@D);R -e 'library(knitr);knit2html("../../$<")'
	rm $(@D)/developer.md

assets/vignettes/input_files.html: ../qtl2/vignettes/input_files.Rmd
	cd $(@D);R -e 'library(knitr);knit2html("../../$<")'
	rm $(@D)/input_files.md


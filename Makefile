all: doc vignettes data
.PHONY: doc vignettes data

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

# build package documentation
doc:
	R -e 'library(devtools);document()'

vignettes: inst/doc/developer.html

inst/doc/%.html: vignettes/%.Rmd
	cd $(@D);R ${R_OPTS} -e 'library(knitr);knit2html("../../$<", "$(@F)")'
#	rm $(@D)/$*.md #<- if .md file created, might want to delete it

data: inst/sampledata/grav2/grav2.yaml inst/sampledata/iron/iron.yaml

inst/sampledata/grav2/grav2.yaml: inst/scripts/grav2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} grav2cross2.R

inst/sampledata/iron/iron.yaml: inst/scripts/iron2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} iron2cross2.R

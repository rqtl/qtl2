all: doc vignettes
.PHONY: doc vignettes

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

# build package documentation
doc:
	R -e 'library(devtools);document()'

vignettes: inst/doc/developer_guide.html inst/doc/input_files.html inst/doc/user_guide.html

inst/doc/%.html: vignettes/%.Rmd
	cd $(@D);R ${R_OPTS} -e 'library(knitr);knit2html("../../$<", "$(@F)")'

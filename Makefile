all: vignettes data external_vignettes extdata
.PHONY: vignettes data external_vignettes extdata

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

VIGNETTES = assets/vignettes/linreg_benchmarks.html assets/vignettes/hmm_benchmarks.html
vignettes: ${VIGNETTES}

EXTERNAL_VIGNETTES = assets/vignettes/developer_guide.html assets/vignettes/input_files.html assets/vignettes/user_guide.html
external_vignettes: ${EXTERNAL_VIGNETTES}

assets/vignettes/linreg_benchmarks.html: assets/vignettes/linreg_benchmarks.Rmd
	cd $(<D);R $(R_OPTS) -e 'rmarkdown::render("$(<F)")'

assets/vignettes/hmm_benchmarks.html: assets/vignettes/hmm_benchmarks.Rmd
	cd $(<D);R $(R_OPTS) -e 'rmarkdown::render("$(<F)")'

data: assets/sampledata/grav2/grav2.yaml assets/sampledata/iron/iron.yaml

assets/sampledata/grav2/grav2.yaml: assets/sampledata/scripts/grav2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} grav2cross2.R

assets/sampledata/iron/iron.yaml: assets/sampledata/scripts/iron2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} iron2cross2.R

assets/vignettes/developer_guide.html: ../qtl2geno/vignettes/developer_guide.Rmd
	cd $(<D); \
	R $(R_OPTS) -e 'rmarkdown::render("$(<F)")'; \
	mv $(@F) ../../Web/$(@D)

assets/vignettes/input_files.html: ../qtl2geno/vignettes/input_files.Rmd
	cd $(<D); \
	R $(R_OPTS) -e 'rmarkdown::render("$(<F)")'; \
	mv $(@F) ../../Web/$(@D)

assets/vignettes/user_guide.html: ../qtl2geno/vignettes/user_guide.Rmd
	cd $(<D); \
	R $(R_OPTS) -e 'rmarkdown::render("$(<F)")'; \
	mv $(@F) ../../Web/$(@D)

EXTDATA = ../qtl2geno/inst/extdata/grav2.zip ../qtl2geno/inst/extdata/iron.zip
extdata: ${EXTDATA}

../qtl2geno/inst/extdata/grav2.zip: assets/sampledata/grav2/grav2.yaml
	cp $(<D)/$(@F) $@

../qtl2geno/inst/extdata/iron.zip: assets/sampledata/iron/iron.yaml
	cp $(<D)/$(@F) $@

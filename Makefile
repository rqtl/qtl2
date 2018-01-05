all: vignettes data external_vignettes extdata
.PHONY: vignettes data external_vignettes extdata

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

VIGNETTES = assets/vignettes/hmm_benchmarks.html assets/vignettes/rqtl_diff.html assets/vignettes/input_files.html assets/vignettes/developer_guide.html assets/vignettes/user_guide.html
vignettes: ${VIGNETTES}

assets/vignettes/%.html: assets/vignettes/%.Rmd ruby/add_navbar.rb ruby/vignette_head.html ruby/vignette_navbar.html
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"
	ruby/add_navbar.rb $@

data: assets/sampledata/grav2/grav2.yaml assets/sampledata/iron/iron.yaml

assets/sampledata/grav2/grav2.yaml: assets/sampledata/scripts/grav2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} grav2cross2.R

assets/sampledata/iron/iron.yaml: assets/sampledata/scripts/iron2cross2.R
	cd $(<D);R CMD BATCH ${R_OPTS} iron2cross2.R

EXTDATA = ../qtl2/inst/extdata/grav2.zip ../qtl2/inst/extdata/iron.zip
extdata: ${EXTDATA}

../qtl2/inst/extdata/grav2.zip: assets/sampledata/grav2/grav2.yaml
	cp $(<D)/$(@F) $@

../qtl2/inst/extdata/iron.zip: assets/sampledata/iron/iron.yaml
	cp $(<D)/$(@F) $@

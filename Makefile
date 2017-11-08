.PHONY: doc

# build package documentation
doc:
	R -e 'devtools::document()'

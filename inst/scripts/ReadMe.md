## Scripts for creating databases of CC variants and mouse genes

This directory contains four R scripts, for creating
[SQLite](https://www.sqlite.org) databases for use with the
Collaborative Cross (CC) and Diversity Outbred (DO) mouse populations:

- `create_ccvariants.R` creates the database `cc_variants.sqlite` which
  combines data on SNPs, indels, and structural variants in the eight
  CC founders.

- `create_mousegenes.R` creates the database `mouse_genes.sqlite`
  which mouse gene locations from the
  [Mouse Genome Informatics (MGI)](http://www.informatics.jax.org/)
  database. Key fields in the resulting database include

  - `source` (e.g., `"MGI"`)
  - `type`
  - `start` (in basepairs)
  - `stop` (in basepairs)
  - `strand`
  - `Name`

- `create_ccvariants_small.R` creates
  `../extdata/cc_variants_small.sqlite`, a small version of
  `cc_variants.sqlite` for use in tests. This contains the variants within
  two small regions (one on chr 2 and one on chr 3).

- `create_mousegenes_small.R` creates
  `../extdata/mouse_genes_small.sqlite`, a small version of
  `mouse_genes.sqlite` for use in tests. This contains the genes with
  `source=="MGI"` that overlap two small regions (one on chr 2 and one
  on chr 3).

- `create_mousegenes_mgi.R` creates `mouse_genes_mgi.sqlite`, a version
  with just the records with `source=="MGI"`.

Using these scripts to constructing these databases requires the
following R packages:

- [data.table](https://github.com/Rdatatable/data.table/wiki)
- [RSQLite](https://rstats-db.github.io/RSQLite/)
- [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
- [R/qtl2scan](https://github.com/rqtl/qtl2scan)

The larger database files created by these scripts are available for
direct download from [figshare](https://figshare.com):

- [`cc_variants.sqlite` doi:10.6084/m9.figshare.5280229.v1](https://doi.org/10.6084/m9.figshare.5280229.v1)
- [`mouse_genes.sqlite` doi:10.6084/m9.figshare.5280238.v3](https://doi.org/10.6084/m9.figshare.5280238.v3)

A smaller version of the mouse genes database, with just the records
with `source=="MGI"`, is also available at
[figshare](https://figshare.com):

- [`mouse_genes_mgi.sqlite` doi:10.6084/m9.figshare.5286019.v4](https://doi.org/10.6084/m9.figshare.5286019.v4)

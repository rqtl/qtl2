## qtl2db scripts

This directory contains four R scripts, for creating
[SQLite](https://www.sqlite.org) databases for use with the
Collaborative Cross (CC) and Diversity Outbred (DO) mouse populations:

- `create_ccsnps.R` creates the database `cc_snps.sqlite` which
  combines data on SNPs, indels, and structural variants in the eight
  CC founders.

- `create_mousegenes.R` creates the database `mouse_genes.sqlite`
  which mouse gene locations from the
  [Mouse Genome Informatics (MGI)](http://www.informatics.jax.org/)
  database.

- `create_ccsnps_small.R` creates
  `../extdata/cc_snps_small.sqlite`, a small version of
  `cc_snps.sqlite` for use in tests. This contains the variants within
  two small regions (one on chr 2 and one on chr 3).

- `create_mousegenes_small.R` creates
  `../extdata/mouse_genes_small.sqlite`, a small version of
  `mouse_genes.sqlite` for use in tests. This contains the genes with
  `source=="MGI"` that overlap two small regions (one on chr 2 and one
  on chr 3).

Using these scripts to constructing these databases requires the
following R packages:

- data.table
- RSQLite
- ...

The larger database files created by these scripts are available for
direct download from [figshare](https://figshare.com):

- [`cc_snps.sqlite`]()
- [`mouse_genes.sqlite`]()

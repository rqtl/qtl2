## Sample input data files (zipped) + small example variant and gene databases

### Sample input files

- [`grav2.zip`](grav2.zip), data from
  [Moore et al. (2013) Genetics 195:1077-1086](http://www.genetics.org/content/195/3/1077.abstract)
  (the second replicate of the RILs)

- [`iron.zip`](iron.zip), data from
  [Grant et al. (2006) Hepatology 44:174-185](https://www.ncbi.nlm.nih.gov/pubmed/16799992)

See the contents at <https://kbroman.org/qtl2/pages/sampledata.html>.

These can be read into R with `qtl2::read_cross2()`.

### Variant and gene databases

- [`cc_variants_small.sqlite`](cc_variants_small.sqlite) - SQLite
  database of variants in the Collaborative Cross founders in two
  small regions

- [`mouse_genes_small.sqlite`](mouse_genes_small.sqlite) - SQLite
  database of mosue genes in two small regions

The larger versions of these database files are available for
direct download from [figshare](https://figshare.com):

- [`cc_variants.sqlite` doi:10.6084/m9.figshare.5280229.v1](https://doi.org/10.6084/m9.figshare.5280229.v2)
- [`mouse_genes.sqlite` doi:10.6084/m9.figshare.5280238.v6](https://doi.org/10.6084/m9.figshare.5280238.v6)

A smaller version of the mouse genes database, with just the records
with `source=="MGI"`, is also available at
[figshare](https://figshare.com):

- [`mouse_genes_mgi.sqlite` doi:10.6084/m9.figshare.5286019.v7](https://doi.org/10.6084/m9.figshare.5286019.v7)

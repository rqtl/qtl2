---
layout: page
title: Resources for mouse genome build 39
description: "Resources for moving to mouse genome build GRCm39 (aka mm11), with a particular focus on mouse genetics projects related to the Collaborative Cross and Diversity Outbred"
---


The mouse genome build [GRCm39 was released on
2020-06-24](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/).
Here we describe available resources for using build 39 with R/qtl2,
and how to move your projects from build 38 to build 39.

We've focused mostly on projects related to the Collaborative Cross or
Diversity Outbred (DO) mice, and on the use of the [MUGA SNP
arrays](https://doi.org/10.1534/g3.115.022087), such as the MegaMUGA
and
[GigaMUGA](https://www.neogen.com/categories/genotyping-arrays/gigamuga/)
arrays.

- The mouse genetic map of [Cox et al
  (2009)](https://doi.org/10.1534/genetics.109.105486) has been
  revised for the new genome build. The main changes concern
  inversions at the centromeres of chromosomes 10 and 14. See [the
  github repository CoxMapV3](https://github.com/kbroman/CoxMapV3).

- The MUGA array annotations have been revised to use positions from
  the GRCm39 genome build. See the [github repository
  MUGAarrays](https://github.com/kbroman/MUGAarrays).

- The data files with founder genotypes, organized for use with data
  in R/qtl2 format, have been revised for build GRCm39.

  - [`MM_processed_files_build39.zip`](https://doi.org/10.6084/m9.figshare.22666336.v1)
  - [`GM_processed_files_build39.zip`](https://doi.org/10.6084/m9.figshare.22666504.v1)
  - [`MMnGM_processed_files_build39.zip`](https://doi.org/10.6084/m9.figshare.22666510.v1)

- We've created an R package, [mmconvert](https://github.com/rqtl/mmconvert), which is
  [available on CRAN](https://cran.r-project.org/package=mmconvert),
  to convert map positions between build GRCm39 and the
  revised Cox genetic map. It is intended to serve the role of the
  "mouse map converter" web service from Gary Churchill's group, which
  is no longer available.

  The mmconvert package includes a function
  `cross2_to_grcm39()` for converting a `cross2` object (created with
  `read_cross2()`, and for a cross using the MegaMUGA and/or GigaMUGA
  arrays) to build GRCm39.

  Here is an example of its use with
  [DO data from Karen Svenson and colleagues](https://github.com/rqtl/qtl2data/tree/main/DO_Svenson291)

  ```r
  library(qtl2)
  file <- paste0("https://raw.githubusercontent.com/rqtl/",
                 "qtl2data/main/DO_Svenson291/svenson.zip")
  do <- read_cross2(file, quiet=FALSE)

  library(mmconvert)
  do_grcm39 <- cross2_to_grcm39(do)
  ```

- A new SQLite database with CC/DO founder variants from Sanger along
  with ensembl genes is [available on
  figshare](https://doi.org/10.6084/m9.figshare.22630642.v1). (Created
  by [Matt Vincent](https://www.jax.org/people/matt-vincent) at the
  Jackson Lab.)

  Download this ~10.2 GB database as `fv.2021.snps.db3`:

  ```r
  download.file("https://figshare.com/ndownloader/files/40157572", "fv.2021.snps.db3")
  ```

  You can then use `create_variant_query_func()` as before, though you
  need to use the `id_field` argument, as follows:

  ```r
  qvf <- create_variant_query_func("fv.2021.snps.db3", id_field="variants_id")
  ```

  The genes database uses different names for several fields, and so use
  `create_gene_query_func()` as follows:

  ```r
  qgf <- create_gene_query_func("fv.2021.snps.db3", chr_field="chromosome", name_field="symbol",
                                start_field="start_position", stop_field="end_position")
  ```

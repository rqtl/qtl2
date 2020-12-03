---
layout: page
title: Preparing Diversity Outbred (DO) mouse data for R/qtl2
description: "How to prepare Diversity Outbred (DO) mouse data, including SNP genotypes from the MegaMUGA/GigaMUGA arrays, for use with R/qtl2."
---

Studies using [Diversity Outbred (DO) mice](https://www.jax.org/strain/009376)
often involve SNP genotype data from the
[MegaMUGA or GigaMUGA arrays](https://csbio.unc.edu/CCstatus/index.py?run=Genotype)
from [GeneSeek](https://genomics.neogen.com/en/mouse-universal-genotyping-array).
To analyze such data with [R/qtl2](https://kbroman.org/qtl2), you first
need to do some data reorganization and pre-processing. In this
document, we will attempt to explain what to do.

### Founder genotypes and SNP maps

With DO data, an important component of the
[R/qtl2 input files](../assets/vignettes/input_files.html) is the SNP
genotypes for the eight founder lines. We also need the genetic and
physical positions of the SNP markers.

[Dan Gatti](https://www.jax.org/research-and-faculty/faculty/research-scientists/daniel-gatti)
has prepared such data for the MegaMUGA and GigaMUGA arrays and made
them available at [ftp.jax.org/MUGA](ftp://ftp.jax.org/MUGA). The files with names
starting `MM_` are for the MegaMUGA array; the files with names
starting `GM_` are for the GigaMUGA array.

We've also placed these files at
[figshare](https://figshare.com):
[doi:10.6084/m9.figshare.c.3879694](https://doi.org/10.6084/m9.figshare.c.3879694).

- [`MM_primary_files.zip` (doi:10.6084/m9.figshare.5404717.v1)](https://doi.org/10.6084/m9.figshare.5404717.v1),
  containing all of the files for the MegaMUGA array
- [`GM_primary_files.zip` (doi:10.6084/m9.figshare.5404675.v1)](https://doi.org/10.6084/m9.figshare.5404675.v1),
  containing all of the files for the GigaMUGA array

The files at [ftp.jax.org/MUGA](ftp://ftp.jax.org/MUGA) and in these
zip files are all in `.Rdata` format, for
[R](https://www.r-project.org). For the GigaMUGA data, the key files
include `GM_snps.Rdata`, containing SNP locations, and
`GM_geno.Rdata`, containing SNP genotypes for the Collaborative Cross
(CC) founder strains (which are also the founders of the DO), as
single-letter codes `N`/`C`/`T`/`A`/`G`/`H`.

To use these data with R/qtl2, we need to reformat the files and
recode the SNP genotypes (for example, as `A`/`H`/`B`/`-`). We've
placed an R script to do this at the figshare site,
[`parse_muga.R` (doi:10.6084/m9.figshare.5405260)](https://doi.org/10.6084/m9.figshare.5405260).
This uses the [R/qtl2convert](https://github.com/rqtl/qtl2convert)
package, in particular `find_consensus_geno()` to determine consensus
genotypes for each founder strain at each SNP from the multiple
individuals that were genotyped, `find_unique_geno()` to determine the
SNP alleles, and `encode_geno()` to encode the SNP genotypes as
`A`/`H`/`B`/`-`.

The processed files are also at figshare:

- [`MM_processed_files.zip` (doi:10.6084/m9.figshare.5404750)](https://doi.org/10.6084/m9.figshare.5404750)
- [`GM_processed_files.zip` (doi:10.6084/m9.figshare.5404759)](https://doi.org/10.6084/m9.figshare.5404759)
- [`MMnGM_processed_files.zip` (doi:10.6084/m9.figshare.5404762)](https://doi.org/10.6084/m9.figshare.5404762)

The `MM_` and `GM_` files are for the MegaMUGA and GigaMUGA arrays,
respectively. The `MMnGM_` file (produced with the R script
[`combine_MMnGM.R` (doi:10.6084/m9.figshare.5405281.v1)](https://doi.org/10.6084/m9.figshare.5405281.v1))
merges the SNPs for the two arrays, for DO projects that have some
mice genotyped with the MegaMUGA array and others genotyped with the
GigaMUGA array.

Each of these zip files contain a series of CSV files. For example,
for the GigaMUGA data, we have:

- `GM_allelecodes.csv`, the allele encodings for each SNP
- `GM_foundergeno*.csv`, the CC founder genotypes, encoded as
  `A`/`H`/`B`/`-`, with one file for each chromosome
- `GM_gmap*.csv`, the genetic map locations for the SNPS (in cM), with
  one file for each chromosome
- `GM_pmap*.csv`, the physical map locations for the SNPS (in Mbp), with
  one file for each chromosome
- `GM_info.csv`, information for all SNPs, including the genetic and
  physical locations and `tier` (an indication of SNP quality)

If your DO data includes just the GigaMUGA array, you'll
need the
[`GM_processed_files.zip`](https://doi.org/10.6084/m9.figshare.5404759) file.
If you've used just the MegaMUGA array, you'll need the
[`MM_processed_files.zip`](https://doi.org/10.6084/m9.figshare.5404750) file.
If your project includes some MegaMUGA and some GigaMUGA arrays,
you'll need the
[`MMnGM_processed_files.zip`](https://doi.org/10.6084/m9.figshare.5404762) file.

Download the appropriate file and unzip it somewhere.



### Encoding the DO genotypes

Now to the actual DO genotypes. GeneSeek will provide a zip file that
contains a series of files. We'll use only the `FinalReport.txt` file,
which is the biggest of them, with one line for each SNP and for each
sample. For 200 mice, it'll be about 1 GB, and maybe a quarter of that
size when compressed as a zip file.

We provide an example R script,
[`geneseek2qtl2.R`](../assets/geneseek2qtl2.R), to convert this file
into what's needed for R/qtl2. It does two things:

- Encodes the genotypes as `A/H/B/-` and writes them to a series of
  CSV files, one per chromosome
- Extracts the SNP array intensities for the SNPs on the X and Y
  chromosome, which we find useful for verifying the sex of the mice.

To use this script, you'll need to edit three lines near the top:

- `codefile` defines the path to the `GM_allelecodes.csv` file (or
  `MM_` or `MMnGM_` version, if you're using MegaMUGA or both arrays
  in your project)

- `ifiles` defines the path to your `FinalReport.txt` file. This can
  be a vector of such file paths, if your genotyping was performed in
  batches

- `ostem` defines the path and file "stem" for the output files. For
  example, if you use `ostem <- "Data/forqtl2"`, the output files will
  be placed in the `Data` subdirectory and will have names like
  `forqtl2_geno1.csv`.

One potential issue to contend with, in reading in the genotype files
from GeneSeek, is that these files may have varying headers at the top
that you may need to skip over.

Another issue regards possible
recoding of the sample identifiers. For example, you may have one batch
where the samples are labeled like `DO-146` and another where they're
labeled simply `146` and another where they're labeled `AA-DO-146`.
Search for the following line, and do some reformatting of the sample
IDs there.

```
# Note: may need to revise the IDs in the second column
```

You'll need to have the
[R/qtl2convert](https://github.com/rqtl/qtl2convert) package
installed. And note that, since this script reads in the full data
into memory, you'll need a computer with appreciable RAM.


The result of the [`geneseek2qtl2.R`](../assets/geneseek2qtl2.R)
script will be a series of CSV files containing the re-coded
genotypes, with one file per chromosome, plus two files containing
the SNP array intensities for the X and Y chromosomes (one file per
allele per chromosome).



### Preparing the phenotype and covariate data

You next need to prepare files with the phenotype data (with strictly
numeric phenotypes) and covariate data (which can include non-numeric
variables). Each should be arranged with samples as rows and variables
as columns, and with the first column being the sample IDs, exactly as
used in the genotype data.

For examples of these files, see the
[qtl2data repository](https://github.com/rqtl/qtl2data), in particular
[`do_pheno.csv`](https://github.com/rqtl/qtl2data/blob/master/DO_Gatti2014/do_pheno.csv)
and
[`do_covar.csv`](https://github.com/rqtl/qtl2data/blob/master/DO_Gatti2014/do_covar.csv),
for the [Gatti et al. (2014)](https://www.g3journal.org/content/ggg/4/9/1623.full.pdf) data.

The covariate data needs to contain a `sex` column, plus a column
giving the generation number (e.g., `ngen`) for each DO animal.

For help in determining the generation numbers, see
the file
[`DO_generation_dates.csv`](../assets/DO_generation_dates.csv) which
lists the generation number for different dates of distribution of
mice from the Jackson Lab.


### Preparing the control file

Finally, you need to prepare a control file which includes details of
the file names and the various encodings and variable names.

We'll use `write_control_file()` to make the control
file. It's a bit complicated, because there's a lot of stuff to
specify, and some of the pieces are sort of confusing.

```r
library(qtl2)
chr <- c(1:19, "X")
write_control_file("forqtl2.json",
                   crosstype="do",
                   description="My awesome DO project",
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("forqtl2_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="forqtl2_pheno.csv",
                   covar_file="forqtl2_covar.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen")
```

The first three arguments are the name of the file to be created (with
extension either `.json` or `.yaml`), the cross type (`"do"` for DO
mice), and a description.

The argument `founder_geno_file` gives the name of the file that contains the founder
genotypes, here a vector of file names since we have the genotypes
split by chromosome. Then `founder_geno_transposed=TRUE` indicates that we
have the founder strains as columns and the markers as rows, rather
than the default orientation with strains as the rows.

The arguments `gmap_file` and `pmap_file` indicate the names of the
files that contain the genetic (`gmap`) and physical (`pmap`) marker
positions, respectively. These are again each a vector of file names, since the
files are split by chromosome.

The argument `geno_file` indicates the name of the file with the DO genotype data (again, a
vector of file names as they're split by chromosome), and we use
`geno_tranposed=TRUE` because the
[`geneseek2qtl2.R`](../assets/geneseek2qtl2.R) script discussed above
writes the CSV files with individuals as columns and markers as rows.

The argument `geno_codes` indicates the encodings for both the DO genotypes
and the founder genotypes. This may be a bit confusing, but we want to
make a list with the names being the codes used in the data and the
values being the integers 1, 2, and 3. So `geno_codes=list(A=1, H=2,
B=3)`.

The argument `xchr` indicates the chromosome identifier for the X
chromosome.

The `pheno_file` and `covar_file` arguments indicate the names of the
files with the phenotype and covariate data, respectively. We then
provide some details about specific covariates that are important: sex
(via `sex_covar`) plus `sex_codes` to explain the encodings, with
the names `F` and `M` being the codes used in the data and the values
`"Female"` and `"Male"` indicating which sexes the codes correspond
to. And finally `crossinfo_covar` indicates the name of the covariate
column that contains the cross information (`"ngen"` for number of
generations, for cross type `"do"`).

It's sometimes useful to create chromosome-specific control
files, for the case that you want to just load one chromosome worth of
data. You can do this by first defining `chr` to be a single value,
e.g. `chr <- "1"`. You can then use exactly the same
`write_control_file()` code, though perhaps changing the first line,
with the name of the file to be produced. For example:

```r
library(qtl2)
chr <- "1"
write_control_file("forqtl2_chr1.json",
                   crosstype="do",
                   description="My awesome DO project",
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("forqtl2_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="forqtl2_pheno.csv",
                   covar_file="forqtl2_covar.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen")
```

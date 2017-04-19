# reformat data from Grant et al. (2006) Hepatology 44:174-185
# in R/qtl2 file format
#
# Abstract of paper at PubMed: https://www.ncbi.nlm.nih.gov/pubmed/16799992
# Data taken from R/qtlbook package, https://github.com/kbroman/qtlbook

library(qtl)
library(qtlbook)
data(iron)

alleles <- attr(iron, "alleles")

odir <- "../iron"

# grab genotypes
g <- pull.geno(iron, chr="-X")
gX <- pull.geno(iron, chr="X")
# convert X genotypes
gX <- reviseXdata("f2", "simple", getsex(iron), geno=gX,
                  cross.attr=attributes(iron))
g <- cbind(g, gX)

# recode genotypes
storage.mode(g) <- "character"
g[is.na(g)] <- "-"
g[g=="1"] <- paste(rep(alleles[1], 2), collapse="")
g[g=="2"] <- paste(alleles[1:2], collapse="")
g[g=="3"] <- paste(rep(alleles[2], 2), collapse="")

g <- cbind(id=as.character(1:nrow(g)), g)
# write genotypes
write.table(g, file.path(odir, "iron_geno.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# write genetic map
map <- pull.map(iron, as.table=TRUE)
map <- cbind(marker=rownames(map), map)
map <- apply(map, 2, function(a) gsub(" ", "", as.character(a)))
write.table(map, file.path(odir, "iron_gmap.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# write phenotypes
phe <- as.matrix(iron$pheno[,1:2])
storage.mode(phe) <- "character"
phe <- cbind(as.character(1:nrow(phe)), phe)
colnames(phe)[1] <- "id"
write.table(phe, file.path(odir, "iron_pheno.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# covariates
covar <- iron$pheno[,3:4]
covar <- sapply(covar, as.character)
# convert pgm to "direction", with more complicated code
covar[covar[,2]=="0",2] <- "(SxB)x(SxB)"
covar[covar[,2]=="1",2] <- "(BxS)x(BxS)"
colnames(covar)[2] <- "cross_direction"
covar <- cbind(as.character(1:nrow(covar)), covar)
colnames(covar)[1] <- "id"
write.table(covar, file.path(odir, "iron_covar.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# phenotype covariates
phecovar <- cbind(pheno=colnames(phe)[-1], tissue=colnames(phe)[-1])
write.table(phecovar, file.path(odir, "iron_phenocovar.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# control info
genotypes <- list(SS=1L, SB=2L, BB=3L)
iron_info <- list(crosstype = "f2",
                  geno = "iron_geno.csv",
                  pheno = "iron_pheno.csv",
                  phenocovar = "iron_phenocovar.csv",
                  covar = "iron_covar.csv",
                  gmap = "iron_gmap.csv",
                  pmap = "iron_pmap.csv",
                  alleles = alleles,
                  genotypes = genotypes,
                  sex = list(covar="sex", f='female', m='male'),
                  cross_info=list(covar="cross_direction", "(SxB)x(SxB)"=0L, "(BxS)x(BxS)"=1L),
                  x_chr = "X",
                  na.strings = c("-", "NA"))

library(yaml)
yaml_file <- file.path(odir, "iron.yaml")
cat("# Data from Grant et al. (2006) Hepatology 44:174-185",
    "# Abstract of paper at PubMed: https://www.ncbi.nlm.nih.gov/pubmed/16799992",
    "# Available as part of R/qtl book package, https://github.com/kbroman/qtlbook",
    file=yaml_file, sep="\n")
cat(as.yaml(iron_info), file=yaml_file, append=TRUE)

# create a version as a zip file
library(qtl2geno)
zip_datafiles(yaml_file)

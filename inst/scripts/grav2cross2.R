# reformat data from Moore et al. (2003) Genetics 195:1077-1086
# in R/qtl2 file format (this is the second replicate)
#
# See QTL Archive, http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013c

url <- "http://qtlarchive.org/grpdoc/moore_2013c/RIL2_GraviInput.csv"
file <- "grav2.csv"
utils::download.file(url, file)

library(qtl)

suppressWarnings(
    grav2 <- read.cross("csv", "", file, na.strings="*", crosstype="riself") )

# remove the file
unlink(file)

# omit blanks from the markernames
for(i in chrnames(grav2)) {
    mn <- names(grav2$geno[[i]]$map)
    mn <- gsub("\\s$", "", gsub("^\\s", "", mn))
    names(grav2$geno[[i]]$map) <- colnames(grav2$geno[[i]]$data) <- mn
}
# avoid superposed markers
grav2 <- jittermap(grav2)

alleles <- c("L", "C")

odir <- "../sampledata/"

# write genotypes
g <- pull.geno(grav2)
storage.mode(g) <- "character"
g[is.na(g)] <- "-"
g[g=="1"] <- paste(alleles[1], collapse="")
g[g=="2"] <- paste(alleles[2], collapse="")
g <- cbind(id=as.character(1:nrow(g)), g)
write.table(g, paste0(odir, "grav2_geno.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# write genetic map
map <- pull.map(grav2, as.table=TRUE)
map <- cbind(marker=rownames(map), map)
map <- apply(map, 2, function(a) gsub(" ", "", as.character(a)))
write.table(map, paste0(odir, "grav2_gmap.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# write phenotypes
phe <- as.matrix(grav2$pheno)
storage.mode(phe) <- "character"
phe <- cbind(as.character(1:nrow(phe)), phe)
colnames(phe)[1] <- "id"
write.table(phe, paste0(odir, "grav2_pheno.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# phenotype covariates
times <- as.numeric(substr(colnames(phe)[-1], 2, nchar(colnames(phe)[-1])))/60
phecovar <- cbind(pheno=colnames(phe)[-1], "time (hrs)"=as.character(times))
write.table(phecovar, paste0(odir, "grav2_phenocovar.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# control info
genotypes <- list(L=1L, C=2L)
grav2_info <- list(crosstype = "riself",
                  geno = "grav2_geno.csv",
                  pheno = "grav2_pheno.csv",
                  phenocovar = "grav2_phenocovar.csv",
                  gmap = "grav2_gmap.csv",
                  alleles = c("L", "C"),
                  genotypes = genotypes,
                  na.strings = c("-", "NA"))

library(yaml)
ofile <- paste0(odir, "grav2.yaml")
cat("# Data from Moore et al. (2013) Genetics 195:1077-1086",
    "# Available at QTL Archive, http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013b",
    file=ofile, sep="\n")
cat(as.yaml(grav2_info), file=ofile, append=TRUE)

# now read it in and save it in the data/ directory
### FIX ME ### need to get read_cross2 working properly first

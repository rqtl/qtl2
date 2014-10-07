# reformat data from Moore et al. (2003) Genetics 195:1077-1086
# in R/qtl2 file format

# data in R/qtlcharts package (http://kbroman.org/qtlcharts)
# also available at QTL Archive, http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013b

library(qtlcharts)
data(grav)

alleles <- attr(grav, "alleles")

odir <- "../sampledata/"

# write genotypes
g <- pull.geno(grav)
storage.mode(g) <- "character"
g[is.na(g)] <- "-"
g[g=="1"] <- paste(alleles[1], collapse="")
g[g=="2"] <- paste(alleles[2], collapse="")
g <- cbind(id=as.character(1:nrow(g)), g)
write.table(g, paste0(odir, "grav_geno.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# write genetic map
map <- pull.map(grav, as.table=TRUE)
map <- cbind(marker=rownames(map), map)
map <- apply(map, 2, function(a) gsub(" ", "", as.character(a)))
write.table(map, paste0(odir, "grav_gmap.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# write phenotypes
phe <- as.matrix(grav$pheno)
storage.mode(phe) <- "character"
phe <- cbind(as.character(1:nrow(phe)), phe)
colnames(phe)[1] <- "id"
write.table(phe, paste0(odir, "grav_pheno.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# phenotype covariates
times <- as.numeric(substr(colnames(phe)[-1], 2, nchar(colnames(phe)[-1])))/60
phecovar <- cbind(pheno=colnames(phe)[-1], "time (hrs)"=as.character(times))
write.table(phecovar, paste0(odir, "grav_phenocovar.csv"), sep=",",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

# control info
grav_info <- list(crosstype = "riself",
                  geno = "grav_geno.csv",
                  pheno = "grav_pheno.csv",
                  phenocovar = "grav_phenocovar.csv",
                  gmap = "grav_gmap.csv",
                  alleles = c("L", "C"),
                  genotypes = c("LL", "CC"),
                  na.strings = "NA")

library(yaml)
ofile <- paste0(odir, "grav.yaml")
cat("# Data from Moore et al. (2013) Genetics 195:1077-1086",
    "# Available at QTL Archive, http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013b",
    file=ofile, sep="\n")
cat(as.yaml(grav_info), file=ofile, append=TRUE)

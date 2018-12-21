# grab intensities from GeneSeek FinalReport.txt files
# convert them to a single big data frame, and save for use with fst

library(fst)

# simple version of data.table::fread()
myfread <- function(filename) data.table::fread(filename, data.table=FALSE, skip=9)

# data at https://doi.org/10.6084/m9.figshare.7359542.v1
#     also see https://github.com/rqtl/qtl2data/tree/master/DO_Svenson291
zip_files <- c("https://ndownloader.figshare.com/files/13599554",
               "https://ndownloader.figshare.com/files/13599572")

# download, unzip, and read the data
dat <- vector("list", length(zip_files))
for(i in seq_along(zip_files)) {
    zipfile <- tempfile()
    download.file(zip_files[i], zipfile)
    unzipped <- unzip(zipfile)
    dat[[i]] <- myfread(unzipped)
    unlink(zipfile)
}

# rbind the results together, saving selected columns
dat <- do.call("rbind", dat)[,c("SNP Name", "Sample ID", "X", "Y")]

# create matrices that are snps x samples
snps <- unique(dat[,"SNP Name"])
samples <- unique(dat[,"Sample ID"])
X <- Y <- matrix(ncol=length(samples), nrow=length(snps))
dimnames(X) <- dimnames(Y) <- list(snps, samples)
for(i in seq(along=samples)) {
    message(i, " of ", length(samples))
    tmp <- dat[dat[,"Sample ID"]==samples[i],]
    X[,samples[i]] <- tmp[,"X"]
    Y[,samples[i]] <- tmp[,"Y"]
}

# bring together in one matrix
result <- cbind(snp=rep(snps, 2),
                channel=rep(c("X", "Y"), each=length(snps)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps), seq_along(snps)+length(snps)))),]
rownames(result) <- 1:nrow(result)

# write to fst file, maximally compressed
write.fst(result, "intensities.fst", compress=100)

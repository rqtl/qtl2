# convert GeneSeek FinalReport files to format for R/qtl2
#
# - creates one genotype CSV file for each chromosome
#
# - also creates 4 files containing the two channels of SNP intensities for markers on the X and Y chr
#   (these are useful for verifying the sex of the mice)

# file containing allele codes for GigaMUGA data
#   - from GM_processed_files.zip, https://doi.org/10.6084/m9.figshare.5404759
codefile <- "GM/GM_allelecodes.csv"

# input files with GigaMUGA genotypes
#  - can be a single file or a vector of multiple files
#  - if samples appear in multiple files, the genotypes in later files
#    will be used in place of genotypes in earlier files
#  - files can be gzipped (".gz" extension)
ifiles <- c("batch1_FinalReport.txt",
            "batch2_FinalReport.txt")

# file "stem" for output files
# output files will be like "gm4qtl2_geno19.csv"
ostem <- "gm4qtl2"

##############################
# define a couple of functions
##############################
# simple version of data.table::fread()
myfread <- function(filename) data.table::fread(filename, data.table=FALSE)

# cbind, replacing matching columns with second set and adding unique ones
cbind_smother <-
  function(mat1, mat2)
  {
      cn1 <- colnames(mat1)
      cn2 <- colnames(mat2)
      m <- (cn2 %in% cn1)
      if(any(m)) {
          mat1[,cn2[m]] <- mat2[,cn2[m],drop=FALSE]
          if(any(!m)) {
              mat1 <- cbind(mat1, mat2[,cn2[!m],drop=FALSE])
          }
      }
      else {
          mat1 <- cbind(mat1, mat2)
      }

      mat1
  }
##############################



# read genotype codes
codes <- myfread(codefile)

full_geno <- NULL
cXint <- cYint <- NULL

for(ifile in ifiles) {
    cat(" -File:", ifile, "\n")
    rezip <- FALSE
    if(!file.exists(ifile)) {
        cat(" -Unzipping file\n")
        system(paste("gunzip", ifile))
        rezip <- TRUE
    }

    cat(" -Reading data\n")
    g <- myfread(ifile)
    # subset to the markers in the codes object
    g <- g[g[,"SNP Name"] %in% codes[,"marker"],]

    # NOTE: may need to revise the IDs in the 2nd column
    # We're assuming the IDs are numbers and creating revised IDs as DO-###
    samples <- unique(g[,"Sample ID"])

    # matrix to contain the genotypes
    geno <- matrix(nrow=nrow(codes), ncol=length(samples))
    dimnames(geno) <- list(codes[,"marker"], samples)

    # fill in matrix
    cat(" -Reorganizing data\n")
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (g[,"Sample ID"]==samples[i])
        geno[g[wh,"SNP Name"],i] <- paste0(g[wh,"Allele1 - Forward"], g[wh,"Allele2 - Forward"])
    }

    cat(" -Encode genotypes\n")
    geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))

    if(is.null(full_geno)) {
        full_geno <- geno
    } else {
        # if any columns in both, use those from second set
        full_geno <- cbind_smother(full_geno, geno)
    }

    # grab X and Y intensities
    cat(" -Grab X and Y intensities\n")
    gX <- g[g[,"SNP Name"] %in% codes[codes$chr=="X","marker"],]
    gY <- g[g[,"SNP Name"] %in% codes[codes$chr=="Y","marker"],]
    cX <- matrix(nrow=sum(codes$chr=="X"),
                 ncol=length(samples))
    dimnames(cX) <- list(codes[codes$chr=="X","marker"], samples)
    cY <- matrix(nrow=sum(codes$chr=="Y"),
                 ncol=length(samples))
    dimnames(cY) <- list(codes[codes$chr=="Y","marker"], samples)
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (gX[,"Sample ID"]==samples[i])
        cX[gX[wh,"SNP Name"],i] <- (gX$X[wh] + gX$Y[wh])/2

        wh <- (gY[,"Sample ID"]==samples[i])
        cY[gY[wh,"SNP Name"],i] <- (gY$X[wh] + gY$Y[wh])/2
    }
    if(is.null(cXint)) {
        cXint <- cX
        cYint <- cY
    } else {
        # if any columns in both, use those from second set
        cXint <- cbind_smother(cXint, cX)
        cYint <- cbind_smother(cYint, cY)
    }

    if(rezip) {
        cat(" -Rezipping file\n")
        system(paste("gzip", ifile))
    }
}

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(cbind(marker=rownames(cXint), cXint),
                       paste0(ostem, "_chrXint.csv"),
                       paste(ostem, "X chr intensities"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint), cYint),
                       paste0(ostem, "_chrYint.csv"),
                       paste(ostem, "Y chr intensities"),
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chr in c(1:19,"X","Y","M")) {
    mar <- codes[codes$chr==chr,"marker"]
    g <- full_geno[mar,]
    qtl2convert::write2csv(cbind(marker=rownames(g), g),
                           paste0(ostem, "_geno", chr, ".csv"),
                           paste0(ostem, " genotypes for chr ", chr),
                           overwrite=TRUE)
}

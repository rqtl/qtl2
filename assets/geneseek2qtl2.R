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
cXint1 <- cXint2 <- NULL
cYint1 <- cYint2 <- NULL

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
    g <- g[g[,1] %in% codes[,1],]

    # NOTE: may need to revise the IDs in the 2nd column
    # We're assuming the IDs are numbers and creating revised IDs as DO-###
    samples <- unique(g[,2])

    # matrix to contain the genotypes
    geno <- matrix(nrow=nrow(codes), ncol=length(samples))
    dimnames(geno) <- list(codes[,1], samples)

    # fill in matrix
    cat(" -Reorganizing data\n")
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (g[,2]==samples[i])
        geno[g[wh,1],i] <- paste0(g[wh,3], g[wh,4])
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
    gX <- g[g[,1] %in% codes[codes$chr=="X",1],]
    gY <- g[g[,1] %in% codes[codes$chr=="Y",1],]
    cX1 <- matrix(nrow=sum(codes$chr=="X"),
                  ncol=length(samples))
    dimnames(cX1) <- list(codes[codes$chr=="X",1], samples)
    cX2 <- cX1
    cY1 <- matrix(nrow=sum(codes$chr=="Y"),
                  ncol=length(samples))
    dimnames(cY1) <- list(codes[codes$chr=="Y",1], samples)
    cY2 <- cY1
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (gX[,2]==samples[i])
        cX1[gX[wh,1],i] <- gX$X[wh]
        cX2[gX[wh,1],i] <- gX$Y[wh]

        wh <- (gY[,2]==samples[i])
        cY1[gY[wh,1],i] <- gY$X[wh]
        cY2[gY[wh,1],i] <- gY$Y[wh]
    }
    if(is.null(cXint1)) {
        cXint1 <- cX1
        cXint2 <- cX2
        cYint1 <- cY1
        cYint2 <- cY2
    } else {
        # if any columns in both, use those from second set
        cXint1 <- cbind_smother(cXint1, cX1)
        cXint2 <- cbind_smother(cXint2, cX2)
        cYint1 <- cbind_smother(cYint1, cY1)
        cYint2 <- cbind_smother(cYint2, cY2)
    }

    if(rezip) {
        cat(" -Rezipping file\n")
        system(paste("gzip", ifile))
    }
}

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(cbind(marker=rownames(cXint1), cXint1),
                       paste0(ostem, "_chrXint1.csv"),
                       paste(ostem, "X chr intensities, channel 1"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cXint2), cXint2),
                       paste0(ostem, "_chrXint2.csv"),
                       paste(ostem, "X chr intensities, channel 2"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint1), cYint1),
                       paste0(ostem, "_chrYint1.csv"),
                       paste(ostem, "Y chr intensities, channel 1"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint2), cYint2),
                       paste0(ostem, "_chrYint2.csv"),
                       paste(ostem, "Y chr intensities, channel 2"),
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chr in c(1:19,"X","Y","M")) {
    mar <- codes[codes$chr==chr,1]
    g <- full_geno[mar,]
    qtl2convert::write2csv(cbind(marker=rownames(g), g),
                           paste0(ostem, "_geno", chr, ".csv"),
                           paste0(ostem, " genotypes for chr ", chr),
                           overwrite=TRUE)
}

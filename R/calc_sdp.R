#' Calculate strain distribution pattern from SNP genotypes
#'
#' Calculate the strain distribution patterns (SDPs) from the strain
#' genotypes at a set of SNPs.
#'
#' @param geno Matrix of SNP genotypes, markers x strains, coded as 1
#' (AA) and 3 (BB). Markers with values other than 1 or 3 are omitted,
#' and monomorphic markers, are omitted.
#'
#' @return A vector of strain distribution patterns: integers between
#' 1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#' binary representation indicates the strain genotypes.
#'
#' @seealso [invert_sdp()], [sdp2char()]
#'
#' @export
#' @examples
#' x <- rbind(m1=c(3, 1, 1, 1, 1, 1, 1, 1),
#'            m2=c(1, 3, 3, 1, 1, 1, 1, 1),
#'            m3=c(1, 1, 1, 1, 3, 3, 3, 3))
#' calc_sdp(x)
calc_sdp <-
    function(geno)
{
    # tolerate data frames, but convert to matrix
    if(!is.matrix(geno)) {
        if(is.data.frame(geno)) {
            geno <- as.matrix(geno)
        } else {
            geno <- rbind(geno)
            dimnames(geno) <- NULL
        }
    }

    n_str <- ncol(geno)

    n_missing <- rowSums(is.na(geno) | (geno != 1 & geno != 3))
    if(any(n_missing > 0)) {
        geno <- geno[n_missing==0,,drop=FALSE]
        warning("Omitting ", sum(n_missing>0), " markers with missing or invalid genotypes")
        if(nrow(geno) == 0)
            stop("No markers with valid data.")
    }

    n_AA <- rowSums(geno==1)
    monomorphic <- (n_AA==0 | n_AA==n_str)
    if(any(monomorphic)) {
        geno <- geno[!monomorphic,,drop=FALSE]
        warning("Omitting ", sum(monomorphic), " monomorphic markers")
    }

    if(nrow(geno) == 0)
        stop("No polymorphic markers.")

    result <- .calc_sdp((geno- 1)/2) # convert matrix to 0/1
    names(result) <- rownames(geno)
    result
}


#' Calculate SNP genotype matrix from strain distribution patterns
#'
#' Calculate the matrix of SNP genotypes from a vector of strain distribution patterns (SDPs).
#'
#' @param sdp Vector of strain distribution patterns (integers between
#' 1 and \eqn{2^n-2} where \eqn{n} is the number of strains.
#' @param n_strains Number of strains
#'
#' @return Matrix of SNP genotypes, markers x strains, coded as 1
#' (AA) and 3 (BB). Markers with values other than 1 or 3 are omitted,
#' and monomorphic markers, are omitted.
#'
#' @seealso [sdp2char()], [calc_sdp()]
#'
#' @export
#' @examples
#' sdp <- c(m1=1, m2=12, m3=240)
#' invert_sdp(sdp, 8)
invert_sdp <-
    function(sdp, n_strains)
{
    if(!is.numeric(sdp)) stop("sdp should be a vector of integers")
    if(!is_pos_number(n_strains)) stop("n_strains should be a single positive integer")

    geno <- .invert_sdp(sdp, n_strains)
    rownames(geno) <- names(sdp)
    geno <- geno*2+1 # convert 0/1 to 1/3

    geno
}

#' Convert strain distribution patterns to character strings
#'
#' Convert a vector of numeric codes for strain distribution patterns to character strings.
#'
#' @param sdp Vector of strain distribution patterns (integers between
#' 1 and \eqn{2^n-2} where \eqn{n} is the number of strains.
#' @param n_strains Number of founder strains (if missing but
#' `strains` is provided, we use the length of `strains`)
#' @param strains Vector of single-letter codes for the strains
#'
#' @return Vector of character strings with the two groups of alleles
#'     separated by a vertical bar (|).
#'
#' @seealso [invert_sdp()], [calc_sdp()]
#'
#' @export
#' @examples
#' sdp <- c(m1=1, m2=12, m3=240)
#' sdp2char(sdp, 8)
#' sdp2char(sdp, strains=c("A", "B", "1", "D", "Z", "C", "P", "W"))
sdp2char <-
    function(sdp, n_strains=NULL, strains=NULL)
{
    if(is.null(n_strains) && !is.null(strains)) {
        n_strains <- length(strains)
    }
    if(!is.null(n_strains) && is.null(strains)) {
        strains <- LETTERS[seq_len(n_strains)]
    }
    if(is.null(n_strains) && is.null(strains)) {
        stop("Provide at least one of n_strains or strains")
    }

    # convert to matrix
    pat <- invert_sdp(sdp, n_strains)

    # convert to character strings
    result <- apply(pat, 1, function(p)
        paste(paste(strains[p==3], collapse=""),
              paste(strains[p==1], collapse=""), sep="|") )

    names(result) <- names(sdp)

    result
}

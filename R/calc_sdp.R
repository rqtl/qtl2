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
#' @seealso \code{\link{invert_sdp}}
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
#' @seealso \code{\link{calc_sdp}}
#'
#' @export
#' @examples
#' sdp <- c(m1=1, m2=12, m3=240)
#' invert_sdp(sdp, 8)
invert_sdp <-
    function(sdp, n_strains)
{
    geno <- .invert_sdp(sdp, n_strains)
    rownames(geno) <- names(sdp)
    geno <- geno*2+1 # convert 0/1 to 1/3

    geno
}

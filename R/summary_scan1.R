# functions to summarize scan1 results

#' Find position with maximum LOD score
#'
#' Return data frame with the positions having maximum LOD score for a
#' particular LOD score column
#'
#' @param scan1_output An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}} or \code{\link{scan1_lmm}}.
#' @param lodcolumn An integer or character string indicating the LOD
#' score column, either as a numeric index or column name.
#' @param chr Option vector of chromosomes to consider.
#' @param na.rm Ignored (take to be TRUE)
#' @param ... Ignored
#'
#' @export
#'
#' @return A data.frame with three columns: chr, pos, and lod score.
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, covar, Xcovar)
#'
#' # maximum of first column
#' max(out)
#'
#' # maximum of spleen column
#' max(out, lodcolumn="spleen")
#'
#' # maximum of first column on chr 2
#' max(out, chr="2")
max_scan1 <-
    function(scan1_output, lodcolumn=1, chr, na.rm=TRUE, ...)
{
    thechr <- chr_scan1(scan1_output)
    thepos <- pos_scan1(scan1_output)
    map <- scan1_output$map

    # to handle output of either scan1() or scan1coef()
    # for coef(), look at the sign
    if("lod" %in% names(scan1_output)) {
        lod <- scan1_output$lod
        sign <- (lod >= 0)*2-1 # sign +1
    }
    else if("coef" %in% names(scan1_output)) {
        lod <- scan1_output$coef
        sign <- (lod >= 0)*2-1 # sign +1/-1
        lod <- abs(lod)
    }
    else stop("Neither lod nor coef found.")

    coln <- colnames(lod)

    mnames <- rownames(lod)

    if(length(lodcolumn) > 1) { # If length > 1, take first value
        warning("lodcolumn should have length 1; one first element used.")
        lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, coln)
        if(is.na(tmp))
            stop('column "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(lod))
        stop("column [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(lod), ")")
    lod <- lod[,lodcolumn]
    sign <- sign[,lodcolumn]

    # subset chromosomes
    if(!missing(chr) && !is.null(chr)) {
        keep <- (thechr %in% chr)
        thechr <- thechr[keep]
        thepos <- thepos[keep]
        mnames <- mnames[keep]
        lod <- lod[keep]*sign[keep] # sign is a kludge to handle scan1coef output
    }

    themax <- which.max(lod)

    result <- data.frame(chr=thechr[themax],
                         pos=thepos[themax],
                         lod=lod[themax],
                         stringsAsFactors=FALSE)
    rownames(result) <- mnames[themax]
    names(result)[3] <- coln[lodcolumn]

    result
}

#' @export
#' @rdname max_scan1
max.scan1 <- max_scan1

#' Overall maximum LOD score
#'
#' Find overall maximum LOD score in genome scan results.
#'
#' @param scan1_output An object of class \code{"scan1"} as returned by
#' \code{\link{scan1}} or \code{\link{scan1_lmm}}.
#' @param chr Option vector of chromosomes to consider.
#'
#' @export
#' @return A single number: the maximum LOD score across all columns and positions for
#' the selected chromosomes.
#'
#' @examples
#' \dontrun{
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, covar, Xcovar)
#'
#' # overall maximum
#' maxlod(out)
#'
#' # maximum on chromosome 2
#' max(out, "2")
#' }
maxlod <-
    function(scan1_output, chr)
{
    # subset by chromosome
    scan1_output <- subset(scan1_output, chr=chr)

    # to handle output of either scan1() or scan1coef()
    # for coef(), look at the sign
    if("lod" %in% names(scan1_output)) {
        lod <- scan1_output$lod
        return(max(lod))
    }
    else if("coef" %in% names(scan1_output)) {
        coef <- scan1_output$coef
        sign <- (coef >= 0)*2-1 # sign +1/-1
        coef <- abs(coef)
        maxcoef <- max(coef)

        # deal with ties
        wh <- which(coef == maxcoef)
        if(length(wh)>1) wh <- sample(wh, 1)

        return(maxcoef * sign[wh])
    }
    else stop("Neither lod nor coef found.")
}

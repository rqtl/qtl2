# summary functions

#' Basic summaries of a cross2 object
#'
#' Basic summaries of a cross2 object.
#'
#' @param cross2 An object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#'
#' @name basic_summaries
#' @keywords utilities
#' @seealso \code{\link{summary.cross2}}
NULL

#' @describeIn basic_summaries Number of phenotyped individuals
#' @export
n_ind <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$pheno)) {
        return(0)
    }

    nrow(cross2$pheno)
}

#' @describeIn basic_summaries IDs of phenotyped individuals
#' @export
ind_ids <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$pheno)) {
        return(NULL)
    }

    rownames(cross2$pheno)
}

#' @export
#' @describeIn basic_summaries Number of genotyped lines
n_lines <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) {
        warning("No genotypes present")
        return(0)
    }

    nrow(cross2$geno[[1]])
}

#' @describeIn basic_summaries IDs of genotyped lines
#' @export
line_ids <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) {
        warning("No genotypes present")
        return(NULL)
    }

    rownames(cross2$geno[[1]])
}

#' @describeIn basic_summaries Number of chromosomes
#' @export
n_chr <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) {
        warning("No genotypes present")
        return(0)
    }

    length(cross2$geno)
}

#' @describeIn basic_summaries Chromosome names
#' @export
chr_names <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) {
        warning("No genotypes present")
        return(NULL)
    }

    names(cross2$geno)
}

#' @describeIn basic_summaries Total number of markers
#' @export
tot_mar <-
function(cross2)
{
    sum(n_mar(cross2))
}

#' @describeIn basic_summaries Number of markers on each chromosome
#' @export
n_mar <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) {
        warning("No genotypes present")
        return(0)
    }

    vapply(cross2$geno, ncol, 1)
}

#' @describeIn basic_summaries Marker names
#' @export
marker_names <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) {
        warning("No genotypes present")
        return(NULL)
    }

    nam <- unlist(lapply(cross2$geno, colnames))
    names(nam) <- NULL
    nam
}


#' @describeIn basic_summaries Number of phenotypes
#' @export
n_pheno <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$pheno)) {
        return(NULL)
    }

    ncol(cross2$pheno)
}

#' @describeIn basic_summaries Phenotype names
#' @export
pheno_names <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$pheno)) {
        return(NULL)
    }

    colnames(cross2$pheno)
}

#' @describeIn basic_summaries Number of covariates
#' @export
n_covar <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$covar)) {
        return(0)
    }

    ncol(cross2$covar)
}

#' @describeIn basic_summaries Covariate names
#' @export
covar_names <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$covar)) {
        return(NULL)
    }

    colnames(cross2$covar)
}
#' @describeIn basic_summaries Number of phenotype covariates
#' @export
n_phenocovar <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$phenocovar)) {
        return(0)
    }

    ncol(cross2$phenocovar)
}

#' @describeIn basic_summaries Phenotype covariate names
#' @export
phenocovar_names <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$phenocovar)) {
        return(NULL)
    }

    colnames(cross2$phenocovar)
}


# print cross
#' Print a cross2 object
#'
#' Print a summary of a cross2 object
#'
#' @param x An object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param ... Ignored.
#'
#' @return None.
#'
#' @export
#' @keywords utilities
print.cross2 <-
function(x, ...)
{
    print(summary(x))

    invisible(NULL)
}

# summary_cross
#' Summary of cross2 object
#'
#' Summarize a cross2 object
#'
#' @param object An object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param ... Ignored.
#'
#' @return None.
#'
#' @export
#' @keywords utilities
#' @seealso \link[=basic_summaries]{basic summaries}
summary.cross2 <-
function(object, ...)
{
    # run checks
    check <- check_cross2(object)

    result <- list(crosstype=object$crosstype,
                   nlines=n_lines(object),
                   nind=n_ind(object),
                   nchr=n_chr(object),
                   nmar=n_mar(object),
                   npheno=n_pheno(object),
                   ncovar=n_covar(object),
                   nphenocovar=n_phenocovar(object),
                   totmar=tot_mar(object))

    class(result) <- c("summary.cross2", "list")
    result
}

#' @export
print.summary.cross2 <-
function(x, ...)
{
    cat('Object of class cross2 (crosstype "', x$crosstype, '")\n\n', sep='')

    toprint <- c("\bNo. genotyped lines"=        x$nlines,
                 "\bNo. phenotyped individuals"= x$nind,
                 "\nNo. phenotypes"=             x$npheno,
                 "\bNo. covariates"=             x$ncovar,
                 "\bNo. phenotype covariates"=   x$nphenocovar,
                 "\nNo. chromosomes"=            x$nchr,
                 "\bTotal markers"=              x$totmar)

    print_aligned(toprint)
    cat("No. markers by chr:\n")
    print(x$nmar)

    invisible(x)
}

print_aligned <-
function(x)
{
    ndig <- ceiling(log10(max(x)))
    n.char <- max(nchar(names(x)))

    format <- paste0("%-", n.char+1, "s  %", ndig, "d\n")

    for(i in seq(along=x))
        cat(sprintf(format, names(x)[i], x[i]))
}

# is an object a cross2 object?
# look for "cross2" within class()
is.cross2 <-
    function(x)
    "cross2" %in% class(x)

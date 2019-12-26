# summary functions

#' Basic summaries of a cross2 object
#'
#' Basic summaries of a cross2 object.
#'
#' @param cross2 An object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#'
#' @name basic_summaries
#' @keywords utilities
#' @seealso [summary.cross2()]
NULL

#' @describeIn basic_summaries Number of individuals (either genotyped or phenotyped)
#' @export
n_ind <-
function(cross2)
{
    length( ind_ids(cross2) )
}

#' @describeIn basic_summaries Number of genotyped individuals
#' @export
n_ind_geno <-
function(cross2)
{
    length( ind_ids_geno(cross2) )
}

#' @describeIn basic_summaries Number of phenotyped individuals
#' @export
n_ind_pheno <-
function(cross2)
{
    length( ind_ids_pheno(cross2) )
}

#' @describeIn basic_summaries Number of individuals with covariate data
#' @export
n_ind_covar <-
function(cross2)
{
    length( ind_ids_covar(cross2) )
}

#' @describeIn basic_summaries Number of individuals with both genotype and phenotype data
#' @export
n_ind_gnp <-
function(cross2)
{
    length( ind_ids_gnp(cross2) )
}

#' @describeIn basic_summaries IDs of individuals (either genotyped or phenotyped)
#' @export
ind_ids <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input cross must have class "cross2"')
    unique( c(ind_ids_geno(cross2), ind_ids_pheno(cross2), ind_ids_covar(cross2)) )
}

#' @describeIn basic_summaries IDs of genotyped individuals
#' @export
ind_ids_geno <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$geno)) return(NULL)

    rownames(cross2$geno[[1]])
}

#' @describeIn basic_summaries IDs of phenotyped individuals
#' @export
ind_ids_pheno <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$pheno)) return(NULL)

    rownames(cross2$pheno)
}

#' @describeIn basic_summaries IDs of individuals with covariate data
#' @export
ind_ids_covar <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input must have class "cross2"')

    if(is.null(cross2$covar)) return(NULL)

    rownames(cross2$covar)
}

#' @describeIn basic_summaries IDs of individuals with both genotype and phenotype data
#' @export
ind_ids_gnp <-
function(cross2)
{
    pheID <- ind_ids_pheno(cross2)
    genID <- ind_ids_geno(cross2)

    # if either is empty, return NULL
    if(is.null(pheID) || is.null(genID)) return(NULL)

    id <- find_common_ids(genID, pheID)

    # if no individuals in common, return NULL
    if(length(id) == 0) return(NULL)

    id
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

#' @describeIn basic_summaries Number of founder strains
#' @export
n_founders <-
    function(cross2)
{
    if("founder_geno" %in% names(cross2)) {
        return(nrow(cross2$founder_geno[[1]]))
    }

    return(2) # if no founder_geno, must be just two founders, right?
}

#' @describeIn basic_summaries Names of founder strains
#' @export
founders <-
    function(cross2)
{
    if("founder_geno" %in% names(cross2)) {
        return(rownames(cross2$founder_geno[[1]]))
    }

    return(cross2$alleles)
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
#' @param x An object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
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
#' @param object An object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param ... Ignored.
#'
#' @return None.
#'
#' @export
#' @keywords utilities
#' @seealso [basic_summaries]
summary.cross2 <-
function(object, ...)
{
    # run checks
    check <- check_cross2(object)

    result <- list(crosstype=object$crosstype,
                   nind=n_ind(object),
                   nind_geno=n_ind_geno(object),
                   nind_pheno=n_ind_pheno(object),
                   nind_gnp=n_ind_gnp(object),
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

    toprint <- c("Total individuals"=           x$nind,
                 "No. genotyped individuals"=   x$nind_geno,
                 "No. phenotyped individuals"=  x$nind_pheno,
                 "No. with both geno & pheno"=  x$nind_gnp,
                 "\nNo. phenotypes"=            x$npheno,
                 "No. covariates"=              x$ncovar,
                 "No. phenotype covariates"=    x$nphenocovar,
                 "\nNo. chromosomes"=           x$nchr,
                 "Total markers"=               x$totmar)

    print_aligned(toprint)

    cat("\nNo. markers by chr:\n")
    print(x$nmar)

    invisible(x)
}

print_aligned <-
function(x)
{
    newline_before <- grepl("\\n", names(x))
    names(x) <- sub("^\\n", "", names(x))

    ndig <- ceiling(log10(max(x, na.rm=TRUE)))
    n.char <- max(nchar(names(x)))

    format <- paste0("%-", n.char+1, "s  %", ndig, "d\n")

    for(i in seq(along=x)) {
        if(newline_before[i]) cat("\n")
        cat(sprintf(format, names(x)[i], x[i]))
    }
}

# is an object a cross2 object?
# look for "cross2" within class()
is.cross2 <-
    function(x)
    inherits(x, "cross2")

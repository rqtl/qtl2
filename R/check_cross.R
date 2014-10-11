# functions to check the integrity of QTL cross data

# check if a crosstype is supported
#   if give_error=TRUE, problem results in error
#   if give_error=FALSE, returns TRUE if supported; FALSE if not
check_crosstype <-
function(crosstype, give_error=TRUE)
{
    z <- try(.crosstype_supported(crosstype), silent=TRUE)

    if("try-error" %in% class(z) || !z) {
        if(give_error)
            stop("Cross type ", crosstype, " not yet supported.")
        return(invisible(FALSE))
    }

    invisible(TRUE) # if no error, return TRUE
}


# count the number of invalid genotypes
# returns a matrix individuals x chromosomes
count_invalid_genotypes <-
function(cross2)
{
    result <- matrix(nrow=nrow(cross2$geno[[1]]),
                     ncol=length(cross2$geno))
    dimnames(result) <- list(rownames(cross2$geno[[1]]),
                             names(cross2$geno))

    cross_info <- t(cross2$cross_info)

    for(i in seq(along=cross2$geno))
        result[,i] <- .count_invalid_genotypes(cross2$crosstype,
                                               t(cross2$geno[[i]]),
                                               cross2$is_x_chr[i],
                                               cross2$is_female,
                                               cross_info)

    result
}

# check_cross2
#' Check a cross2 object
#'
#' Check the integrity of the data within a cross2 object.
#'
#' @param cross2 An object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#'
#' @return If everything is correct, returns \code{TRUE}; otherwise \code{FALSE},
#' with attributes that give the problems.
#'
#' @details Checks whether a cross2 object meets the
#' specifications. Problems are issued as warnings.
#'
#' @export
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' check_cross2(grav2)

check_cross2 <-
function(cross2)
{
    result <- TRUE

    # required pieces
    #     crosstype
    #     pheno
    #     geno
    #     gmap
    #     is_female
    #     is_x_chr
    #     cross_info

    # optional pieces
    #     covar
    #     phenocovar
    #     pmap
    #     linemap (but required if nrow(pheno) != nrow(geno[[1]])
    #     foundergeno (required for many crosstypes...add need_foundergeno function?)

    crosstype <- cross2$crosstype
    if(is.null(crosstype)) {
        result <- FALSE
        warning("crosstype is missing")
    }
    else if(!check_crosstype(cross2$crosstype)) {
        result <- FALSE
        warning("Crosstype ", cross2$crosstype, " is not supported")
    }

    # pheno
    pheno <- cross2$pheno
    if(is.null(pheno)) { # pheno required
        result <- FALSE
        warning("pheno is missing")
    }
    else {
        if(is.null(rownames(pheno))) {
            result <- FALSE
            warning("rownames(pheno) is NULL")
        }
        if(is.null(colnames(pheno))) {
            result <- FALSE
            warning("colnames(pheno) is NULL")
        }
        if(!is.matrix(pheno)) {
            result <- FALSE
            warning("pheno is not a matrix")
        }
        if(storage.mode(pheno) != "double") {
            result <- FALSE
            warning("pheno is not stored as doubles but rather ", storage.mode(pheno))
        }

        # covar
        covar <- cross2$covar
        if(!is.null(covar)) { # covar is optional
            if(nrow(pheno) != nrow(covar)) {
                result <- FALSE
                warning("nrow(pheno) (", nrow(pheno), ") != nrow(covar) (", nrow(covar), ")")
            }
            else if(any(rownames(pheno) != rownames(covar))) {
                result <- FALSE
                warning("rownames(pheno) != rownames(covar)")
            }
        }

        # phenocovar
        phenocovar <- cross2$phenocovar
        if(!is.null(phenocovar)) { # phenocovar is optional
            if(ncol(pheno) != nrow(phenocovar)) {
                result <- FALSE
                warning("ncol(pheno) (", ncol(pheno), ") != nrow(phenocovar) (", nrow(phenocovar), ")")
            }
            else if(any(colnames(pheno) != rownames(phenocovar))) {
                result <- FALSE
                warning("colnames(pheno) != rownames(phenocovar)")
            }
        }
    }

    geno <- cross2$geno
    gmap <- cross2$gmap
    if(is.null(geno)) {
        result <- FALSE
        warning("geno is missing")
    }
    if(is.null(gmap)) {
        result <- FALSE
        warning("gmap is missing")
    }
    if(is.null(geno)) {
        # compare geno to gmap
        if(!is.null(gmap)) {
            if(length(geno) != length(gmap)) {
                result <- FALSE
                warning("length(geno) (", length(geno), ") != length(gmap) (", length(gmap), ")")
            }
            else {
                if(any(names(geno) != names(gmap))) {
                    result <- FALSE
                    warning("names(geno) != names(gmap)")
                }
                for(i in seq(along=geno)) {
                    if(length(geno[[i]]) != length(gmap[[i]])) {
                        result <- FALSE
                        warning("Mismatch between geno and gmap in no. markers on chr ", names(geno)[i])
                    }
                    else if(any(colnames(geno[[i]]) != names(gmap[[i]]))) {
                        result <- FALSE
                        warning("Mismatch in marker names between geno and gmap on chr ", names(geno)[i])
                    }
                }
            }
        }

        # is_female
        is_female <- cross2$is_female
        if(is.null(is_female)) {
            result <- FALSE
            warning("is_female is missing")
        }
        else if(nrow(geno) != length(is_female)) {
            result <- FALSE
            warning("nrow(geno) (", nrow(geno), ") != length(is_female) (", length(is_female), ")")
        }
        else if(any(rownames(geno) != names(is_female))) {
            result <- FALSE
            warning("rownames(geno) != names(is_female)")
        }
        if(!is.logical(is_female)) {
            result <- FALSE
            warning("is_female is not logical")
        }

        # cross_info
        cross_info <- cross2$cross_info
        if(is.null(cross_info)) {
            result <- FALSE
            warning("cross_info is missing")
        }
        else {
            if(nrow(geno) != nrow(cross_info)) {
                result <- FALSE
                warning("nrow(geno) (", nrow(geno), ") != nrow(cross_info) (", nrow(cross_info), ")")
            }
            else if(any(rownames(geno) != nrownames(cross_info))) {
                result <- FALSE
                warning("rownames(geno) != nrownames(cross_info)")
            }
            if(!is.matrix(cross_info)) {
                result <- FALSE
                warning("cross_info is not a matrix")
            }
            if(storage.mode(cross_info) != "integer") {
                result <- FALSE
                warning("cross_info is not stored as integers but rather ", storage.mode(cross_info))
            }
        }

        linemap <- cross2$linemap
        if(!is.null(linemap) && !is.null(pheno)) {
            if(length(linemap) != nrow(pheno)) {
                result <- FALSE
                warning("length(linemap) (", length(linemap), ") != nrow(pheno) (", nrow(pheno), ")")
            }
            else if(any(names(linemap) != rownames(pheno))) {
                result <- FALSE
                warning("names(linemap) != rownames(pheno)")
            }

            if(any(is.na(match(linemap, rownames(geno))))) {
                result <- FALSE
                warning("Some lines in linemap are not in rownames(geno)")
            }
        }

        if(is.null(linemap) && nrow(pheno) != nrow(geno)) {
            result <- FALSE
            warning("linemap missing but nrow(pheno) != nrow(geno)")
        }

        foundergeno <- cross2$foundergeno
        if(!is.null(foundergeno)) { # foundergeno is optional

            if(length(geno) != length(foundergeno)) {
                result <- FALSE
                warning("length(geno) (", length(geno), ") != length(foundergeno) (", length(foundergeno), ")")
            }
            else {
                if(any(names(geno) != names(foundergeno))) {
                    result <- FALSE
                    warning("names(geno) != names(foundergeno)")
                }
                for(i in seq(along=geno)) {
                    if(length(geno[[i]]) != length(foundergeno[[i]])) {
                        result <- FALSE
                        warning("Mismatch between geno and foundergeno in no. markers on chr ", names(geno)[i])
                    }
                    else if(any(colnames(geno[[i]]) != names(foundergeno[[i]]))) {
                        result <- FALSE
                        warning("Mismatch in marker names between geno and foundergeno on chr ", names(geno)[i])
                    }
                }
            }
        }

    }
    else { # geno is absent
        if(is.null(cross2$is_female))
            warning("is_female is missing")

        if(is.null(cross2$cross_info))
            warning("cross_info is missing")

        if(is.null(cross2$is_x_chr))
            warning("is_x_chr is missing")
    }


    pmap <- cross2$pmap
    if(!is.null(pmap) && !is.null(gmap)) { # pmap is optional
        if(length(gmap) != length(pmap)) {
            result <- FALSE
            warning("length(gmap) (", length(gmap), ") != length(pmap) (", length(pmap), ")")
        }
        else {
            if(any(names(gmap) != names(pmap))) {
                result <- FALSE
                warning("names(gmap) != names(pmap)")
            }
            for(i in seq(along=gmap)) {
                if(length(gmap[[i]]) != length(pmap[[i]])) {
                    result <- FALSE
                    warning("Mismatch between gmap and pmap in no. markers on chr ", names(gmap)[i])
                }
                else if(any(colnames(gmap[[i]]) != names(pmap[[i]]))) {
                    result <- FALSE
                    warning("Mismatch in marker names between gmap and pmap on chr ", names(gmap)[i])
                }
            }
        }
    }

    if(!result) { # check genotypes only if everything else is okay
        n_invalid <- count_invalid_genotypes(cross2)
        if(sum(n_invalid)>0)
            warning(sum(n_invalid), " genotypes in cross")
    }

    result
}


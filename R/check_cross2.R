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
    if(!is.cross2(cross2))
        stop('Input cross must have class "cross2"')

    result <- matrix(nrow=nrow(cross2$geno[[1]]),
                     ncol=length(cross2$geno))
    dimnames(result) <- list(rownames(cross2$geno[[1]]),
                             names(cross2$geno))

    ind <- rownames(cross2$geno[[1]])

    # handle case of missing cross_info or is_female
    cross_info <- handle_null_crossinfo(cross2$cross_info, ind)
    is_female <- handle_null_isfemale(cross2$is_female, ind)
    is_x_chr <- handle_null_isxchr(cross2$is_x_chr, names(cross2$geno))

    cross_info <- t(cross_info)

    for(i in seq(along=cross2$geno))
        result[,i] <- .count_invalid_genotypes(cross2$crosstype,
                                               t(cross2$geno[[i]]),
                                               is_x_chr[i],
                                               is_female,
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
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' check_cross2(grav2)

check_cross2 <-
function(cross2)
{
    if(!is.cross2(cross2))
        stop('Input cross must have class "cross2"')

    result <- TRUE

    # required pieces
    #     crosstype
    #     geno
    #     gmap
    #     is_x_chr

    # sometimes required
    #     is_female
    #     cross_info
    #     linemap (required if nrow(pheno) != nrow(geno[[1]])
    #     founder_geno (required for some crosstypes)

    # optional pieces
    #     pheno
    #     covar
    #     phenocovar
    #     pmap

    crosstype <- cross2$crosstype
    if(is.null(crosstype)) {
        result <- FALSE
        warning("crosstype is missing")
    }
    else if(!check_crosstype(cross2$crosstype)) {
        result <- FALSE
        warning("Crosstype ", cross2$crosstype, " is not supported")
    }

    # check for required pieces
    geno <- cross2$geno
    if(is.null(geno)) {
        result <- FALSE
        warning("geno is missing")
    }

    gmap <- cross2$gmap
    if(is.null(gmap)) {
        result <- FALSE
        warning("gmap is missing")
    }

    # if either of those pieces is missing, we'll just return
    if(!result) return(result)

    # generally required pieces, but not worth stopping
    is_x_chr <- cross2$is_x_chr
    if(is.null(is_x_chr)) {
        result <- FALSE
        warning("is_x_chr is missing")
        is_x_chr <- rep(FALSE, length(geno))
        names(is_x_chr) <- names(geno)
        cross2$is_x_chr <- is_x_chr
    }
    else {
        chr <- names(cross2$geno)
        # double-check X chr in case that one of the chromosomes is called "X" or "x"
        if(any(chr=="X" | chr=="x")) {
          if(!all(is_x_chr[chr=="X" | chr=="x"]) ||
             any(is_x_chr[chr!="X" & chr!="x"]))
              warning("Suspicious indication of X chromosome:\n",
                      "    autosomes: ", paste(chr[!is_x_chr], collapse=" "), "\n",
                      "    X chr: ", paste(chr[is_x_chr], collapse=" "))
      }
    }

    if(!check_handle_x_chr(crosstype, any(is_x_chr))) {
        result <- FALSE
        warning("X chr not handled for cross type ", crosstype)
    }

    is_female <- cross2$is_female
    to_pass <- is_female
    if(is.null(to_pass)) to_pass <- logical(0)
    if(!check_is_female_vector(crosstype, to_pass, any(is_x_chr))) {
        result <- FALSE
        if(is.null(is_female)) warning("is_female is missing")
        else warning("is_female is misspecified")
    }

    cross_info <- cross2$cross_info
    to_pass <- cross_info
    if(is.null(to_pass)) to_pass <- matrix(0L,0,0)
    if(!check_crossinfo(crosstype, to_pass, any(is_x_chr))) {
        result <- FALSE
        if(is.null(cross_info)) warning("cross_info is missing")
        else warning("cross_info is misspecified")
    }

    # pheno
    pheno <- cross2$pheno
    if(!is.null(pheno)) { # pheno optional
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

    # check dimnames of geno
    for(i in seq(along=geno)) {
        if(is.null(rownames(geno[[i]]))) {
            result <- FALSE
            warning("missing geno rownames for chr ", names(geno)[i])
        }
        if(is.null(colnames(geno[[i]]))) {
            result <- FALSE
            warning("missing geno colnames for chr ", names(geno)[i])
        }
        if(nrow(geno[[i]]) != nrow(geno[[1]])) {
            result <- FALSE
            warning("nrow in geno for chr ", names(geno)[i],
                    " (", nrow(geno[[i]]), ") doesn't match that for chr ",
                    names(geno)[1], " (", nrow(geno[[1]]), ")")
        }
    }

    # compare geno to gmap
    if(length(geno) != length(gmap)) {
        result <- FALSE
        warning("length(geno) (", length(geno), ") != length(gmap) (", length(gmap), ")")
    }
    else {
        if(is.null(names(geno))) {
            result <- FALSE
            warning("names(geno) is missing")
        }
        if(any(names(geno) != names(gmap))) {
            result <- FALSE
            warning("names(geno) != names(gmap)")
        }
        for(i in seq(along=geno)) {
            if(ncol(geno[[i]]) != length(gmap[[i]])) {
                result <- FALSE
                warning("Mismatch between geno and gmap in no. markers on chr ", names(geno)[i])
            }
            else if(any(colnames(geno[[i]]) != names(gmap[[i]]))) {
                result <- FALSE
                warning("Mismatch in marker names between geno and gmap on chr ", names(geno)[i])
            }
        }
    }

    # Markers in gmap in non-decreasing order
    d <- vapply(gmap, function(x) ifelse(length(x)==1, 0, min(diff(x))), 0)
    if(any(d < 0)) {
        result <- FALSE
        warning("Markers not in order in gmap on chr ", paste(names(gmap)[d < 0], collapse=", "))
    }

    # compare geno to is_female
    is_female <- cross2$is_female
    if(!is.null(is_female)) {
        if(nrow(geno[[1]]) != length(is_female)) {
            result <- FALSE
            warning("length(is_female) (", length(is_female), ") != nrow(geno[[1]]) (", nrow(geno[[1]]), ")")
        }
        else if(any(names(is_female) != rownames(geno[[1]]))) {
            result <- FALSE
            warning("names(is_female) != rownames(geno[[1]])")
        }
        if(!is.logical(is_female)) {
            result <- FALSE
            warning("is_female is not logical")
        }
        if(any(is.na(is_female))) {
            result <- FALSE
            warning("is_female has ", sum(is.na(is_female)), " missing values")
        }
    }

    # compare geno to cross_info
    if(!is.null(cross_info)) {
        if(nrow(geno[[1]]) != nrow(cross_info)) {
            result <- FALSE
            warning("nrow(cross_info) (", nrow(cross_info), ") != nrow(geno[[1]]) (", nrow(geno[[1]]), ")")
        }
        else if(any(rownames(cross_info) != rownames(geno[[1]]))) {
            result <- FALSE
            warning("rownames(cross_info) != rownames(geno[[1]])")
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

    # compare individuals in geno and pheno
    if(!is.null(pheno)) {
        id <- find_common_ids(rownames(pheno), rownames(geno[[1]]))
        if(length(id) == 0)
            warning("No individuals in common between geno and pheno")
    }

    # check is_x_chr
    if(length(is_x_chr) != length(geno)) {
        result <- FALSE
        warning("length(is_x_chr) (", length(is_x_chr), ") != length(geno) (", length(geno), ")")
    }
    else if(any(names(is_x_chr) != names(geno))) {
        result <- FALSE
        warning("names(is_x_chr) != names(geno)")
    }

    # founder_geno
    founder_geno <- cross2$founder_geno
    if(need_founder_geno(crosstype)) {
        if(is.null(founder_geno)) {
            result <- FALSE
            warning("founder_geno not provided but needed.")
        }
        else {

            if(length(geno) != length(founder_geno)) {
                result <- FALSE
                warning("length(geno) (", length(geno), ") != length(founder_geno) (", length(founder_geno), ")")
            }
            else {
                if(any(names(geno) != names(founder_geno))) {
                    result <- FALSE
                    warning("names(geno) != names(founder_geno)")
                }
                for(i in seq(along=geno)) {
                    if(ncol(geno[[i]]) != ncol(founder_geno[[i]])) {
                        result <- FALSE
                        warning("Mismatch between geno and founder_geno in no. markers on chr ", names(geno)[i])
                    }
                    else if(any(colnames(geno[[i]]) != colnames(founder_geno[[i]]))) {
                        result <- FALSE
                        warning("Mismatch in marker names between geno and founder_geno on chr ", names(geno)[i])
                    }

                    # check values
                    if(!test_founder_geno_values(crosstype, founder_geno[[i]])) {
                        result <- FALSE
                    }
                }
            }
        }
    }

    # pmap
    pmap <- cross2$pmap
    if(!is.null(pmap)) { # pmap is optional
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

        # markers in order?
        d <- vapply(pmap, function(x) ifelse(length(x)==1, 0, min(diff(x))), 0)
        if(any(d < 0)) {
            result <- FALSE
            warning("Markers not in order in pmap on chr ", paste(names(pmap)[d < 0], collapse=", "))
        }
    }

    if(!result) { # check genotypes only if everything else is okay
        n_invalid <- count_invalid_genotypes(cross2)
        if(sum(n_invalid)>0)
            warning(sum(n_invalid), " genotypes in cross")
    }

    result
}

handle_null_crossinfo <-
    function(crossinfo, ids)
{
    n_ind <- length(ids)
    if(is.null(crossinfo)) {
        crossinfo <- matrix(0L, nrow=n_ind, ncol=0)
        rownames(crossinfo) <- ids
    }
    crossinfo
}

handle_null_isfemale <-
    function(isfemale, ids)
{
    n_ind <- length(ids)
    if(is.null(isfemale)) {
        isfemale <- rep(FALSE, n_ind)
        names(isfemale) <- ids
    }
    isfemale
}

handle_null_isxchr <-
    function(is_x_chr, chrnames)
{
    if(is.null(is_x_chr)) {
        is_x_chr <- rep(FALSE, length(chrnames))
        names(is_x_chr) <- chrnames
    }
    is_x_chr
}

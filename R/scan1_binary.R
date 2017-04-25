# stuff for genome scan with binary trait

# scan1 function taking nicely aligned data with no missing values
#
# Here genoprobs is a plain 3d array
scan1_binary_clean <-
    function(genoprobs, pheno, addcovar, intcovar,
             weights, add_intercept=TRUE, maxit, tol, qr_tol,
             intcovar_method)
{
    n <- nrow(pheno)
    if(add_intercept)
        addcovar <- cbind(rep(1,n), addcovar) # add intercept

    if(is.null(intcovar)) { # no interactive covariates

        if(is.null(weights)) { # no weights
            return( scan_binary_onechr(genoprobs, pheno, addcovar, maxit, tol, qr_tol) )
        } else { # weights included
            return( scan_binary_onechr_weighted(genoprobs, pheno, addcovar, weights, maxit, tol, qr_tol) )
        }

    } else { # interactive covariates
        # high- and low-memory versions of functions
        if(intcovar_method=="highmem")
            scanf <- c(scan_binary_onechr_intcovar_highmem,
                       scan_binary_onechr_intcovar_weighted_highmem)
        else
            scanf <- c(scan_binary_onechr_intcovar_lowmem,
                       scan_binary_onechr_intcovar_weighted_lowmem)

        if(is.null(weights)) { # no weights
            return( scanf[[1]](genoprobs, pheno, addcovar, intcovar, maxit, tol, qr_tol) )
        } else { # weights included
            return( scanf[[2]](genoprobs, pheno, addcovar, intcovar, weights, maxit, tol, qr_tol) )
        }

    }
}

# calculate null log likelihood, with nicely aligned data with no missing values
null_binary_clean <-
    function(pheno, addcovar, weights, add_intercept=TRUE, maxit, tol, qr_tol)
{
    n <- nrow(pheno)
    if(add_intercept)
        addcovar <- cbind(rep(1,n), addcovar) # add intercept

    if(is.null(weights) || length(weights)==0) { # no weights
        result <- apply(pheno, 2, function(phe) calc_ll_binreg(addcovar, phe, maxit, tol, qr_tol))
    } else { # weights included
        result <- apply(pheno, 2, function(phe) calc_ll_binreg_weighted(addcovar, phe, weights, maxit, tol, qr_tol))
    }

    as.numeric(result)
}

# check phenotype columns are binary (or at least between 0 and 1)
# if not, omit the non-binary columns
check_binary_pheno <-
    function(pheno, tol=1e-8)
{
    binary_col <- apply(pheno, 2, function(a) all(is.na(a) | (a >= 0 & a <= 1)))

    if(!all(binary_col)) {
        if(!any(binary_col))
            stop("Phenotypes are not binary (values in [0, 1])")
        if(!all(binary_col)) {
            warning("Not all phenotypes are binary. Dropping ",
                    sum(!binary_col), " columns")
            pheno <- pheno[,binary_col,drop=FALSE]
        }
    }

    round(pheno)
}

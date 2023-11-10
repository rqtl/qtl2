# find_dup_markers
#
#' Find markers with identical genotype data
#'
#' Identify sets of markers with identical genotype data.
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#'
#' @param chr Optional vector specifying which chromosomes to consider.
#' This may be a logical, numeric, or character string vector.
#'
#' @param exact_only If TRUE, look only for markers that have matching
#' genotypes and the same pattern of missing data; if FALSE, also look for
#' cases where the observed genotypes at one marker match those at
#' another, and where the first marker has missing genotype whenever the
#' genotype for the second marker is missing.
#'
#' @param adjacent_only If TRUE, look only for sets of markers that are
#' adjacent to each other.
#'
#' @return A list of marker names; each component is a set of markers whose
#' genotypes match one other marker, and the name of the component is the
#'  name of the marker that they match.
#'
#' @details
#'  If `exact.only=TRUE`, we look only for groups of markers whose
#'  pattern of missing data and observed genotypes match exactly.  One
#'  marker (chosen at random) is selected as the name of the group (in the
#'  output of the function).
#'
#'  If `exact.only=FALSE`, we look also for markers whose observed genotypes
#'  are contained in the observed genotypes of another marker.  We use a
#'  pair of nested loops, working from the markers with the most observed
#'  genotypes to the markers with the fewest observed genotypes.
#'
#' @export
#' @keywords utilities
#'
#' @seealso [drop_markers()], [drop_nullmarkers()], [reduce_markers()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' dup <- find_dup_markers(grav2)
#' grav2_nodup <- drop_markers(grav2, unlist(dup))

find_dup_markers <-
    function(cross, chr, exact_only=TRUE, adjacent_only=FALSE)
{
    if(!missing(chr)) cross <- subset(cross, chr=chr)

    if(!is.cross2(cross))
        stop('Input cross should be a "cross2" object.')

    # genotypes as one matrix; include founder_geno

    g <- do.call("cbind", cross$geno)
    if(!is.null(cross$founder_geno)) {
        g <- rbind(g, do.call("cbind", cross$founder_geno))
    }

    markers <- colnames(g)
    markerloc <- lapply(n_mar(cross), function(a) 1:a)
    if(length(markerloc) > 1) {
        for(j in 2:length(markerloc))
            markerloc[[j]] <- markerloc[[j]] + max(markerloc[[j-1]]) + 10
    }
    markerloc <- unlist(markerloc)

    if(exact_only) {
        g[is.na(g)] <- 0

        # genotype data patterns
        pat <- apply(g, 2, paste, collapse="")

        # table of unique values
        tab <- table(pat)

        # no duplicates; return
        if(all(tab == 1)) return(NULL)

        wh <- which(tab > 1)
        theloc <- themar <- vector("list", length(wh))
        for(i in seq(along=wh)) {
            themar[[i]] <- names(pat)[pat==names(tab)[wh[i]]]
            theloc[[i]] <- markerloc[pat==names(tab)[wh[i]]]
        }

        if(adjacent_only) {
            extraloc <- list()
            extramar <- list()
            for(i in seq(along=theloc)) {
                d <- diff(theloc[[i]]) # look for adjacency
                if(any(d>1)) { # split into adjacent groups
                    temp <- which(d>1)
                    first <- c(1, temp+1)
                    last <- c(temp, length(theloc[[i]]))
                    for(j in 2:length(first)) {
                        extraloc[[length(extraloc)+1]] <- theloc[[i]][first[j]:last[j]]
                        extramar[[length(extramar)+1]] <- themar[[i]][first[j]:last[j]]
                    }
                    themar[[i]] <- themar[[i]][first[1]:last[1]]
                    theloc[[i]] <- theloc[[i]][first[1]:last[1]]
                }
            }
            themar <- c(themar, extramar)
            theloc <- c(theloc, extraloc)

            nm <- sapply(themar, length)
            if(all(nm==1)) return(NULL) # nothing left
            themar <- themar[nm>1]
            theloc <- theloc[nm>1]
        }

        # order by first locus
        o <- order(sapply(theloc, min))
        themar <- themar[o]

        randompics <- sapply(themar, function(a) sample(length(a), 1))
        for(i in seq(along=themar)) {
            names(themar)[i] <- themar[[i]][randompics[i]]
            themar[[i]] <- themar[[i]][-randompics[i]]
        }

    }
    else {
        themar <- NULL

        ntyp <- n_typed(cross, "marker")
        o <- order(ntyp, decreasing=TRUE)

        g[is.na(g)] <- 0
        result <- .find_dup_markers_notexact(g, o, markerloc, adjacent_only)

        if(all(result==0)) return(NULL)
        u <- unique(result[result != 0])
        themar <- vector("list", length(u))
        names(themar) <- colnames(g)[u]
        for(i in seq(along=themar))
            themar[[i]] <- colnames(g)[result==u[i]]
    }

    themar
}

#' Compare two sets of genotype probabilities
#'
#' Compare two sets of genotype probabilities for one individual on a single chromosome.
#'
#' @md
#'
#' @param probs1 Genotype probabilities (as produced by [qtl2geno::calc_genoprob()])
#' or allele dosages (as produced by [qtl2geno::genoprob_to_alleleprob()]).
#' @param probs2 A second set of genotype probabilities, just like `probs1`.
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosome; a single character string.
#' @param minprob Minimum probability for inferring genotypes (passed to [maxmarg()]).
#' @param minmarkers Minimum number of markers in results.
#' @param minwidth Minimum width in results.
#' @param annotate If TRUE, add some annotations to the `geno1` and
#' `geno2` columns to indicate, where they differ, which one
#' matches what appears to be the best genotype. (`*` = matches
#' the best genotype; `-` = lower match).
#'
#' @details
#' The function does the following:
#' - Reduce the probabilities to a set of common locations that also appear in `cross`.
#' - Use [maxmarg()] to infer the genotype at every position using each set of probabilities.
#' - Identify intervals where the two inferred genotypes are constant.
#' - Within each segment, compare the observed SNP genotypes to the founders' genotypes.
#'
#' @return
#' A data frame with each row corresponding to an interval over which
#' `probs1` and `probs2` each have a fixed inferred genotype. Columns
#' include the two inferred genotypes, the start and end points and
#' width of the interval, and when founder genotypes are in `cross`,
#' the proportions of SNPs where the individual matches each possible
#' genotypes.
#'
#' @seealso `plot_genoprobcomp()` in [R/qtl2plot](https://github.com/rqtl/qtl2plot)
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' iron <- iron[1,"2"]   # subset to first individual on chr 2
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # in presence of a genotyping error, how much does error_prob matter?
#' iron$geno[[1]][1,3] <- 3
#' pr_e <- calc_genoprob(iron, map, error_prob=0.002)
#' pr_ne <- calc_genoprob(iron, map, error_prob=1e-15)
#'
#' compare_genoprob(pr_e, pr_ne, iron, minmarkers=1, minprob=0.5)
#'
#' @export
compare_genoprob <-
    function(probs1, probs2, cross, ind=1, chr=NULL,
             minprob=0.95, minmarkers=10, minwidth=0,
             annotate=FALSE)
{
    # check inputs
    if(is.null(chr)) chr <- names(probs1)[1]
    if(length(chr) > 1) {
        warning("chr should have length 1; using the first value")
        chr <- chr[1]
    }
    if(!(chr %in% names(probs1) && chr %in% names(probs2)))
        stop("chr ", chr, " not found in both probs")
    if(!(chr %in% chr_names(cross))) stop("chr ", chr, " not found in cross")
    if(length(ind) > 1) {
        warning("ind should have length 1; using the first value")
        ind <- ind[1]
    }

    # pull out individual's probs; make it positions x probs
    if(is.character(ind) && !(ind %in% rownames(probs1[[chr]]) && ind %in% rownames(probs2[[chr]])))
        stop("ind ", ind, " not found in both sets of probs")
    if(is.numeric(ind)) {
        if(nrow(probs1[[chr]]) != nrow(probs2[[chr]]))
            stop("If ind is numeric, must have nrow(probs1) == nrow(probs2)")
        if(ind < 1 || ind > nrow(probs1[[chr]]))
            stop("ind ", ind, " should be in the range [1, ", nrow(probs1[[chr]]), "]")
        ind <- rownames(probs1[[1]])[ind]
    }

    # make sure they have the same genotypes
    if(ncol(probs1[[chr]]) != ncol(probs2[[chr]]) || !all(colnames(probs1[[chr]]) == colnames(probs2[[chr]])))
        stop("Need ncol(probs1) == ncol(probs2) and same colnames")
    geno_names <- colnames(probs1[[chr]])

    # infer genotypes
    ginf1 <- maxmarg(probs1[ind,chr], minprob=minprob)[[1]]
    ginf2 <- maxmarg(probs2[ind,chr], minprob=minprob)[[1]]

    if("pmap" %in% names(cross)) {
        map <- cross$pmap
        map_units <- "Mbp"
    } else {
        map <- cross$gmap
        map_units <- "cM"
    }
    gmap <- cross$gmap
    pmap <- cross$pmap

    # pull out selected chromosomes
    probs1 <- probs1[[chr]]
    probs2 <- probs2[[chr]]
    map <- map[[chr]]

    if(!is.null(gmap)) gmap <- gmap[[chr]]
    if(!is.null(pmap)) pmap <- pmap[[chr]]
    geno <- cross$geno[[chr]]

    # find common markers
    mar1 <- dimnames(probs1)[[3]]
    mar2 <- dimnames(probs2)[[3]]
    marm <- names(map)
    mar <- mar1[mar1 %in% mar2 & mar1 %in% marm]
    if(length(mar) < 2) stop("<2 markers in common")
    map <- map[mar]
    if(!is.null(gmap)) gmap <- gmap[mar]
    if(!is.null(pmap)) pmap <- pmap[mar]
    geno <- geno[,mar,drop=FALSE]
    ginf1 <- ginf1[,mar]
    ginf2 <- ginf2[,mar]

    if(!is.null(cross$founder_geno)) {
        founder_geno <- cross$founder_geno[[chr]]
        founder_geno <- founder_geno[,mar,drop=FALSE]
    }

    # now pull out individual from the genotypes
    geno <- geno[ind,]

    # find segments with constant (and non-missing calls)
    ginf <- paste(geno_names[ginf1], geno_names[ginf2], sep=":")
    ginf_num <- (ginf1-1) * ncol(probs2) + ginf2
    ginf_num[is.na(ginf_num)] <- -1
    jump <- which(diff(ginf_num) != 0)

    value_num <- c(ginf_num[1], ginf_num[jump+1])
    value <- c(ginf[1], ginf[jump+1])
    left_i <- c(1, jump+1)
    right_i <- c(jump, length(map))
    n_markers=right_i - left_i + 1
    left <- c(map[1], map[jump+1])
    right <- c(map[jump], map[length(map)])
    width <- right - left
    if(!is.null(gmap)) gwidth <- c(gmap[jump], gmap[length(gmap)]) - c(gmap[1], gmap[jump+1])
    if(!is.null(pmap)) pwidth <- c(pmap[jump], pmap[length(pmap)]) - c(pmap[1], pmap[jump+1])

    result <- data.frame(geno1=geno_names[c(ginf1[1], ginf1[jump+1])],
                         geno2=geno_names[c(ginf2[1], ginf2[jump+1])],
                         left_index = left_i,
                         right_index = right_i,
                         n_markers = n_markers,
                         left = left,
                         right = right,
                         width_Mbp = width,
                         width_cM = width,
                         stringsAsFactors=FALSE)

    # add units to left/right
    colnames(result)[6:7] <- paste(colnames(result)[6:7], map_units, sep="_")

    col2drop <- NULL
    if(!is.null(gmap)) result$width_cM <- gwidth
    else col2drop <- 9
    if(!is.null(pmap)) result$width_Mbp <- pwidth
    else col2drop <- c(col2drop, 8)
    if(!is.null(col2drop)) result <- result[,-col2drop,drop=FALSE]

    result <- result[value_num != -1 & n_markers >= minmarkers & width >= minwidth, , drop=FALSE]
    if(nrow(result) > 0) rownames(result) <- 1:nrow(result)

    ### if no founder_geno in cross, return what we've got
    if(is.null(cross$founder_geno)) return(result)

    # construct F1-type individuals
    f1 <- infer_f1_geno(founder_geno)

    # adjust genotype labels
    if("is_x_chr" %in% names(cross) && cross$is_x_chr[chr] &&
       "is_female" %in% names(cross) && !cross$is_female[ind])
        rownames(founder_geno) <- paste0(rownames(founder_geno), "Y")
    else
        rownames(founder_geno) <- paste0(rownames(founder_geno), rownames(founder_geno))

    # add F1s to founders
    founder_geno <- rbind(founder_geno, f1)

    # calculate proportion SNPs that match each possible genotype
    p_match <- matrix(nrow=nrow(result), ncol=nrow(founder_geno))
    colnames(p_match) <- rownames(founder_geno)
    for(k in 1:nrow(result)) {
        gsub <- geno[result$left_index[k]:result$right_index[k]]
        p_match[k,] <- apply(founder_geno[,result$left_index[k]:result$right_index[k],drop=FALSE], 1,
                             function(a) sum(a != 0 & gsub != 0 & gsub == a) / sum(a != 0 & gsub != 0))
    }

    # paste detailed results onto the end
    result <- cbind(result,
                    match1=rep(0, nrow(result)),
                    match2=rep(0, nrow(result)),
                    match_best=rep(0, nrow(result)),
                    what_best=rep("", nrow(result)),
                    match_next=rep(0, nrow(result)),
                    what_next=rep("", nrow(result)),
                    stringsAsFactors=FALSE)

    for(i in 1:nrow(result)) {
        result$match1[i] <- p_match[i,result$geno1[i]]
        result$match2[i] <- p_match[i,result$geno2[i]]
        z <- sort(p_match[i,], decreasing=TRUE)
        result$match_best[i] <- z[1]
        result$match_next[i] <- z[2]
        result$what_best[i] <- names(z)[1]
        result$what_next[i] <- names(z)[2]
    }

    if(annotate) result <- compare_genoprob_add_annotation(result)

    # drop the index columns
    result <- result[,-grep("index", colnames(result))]

    # add matrix of proportion of matches for each possible genotype, as an attribute
    rownames(p_match) <- rownames(result)
    attr(result, "prop_match") <- p_match

    result
}

# from founders, infer all F1 SNP genotypes
infer_f1_geno <-
    function(founders)
{
    n <- nrow(founders)
    n_f1 <- choose(n, 2)
    result <- matrix(0, nrow=n_f1, ncol=ncol(founders))
    dimnames(result) <- list(1:n_f1, colnames(founders))

    f <- rownames(founders)
    k <- 1
    for(i in 1:(n-1)) {
        for(j in (i+1):n) {
            rownames(result)[k] <- paste0(f[i], f[j])

            result[k,founders[i,]==1 & founders[j,]==1] <- 1
            result[k,founders[i,]==3 & founders[j,]==3] <- 3
            result[k,founders[i,]==1 & founders[j,]==3] <- 2
            result[k,founders[i,]==3 & founders[j,]==1] <- 2

            k <- k + 1
        }
    }

    result
}

compare_genoprob_add_annotation <-
    function(x)
{
    if("match_best" %in% names(x)) { # add some annotations
        for(i in 1:nrow(x)) {
            if(x[i,"geno1"]==x[i,"geno2"]) {
                x[i,"geno1"] <- paste0(x[i,"geno1"], " ")
                x[i,"geno2"] <- paste0(x[i,"geno2"], " ")
            } else {
                if(x$match1[i] >= x$match_best[i])
                    x[i,"geno1"] <- paste0(x[i,"geno1"], "*")
                else
                    x[i,"geno1"] <- paste0(x[i,"geno1"], "-")
                if(x$match2[i] >= x$match_best[i])
                    x[i,"geno2"] <- paste0(x[i,"geno2"], "*")
                else
                    x[i,"geno2"] <- paste0(x[i,"geno2"], "-")
            }
        }
    }

    x
}

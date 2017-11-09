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
#' @param map_type Whether to use the physical or genetic map for
#' widths and positions (if selected map is not available, use
#' whichever one is available.
#' @param minprob Minimum probability for inferring genotypes (passed to [maxmarg()]).
#' @param minmarkers Minimum number of markers in results.
#' @param minwidth Minimum width in results.
#'
#' @details
#' The function does the following:
#' - Reduce the probabilities to a set of common locations that also appear in `cross`.
#' - Use [maxmarg()] to infer the genotype at every position using each set of probabilities.
#' - Identify intervals where the two inferred genotypes are constant.
#' - Within each segment, compare the observed SNP genotypes to the founders' genotypes.
#'
#' @return
#' An object of class `"compare_genoprob`", which is a data frame with
#' each row corresponding to an interval over which `probs1` and
#' `probs2` each have a fixed inferred genotype. Columns include the
#' two inferred genotypes, the start and end points and width of the
#' interval, and when founder genotypes are in `cross`, the
#' proportions of SNPs where the individual matches each possible
#' genotypes.
#'
#' @seealso [qtl2plot::plot_genoprobcomp()]
#'
#' @export
compare_genoprob <-
    function(probs1, probs2, cross, ind=1, chr=NULL, minprob=0.95,
             minmarkers=10)
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
    }

    # make sure they have the same genotypes
    if(ncol(probs1[[chr]]) != ncol(probs2[[chr]]) || !all(colnames(probs1[[chr]]) == colnames(probs2[[chr]])))
        stop("Need ncol(probs1) == ncol(probs2) and same colnames")
    geno_names <- colnames(probs1[[chr]])

    # infer genotypes
    ginf1 <- maxmarg(probs1[ind,chr], minprob=minprob)[[1]]
    ginf2 <- maxmarg(probs2[ind,chr], minprob=minprob)[[1]]

    # pull out selected chromosomes
    probs1 <- probs1[[chr]]
    probs2 <- probs2[[chr]]
    pmap <- cross$pmap[[chr]]
    gmap <- cross$gmap[[chr]]
    geno <- cross$geno[[chr]]
    founder_geno <- cross$founder_geno[[chr]]

    # find common markers
    mar1 <- dimnames(probs1)[[3]]
    mar2 <- dimnames(probs2)[[3]]
    marm <- names(pmap)
    mar <- mar1[mar1 %in% mar2 & mar1 %in% marm]
    if(length(mar) < 2) stop("<2 markers in common")
    pmap <- pmap[mar]
    gmap <- gmap[mar]
    geno <- geno[,mar,drop=FALSE]
    founder_geno <- founder_geno[,mar,drop=FALSE]
    ginf1 <- ginf1[,mar]
    ginf2 <- ginf2[,mar]

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
    right_i <- c(jump, length(pmap))
    n_markers=right_i - left_i + 1
    left_Mbp <- c(pmap[1], pmap[jump+1])
    right_Mbp <- c(pmap[jump], pmap[length(pmap)])
    width_Mbp <- right_Mbp - left_Mbp
    width_cM <- c(gmap[jump], gmap[length(gmap)]) - c(gmap[1], gmap[jump+1])

    result <- data.frame(geno1=geno_names[c(ginf1[1], ginf1[jump+1])],
                         geno2=geno_names[c(ginf2[1], ginf2[jump+1])],
                         left_index = left_i,
                         right_index = right_i,
                         n_markers = n_markers,
                         left_Mbp = left_Mbp,
                         right_Mbp = right_Mbp,
                         width_Mbp = width_Mbp,
                         width_cM = width_cM,
                         stringsAsFactors=FALSE)

    result <- result[value_num != -1,]

    # look at corresponding markers
    f1 <- infer_f1_geno(founder_geno)
    if(cross$is_x_chr[chr] && !cross$is_female[ind])
        rownames(founder_geno) <- paste0(rownames(founder_geno), "Y")
    else
        rownames(founder_geno) <- paste0(rownames(founder_geno), rownames(founder_geno))
    founder_geno <- rbind(founder_geno, f1)


    p_match <- matrix(nrow=nrow(result), ncol=nrow(founder_geno))
    colnames(p_match) <- rownames(founder_geno)

    for(k in 1:nrow(result)) {
        gsub <- geno[result$left_index[k]:result$right_index[k]]
        p_match[k,] <- apply(founder_geno[,result$left_index[k]:result$right_index[k],drop=FALSE], 1,
                             function(a) sum(a != 0 & gsub != 0 & gsub == a) / sum(a != 0 & gsub != 0))
    }

    result <- cbind(result,
                    match1=rep(0, nrow(result)),
                    match2=rep(0, nrow(result)),
                    match_best=rep(0, nrow(result)),
                    what_best=rep("", nrow(result)),
                    match_next=rep(0, nrow(result)),
                    what_next=rep("", nrow(result)),
                    stringsAsFactors=FALSE)
    rownames(result) <- 1:nrow(result)

    for(i in 1:nrow(result)) {
        result$match1[i] <- p_match[i,result$geno1[i]]
        result$match2[i] <- p_match[i,result$geno2[i]]
        z <- sort(p_match[i,], decreasing=TRUE)
        result$match_best[i] <- z[1]
        result$match_next[i] <- z[2]
        result$what_best[i] <- names(z)[1]
        result$what_next[i] <- names(z)[2]

        if(result[i,1]==result[i,2]) {
            result[i,1] <- paste0(result[i,1], " ")
            result[i,2] <- paste0(result[i,2], " ")
        } else {
            if(result$match1[i] >= result$match_best[i])
                result[i,1] <- paste0(result[i,1], "*")
            else
                result[i,1] <- paste0(result[i,1], "-")
            if(result$match2[i] >= result$match_best[i])
                result[i,2] <- paste0(result[i,2], "*")
            else
                result[i,2] <- paste0(result[i,2], "-")
        }
    }

    result[,c(1:2,5:ncol(result))]
}


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

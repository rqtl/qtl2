# snpprob_from_cross: genoprob_to_snpprob where snpinfo is a cross2 object
# ...construct snpinfo for all markers with no missing founder genotypes
# ...apply genoprob_to_snpprob with that

snpprob_from_cross <-
    function(genoprobs, cross)
{
    if(!("founder_geno" %in% names(cross)))
        stop("cross doesn't contain founder_geno")

    ### calculate snp info from cross
    # drop markers that are missing any founder genotypes
    cross <- drop_markers(cross, unlist(lapply(cross$founder_geno, function(a) colnames(a)[colSums(a==0)>0])))

    # get physical map, if available; otherwise genetic map
    map <- cross$pmap
    if(is.null(map)) map <- cross$gmap

    sdp <- calc_sdp(t(do.call("cbind", cross$founder_geno)))

    snpinfo <- data.frame(chr=rep(chr_names(cross), n_mar(cross)),
                          pos=unlist(map),
                          sdp=sdp,
                          snp=names(sdp))

    for(i in seq_along(genoprobs)) {
        genoprobs[[i]] <- genoprobs[[i]][,,colnames(cross$geno[[i]]),drop=FALSE]
    }

    snpinfo <- index_snps(map, snpinfo)

    genoprob_to_snpprob(genoprobs, snpinfo)
}

# reduce snpinfo to the indexed snps
reduce_to_indexed_snps <- function(snpinfo)
{
    if(length(unique(snpinfo$chr)) > 1) { # more than one chromosome
        # rows for each chromosome
        rows_spl <- split(1:nrow(snpinfo), factor(snpinfo$chr, unique(snpinfo$chr)))

        # apply this function to one chromosome at a time
        result <- vector("list", length(rows_spl))
        for(chr in seq_along(rows_spl)) {
            result[[chr]] <- reduce_to_indexed_snps(snpinfo[rows_spl[[chr]],,drop=FALSE])
        }

        # re-combine by rows
        return( do.call("rbind", result) )
    }

    # now working with just one chromosome
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo))) {
        stop("index seems messed up for chromosome ", snpinfo$chr[1])
    }

    result <- snpinfo[unique(snpinfo$index), , drop=FALSE]
    result$index <- seq_len(nrow(result))

    result
}

# function for subsetting by chromosome
#
# chr = logical, numeric, or character strings
#       (numeric converted to character strings)
# all_chr = vector of character strings with all chromosome names
#
# result = vector of character strings (the selected ones)
subset_chr <-
    function(chr, all_chr)
{
    if(is.logical(chr)) {
        if(length(chr) != length(all_chr))
            stop("chr is logical but length [", length(chr), "] != n_chr is x [",
                 length(all_chr), "]")
        chr <- all_chr[chr]
    } else {
        chr <- as.character(chr)

        # look for negatives; turn to positives
        if(any(grepl("^\\-", chr))) {
            if(!all(grepl("^\\-", chr)))
                stop("Can't mix negative and positive chr subscripts")
            chr <- sub("^\\-", "", chr)
            if(!all(chr %in% all_chr)) {
                if(!any(chr %in% all_chr))
                    stop("None of the chr found in the cross object")
                warning("Some chr not found: ", paste(chr[!(chr %in% all_chr)], collapse=", "))
                chr <- chr[chr %in% all_chr]
            }
            chr <- all_chr[!(all_chr %in% chr)]
        }

        if(!all(chr %in% all_chr)) {
            if(!any(chr %in% all_chr)) stop("None of the chromosomes in cross")
            warning("Some chr not in cross: ", paste(chr[!(chr %in% all_chr)], collapse=", "))
            chr <- chr[chr %in% all_chr]
        }
    }

    chr
}

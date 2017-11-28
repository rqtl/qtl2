# assign some allele codes based on how many you need
assign_allele_codes <-
    function(num_alleles, genotypes=NULL)
{
    if(!is.null(genotypes)) { # genotypes provided; find unique letters
        alleles <- sort(unique(unlist(strsplit(genotypes, ""))))
        if(num_alleles <= length(alleles))
            return(alleles[1:num_alleles])
    }

    # otherwise, use capital letters
    if(num_alleles <= length(LETTERS))
        return(LETTERS[1:num_alleles])

    # ack! need more than 26. Use 2-letter codes.
    alleles <- paste0(LETTERS, rep(LETTERS, each=length(LETTERS)))
    if(num_alleles <= length(alleles))
        return(alleles[1:num_alleles])

    # ack! need more than 26^2. Just use numbers
    as.character(1:num_alleles)
}

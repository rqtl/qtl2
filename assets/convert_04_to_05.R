# R/qtl2 version 0.5 includes major changes in data structures
#
# These functions are for converting objects from the formats
# for R/qtl2 version0.4 to that for R/qtl2 version 0.5.

# convert output of calc_genoprob
convert_probs <-
    function(probs)
{
    attr(probs$probs, "is_x_chr") <- probs$is_x_chr
    attr(probs$probs, "crosstype") <- probs$crosstype
    attr(probs$probs, "alleles") <- probs$alleles
    attr(probs$probs, "alleleprobs") <- probs$alleleprobs
    class(probs$probs) <- class(probs)

    probs$probs
}


# convert output of scan1/scan1coef
convert_scan1 <-
    function(scan1_output)
{
    if("lod" %in% names(scan1_output)) {
        result <- scan1_output$lod
    } else if("coef" %in% names(scan1_output)) {
        result <- scan1_output$coef
    } else {
        stop("Can't find lod or coef in input")
    }

    attr(result, "sample_size") <- scan1_output$sample_size
    attr(result, "hsq") <- scan1_output$hsq
    attr(result, "SE") <- scan1_output$SE
    class(result) <- class(scan1_output)

    result
}

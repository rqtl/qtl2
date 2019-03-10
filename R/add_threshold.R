#' Add thresholds to genome scan plot
#'
#' Draw line segments at significance thresholds for a genome scan plot
#'
#' @param map Marker map used in the genome scan plot
#' @param thresholdA Autosomal threshold. Numeric or a list. (If a
#' list, the `"A"` component is taken to be `thresholdA` and the
#' `"X"` component is taken to be `thresholdX`.)
#' @param thresholdX X chromosome threshold (if missing, assumed to be the same as `thresholdA`)
#' @param chr Chromosomes that were included in the plot
#' @param gap Gap between chromosomes in the plot. Default is 1\% of the total genome length.
#' @param ... Additional arguments passed to [graphics::segments()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,c(2,16,"X")]}
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#' Xcovar <- get_x_covar(iron)
#' out <- scan1(probs, iron$pheno[,1], Xcovar=Xcovar)
#' \dontrun{operm <- scan1perm(probs, iron$pheno[,1], addcovar=Xcovar,
#'                    n_perm=1000, perm_Xsp=TRUE, chr_lengths=chr_lengths(map))}
#' \dontshow{operm <- scan1perm(probs, iron$pheno[,1], addcovar=Xcovar,
#'                    n_perm=10, perm_Xsp=TRUE, chr_lengths=chr_lengths(map))}
#'
#' plot(out, map)
#' add_threshold(map, summary(operm), col="violetred", lty=2)
#' @importFrom graphics segments
#' @importFrom stats setNames
#' @export
add_threshold <-
    function(map, thresholdA, thresholdX=NULL, chr=NULL, gap=NULL, ...)
{
    if(is.null(gap)) gap <- sum(chr_lengths(map))/100

    if(is.list(thresholdA)) {
        thresholdX <- thresholdA$X
        thresholdA <- thresholdA$A
    }
    if(is.null(thresholdX)) thresholdX <- thresholdA

    # which chromosomes are X chr?
    is_x_chr <- attr(map, "is_x_chr")
    if(is.null(is_x_chr)) {
        is_x_chr <- setNames(names(map) %in% c("X", "x"), names(map))
    }

    # subset map and is_x_chr
    if(!is.null(chr)) {
        map <- map[chr]
        is_x_chr <- is_x_chr[chr]
    }

    if(length(map) == 1) { # one chromosome
        abline(h=ifelse(is_x_chr, thresholdX, thresholdA), ...)
    }
    else {
        start <- vapply(map, min, na.rm=TRUE, 0) - gap/2
        end <- vapply(map, max, na.rm=TRUE, 0) + gap/2

        for(chr in names(map)) {
            h <- ifelse(is_x_chr[chr], thresholdX, thresholdA)
            segments(xpos_scan1(map, gap=gap, thechr=chr, thepos=start[chr]), h,
                     xpos_scan1(map, gap=gap, thechr=chr, thepos=end[chr]), h, ...)
        }
    }
}

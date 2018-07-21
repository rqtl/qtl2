#' Clean scan1 output
#'
#' Clean scan1 output by replacing negative values with NA and remove
#' rows where all values are NA.
#'
#' @md
#'
#' @param object Output of [scan1()].
#' @param ... Ignored at present
#'
#' @return The input object with negative values replaced with NAs and then rows with all NAs removed.
#'
#' @export
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,"2"]}
#' pr <- calc_genoprob(iron)
#' out <- scan1(pr, iron$pheno)
#'
#' out <- clean(out)

clean_scan1 <-
    function(object, ...)
{
    cl <- class(object)
    at <- attributes(object)
    if(!("scan1" %in% cl)) {
        stop("object is not of class scan1")
    }

    object[object < 0] <- NA

    omit <- apply(object, 1, function(a) all(is.na(a)))
    if(any(omit)) {
        object <- object[!omit, , drop=FALSE]
        class(object) <- cl

        # retain attributes
        for(a in names(at)) {
            if(!(a %in% c("dim", "dimnames", "class"))) {
                attr(object, a) <- at[[a]]
            }
        }
    }

    object
}


#' @rdname clean_scan1
#' @export
clean.scan1 <- clean_scan1

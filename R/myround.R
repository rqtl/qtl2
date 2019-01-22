# Round a number, preserving extra 0's. (taken from R/broman, https://github.com/kbroman/broman)
myround <-
    function(x, digits=1)
{
    if(digits < 1)
        stop("This is intended for the case digits >= 1.")

    if(length(digits) > 1) {
        digits <- digits[1]
        warning("Using only digits[1]")
    }

    tmp <- sprintf(paste("%.", digits, "f", sep=""), x)

    # deal with "-0.00" case
    zero <- paste0("0.", paste(rep("0", digits), collapse=""))
    tmp[tmp == paste0("-", zero)] <- zero

    tmp
}

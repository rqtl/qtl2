# check if a crosstype is supported
#   if give_error=TRUE, problem results in error
#   if give_error=FALSE, returns TRUE if supported; FALSE if not
check_crosstype <-
function(crosstype, give_error=TRUE)
{
    z <- try(.crosstype_supported(crosstype), silent=TRUE)

    if("try-error" %in% class(z) || !z) {
        if(give_error)
            stop("Cross type ", crosstype, " not yet supported.")
        return(invisible(FALSE))
    }

    invisible(TRUE) # if no error, return TRUE
}

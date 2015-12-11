# utility for grabbing "..." args
#
# dotargs = list(...) from function call
# argname = character string of the argument to grab
# default = default value for argument
# values  = optional vector of character strings with possible values
grab_dots <-
    function(dotargs, argname, default, values)
{
    if(argname %in% names(dotargs)) {
        arg <- dotargs[[argname]]
        if(!missing(values) && !is.null(values) && !(arg %in% values)) {
            warning(argname, ' "', arg, '" not valid; using "',
                    default, '".')
            arg <- default
        }
    }
    else arg <- default
    arg
}

# give warning if there are extra named arguments
#
# dotargs = list(...) from function call
# expected = vector of character strings of the names of possible args anticipated
check_extra_dots <-
    function(dotargs, expected)
{
    args <- names(dotargs)
    if(!all(args %in% expected)) {
        extra <- args[!(args %in% expected)]
        if(length(extra) == 1)
            warning("Extra argument ignored: ", extra)
        else
            warning("Extra arguments ignored: ",
                    paste(extra, collapse=" "))
    }
}

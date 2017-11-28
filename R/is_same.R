# check that two objects are the same
is_same <-
    function(a, b)
{
    if(is.null(a) && is.null(b)) return(TRUE)
    if(is.null(a) || is.null(b)) return(FALSE)

    if(is.list(a) && is.list(b)) {
        if(length(a) != length(b)) return(FALSE)
        for(i in seq_along(a))
            if(!is_same(a[[i]], b[[i]])) return(FALSE)
        return(TRUE)
    }
    if(is.list(a) || is.list(b)) return(FALSE)

    attra <- attributes(a)
    attrb <- attributes(b)

    if(length(a) != length(b)) return(FALSE)
    if(!all((is.na(a) & is.na(b)) |
            !is.na(a) & !is.na(b) & a==b)) return(FALSE)

    if(is.null(names(attra))) {
        if(!is.null(names(attrb))) return(FALSE)
    }
    else {
        if(length(attra) != length(attrb)) return(FALSE)
        attra <- attrb[order(names(attra))]
        attrb <- attrb[order(names(attrb))]
        if(!all(names(attra) == names(attrb))) return(FALSE)
        nam <- names(attra)
        nam <- nam[nam != "names"]
        for(i in names(attra)) {
            if(!is_same(attra[[i]], attrb[[i]])) return(FALSE)
        }
    }

    TRUE
}

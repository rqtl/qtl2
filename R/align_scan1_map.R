# align scan1 output with a map
align_scan1_map <-
    function(scan1_output, map)
{
    if(!is.list(map)) stop("map should be a list")

    scan1_names <- rownames(scan1_output)
    map_names <- map2markernames(map)
    map_chr <- map2chr(map)

    # perfectly fine
    if(length(scan1_names) == length(map_names) &&
       all(scan1_names == map_names))
        return(list(scan1_output=scan1_output, map=map))

    if(!any(scan1_names %in% map_names))
        stop("scan1 output and map have no markers in common")

    # subset scan1_output to markers in the map
    #    and use order as in map
    scan1_attr <- attributes(scan1_output)
    names_ordered <- map_names[map_names %in% scan1_names]
    scan1_output <- scan1_output[names_ordered,,drop=FALSE]
    for(a in c("SE", "sample_size", "hsq", "class"))
        attr(scan1_output, a) <- scan1_attr[[a]]
    scan1_names <- rownames(scan1_output)

    # subset map to markers in scan1_output
    if(any(!(map_names %in% scan1_names))) {
        map_attr <- attributes(map)
        keep <- map_names %in% scan1_names
        uchr <- unique(map_chr[keep])
        map <- map[uchr]
        if("is_x_chr" %in% names(map_attr))
            map_attr$is_x_chr <- map_attr$is_x_chr[uchr]
        for(i in seq_along(map))
            map[[i]] <- map[[i]][names(map[[i]]) %in% scan1_names]
        for(a in c("is_x_chr", "class"))
            attr(map, a) <- map_attr[[a]]
    }

    list(scan1_output=scan1_output, map=map)
}

# grab marker names as a vector
map2markernames <-
    function(map)
{
    nam <- unlist(lapply(map, names))
    names(nam) <- NULL
    nam
}

# grab chromosome IDs as a vector
map2chr <-
    function(map)
{
    chr <- rep(names(map), vapply(map, length, 0))
    names(chr) <- map2markernames(map)
    chr
}

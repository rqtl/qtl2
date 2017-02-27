# create_marker_index
#
# find marker names in map and create index of locations
# (indexes start at 0)

create_marker_index <-
    function(marker_names, map)
{
    if(length(marker_names) != length(map))
        stop("Different numbers of chromosomes")
    if(!all(names(marker_names) == names(map)))
        stop("Different chromosome names")

    index <- vector("list", length(map))
    for(i in seq(along=marker_names)) {
        if(!all(marker_names[[i]] %in% names(map[[i]])))
            stop("Some markers not found in map")
        index[[i]] <- match(names(map[[i]]), marker_names[[i]]) - 1
        index[[i]][is.na(index[[i]])] <- -1
    }
    names(index) <- names(map)

    index
}

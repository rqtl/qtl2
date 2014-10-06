# create_pseudomarker_map

# pseudomarker map as grid, ignoring markers
#
# map = vector of marker positions; names = marker names
# step = distance between pseudomarkers
# off_end = amount to go beyond the terminal markers
# tol = tolerance for determining whether a marker hits a pseudomarker
# pmar_stem = leading part of name of pseudomarkers
create_pseudomarker_map_grid <-
function(map, step, off_end=0, tol=0.01, pmar_stem="loc")
{
    if(any(is.na(map))) stop("map values can't be missing")
    if(step==0) {
        attr(map, "index") <- seq(along=map)
        attr(map, "stepwidth") <- "grid"
        attr(map, "step") <- step
        attr(map, "off_end") <- off_end
        return(map)
    }
    if(step < 0) stop("step should be >= 0")
    if(tol >= step) stop("tol should be << step")
    if(off_end < 0) stop("off_end should be >= 0")
    if(tol < 0) stop("tol should be >= 0")

    # locations of pseudomarkers
    pmar <- seq(min(map)-off_end, max(map)+off_end, by=step)

    # for each pseudomarker, distance to nearest marker
    d <- abs(outer(map, pmar, "-"))
    mind <- apply(d, 2, min)

    # omit pseudomarkers that hit a marker
    to_omit <- (mind < tol)
    pmar <- pmar[!to_omit]
    if(length(pmar) == 0) {
        attr(map, "index") <- seq(along=map)
        attr(map, "grid") <- rep(TRUE, length(map))
        attr(map, "stepwidth") <- "grid"
        attr(map, "step") <- step
        attr(map, "off_end") <- off_end
        return(map)
    }

    # markers that are on pseudomarker grid
    grid <- rep(FALSE, length(map))
    if(any(to_omit))
        grid[apply(d[,to_omit,drop=FALSE], 2, which.min)] <- TRUE

    # no. digits to use in name
    digits <- ceiling(-log10(step))
    digits <- ifelse(digits < 0, 0, digits)

    names(pmar) <- paste0(pmar_stem, round(pmar, digits))

    # index of markers vs pseudomarkers
    index <- c(seq(along=map), rep(0, length(pmar)))
    # logical vector indicating grid
    grid <- c(grid, rep(TRUE, length(pmar)))

    # add pseudomarkers to markers
    map <- c(map, pmar)

    # sort the map
    o <- order(map)
    map <- map[o]
    attr(map, "stepwidth") <- "grid"
    attr(map, "index") <- index[o]
    attr(map, "grid") <- grid[o]
    attr(map, "step") <- step
    attr(map, "off_end") <- off_end

    map
}

# pseudomarker map, minimal number of pseudomarkers to add
#
# map = vector of marker positions; names = marker names
# step = min distance between markers or pseudomarkers
# off_end = amount to go beyond the terminal markers
# tol = tolerance for determining whether a marker hits a pseudomarker
# pmar_stem = leading part of name of pseudomarkers
create_pseudomarker_map_minimal <-
function(map, step, off_end=0, tol=0.01, pmar_stem="loc")
{
    if(any(is.na(map))) stop("map values can't be missing")
    if(step==0) {
        attr(map, "index") <- seq(along=map)
        attr(map, "stepwidth") <- "max"
        attr(map, "step") <- step
        attr(map, "off_end") <- off_end
        return(map)
    }
    if(step < 0) stop("step should be >= 0")
    if(tol >= step) stop("tol should be << step")
    if(off_end < 0) stop("off_end should be >= 0")
    if(tol < 0) stop("tol should be >= 0")

    # distance between markers
    d <- diff(map)
    if(any(d < 0))
        stop("marker positions should be sorted")

    pmar <- c(seq(min(map), min(map)-off_end, by=-step)[-1],
              seq(max(map), max(map)+off_end, by=step)[-1])

    if(length(map) == 1) { # just one marker
        index <- c(1, rep(0, length(pmar)))
        map <- c(map, pmar)
        o <- order(map)
        map <- map[o]
        attr(map, "index") <- index[o]
        attr(map, "stepwidth") <- "max"
        attr(map, "step") <- step
        attr(map, "off_end") <- off_end
        return(map)
    }

    # intervals needing pseudomarkers (d = distances between adjacent markers
    int2fill <- which(d > step)

    pmar <- c(pmar, unlist(lapply(int2fill, function(a) {
        n <- ceiling((map[a+1]-map[a])/step)+1
        seq(map[a], map[a+1], length.out=n)[-c(1, n)] })))

    if(length(pmar) == 0) {
        attr(map, "index") <- seq(along=map)
        attr(map, "stepwidth") <- "max"
        attr(map, "step") <- step
        attr(map, "off_end") <- off_end
        return(map)
    }

    # no. digits to use in name
    digits <- ceiling(-log10(step))
    digits <- ifelse(digits < 0, 0, digits)

    names(pmar) <- paste0(pmar_stem, round(pmar, digits))

    # index of markers vs pseudomarkers
    index <- c(seq(along=map), rep(0, length(pmar)))

    # add pseudomarkers to markers
    map <- c(map, pmar)

    # sort the map
    o <- order(map)
    map <- map[o]
    attr(map, "stepwidth") <- "max"
    attr(map, "step") <- step
    attr(map, "off_end") <- off_end
    attr(map, "index") <- index[o]

    map
}

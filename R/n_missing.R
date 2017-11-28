# n_missing
#' Count missing genotypes
#'
#' Number (or proportion) of missing (or non-missing) genotypes by individual or marker
#'
#' @md
#'
#' @param cross An object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param by Whether to summarize by individual or marker
#' @param summary Whether to take count or proportion
#'
#' @return Vector of counts (or proportions) of missing (or non-missing) genotypes.
#'
#' @describeIn n_missing Count missing genotypes
#' @export
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' nmis_ind <- n_missing(iron)
#' pmis_mar <- n_typed(iron, "mar", "proportion")
#' plot(nmis_ind, xlab="Individual", ylab="No. missing genotypes")
#' plot(pmis_mar, xlab="Markers", ylab="Prop. genotyped")
n_missing <-
    function(cross, by=c("individual", "marker"), summary=c("count", "proportion"))
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    by <- match.arg(by)
    summary <- match.arg(summary)

    result <- NULL
    denom <- 0
    for(i in seq(along=cross$geno)) {
        if(by=="individual") {
            this_result <- rowSums(cross$geno[[i]]==0)
            denom <- denom + ncol(cross$geno[[i]])
        }
        else
            this_result <- colSums(cross$geno[[i]]==0)


        if(is.null(result)) result <- this_result
        else {
            if(by=="individual")
                result <- result + this_result
            else
                result <- c(result, this_result)
        }
    }

    if(summary == "proportion") {
        if(by=="marker")
            denom <- nrow(cross$geno[[1]])
        result <- result/denom
    }

    result
}

#' @describeIn n_missing Count genotypes
#' @export
n_typed <-
    function(cross, by=c("individual", "marker"), summary=c("count", "proportion"))
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    by <- match.arg(by)
    summary <- match.arg(summary)

    opp <- n_missing(cross, by, summary)
    if(summary=="proportion") return(1-opp)
    else if(by=="individual") return(tot_mar(cross) - opp)
    else return(n_ind(cross) - opp)
}

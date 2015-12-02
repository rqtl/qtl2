# get_x_covar
#' Get X chromosome covariates
#'
#' Get the matrix of covariates to be used for the null hypothesis when
#' performing QTL analysis with the X chromosome.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#'
#' @return A matrix of size individuals x no. covariates.
#'
#' @details For most crosses, the result is either a matrix with no
#' columns (indicating no additional covariates are needed) or a
#' matrix with a single column containing sex indicators (1 for males
#' and 0 for females).
#'
#'  For an intercross, we also consider cross direction. There are
#' four cases: (1) All male or all female but just one direction: no
#' covariate; (2) All female but both directions: covariate indicating
#' cross direction; (3) Both sexes, one direction: covariate
#' indicating sex; (4) Both sexes, both directions: a covariate
#' indicating sex and a covariate that is 1 for females from the
#' reverse direction and 0 otherwise.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' xcovar <- get_x_covar(iron)

get_x_covar <-
function(cross)
{
    # check input
    if(class(cross) != "cross2")
        stop('Input cross must have class "cross2"')

    x <- .get_x_covar(cross$crosstype, cross$is_female, cross$cross_info)

    if(ncol(x)==0) return(NULL)

    # add rownames
    nam <- names(cross$is_female)
    if(!is.null(nam)) rownames(x) <- nam

    x
}

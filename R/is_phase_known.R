#' Is cross phase-known
#'
#' Determine if a cross is of phase-known type
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#'
#' @return TRUE/FALSE indicating whether the cross is phase-known or not.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' is_phase_known(grav2)

is_phase_known <-
    function(cross)
{
    .is_phase_known(cross$crosstype)
}

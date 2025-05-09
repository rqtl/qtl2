#' @name CCcolors
#' @aliases CCcolors
#'
#' @title Collaborative Cross colors
#'
#' @description A vector of 8 colors for use with the mouse
#' Collaborative Cross and Diversity Outbreds.
#'
#' @details
#' `CCorigcolors` are the original eight colors for the Collaborative Cross founder strains.
#' `CCaltcolors` are a slightly modified version, but still not color-blind friendly.
#' `CCcolors` are derived from the the Okabe-Ito color blind friendly palette in
#'  Wong (2011) Nature Methods \doi{10.1038/nmeth.1618}.
#'
#' @source <https://web.archive.org/web/20250215070655/https://csbio.unc.edu/CCstatus/index.py?run=AvailableLines.information>
#'
#' @keywords datasets
#'
#' @examples
#' data(CCcolors)
#' data(CCaltcolors)
#' data(CCorigcolors)
#  c(AJ = "#F0E442", B6 = "#555555", `129` = "#E69F00", NOD = "#0072B2",
#    NZO = "#56B4E9", CAST = "#009E73", PWK = "#D55E00", WSB = "#CC79A7")
NULL

#' @rdname CCcolors
#' @name CCorigcolors
#' @aliases CCorigcolors
# c(AJ = "#FFFF00", B6 = "#888888", `129` = "#FF8888", NOD = "#1111FF",
#   NZO = "#00CCFF", CAST = "#00AA00", PWK = "#FF0000", WSB = "#9900EE")
NULL

#' @rdname CCcolors
#' @name CCaltcolors
#' @aliases CCaltcolors
# c(AJ = "#FFDC00", B6 = "#888888", `129` = "#F08080", NOD = "#0064C9",
#   NZO = "#7FDBFF", CAST = "#2ECC40", PWK = "#FF4136", WSB = "#B10DC9")
NULL

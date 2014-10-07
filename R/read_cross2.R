# read_cross2
#' Read QTL data from files
#'
#' Read QTL data from a set of files
#'
#' @param yaml_file Character string with path to the yaml file
#' containing all of the control information.
#'
#' @return Object of class \code{"cross2"}. For details, see the
#' R/qtl2 Developer Guide. **FIX ME (put the url here)**
#'
#' @details A control file in YAML format contains information above
#' basic parameters as well as the names of the series of data files
#' to be read. Sample data files at **FIX ME (put url here)**.
#'
#' @export
#' @keywords IO
#' @examples
#' yaml_file <- system.file("sampledata", "grav.yaml", package="qtl2")
#' grav <- read_cross2(yaml_file)
read_cross2 <-
function(yaml_file)
{
    # directory containing the data
    dir <- dirname(yaml_file)

    # load the control file
    control <-  yaml::yaml.load_file(yaml_file)

    # grab cross type
    if(!("crosstype" %in% names(control)))
        stop("crosstype not found")
    output <- list(crosstype=control$crosstype)

    # grab the major bits
    major_files <- c("geno", "gmap", "pmap", "pheno", "covar", "phenocovar")
    for(i in major_files) {
        if(i %in% names(control)) {
            tmp <- data.table::fread(file.path(dir, control[[i]]),
                                     na.strings=control$na.strings,
                                     sep=control$sep)
            class(tmp) <- "data.frame"
            output[[i]] <- tmp
        }
    }

    output
}

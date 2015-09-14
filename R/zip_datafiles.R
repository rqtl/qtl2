# zip_datafiles
#' Zip a set of data files
#'
#' Zip a set of data files (in format read by \code{\link{read_cross2}}).
#'
#' @param control_file Character string with path to the control file
#' (\href{http://www.yaml.org}{YAML} or \href{http://www.json.org/}{JSON})
#' containing all of the control information.
#' @param zip_file Name of zip file to use. If omitted, we use the
#' stem of \code{control_file} but with a \code{.zip} extension.
#'
#' @return Character string with the file name of the zip file that
#' was created.
#'
#' @details The input \code{control_file} is the control file (in
#' \href{http://www.yaml.org}{YAML} or \href{http://www.json.org/}{JSON} format)
#' to be read by \code{\link{read_cross2}}.  (See the
#' \href{http://kbroman.org/qtl2/pages/sampledata.html}{sample data files} and the
#' \href{http://kbroman.org/qtl2/assets/vignettes/input_files.html}{vignette describing the input file format}.)
#'
#' The \code{\link[utils]{zip}} function is used to do the zipping.
#'
#' @export
#' @keywords IO
#' @seealso \code{\link{read_cross2}}, sample data files at \url{http://kbroman.org/qtl2/pages/sampledata.html}
#' @examples
#' \dontrun{
#' control_file <- "~/grav2_data/grav2.yaml"
#' zip_datafiles(control_file, "grav2.zip")
#' }
zip_datafiles <-
function(control_file, zip_file)
{
    control_file <- path.expand(control_file)
    if(!(file.exists(control_file)))
        stop("The control file (", control_file, ") doesn't exist.")

    dir <- dirname(control_file)

    if(missing(zip_file) || is.null(zip_file))
        zip_file <- gsub("\\.[a-z]+$", ".zip", control_file)

    # read control file
    control <-  read_control_file(control_file)

    # get all of the file names
    sections <- c("geno", "gmap", "pmap", "pheno", "covar", "phenocovar", "founder_geno")
    files <- basename(control_file)
    for(section in sections) {
        if(section %in% names(control))
            files <- c(files, control[[section]])
    }

    # sex and cross_info as files?
    sections <- c("sex", "cross_info")
    for(section in sections) {
        if(section %in% names(control)) {
            if("file" %in% names(control[[section]]))
                files <- c(files, control[[section]][["file"]])
        }
    }

    # linemap as a file?
    if("linemap" %in% names(control)) {
        filename <- file.path(dir, control[["linemap"]])
        if(file.exists(filename))
            files <- c(files, control[["linemap"]])
    }

    # zip the files (with the -j flag, directory info not included in zip file)
    utils::zip(zip_file, file.path(dir, files), flags="-j")

    invisible(zip_file)
}

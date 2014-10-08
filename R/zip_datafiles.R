# zip_datafiles
#' Zip a set of data files
#'
#' Zip a set of data files (in format read by \code{\link{read_cross2}}).
#'
#' @param yaml_file Character string with path to the yaml file
#' containing all of the control information.
#' @param zip_file Name of zip file to use. If omitted, we use the
#' stem of \code{yaml_file} but with a \code{.zip} extension.
#'
#' @return Character string with the file name of the zip file that
#' was created.
#'
#' @details The input \code{yaml_file} is the control file (in
#' \href{http://www.yaml.org}{YAML} format for data to in the format
#' to be read by \code{\link{read_cross2}}.  (See the
#' \href{http://kbroman.org/qtl2/pages/sampledata.html}{sample data files} and the
#' \href{http://kbroman.org/qtl2/assets/vignettes/input_files.html}{vignette describing the input file format}.)
#'
#' The \code{\link[utils]{zip}} function is used to do the zipping.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @keywords IO
#' @seealso \code{\link{read_cross2}}, sample data files at \url{http://kbroman.org/qtl2/pages/sampledata.html}
#' @examples
#' \dontrun{
#' yaml_file <- "~/grav2_data/grav2.yaml"
#' zip_datafiles(yaml_file, "grav2.zip")
#' }
zip_datafiles <-
function(yaml_file, zip_file)
{
    yaml_file <- path.expand(yaml_file)
    if(!(file.exists(yaml_file)))
        stop("The yaml file (", yaml_file, ") doesn't exist.")

    dir <- dirname(yaml_file)

    if(missing(zip_file) || is.null(zip_file))
        zip_file <- gsub("\\.yaml$", ".zip", yaml_file)

    # read yaml file
    control <-  yaml::yaml.load_file(yaml_file)

    # get all of the file names
    sections <- c("geno", "gmap", "pmap", "pheno", "covar", "phenocovar", "founder_geno")
    files <- basename(yaml_file)
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
        if(file.exists(as_filename))
            files <- c(files, control[["linemap"]])
    }

    # zip the files (with the -j flag, directory info not included in zip file)
    utils::zip(zip_file, file.path(dir, files), flags="-j")

    invisible(zip_file)
}

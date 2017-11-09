# zip_datafiles
#' Zip a set of data files
#'
#' Zip a set of data files (in format read by \code{\link{read_cross2}}).
#'
#' @param control_file Character string with path to the control file
#' (\href{http://www.yaml.org}{YAML} or \href{http://www.json.org/}{JSON})
#' containing all of the control information.
#' @param zip_file Name of zip file to use. If NULL, we use the
#' stem of \code{control_file} but with a \code{.zip} extension.
#' @param quiet If \code{FALSE}, print progress messages.
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
#' The files should all be contained within the directory where the
#' \code{control_file} sits, or in a subdirectory of that directory.
#' If file paths use \code{..}, these get stripped by zip, and so the
#' resulting zip file may not work with \code{\link{read_cross2}}.
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
function(control_file, zip_file=NULL, quiet=TRUE)
{
    control_file <- path.expand(control_file)
    if(!(file.exists(control_file)))
        stop("The control file (", control_file, ") doesn't exist.")

    dir <- dirname(control_file)

    if(is.null(zip_file))
        zip_file <- sub("\\.[a-z]+$", ".zip", control_file)

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

    # flag for quiet
    zip_flags <- ifelse(quiet, "-q", "")

    # move to the directory with the files
    cwd <- getwd()
    on.exit(setwd(cwd)) # move back on exit
    setwd(dir)

    # check for ".." in file paths
    patterns <- c("^\\.\\.\\/", "\\/\\.\\.\\/")
    patterns <- gsub("/", .Platform$file.sep, patterns, fixed=TRUE)
    if(any(grepl(patterns[1], files) | grepl(patterns[2], files))) {
        warning('zip strips ".." from file paths, so ', zip_file,
                ' may not work with read_cross2().')
    }

    # do the zipping
    utils::zip(zip_file, files, flags=zip_flags)

    invisible(zip_file)
}

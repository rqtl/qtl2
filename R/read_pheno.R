#' Read phenotype data
#'
#' Read phenotype data from a CSV file (and, optionally, phenotype
#' covariate data from a separate zip file). The CSV files may be
#' contained in zip files, separately or togther.
#'
#' @param file Character string with path to the phenotype data file
#' (or a zip file containing both the phenotype and phenotype
#' covariate files).
#' @param phenocovarfile Character string with path to the phenotype
#' covariate file. This can be a separate CSV or zip file; if a zip
#' file, it must contain exactly one CSV file. Alternatively, if the
#' \code{file} argument indicates a zip file that contains two files
#' (phenotypes and phenotype covariates), then this
#' \code{phenocovarfile} argument must indicate the base name for the
#' phenotype covariate file.
#' @param sep the field separator character
#' @param na.strings a character vector of strings which are to be
#' interpreted as \code{NA} values.
#' @param comment.char A character vector of length one containing a
#' single character to denote comments within the CSV files.
#' @param transpose If TRUE, the phenotype data will be transposed. The
#' phenotype covariate information is \bold{never} transposed.
#' @param quiet If \code{FALSE}, print progress messages.
#'
#' @return Either a matrix of phenotype data, or a list containing
#' \code{pheno} (phenotype matrix) and \code{phenocovar} (phenotype
#' covariate matrix).
#'
#' @export
#' @keywords IO
#' @seealso \code{\link{read_cross2}},
#' sample data files at \url{http://kbroman.org/qtl2/pages/sampledata.html}
#' and \url{https://github.com/rqtl/qtl2data}
#'
#' @examples
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/Gough/gough_pheno.csv")
#' phe <- read_pheno(file)
#'
#' phecovfile <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                      "qtl2data/master/Gough/gough_phenocovar.csv")
#' phe_list <- read_pheno(file, phecovfile)
#' }
read_pheno <-
    function(file, phenocovarfile=NULL, sep=",", na.strings=c("-", "NA"),
             comment.char="#", transpose=FALSE, quiet=TRUE)
{
    # handle zip file, which may contain phenocovarfile
    if(grepl("\\.zip$", file)) { # zip file
        dir <- tempdir()
        if(is_web_file(file)) {
            tmpfile <- tempfile()
            if(!quiet) message(" - downloading ", file, "\n       to ", tmpfile)
            utils::download.file(file, tmpfile, quiet=TRUE)
            file <- tmpfile
            on.exit(unlink(tmpfile))
        }

        if(!quiet) message(" - unzipping ", file, "\n       to ", dir)
        file <- path.expand(file)
        stop_if_no_file(file)
        unzipped_files <- utils::unzip(file, exdir=dir)
        basenames <- basename(unzipped_files)

        on.exit({ # clean up when done
            if(!quiet) message(" - cleaning up")
            unlink(unzipped_files)
        }, add=TRUE)

        other_files <- unzipped_files
        if(!is.null(phenocovarfile) && any(basenames == phenocovarfile)) {
            # look for phenocovarfile in the unzipped files
            if(sum(basenames==phenocovarfile) > 1)
                stop("Multiple copies of ", phenocovarfile, " in zip file")
            other_files <- unzipped_files[basenames != phenocovarfile]
            phenocovarfile <- unzipped_files[basenames==phenocovarfile]
        }

        # halt if more than one file left
        if(length(other_files) > 1) {
            stop("Unclear which file to read: ",
                 paste(basename(other_files), collapse=" "))
        }

        # redefine file with new name
        file <- other_files
    }

    # directory containing the data
    file <- path.expand(file)

    # read file
    stop_if_no_file(file)
    if(!quiet) message(" - Reading ", basename(file))
    pheno <- read_csv_numer(file, sep=sep, na.strings=na.strings,
                            comment.char=comment.char, transpose=transpose)

    # read phenocovarfile
    if(!is.null(phenocovarfile)) {
        if(!quiet) message(" - Reading ", basename(phenocovarfile))
        phenocovar <- read_phenocovar(phenocovarfile, sep=sep, na.strings=na.strings,
                                      comment.char=comment.char, quiet=quiet)
        return(list(pheno=pheno, phenocovar=phenocovar))
    }

    pheno
}

# read phenocovar file (to handle case where it's also zipped)
read_phenocovar <-
    function(file, phenocovarfile=NULL, sep=",", na.strings=c("-", "NA"),
             comment.char="#", quiet=TRUE)
{
    # handle zip file, which may contain phenocovarfile
    if(grepl("\\.zip$", file)) { # zip file
        dir <- tempdir()
        if(is_web_file(file)) {
            tmpfile <- tempfile()
            if(!quiet) message(" - downloading ", file, "\n       to ", tmpfile)
            utils::download.file(file, tmpfile, quiet=TRUE)
            file <- tmpfile
            on.exit(unlink(tmpfile))
        }

        if(!quiet) message(" - unzipping ", file, "\n       to ", dir)
        file <- path.expand(file)
        stop_if_no_file(file)
        unzipped_files <- utils::unzip(file, exdir=dir)

        on.exit({ # clean up when done
            if(!quiet) message(" - cleaning up")
            unlink(unzipped_files)
        }, add=TRUE)

        if(length(unzipped_files) > 1) {
            stop("Unclear which file to read: ",
                 paste(basename(unzipped_files), collapse=" "))
        }

        # redefine file with new name
        file <- unzipped_files
    }

    # directory containing the data
    file <- path.expand(file)

    # read file
    stop_if_no_file(file)
    if(!quiet) message(" - Reading ", basename(file))
    read_csv(file, sep=sep, na.strings=na.strings,
             comment.char=comment.char, transpose=FALSE)
}

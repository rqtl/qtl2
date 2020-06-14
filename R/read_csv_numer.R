#' Read a csv file that has numeric columns
#'
#' Read a csv file via [data.table::fread()] using a
#' particular set of options, including the ability to transpose the
#' result. This version assumes that the contents other than the first
#' column and the header row are strictly numeric.
#'
#' @param filename Name of input file
#' @param sep Field separator
#' @param na.strings Missing value codes
#' @param comment.char Comment character; rest of line after this character is ignored
#' @param transpose If TRUE, transpose the result
#' @param rownames_included If TRUE, the first column is taken to be row names.
#'
#' @return Data frame
#'
#' @details Initial two lines can contain comments with number of rows
#' and columns. Number of columns includes an ID column; number of
#' rows does not include the header row.
#'
#' The first column is taken to be a set of row names
#'
#' @seealso [read_csv()]
#'
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' \dontrun{mydata <- read_csv_numer("myfile.csv", transpose=TRUE)}
read_csv_numer <-
    function(filename, sep=",", na.strings=c("NA", "-"), comment.char="#", transpose=FALSE,
             rownames_included=TRUE)
{
    # read header and extract expected number of rows and columns
    # (number of columns includes ID column; number of rows does *not* include header row)
    header <- read_header(filename, comment.char=comment.char)
    expected_dim <- extract_dim_from_header(header)

    # determine expected columns
    expected_ncol <- expected_dim[2]
    if(is.na(expected_ncol)) { # no header about no. columns
        x <- data.table::fread(filename, na.strings=na.strings, sep=sep, header=TRUE,
                               verbose=FALSE, showProgress=FALSE, data.table=FALSE,
                               nrow=1, skip=length(header))
        expected_ncol <- ncol(x)
    }
    if(expected_ncol == 1) colClasses <- "character"
    else colClasses <- c("character", rep("numeric", expected_ncol-1))

    # read the data
    x <- data.table::fread(filename, na.strings=na.strings, sep=sep, header=TRUE,
                           verbose=FALSE, showProgress=FALSE, data.table=FALSE,
                           colClasses=colClasses, skip=length(header))

    # check that number of rows and columns match expected from header
    for(i in 1:2) {
        labels <- c("rows", "columns")
        if(!is.na(expected_dim[i])) { # nrows given
            if(dim(x)[i] != expected_dim[i])
                stop('In file "', filename, '", no. ', labels[i],
                     ' (', dim(x)[i], ') != expected (', expected_dim[i], ')')
        }
    }

    # move first column to row names
    if(rownames_included)
        x <- firstcol2rownames(x, filename)

    # transpose if requested
    if(transpose)
        x <- as.data.frame(t(x), stringAsFactors=FALSE)

    x
}

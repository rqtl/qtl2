#' Read a csv file
#'
#' Read a csv file via \code{\link[data.table]{fread}} using a
#' particular set of options, including the ability to transpose the
#' result.
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
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' \dontrun{mydata <- read_csv("myfile.csv", transpose=TRUE)}
read_csv <-
    function(filename, sep=",", na.strings=c("NA", "-"), comment.char="#", transpose=FALSE,
             rownames_included=TRUE)
{
    # read header and extract expected number of rows and columns
    # (number of columns includes ID column; number of rows does *not* include header row)
    header <- read_header(filename, comment.char=comment.char)
    expected_dim <- extract_dim_from_header(header)

    # read the data
    x <- data.table::fread(filename, na.strings=na.strings, sep=sep, header=TRUE,
                           verbose=FALSE, showProgress=FALSE, data.table=FALSE,
                           colClasses="character", skip=length(header))

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

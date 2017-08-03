#' Create a function to query variants
#'
#' Create a function that will connect to a SQLite database of founder
#' variant information and return a data frame with variants for a
#' selected region.
#'
#' @param dbfile Name of database file
#' @param db Optional database connection (provide one of \code{file} and \code{db}).
#' @param table_name Name of table in the database
#' @param chr_field Name of chromosome field
#' @param pos_field Name of position field
#'
#' @return Function with three arguments, \code{chr}, \code{start},
#'     and \code{end}, which returns a data frame with the variants in
#'     that region, with \code{start} and \code{end} being in Mbp. The
#'     output should contain at least the columns \code{chr} and
#'     \code{pos}, the latter being position in Mbp.
#'
#' @details Note that this function assumes that the database has a
#'     \code{pos} field that is in basepairs, but the selection uses
#'     \code{start} and \code{end} positions in Mbp, and the output
#'     data frame should have \code{pos} in Mbp.
#'
#' @export
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery
#'
#' @examples
#' \dontrun{
#' # create query function by connecting to file
#' query_variants <- create_variant_query_func("cc_variants.sqlite")
#' # query_variants will connect and disconnect each time
#' variants <- query_variants("19", 25.1, 26.1)
#'
#' # connect and disconnect separately
#' library(RSQLite)
#' db <- dbConnect(SQLite(), "cc_variants.sqlite")
#' query_variants <- create_variant_query_func(db=db)
#' variants <- query_variants("19", 25.1, 26.1)
#' dbDisconnect(db)
#' }

create_variant_query_func <-
    function(dbfile=NULL, db=NULL, table_name="variants", chr_field="chr", pos_field="pos")
{
    if(!is.null(db)) {
        if(!is.null(dbfile))
            warning("Provide just one of dbfile or db; using db")

        query_func <- function(chr, start, end) {

            # convert start and end to basepairs
            start <- round(start*1e6)
            end <- round(end*1e6)

            # bits of the query
            chrselect <- paste0(chr_field, " == '", chr, "'")
            startselect <- paste0(pos_field, " >= ", start)
            endselect <- paste0(pos_field, " <= ", end)

            # the full query
            query <- paste0("SELECT * FROM ", table_name, " WHERE ",
                            chrselect, " AND ",
                            startselect, " AND ",
                            endselect)

            # do the query
            result <- RSQLite::dbGetQuery(db, query)

            # include pos column in Mbp
            result$pos <- result[,pos_field]/1e6

            result
        }

    }
    else {
        if(is.null(dbfile) || dbfile=="")
            stop("Provide either dbfile or db")

        query_func <- function(chr, start, end)
        {
            db <- RSQLite::dbConnect(RSQLite::SQLite(), dbfile)
            on.exit(RSQLite::dbDisconnect(db)) # disconnect on exit

            # convert start and end to basepairs
            start <- round(start*1e6)
            end <- round(end*1e6)

            # bits of the query
            chrselect <- paste0(chr_field, " == '", chr, "'")
            startselect <- paste0(pos_field, " >= ", start)
            endselect <- paste0(pos_field, " <= ", end)

            # the full query
            query <- paste0("SELECT * FROM ", table_name, " WHERE ",
                            chrselect, " AND ",
                            startselect, " AND ",
                            endselect)

            # do the query
            result <- RSQLite::dbGetQuery(db, query)

            # include pos column in Mbp
            result$pos <- result[,pos_field]/1e6

            result
        }
    }

    query_func
}

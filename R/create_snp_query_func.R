#' Create a function to query SNPs
#'
#' Create a function that will connect to a SQLite database of founder SNP
#' information and return a data frame with SNP information for a
#' selected region.
#'
#' @param dbfile Name of database file
#' @param db Optional database connection (provide one of \code{file} and \code{db}).
#' @param table_name Name of table in the database
#' @param chr_field Name of chromosome field
#' @param pos_field Name of position field
#'
#' @return Function with three arguments, \code{chr}, \code{start}, and \code{end},
#' which returns a data frame with the SNPs in that region. The output should contain
#' at least the columns \code{chr} and \code{pos}, the latter being position in Mbp.
#'
#' @export
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery
#'
#' @examples
#' \dontrun{
#' # create query function by connecting to file
#' query_snps <- create_snp_query_func("ccvariants.sqlite")
#' # query_snps will connect and disconnect each time
#' snps <- query_snps("19", 25.1, 26.1)
#'
#' # connect and disconnect separately
#' db <- dbConnect(SQLite(), "ccvariants.sqlite")
#' query_snps <- create_snp_query_func(db=db)
#' snps <- query_snps("19", 25.1, 26.1)
#' dbDisconnect(db)
#' }

create_snp_query_func <-
    function(dbfile=NULL, db=NULL, table_name="snps", chr_field="chr", pos_field="pos_Mbp")
{
    if(!is.null(db)) {
        if(!is.null(dbfile))
            warning("Provide just one of dbfile or db; using db")

        query_func <- function(chr, start, end) {

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
            RSQLite::dbGetQuery(db, query)
        }

    }
    else {
        if(is.null(dbfile) || dbfile=="")
            stop("Provide either dbfile or db")

        query_func <- function(chr, start, end)
        {
            db <- RSQLite::dbConnect(RSQLite::SQLite(), dbfile)
            on.exit(RSQLite::dbDisconnect(db)) # disconnect on exit

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
            RSQLite::dbGetQuery(db, query)
        }
    }

    query_func
}

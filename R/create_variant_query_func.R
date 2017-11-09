#' Create a function to query variants
#'
#' Create a function that will connect to a SQLite database of founder
#' variant information and return a data frame with variants for a
#' selected region.
#'
#' @md
#'
#' @param dbfile Name of database file
#' @param db Optional database connection (provide one of `file` and `db`).
#' @param table_name Name of table in the database
#' @param chr_field Name of chromosome field
#' @param pos_field Name of position field
#' @param filter Additional SQL filter (as a character string)
#'
#' @return Function with three arguments, `chr`, `start`,
#'     and `end`, which returns a data frame with the variants in
#'     that region, with `start` and `end` being in Mbp. The
#'     output should contain at least the columns `chr` and
#'     `pos`, the latter being position in Mbp.
#'
#' @details Note that this function assumes that the database has a
#'     `pos` field that is in basepairs, but the selection uses
#'     `start` and `end` positions in Mbp, and the output
#'     data frame should have `pos` in Mbp.
#'
#' @export
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery
#'
#' @examples
#' # create query function by connecting to file
#' dbfile <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2db")
#' query_variants <- create_variant_query_func(dbfile)
#' # query_variants will connect and disconnect each time
#' variants <- query_variants("2", 97.0, 98.0)
#'
#' # create query function to just grab SNPs
#' query_snps <- create_variant_query_func(dbfile, filter="type=='snp'")
#' # query_variants will connect and disconnect each time
#' snps <- query_snps("2", 97.0, 98.0)
#'
#' # connect and disconnect separately
#' library(RSQLite)
#' db <- dbConnect(SQLite(), dbfile)
#' query_variants <- create_variant_query_func(db=db)
#' variants <- query_variants("2", 97.0, 98.0)
#' dbDisconnect(db)

create_variant_query_func <-
    function(dbfile=NULL, db=NULL, table_name="variants",
             chr_field="chr", pos_field="pos", filter=NULL)
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

            if(!is.null(filter) && filter != "")
                query <- paste0(query, " AND (", filter, ")")

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
            if(!file.exists(dbfile))
                stop("File ", dbfile, " doesn't exist")

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

            if(!is.null(filter) && filter != "")
                query <- paste0(query, " AND (", filter, ")")

            # do the query
            result <- RSQLite::dbGetQuery(db, query)

            # include pos column in Mbp
            result$pos <- result[,pos_field]/1e6

            result
        }
    }

    query_func
}

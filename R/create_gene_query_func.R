#' Create a function to query genes
#'
#' Create a function that will connect to a SQLite database of gene
#' information and return a data frame with gene information for a
#' selected region.
#'
#' @param dbfile Name of database file
#' @param db Optional database connection (provide one of `file` and `db`).
#' @param table_name Name of table in the database
#' @param chr_field Name of chromosome field
#' @param start_field Name of field with start position (in basepairs)
#' @param stop_field Name of field with stop position (in basepairs)
#' @param filter Additional SQL filter (as a character string).
#'
#' @return Function with three arguments, `chr`, `start`,
#'     and `end`, which returns a data frame with the genes
#'     overlapping that region, with `start` and `end` being
#'     in Mbp. The output should contain at least the columns
#'     `Name`, `chr`, `start`, and `stop`, the
#'     latter two being positions in Mbp.
#'
#' @details Note that this function assumes that the database has
#'     `start` and `stop` fields that are in basepairs, but
#'     the selection uses positions in Mbp, and the output data frame
#'     should have `start` and `stop` columns in Mbp.
#'
#' Also note that a SQLite database of MGI mouse genes
#' is available at figshare:
#' [doi:10.6084/m9.figshare.5286019.v7](https://doi.org/10.6084/m9.figshare.5286019.v7)
#'
#' @export
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery
#'
#' @examples
#' # create query function by connecting to file
#' dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2")
#' query_genes <- create_gene_query_func(dbfile, filter="(source=='MGI')")
#' # query_genes will connect and disconnect each time
#' genes <- query_genes("2", 97.0, 98.0)
#'
#' # connect and disconnect separately
#' library(RSQLite)
#' db <- dbConnect(SQLite(), dbfile)
#' query_genes <- create_gene_query_func(db=db, filter="(source=='MGI')")
#' genes <- query_genes("2", 97.0, 98.0)
#' dbDisconnect(db)

create_gene_query_func <-
    function(dbfile=NULL, db=NULL, table_name="genes",
             chr_field="chr", start_field="start", stop_field="stop",
             filter=NULL)
{
    if(!is.null(db)) {
        if(!is.null(dbfile))
            warning("Provide just one of dbfile or db; using db")

        query_func <- function(chr, start, end) {

            # convert input positions from Mbp to basepairs
            start <- round(start * 1e6)
            end <- round(end * 1e6)

            # bits of the query
            chrselect <- paste0(chr_field, " == '", chr, "'")
            pos_select1 <- paste0("(", start_field, " >= ", start, " AND ", start_field, " <= ", end, ")")
            pos_select2 <- paste0("(", stop_field, " >= ", start, " AND ", stop_field, " <= ", end, ")")
            pos_select3 <- paste0("(", start_field, " <= ", end, " AND ", stop_field, " >= ", start, ")")

            # the full query
            query <- paste0("SELECT * FROM ", table_name, " WHERE ",
                            chrselect, " AND ( ",
                            pos_select1, " OR ", pos_select2, " OR ", pos_select3, " )")
            if(!is.null(filter) && filter != "")
                query <- paste0(query, " AND (", filter, ")")

            # do the query
            result <- RSQLite::dbGetQuery(db, query)

            # include start/stop columns in Mbp
            result$start <- result[,start_field]/1e6
            result$stop <- result[,stop_field]/1e6

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

            # convert input positions from Mbp to basepairs
            start <- round(start * 1e6)
            end <- round(end * 1e6)

            # bits of the query
            chrselect <- paste0(chr_field, " == '", chr, "'")
            pos_select1 <- paste0("(", start_field, " >= ", start, " AND ", start_field, " <= ", end, ")")
            pos_select2 <- paste0("(", stop_field, " >= ", start, " AND ", stop_field, " <= ", end, ")")
            pos_select3 <- paste0("(", start_field, " <= ", end, " AND ", stop_field, " >= ", start, ")")

            # the full query
            query <- paste0("SELECT * FROM ", table_name, " WHERE ",
                            chrselect, " AND ( ",
                            pos_select1, " OR ", pos_select2, " OR ", pos_select3, " )")
            if(!is.null(filter) && filter != "")
                query <- paste0(query, " AND (", filter, ")")

            # do the query
            result <- RSQLite::dbGetQuery(db, query)

            # include start/stop columns in Mbp
            result$start <- result[,start_field]/1e6
            result$stop <- result[,stop_field]/1e6

            result
        }
    }

    query_func
}

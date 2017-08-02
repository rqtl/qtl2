#' Create a function to query genes
#'
#' Create a function that will connect to a SQLite database of gene
#' information and return a data frame with gene information for a
#' selected region.
#'
#' @param dbfile Name of database file
#' @param db Optional database connection (provide one of \code{file} and \code{db}).
#' @param table_name Name of table in the database
#' @param chr_field Name of chromosome field
#' @param start_field Name of field with start position (in Mbp)
#' @param stop_field Name of field with stop position (in Mbp)
#' @param filter Additional SQL filter (as a character string).
#'
#' @return Function with three arguments, \code{chr}, \code{start},
#'     and \code{end}, which returns a data frame with the genes
#'     overlapping that region. The output should contain at least the
#'     columns \code{Name}, \code{chr}, \code{start}, and \code{stop}, the latter
#'     two being positions in Mbp.
#'
#' @export
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery
#'
#' @examples
#' \dontrun{
#' # create query function by connecting to file
#' dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2db")
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
#' }

create_gene_query_func <-
    function(dbfile=NULL, db=NULL, table_name="genes",
             chr_field="chr", start_field="start_Mbp", stop_field="stop_Mbp",
             filter=NULL)
{
    if(!is.null(db)) {
        if(!is.null(dbfile))
            warning("Provide just one of dbfile or db; using db")

        query_func <- function(chr, start, end) {

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

            # clean up start/stop names: make them 'start' and 'stop'
            start_col <- result$start_field
            stop_col <- result$stop_field

            cols <- names(result)
            cols_keep <- cols[is.na(match(cols, c("start", "stop", start_field, stop_field)))]
            result <- result[,cols,drop=FALSE]
            result$start <- start_col
            result$stop <- stop_col

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

            # clean up start/stop names: make them 'start' and 'stop'
            start_col <- result$start_field
            stop_col <- result$stop_field

            cols <- names(result)
            cols_keep <- cols[is.na(match(cols, c("start", "stop", start_field, stop_field)))]
            result <- result[,cols,drop=FALSE]
            result$start <- start_col
            result$stop <- stop_col

            result
        }
    }

    query_func
}

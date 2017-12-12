# used to skip tests when not on *my* computer in interactive mode
# needed for testing mult-core code, which isn't allowed within R CMD check
isnt_karl <-
    function()
{
    !(interactive() &&
      identical(Sys.getenv("KARL_LOCAL"), "true"))
}

# is a number? (from assertthat)
is_number <- function(x) is.numeric(x) && length(x)==1
is_nonneg_number <- function(x) is_number(x) && x >= 0
is_pos_number <- function(x) is_number(x) && x > 0

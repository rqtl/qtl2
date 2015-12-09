# used to skip tests when not on *my* computer in interactive mode
# needed for testing mult-core code, which isn't allowed within R CMD check
isnt_karl <-
    function()
{
    !(interactive() &&
      identical(Sys.getenv("KARL_LOCAL"), "true"))
}

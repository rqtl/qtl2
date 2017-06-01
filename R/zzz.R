# based on code in the tidyverse package
#     https://github.com/tidyverse/tidyverse

.onAttach <- function(...)
{
    needed <- core[!is_attached(core)]
    if(length(needed)==0)
        return()

    qtl2_attach()
}

is_attached <- function(x)
{
    paste0("package:", x) %in% search()
}

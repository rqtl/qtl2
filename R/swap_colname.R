# swap_colname
#
# replace column name in dataframe
# column labeled orig_name becomes new_name
# if there had already been a column labeled new_name,
#     we append "_old"
#
# (used in create_variant_query_func and create_gene_query_func)
swap_colname <-
    function(df, orig_name, new_name)
{
    if(orig_name == new_name) return(df)

    cn <- colnames(df)
    if(!any(cn == orig_name)) {
        stop('field "', orig_name, '" not found')
    }

    if(any(cn==new_name)) {
        cn[cn==new_name] <- paste0(new_name, "_old")
    }
    cn[cn==orig_name] <- new_name

    colnames(df) <- cn
    df
}

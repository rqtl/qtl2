# dimensions of calc_genoprob object
dim.calc_genoprob <-
    function(x)
{
  vapply(x, dim, rep(1,3))
}

# dimnames of calc_genoprob object
dimnames.calc_genoprob <-
    function(x)
{
  dnames <- lapply(x, dimnames)

  list(ind = dnames[[1]][[1]],
       gen = lapply(dnames, '[[', 2),
       mar = lapply(dnames, '[[', 3))
}

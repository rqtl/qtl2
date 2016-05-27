# test function
calc_width <-
    function(pos)
{
    c(pos[2]-pos[1],
      (pos[-(1:2)] - pos[-(length(pos)-c(0,1))])/2,
      pos[length(pos)]-pos[length(pos)-1])
}

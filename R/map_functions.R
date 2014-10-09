# map functions
mf <-
function(d, map.function=c("haldane", "kosambi", "c-f", "morgan"))
{
    map.function <- match.arg(map.function)
    switch(map.function,
           haldane=mf.h(d),
           kosambi=mf.k(d),
           "c-f"=  mf.cf(d),
           morgan= mf.m(d))
}

imf <-
function(r, map.function=c("haldane", "kosambi", "c-f", "morgan"))
{
    map.function <- match.arg(map.function)
    switch(map.function,
           haldane=imf.h(r),
           kosambi=imf.k(r),
           "c-f"=  imf.cf(r),
           morgan= imf.m(r))
}


mf.h <-
function(d)
    0.5*(1-exp(-d/50))
imf.h <-
function(r)
{
    r[r >= 0.5] <- 0.5-1e-14
    -50*log(1-2*r)
}

mf.k <-
function(d)
    0.5*tanh(d/50)

imf.k <-
function(r)
{
    r[r >= 0.5] <- 0.5-1e-14
    50*atanh(2*r)
}

mf.m <-
function(d)
    pmin(d/100, 0.5)

imf.m <-
function(r)
    pmin(r*100, 50)

# carter-falconer
imf.cf <- function(r) {
    r[r >= 0.5] <- 0.5-1e-14
    12.5*(log(1+2*r)-log(1-2*r))+25*atan(2*r)
}

mf.cf <-
function(d)
{
    d[d >= 300] <- 300

    icf <- function(r,d)
        imf.cf(r)-d

    sapply(d,function(a) {
        if(a==0) return(0)
        uniroot(icf, c(0,0.5-1e-14),d=a,tol=1e-12)$root })
}


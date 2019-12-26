# convert2cross2
#' Convert R/qtl cross object to new format
#'
#' Convert a cross object from the R/qtl format to the R/qtl2 format
#'
#' @param cross An object of class `"cross"`; see
#' [qtl::read.cross()] for details.
#'
#' @return Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#'
#' @export
#' @keywords utilities
#' @seealso [read_cross2()]
#' @examples
#' library(qtl)
#' data(hyper)
#' hyper2 <- convert2cross2(hyper)
convert2cross2 <-
function(cross)
{
    if(is.null(cross)) stop("cross is NULL")
    crosstype <- rqtl1_crosstype(cross)
    check_crosstype(crosstype)
    result <- list(crosstype=crosstype)
    n.ind <- rqtl_nind(cross)

    # genetic map, and grab chrtype
    result$gmap <- rqtl_pull_map(cross)
    class(result$gmap) <- "list"
    result$is_x_chr <- vapply(result$gmap, function(a) rqtl1_chrtype(a)=="X", TRUE)
    for(i in seq(along=result$gmap))
        class(result$gmap[[i]]) <- "numeric"

    # ids
    ids <- rqtl_getid(cross)
    if(is.null(ids)) ids <- as.character(1:n.ind)

    # split out genotype data
    result$geno <- lapply(cross$geno, "[[", "data")
    for(i in seq(along=result$geno))
        rownames(result$geno[[i]])  <- ids

    # sex/pgm
    sexpgm <- rqtl_getsex(cross)
    if(is.null(sexpgm$sex))
        result$is_female <- rep(FALSE, n.ind)
    else result$is_female <- (sexpgm$sex == 0)
    names(result$is_female) <- ids

    if(is.null(sexpgm$pgm))
        result$cross_info <- matrix(0L, ncol=1, nrow=n.ind) # if missing, assume they're all 0's
    else result$cross_info <- matrix(as.integer(sexpgm$pgm))
    rownames(result$cross_info) <- ids

    # convert X chr genotypes
    if(any(result$is_x_chr) && crosstype %in% c("bc", "f2", "bcsft")) { # bcsft not really supported yet
        for(i in which(result$is_x_chr)) {
            result$geno[[i]] <- rqtl_reviseXdata(crosstype, "simple",
                                                sexpgm, geno=result$geno[[i]],
                                                cross.attr=attributes(cross), force=TRUE)
        }
    }

    # in genotypes, replace NAs with 0s
    for(i in seq(along=result$geno)) {
        result$geno[[i]][is.na(result$geno[[i]])] <- 0L
        storage.mode(result$geno[[i]]) <- "integer"
    }

    # phenotypes: pull out numeric columns and put the rest in covariates
    phe <- cross$pheno

    # put sex and pgm in covar data even if numeric
    phe.names <- colnames(phe)
    is_sex <- is_pgm <- rep(FALSE, ncol(phe))
    is_sex[grep("^[Ss][Ee][Xx]$", phe.names)] <- TRUE
    is_pgm[grep("^[Pp][Gg][Mm]$", phe.names)] <- TRUE
    numercol <- vapply(phe, is.numeric, TRUE) & (!is_sex) & (!is_pgm)
    covar <- phe[,!numercol,drop=FALSE]
    if(ncol(covar) > 0) {
        result$covar <- covar
        rownames(result$covar) <- ids
    }
    phe <- phe[, numercol, drop=FALSE]
    if(ncol(phe) > 0) {
        result$pheno <- as.matrix(phe)
        rownames(result$pheno) <- ids
        storage.mode(result$pheno) <- "double"
    }

    # alleles
    alleles <- attr(cross, "alleles")
    if(!is.null(alleles)) result$alleles <- alleles
    else result$alleles <- LETTERS[seq_len(nalleles(result$crosstype))]

    class(result) <- "cross2"

    check_cross2(result) # double-check

    result
}


rqtl_nind <-
    function(object)
{
    if(!inherits(object, "cross"))
        stop("Input should have class \"cross\".")

    n.ind1 <- nrow(object$pheno)
    n.ind2 <- sapply(object$geno,function(x) nrow(x$data))
    if(any(n.ind2 != n.ind1))
        stop("Different numbers of individuals in genotypes and phenotypes.")
    n.ind1
}

# R/qtl pull.map
rqtl_pull_map <-
    function(cross)
{
    if(!inherits(cross, "cross"))
        stop("Input should have class \"cross\".")

    result <- lapply(cross$geno,function(a) {
        b <- a$map
        class(b) <- as.character(class(a))
        b })
    class(result) <- "map"
    result
}

# R/qtl getid
rqtl_getid <-
    function(cross)
{
    phe <- cross$pheno
    nam <- names(phe)
    if("id" %in% nam) {
        id <- phe$id
        phenam <- "id"
    }
    else if("ID" %in% nam) {
        id <- phe$ID
        phenam <- "ID"
    }
    else if("Id" %in% nam) {
        id <- phe$Id
        phenam <- "Id"
    }
    else if("iD" %in% nam) {
        id <- phe$iD
        phenam <- "iD"
    }
    else {
        id <- NULL
        phenam <- NULL
    }

    if(is.factor(id))
        id <- as.character(id)

    attr(id, "phenam") <- phenam

    id
}

# R/qtl getsex
rqtl_getsex <-
    function(cross)
{
    type <- rqtl1_crosstype(cross)
    if(type != "bc" && type != "f2" && type != "4way") return(list(sex=NULL, pgm=NULL))

    phe.names <- names(cross$pheno)

    sex.column <- grep("^[Ss][Ee][Xx]$", phe.names)
    pgm.column <- grep("^[Pp][Gg][Mm]$", phe.names)

    if(length(sex.column)==0) { # no sex included
        sex <- NULL
    }
    else {
        if(length(sex.column)>1)
            warning("'sex' included multiple times.  Using the first one.")
        temp <- cross$pheno[,sex.column[1]]
        if(is.numeric(temp)) {
            if(any(!is.na(temp) & temp != 0 & temp != 1)) {
                warning("Sex column should be coded as 0=female 1=male; sex ignored.")
                sex <- NULL
            }
            else sex <- temp
        }
        else {
            if(!is.factor(temp)) temp <- as.factor(temp)

            if(length(levels(temp)) == 1) {
                if(levels(temp) == "F" || levels(temp)=="f" ||
                   toupper(levels(temp)) == "FEMALE") sex <- rep(0,rqtl_nind(cross))
                else if(levels(temp) == "M" || levels(temp)=="m" ||
                        toupper(levels(temp)) == "MALE") sex <- rep(1,rqtl_nind(cross))
                else
                    warning("Sex column should be coded as 0=female 1=male; sex ignored.")
            }
            else if(length(levels(temp)) > 2) {
                warning("Sex column should be coded as a two-level factor; sex ignored.")
                sex <- NULL
            }
            else { # is a factor with two levels
                lev <- levels(temp)
                if(length(grep("^[Ff]",lev))>0 &&
                   length(males <- grep("^[Mm]",lev))>0) {
                    temp <- as.character(temp)
                    sex <- rep(0,length(temp))
                    sex[is.na(temp)] <- NA
                    sex[!is.na(temp) & temp==lev[males]] <- 1
                }
                else
                    warning("Don't understand levels in sex column; sex ignored.")
            }
        }
    }

    if(length(pgm.column)==0 || type=="4way") { # no pgm included
        pgm <- NULL
    }
    else {
        if(length(pgm.column)>1)
            warning("'pgm' included multiple times.  Using the first one.")
        temp <- cross$pheno[,pgm.column[1]]
        if(!is.numeric(temp))
            temp <- as.numeric(temp)-1
        if(any(!is.na(temp) & temp != 0 & temp != 1)) {
            warning("pgm column should be coded as 0/1; pgm ignored.")
            pgm <- NULL
        }
        else pgm <- temp
    }

    if(!is.null(sex) && any(is.na(sex))) {
        if(all(sex[!is.na(sex)]==1)) {
            warning(sum(is.na(sex)), " individuals with missing sex; assuming they're male like the others")
            sex[is.na(sex)] <- 1
        }
        else if(all(sex[!is.na(sex)]==0)) {
            warning(sum(is.na(sex)), " individuals with missing sex; assuming they're female like the others")
            sex[is.na(sex)] <- 0
        }
        else {
            warning(sum(is.na(sex)), " individuals with missing sex; assuming they're female")
            sex[is.na(sex)] <- 0
        }
    }

    if(!is.null(pgm) && any(is.na(pgm))) {
        if(all(pgm[!is.na(pgm)]==1)) {
            warning(sum(is.na(pgm)), " individuals with missing pgm; assuming pgm==1 like the others")
            pgm[is.na(pgm)] <- 1
        }
        else if(all(pgm[!is.na(pgm)]==0)) {
            warning(sum(is.na(pgm)), " individuals with missing pgm; assuming pgm==0 like the others")
            pgm[is.na(pgm)] <- 0
        }
        else {
            warning(sum(is.na(pgm)), " individuals with missing pgm; assuming pgm==0")
            pgm[is.na(pgm)] <- 0
        }
    }

    list(sex=sex,pgm=pgm)
}

# R/qtl getgenonames
#
# get names of genotypes
# used in discan, effectplot, plotPXG, scanone, scantwo, vbscan, reviseXdata
# cross.attr gives the cross attributes
rqtl_getgenonames <-
    function(type=c("f2","bc","riself","risib","4way","dh","haploid","special","bcsft"),
             chrtype=c("A","X"), expandX=c("simple","standard","full"),
             sexpgm, cross.attr=NULL)
{
    type <- match.arg(type)
    chrtype <- match.arg(chrtype)
    expandX <- match.arg(expandX)

    ## Treat bcsft as bc if no intercross generations; otherwise as f2.
    if(type == "bcsft") {
        if(cross.attr$scheme[2] == 0)
            type <- "bc"
        else
            type <- "f2"
    }

    if(chrtype=="X") {
        sex <- sexpgm$sex
        pgm <- sexpgm$pgm
    }

    if(type=="special") return(cross.attr$genotypes)

    if(is.null(cross.attr) || !("alleles" %in% names(cross.attr))) {
        if(type == "4way") alleles <- LETTERS[1:4]
        else alleles <- LETTERS[1:2]
    }
    else
        alleles <- cross.attr$alleles

    tempgn <- c(paste(rep(alleles[1],2),collapse=""),
                paste(alleles,collapse=""),
                paste(rep(alleles[2],2),collapse=""),
                paste(alleles[1],"Y",sep=""),
                paste(alleles[2],"Y",sep=""))

    # get rid of missing sex and pgm values, if there are any
    if(chrtype=="X") {
        if(length(sex)>0) sex <- sex[!is.na(sex)]
        if(length(pgm)>0) pgm <- pgm[!is.na(pgm)]
    }

    if(type=="riself" || type=="risib" || type=="dh")
        gen.names <- tempgn[c(1,3)]

    else if(type=="haploid")
        gen.names <- alleles

    else if(type == "4way") {
        if(chrtype=="A")
            gen.names <- c(paste(alleles[1],alleles[3],sep=""),
                           paste(alleles[2],alleles[3],sep=""),
                           paste(alleles[1],alleles[4],sep=""),
                           paste(alleles[2],alleles[4],sep=""))
        else
            gen.names <- c(paste(alleles[1],alleles[3],sep=""),
                           paste(alleles[2],alleles[3],sep=""),
                           paste(alleles[1],"Y",sep=""),
                           paste(alleles[2],"Y",sep=""))
    }

    else if(type == "bc") {

        if(chrtype=="A") # autosome
            gen.names <- tempgn[1:2] # AA AB

        else { # X chromosome

            #                 simple     standard       full
            #   -both sexes   A-/AB/BY   AA/AB/AY/BY    same as std
            #   -all females  AA/AB      same           same
            #   -all males    AY/BY      same           same

            if(length(sex)==0 || all(sex==0)) # all females
                gen.names <- tempgn[1:2] # AA AB
            else if(all(sex==1)) # all males
                gen.names <- tempgn[4:5] # AY BY
            else { # some of each
                if(expandX == "simple")
                    gen.names <- c(paste(alleles[1], "-", sep=""),
                                   tempgn[c(2,5)]) # A-, AB, BY
                else gen.names <- tempgn[c(1,2,4,5)]  # AA,AB,AY,BY
            }
        }
    }

    else { # intercross
        if(chrtype == "A")  # autosomal
            gen.names <- tempgn[1:3]
        else { # X chromsome

            # both crosses     simple     standard         full
            #   -both sexes   A-/AB/B-    AA/AB/BB/AY/BY   AA/AB1/AB2/BB/AY/BY
            #   -all females  AA/AB/BB    same as simple   AA/AB1/AB2/BB
            #   -all males    AY/BY       same             same
            # forw cross
            #   -both sexes   A-/AB/BY    AA/AB/AY/BY      same as std
            #   -all females  AA/AB       same             same
            #   -all males    AY/BY       same             same
            # backw cross
            #   -both sexes   B-/AB/AY    BB/AB/AY/BY      same as std
            #   -all females  BB/AB       same             same
            #   -all males    AY/BY       same             same

            if(length(sex)==0 || all(sex==0)) { # all females
                if(length(pgm)==0 || all(pgm==0)) # all forw dir
                    gen.names <- tempgn[1:2] # AA AB
                else if(all(pgm==1))  # all backw dir
                    gen.names <- tempgn[3:2] # BB AB
                else { # some of each direction
                    if(expandX=="full")
                        gen.names <- c(tempgn[1],
                                       paste(tempgn[2],c("f","r"), sep=""),
                                       tempgn[3])
                    else gen.names <- tempgn[1:3]
                }
            }
            else if(all(sex==1))  # all males
                gen.names <- tempgn[4:5]
            else { # some of each sex
                if(length(pgm)==0 || all(pgm==0)) { # all forw
                    if(expandX=="simple")
                        gen.names <- c(paste(alleles[1],"-", sep=""),
                                       tempgn[c(2,5)])
                    else gen.names <- tempgn[c(1,2,4,5)]
                }
                else if (all(pgm==1)) { # all backw
                    if(expandX=="simple")
                        gen.names <- c(paste(alleles[2], "-",sep=""),
                                       tempgn[c(2,4)])
                    else gen.names <- tempgn[c(3,2,4,5)]
                }
                else { # some of each dir
                    if(expandX=="simple")
                        gen.names <- c(paste(alleles[1],"-",sep=""),
                                       tempgn[2],
                                       paste(alleles[2],"-",sep=""))
                    else if(expandX=="standard")
                        gen.names <- tempgn
                    else
                        gen.names <- c(tempgn[1],
                                       paste(tempgn[2],c("f","r"),sep=""),
                                       tempgn[3:5])
                }
            }
        }
    }

    gen.names
}

# R/qtl reviseXdata reduced to geno case
rqtl_reviseXdata <-
    function(type=c("f2","bc","bcsft"), expandX=c("simple","standard","full"),
             sexpgm, geno, cross.attr, force=FALSE)
{
    type <- match.arg(type)
    expandX <- match.arg(expandX)

    ## Treat bcsft as bc if no intercross generations; otherwise as f2.
    if(type == "bcsft") {
        if(cross.attr$scheme[2] == 0)
            type <- "bc"
        else
            type <- "f2"
    }

    sex <- sexpgm$sex
    pgm <- sexpgm$pgm

    # get genonames
    genonames <- rqtl_getgenonames(type, "X", expandX, sexpgm, cross.attr)

    if(type == "bc") { # backcross

        if(length(sex)==0 || ((all(sex==0) || all(sex==1)) && !force)) { # all one sex
            # no changes necessary
            return(geno)
        }

        else { # both sexes

            gmale <- geno[sex==1,]
            if(expandX=="simple")
                gmale[!is.na(gmale) & gmale==2] <- 3
            else {
                gmale[!is.na(gmale) & gmale==1] <- 3
                gmale[!is.na(gmale) & gmale==2] <- 4
            }
            geno[sex==1,] <- gmale
            return(geno)
        } # end of "both sexes" / backcross

    } # end of backcross

    else { # intercross

        if(length(sex)==0 || all(sex==0)) { # all females

            if(length(pgm)==0 || ((all(pgm==0) || all(pgm==1)) && !force)) { # one dir, females
                return(geno)
            }

            else { # both dir, females
                gback <- geno[pgm==1,]
                if(expandX!="full") {
                    gback[!is.na(gback) & gback==1] <- 3
                    geno[pgm==1,] <- gback
                }
                else {
                    gback[!is.na(gback) & gback==1] <- 4
                    gback[!is.na(gback) & gback==2] <- 3
                    geno[pgm==1,] <- gback
                }
                return(geno)
            }
        }
        else if(all(sex==1) && !force)  { # all males
            return(geno)
        }

        else { # both sexes

            if(length(pgm)==0 || all(pgm==0)) { # both sexes, forw dir
                gmale <- geno[sex==1,]
                if(expandX=="simple")
                    gmale[!is.na(gmale) & gmale==2] <- 3
                else {
                    gmale[!is.na(gmale) & gmale==1] <- 3
                    gmale[!is.na(gmale) & gmale==2] <- 4
                }
                geno[sex==1,] <- gmale
                return(geno)

            } # both sexes, forw dir

            if(all(pgm==1) && !force) { # both sexes, backw dir
                gmale <- geno[sex==1,]
                if(expandX!="full") {
                    gmale[!is.na(gmale) & gmale==1] <- 3
                    gmale[!is.na(gmale) & gmale==2] <- 1
                }
                else {
                    gmale[!is.na(gmale) & gmale==1] <- 3
                    gmale[!is.na(gmale) & gmale==2] <- 4
                }
                geno[sex==1,] <- gmale
                return(geno)
            } # both sexes, backw dir

            else { # both dir, both sexes
                gmale <- geno[sex==1,]
                gfemaler <- geno[sex==0 & pgm==1,]
                if(expandX=="simple") {
                    gmale[!is.na(gmale) & gmale==2] <- 3
                    gfemaler[!is.na(gfemaler) & gfemaler==1] <- 3
                }
                else if(expandX=="standard") {
                    gmale[!is.na(gmale) & gmale==1] <- 4
                    gmale[!is.na(gmale) & gmale==2] <- 5
                    gfemaler[!is.na(gfemaler) & gfemaler==1] <- 3
                }
                else {
                    gmale[!is.na(gmale) & gmale==1] <- 5
                    gmale[!is.na(gmale) & gmale==2] <- 6
                    gfemaler[!is.na(gfemaler) & gfemaler==1] <- 4
                    gfemaler[!is.na(gfemaler) & gfemaler==2] <- 3
                }
                geno[sex==1,] <- gmale
                geno[sex==0 & pgm==1,] <- gfemaler
                return(geno)
            }
        }

    } # end of intercross

}

# R/qtl1 cross type
rqtl1_crosstype <-
    function(cross)
    {
        type <- class(cross)
        type <- type[type != "cross" & type != "list"]
        if(length(type) > 1) {
            warning("cross has multiple classes")
        }
        type[1]
    }

rqtl1_chrtype <-
    function(object)
    {
        if(inherits(object, "X")) return("X")
        "A"
    }

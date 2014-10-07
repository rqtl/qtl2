# read_cross2
#' Read QTL data from files
#'
#' Read QTL data from a set of files
#'
#' @param yaml_file Character string with path to the yaml file
#' containing all of the control information.
#'
#' @return Object of class \code{"cross2"}. For details, see the
#' R/qtl2 Developer Guide. **FIX ME (put the url here)**
#'
#' @details A control file in YAML format contains information above
#' basic parameters as well as the names of the series of data files
#' to be read. Sample data files at **FIX ME (put url here)**.
#'
#' @export
#' @keywords IO
#' @seealso Sample data files at \url{http://kbroman.org/qtl2/pages/sampledata.html}
#' @examples
#' yaml_file <- "http://kbroman.org/qtl2/assets/sampledata/grav2/grav2.yaml"
#' grav2 <- read_cross2(yaml_file)
read_cross2 <-
function(yaml_file)
{
    # directory containing the data
    dir <- dirname(yaml_file)

    # load the control file
    control <-  yaml::yaml.load_file(yaml_file)

    # grab cross type
    if("crosstype" %in% names(control))
        output <- list(crosstype=control$crosstype)
    else
        stop("crosstype not found")

    # grab genotype encodings or use defaults
    if("genotypes" %in% names(control))
        genotypes <- control$genotypes
    else
        genotypes <- list(A=1, H=2, B=3, D=4, C=5)

    # grab the major bits
    major_files <- c("geno", "gmap", "pmap", "pheno", "covar", "phenocovar", "founder_geno")
    for(section in major_files) {
        if(section %in% names(control)) {
            sheet <- data.table::fread(file.path(dir, control[[section]]),
                                       na.strings=control$na.strings,
                                       sep=control$sep,                     # note that nulls just get ignored
                                       verbose=FALSE, showProgress=FALSE, data.table=FALSE)

            # treat first column as rownames
            sheet <- firstcol2rownames(sheet, section)

            # change genotype codes and convert phenotypes to numeric matrix
            if(section=="geno" || section=="founder_geno")
                sheet <- encode_geno(sheet, genotypes)
            else if(section=="pheno")
                sheet <- pheno2matrix(sheet)

            output[[section]] <- sheet
        }
    }

    # pull out a map
    if("gmap" %in% names(output))
        map <- output$gmap
    else if("pmap" %in% output)
        map <- output$pmap
    else stop("Need a genetic or physical marker map")

    if(!("geno" %in% names(output)))
       stop("No genotype data found.")

    # split genotypes by chromosome
    geno <- c("geno", "founder_geno")
    for(section in geno) {
        if(section %in% names(output))
            output[[section]] <- split_geno(output[[section]], map)
    }

    # split maps by chr
    maps <- c("gmap", "pmap")
    for(section in maps) {
        if(section %in% names(output))
            output[[section]] <- split_map(output[[section]])
    }

    # X chr?
    chr <- names(output$geno)
    output$is_x_chr <- rep(FALSE, length(chr))
    names(output$is_x_chr) <- chr
    if("x_chr" %in% names(control)) {
        x_chr <- control$x_chr # name of X chromosome
        output$is_x_chr[x_chr] <- TRUE
    }

    # sex
    output$is_female <- convert_sex(control$sex, output$covar, control$sep, dir)
    if(is.null(output$is_female)) { # missing; assume all FALSE
        output$is_female <- rep(FALSE, nrow(output$geno[[1]]))
        names(output$is_female) <- rownames(output$geno[[1]])
    }

    # cross_info
    output$cross_info <- convert_cross_info(control$cross_info, output$covar, control$sep, dir)
    if(is.null(output$cross_info)) { # missing; make a 0-column matrix
        output$cross_info <- matrix(ncol=0, nrow=nrow(output$geno[[1]]))
        rownames(output$cross_info) <- rownames(output$geno[[1]])
    }

    # alleles?
    if("alleles" %in% names(control))
        output$alleles <- control$alleles

    class(output) <- "cross2"
    output
}

# convert genotype data, using genotype encodings
# genotypes is a list with names = code in data
#                     and values = new numeric code
encode_geno <-
function(geno, genotypes)
{
    newgeno <- geno <- as.matrix(geno)

    if(any(unlist(genotypes)==0))
        stop("Can't encode genotypes as 0, that's used for missing values.")

    # genotype codes
    gnames <- names(genotypes)
    codes <- genotypes

    # any mismatches?
    gch <- as.character(geno)
    mismatch <- !is.na(gch) & is.na(match(gch,gnames))
    if(any(mismatch)) {
        warning(sum(mismatch), " genotypes treated as missing: ",
                paste0('"', unique(gch[mismatch]), '"', collapse=", "))
        newgeno[mismatch] <- NA
    }

    # re-code
    for(g in gnames)
        newgeno[!is.na(geno) & geno==g] <- genotypes[[g]]

    # turn missing values to 0s
    newgeno[is.na(newgeno)] <- "0"

    #
    storage.mode(newgeno) <- "integer"

    newgeno
}

# convert a phenotype data.frame to a matrix of doubles
pheno2matrix <-
function(pheno)
{
    pheno <- as.matrix(pheno)
    storage.mode(pheno) <- "double"

    pheno
}

# make the first column of a data frame the row names
firstcol2rownames <-
function(covar, name="")
{
    # assume the first column is the ID
    idcolumn <- covar[,1]

    # look for duplicates
    check4duplicates(idcolumn, name)

    # drop that column and make it the row names
    covar <- covar[,-1,drop=FALSE]
    rownames(covar) <- as.character(idcolumn)

    covar
}

# check a list of IDs for duplicates
check4duplicates <-
function(ids, name="")
{
    if(name != "") name <- paste("in", name)
    dup <- duplicated(ids)
    if(any(dup)) {
        stop("Not all ids ", name, " are unique: ",
             paste0('"', unique(ids[dup]), '"', collapse=", "))
    }
    TRUE
}

# split genotype data into list of chromosomes
split_geno <-
function(geno, map)
{
    gmark <- colnames(geno)
    mmark <- rownames(map)

    m <- match(gmark, mmark)
    if(any(is.na(m)))
        stop("Some markers in genotype data not found in map: ",
             paste0('"', gmark[is.na(m)], '"', collapse=", "))

    chr <- map[match(colnames(geno), rownames(map)), 1]
    pos <- map[match(colnames(geno), rownames(map)), 2]

    uchr <- unique(chr)
    newgeno <- vector("list", length(uchr))
    names(newgeno) <- uchr
    for(i in uchr) {
        wh <- which(chr==i)
        o <- order(pos[chr==i])
        newgeno[[i]] <- geno[,wh[o], drop=FALSE]
    }

    newgeno
}

# split a data frame with chr, pos (and marker names as rows)
#   into a list of chromosomes
split_map <-
function(map)
{
    pos <- map[,2]
    names(pos) <- rownames(map)

    lapply(split(pos, map[,1]), sort)
}

# grab sex information
convert_sex <-
function(sex_control, covar, sep, dir)
{
    if("covar" %in% names(sex_control)) { # sex within the covariates
        sex <- covar[,sex_control[["covar"]], drop=FALSE]
    }
    else if("file" %in% names(sex_control)) { # look for file
        sex <- data.table::fread(file.path(dir, sex_control[["file"]]),
                                 verbose=FALSE, showProgress=FALSE, data.table=FALSE)
        sex <- firstcol2rownames(sex)
    }
    else return(NULL)

    # convert to vector
    id <- rownames(sex)
    sex <- sex[,1]
    names(sex) <- id

    # missing values?
    if(any(is.na(sex))) {
        stop(sum(is.na(sex)), " missing sexes (sex can't be missing).")
    }

    # grab the rest of sex_control as conversion codes
    codes <- sex_control[is.na(match(names(sex_control), c("file", "covar")))]
    sexcode <- names(codes)

    if(length(codes)==0) { # no codes, pass over as is
        storage.mode(sex) <- "integer"
        return(sex==0)
    }

    # any mismatches?
    sexch <- as.character(sex)
    mismatch <- !is.na(sexch) & is.na(match(sexch,sexcode))
    if(any(mismatch)) {
        stop(sum(mismatch), " sexes don't match the codes (sex can't be missing): ",
                paste0('"', unique(sexch[mismatch]), '"', collapse=", "))
    }

    # re-code
    newsex <- sex
    for(code in sexcode)
        newsex[sex==code] <- codes[[code]]

    storage.mode(newsex) <- "integer"

    (newsex==0)
}

# grab cross_info
convert_cross_info <-
function(cross_info_control, covar, sep, dir)
{
    if("covar" %in% names(cross_info_control)) { # cross_info within the covariates
        cross_info <- covar[,cross_info_control[["covar"]], drop=FALSE]
    }
    else if("file" %in% names(cross_info_control)) { # look for file
        cross_info <- data.table::fread(file.path(dir, cross_info_control[["file"]]),
                                 verbose=FALSE, showProgress=FALSE, data.table=FALSE)
        if(any(is.na(cross_info)))
            stop(sum(is.na(cross_info)), " missing values in cross_info (cross_info can't be missing.")
        return(firstcol2rownames(cross_info))
    }
    else return(NULL)

    # make it a matrix
    cross_info <- as.matrix(cross_info)

    if(any(is.na(cross_info))) {
        stop(sum(is.na(cross_info)), " missing values in cross_info (cross_info can't be missing).")
    }

    # grab the rest of cross_info_control as conversion codes
    codes <- cross_info_control[is.na(match(names(cross_info_control), c("file", "covar")))]
    cicode <- names(codes)

    if(length(codes)==0) { # no codes, pass over as is
        storage.mode(cross_info) <- "integer"
        return(cross_info)
    }

    # any mismatches?
    cich <- as.character(cross_info)
    mismatch <- !is.na(cich) & is.na(match(cich,cicode))
    if(any(mismatch)) {
        stop(sum(mismatch), " cross_info vals don't match the codes (cross_info can't be missing): ",
                paste0('"', unique(cich[mismatch]), '"', collapse=", "))
    }
    # re-code
    newci <- cross_info
    for(code in cicode)
        newci[cross_info==code] <- codes[[code]]

    storage.mode(newci) <- "integer"

    newci
}


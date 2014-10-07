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
#' @examples
#' yaml_file <- system.file("sampledata", "grav.yaml", package="qtl2")
#' grav <- read_cross2(yaml_file)
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
                                     sep=control$sep) # note that nulls just get ignored
            class(sheet) <- "data.frame"

            # treat first column as rownames
            sheet <- firstcol2rownames(sheet, section)

            if(section=="geno" || section=="founder_geno")
                sheet <- encode_geno(sheet, genotypes)
            else if(section=="pheno")
                sheet <- pheno2matrix(sheet)

            output[[section]] <- sheet
        }
    }

    # alleles?
    if("alleles" %in% names(control))
        output$alleles <- control$alleles

    class(output) <- "cross2"
    output
}

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

pheno2matrix <-
function(pheno)
{
    pheno <- as.matrix(pheno)
    storage.mode(pheno) <- "double"

    pheno
}

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

check4duplicates <-
function(ids, name="")
{
    if(name != "") name <- paste("in", name)
    dup <- duplicated(ids)
    if(any(dup)) {
        stop("Not all ids ", name, " are unique: ",
             paste0('"', unique(id[dup]), '"', collapse=", "))
    }
    TRUE
}

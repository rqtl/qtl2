# read_cross2
#' Read QTL data from files
#'
#' Read QTL data from a set of files
#'
#' @param file Character string with path to the
#' \href{http://www.yaml.org}{YAML} or \href{http://www.json.org/}{JSON} file containing all of the control
#' information. This could instead be a zip file containing all of the
#' data files, in which case the contents are unzipped to a temporary
#' directory and then read.
#' @param quiet If \code{FALSE}, print progress messages.
#'
#' @return Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#'
#' @details
#' A control file in \href{http://www.yaml.org}{YAML} or
#' \href{http://www.json.org/}{JSON} format contains information
#' about basic parameters as well as the names of the series of
#' data files to be read. See the
#' \href{http://kbroman.org/qtl2/pages/sampledata.html}{sample
#' data files} and the
#' \href{http://kbroman.org/qtl2/assets/vignettes/input_files.html}{vignette
#' describing the input file format}.
#'
#' @export
#' @keywords IO
#' @seealso \code{\link{read_pheno}}, \code{\link{write_control_file}},
#' sample data files at \url{http://kbroman.org/qtl2/pages/sampledata.html}
#' and \url{http://github.com/rqtl/qtl2data}
#'
#' @examples
#' \dontrun{
#' yaml_file <- "http://kbroman.org/qtl2/assets/sampledata/grav2/grav2.yaml"
#' grav2 <- read_cross2(yaml_file)
#' }
#' zip_file <- system.file("extdata", "grav2.zip", package="qtl2geno")
#' grav2 <- read_cross2(zip_file)
read_cross2 <-
function(file, quiet=TRUE)
{
    if(length(grep("\\.zip$", file)) > 0) { # zip file
        dir <- qtl2_temp_dir()

        if(is_web_file(file)) {
            tmpfile <- tempfile()
            if(!quiet) message(" - downloading ", file, "\n       to ", tmpfile)
            utils::download.file(file, tmpfile, quiet=TRUE)
            file <- tmpfile
            on.exit(unlink(tmpfile))
        }

        if(!quiet) message(" - unzipping ", file, "\n       to ", dir)
        file <- path.expand(file)
        stop_if_no_file(file)
        unzipped_files <- utils::unzip(file, exdir=dir)

        # remove temporary directory on exit
        on.exit({
            if(!quiet) message(" - cleaning up")
            unlink(dir, recursive=TRUE)
        }, add=TRUE)

        # ignore "__MACOSX/._" files
        unzipped_files <- grep("__MACOSX/._", unzipped_files, fixed=TRUE, invert=TRUE, value=TRUE)

        if(any(grepl("\\.yaml$", unzipped_files))) {
            file <- grep("\\.yaml$", unzipped_files, value=TRUE)
            if(length(file) > 1)
                stop("The zip file contains multiple yaml files")
            if(any(grepl("\\.json$", unzipped_files)))
                warning("The zip file contains both YAML and JSON files; using the YAML file.")
        }
        else if(any(grepl("\\.json$", unzipped_files))) {
            file <- grep("\\.json$", unzipped_files, value=TRUE)
            if(length(file) > 1)
                stop("The zip file contains multiple json files")
        }
        else {
            stop('No ".yaml" or ".json" control file found')
        }
    }

    # directory containing the data
    file <- path.expand(file)
    dir <- dirname(file)

    # load the control file
    stop_if_no_file(file)
    control <-  read_control_file(file)

    # for keeping track of use of control stuff
    used_control <- rep(FALSE, length(control))
    names(used_control) <- names(control)
    used_control["sep"] <- used_control["na.strings"] <-
        used_control["comment.char"] <- used_control["description"] <-
            used_control["comments"] <- TRUE

    # grab cross type
    if("crosstype" %in% names(control)) {
        output <- list(crosstype=control$crosstype)
        used_control["crosstype"] <- TRUE # indicate that we used it
    }
    else
        stop("crosstype not found")

    # grab genotype encodings or use defaults
    if("genotypes" %in% names(control)) {
        genotypes <- control$genotypes
        used_control["genotypes"] <- TRUE # indicate that we used it
    }
    else
        genotypes <- list(A=1, H=2, B=3, D=4, C=5)

    # grab the major bits
    sections <- c("geno", "gmap", "pmap", "pheno", "covar", "phenocovar", "founder_geno")
    for(section in sections) {
        if(section %in% names(control)) {
            used_control[section] <- TRUE # indicate that we used it
            if(!quiet) message(" - reading ", section)

            # transposed?
            tr <- paste0(section, "_transposed")
            used_control[tr] <- TRUE # indicate that we used it
            tr <- tr %in% names(control) && control[[tr]]

            filename <- control[[section]]

            if(length(filename)==1) { # single file
                filename <- file.path(dir, filename)
                stop_if_no_file(filename)

                # read file
                sheet <- read_csv(filename, na.strings=control$na.strings, sep=control$sep,
                                  comment.char=control$comment.char, transpose=tr,
                                  rownames_included=TRUE)
            }
            else { # vector of files
                # add dir to paths
                filenames <- vapply(filename, function(a) file.path(dir, a), "")

                # check that the all exist
                lapply(filenames, stop_if_no_file)

                if(section=="gmap" || section=="pmap") {
                    # rbind the chromosomes together
                    sheet <- NULL
                    for(filename in filenames)
                        sheet <- rbind(sheet,
                                       read_csv(filename, na.strings=control$na.strings,
                                                sep=control$sep, comment.char=control$comment.char,
                                                transpose=tr, rownames_included=TRUE))
                }
                else {
                    # read all of the files and cbind(), matching on row names
                    sheet <- read_mult_csv(filenames, na.strings=control$na.strings,
                                           sep=control$sep, comment.char=control$comment.char,
                                           transpose=tr)
                }
            }
            if(length(unique(colnames(sheet))) != ncol(sheet))
                warning("Duplicate column names in ", section, " data")

            # change genotype codes and convert phenotypes to numeric matrix
            if(section=="geno" || section=="founder_geno") {
                if(!quiet) message(" - encoding ", section)
                sheet <- recode_geno(sheet, genotypes)
            }
            else if(section=="pheno") {
                sheet <- pheno2matrix(sheet)
            }

            output[[section]] <- sheet
        }
    }

    # stop if no genotypes
    if(!("geno" %in% names(output)))
       stop("No genotype data found.")

    # omit markers not in common among geno, gmap, pmap, founder_geno
    output <- drop_incomplete_markers( output )

    # pull out a map; make it numeric
    if("gmap" %in% names(output))
        map <- output$gmap
    else if("pmap" %in% names(output))
        map <- output$pmap
    else stop("Need a genetic or physical marker map")
    map[,2] <- as.numeric(map[,2])
    map <- reorder_map_table(map)

    # split maps by chr
    maps <- c("gmap", "pmap")
    for(section in maps) {
        if(section %in% names(output)) {
            if(!quiet) message(" - splitting up ", section)
            output[[section]] <- split_map(output[[section]], unique(map[,1]))
        }
    }

    # split genotypes by chromosome
    geno <- c("geno", "founder_geno")
    for(section in geno) {
        if(section %in% names(output)) {
            if(!quiet) message(" - splitting up ", section)
            output[[section]] <- split_geno(output[[section]], map)
        }
    }

    # X chr?
    chr <- names(output$geno)
    output$is_x_chr <- rep(FALSE, length(chr))
    names(output$is_x_chr) <- chr
    if("x_chr" %in% names(control)) {
        x_chr <- control$x_chr # name of X chromosome
        used_control["x_chr"] <- TRUE # indicate that we used it

        # deal with case of multiple chromosomes, and with not all present in genotypes
        found <- x_chr[x_chr %in% names(output$is_x_chr)]
        output$is_x_chr[found] <- TRUE
        notfound <- x_chr[!(x_chr %in% names(output$is_x_chr))]
        if(length(notfound) == 1) {
            warning('Chromosome "', notfound, '" not in genotype data but in control file as x_chr')
        } else if(length(notfound) > 1) {
            warning('Chromosomes ', paste0('"', notfound, '"', collapse=", "),
                    ' not in genotype data but in control file as x_chr')
        }
    }

    # sex
    output$is_female <- convert_sex(control$sex, output$covar, control$sep,
                                    control$comment.char, dir, quiet=quiet)
    if(is.null(output$is_female)) { # missing; assume all FALSE
        output$is_female <- rep(FALSE, nrow(output$geno[[1]]))
        names(output$is_female) <- rownames(output$geno[[1]])
    }
    used_control["sex"] <- TRUE # indicate that we used it

    # cross_info
    output$cross_info <- convert_cross_info(control$cross_info, output$covar, control$sep,
                                            control$comment.char, dir, quiet=quiet)
    if(is.null(output$cross_info)) { # missing; make a 0-column matrix
        output$cross_info <- matrix(0L, ncol=0, nrow=nrow(output$geno[[1]]))
        rownames(output$cross_info) <- rownames(output$geno[[1]])
    }
    used_control["cross_info"] <- TRUE # indicate that we used it

    # alleles?
    n_alleles <- nalleles(output$crosstype)
    if("alleles" %in% names(control)) {
        output$alleles <- control$alleles
        used_control["alleles"] <- TRUE # indicate that we used it
        if(n_alleles != length(output$alleles))
            stop("length(alleles) [", length(output$alleles),
                 "] != expected number [", n_alleles, "]")
    } else {
        output$alleles <- LETTERS[1:n_alleles]
    }

    class(output) <- "cross2"

    # force genotypes, is_female, and cross_info to be aligned
    gfc <- align_geno_sex_cross(output$geno,
                                output$is_female,
                                output$cross_info)
    output$geno <- gfc$geno
    output$is_female <- gfc$is_female
    output$cross_info <- gfc$cross_info

    # force same markers in geno, founder_geno, and cross_info
    for(i in seq(along=output$geno)) {
        gmar <- colnames(output$geno[[i]])
        mmar <- names(output$gmap[[i]])
        pmar <- names(output$pmap[[i]])
        if(is.null(pmar)) pmar <- mmar
        if(is.null(mmar)) mmar <- pmar
        fmar <- colnames(output$founder_geno[[i]])
        if(is.null(fmar)) fmar <- gmar

        mar <- find_common_ids(gmar, mmar, pmar, fmar)
        tmar <- length(unique(c(gmar, mmar, pmar, fmar)))
        if(length(mar) < tmar)
            warning("Omitting ", tmar - length(mar),
                    " markers on chr ", names(output$geno)[i],
                    " that ", ifelse(tmar-length(mar)==1, "is", "are"),
                    " missing position, genotypes, or founder genotypes.")

        output$geno[[i]] <- output$geno[[i]][,mar,drop=FALSE]
        output$gmap[[i]] <- output$gmap[[i]][mar]
        if(!is.null(output$pmap))
            output$pmap[[i]] <- output$pmap[[i]][mar]
        if(!is.null(output$founder_geno))
            output$founder_geno[[i]] <- output$founder_geno[[i]][,mar,drop=FALSE]
    }

    # if pheno and geno have same individuals, make same order
    if(!is.null(output$pheno)) {
        ind_g <- rownames(output$geno[[1]])
        ind_p <- rownames(output$pheno)
        if(length(ind_g) == length(ind_p)) {
            if(!all(ind_g == ind_p) && all(sort(ind_g) == sort(ind_p))) {
                output$pheno <- output$pheno[ind_g,,drop=FALSE]
            }
        }
    }
    # same with geno and covar
    if(!is.null(output$covar)) {
        ind_g <- rownames(output$geno[[1]])
        ind_p <- rownames(output$covar)
        if(length(ind_g) == length(ind_p)) {
            if(!all(ind_g == ind_p) && all(sort(ind_g) == sort(ind_p))) {
                output$covar <- output$covar[ind_g,,drop=FALSE]
            }
        }
    }

    check_cross2(output) # run all the checks

    if(any(!used_control))
        warning("Used control information: ",
                paste0('"', names(used_control)[!used_control], '"', collapse=", "))

    # add is_x_chr as attributes to gmap and pmap
    if("is_x_chr" %in% names(output)) {
        for(obj in c("gmap", "pmap")) {
            if(obj %in% names(output))
                attr(output[[obj]], "is_x_chr") <- output$is_x_chr
        }
    }

    output
}

# convert genotype data, using genotype encodings
# genotypes is a list with names = code in data
#                     and values = new numeric code
recode_geno <-
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
    for(g in gnames) {
        if(g != genotypes[[g]]) # if g==genotypes[[g]], don't need to switch
            newgeno[!is.na(geno) & geno==g] <- genotypes[[g]]
    }

    # turn missing values to 0s
    newgeno[is.na(newgeno)] <- "0"

    storage.mode(newgeno) <- "integer"

    newgeno
}

# convert a phenotype data.frame to a matrix of doubles
pheno2matrix <-
function(pheno)
{
    for(i in 1:ncol(pheno))
        pheno[,i] <- as.numeric(pheno[,i])
    pheno <- as.matrix(pheno)
    storage.mode(pheno) <- "double"

    pheno
}

# make the first column of a data frame the row names
firstcol2rownames <-
function(mat, name="")
{
    # assume the first column is the ID
    idcolumn <- mat[,1]

    # look for duplicates
    check4duplicates(idcolumn, name)

    # drop that column and make it the row names
    mat <- mat[,-1,drop=FALSE]
    rownames(mat) <- as.character(idcolumn)

    mat
}

# check a list of IDs for duplicates
check4duplicates <-
function(ids, name="")
{
    if(name != "") name <- paste(" in", name)
    dup <- duplicated(ids)
    if(any(dup)) {
        stop("Not all ids", name, " are unique: ",
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
    map <- map[mmark %in% gmark,,drop=FALSE] # toss from map any markers not in geno

    # sort map
    map <- map[order(factor(map[,1], levels=unique(map[,1])), map[,2], seq_len(nrow(map))),,drop=FALSE]

    geno <- geno[,rownames(map),drop=FALSE]

    chr <- map[,1]
    uchr <- unique(chr)
    newgeno <- vector("list", length(uchr))
    names(newgeno) <- uchr
    for(i in uchr) {
        newgeno[[i]] <- geno[,which(chr==i), drop=FALSE]
    }

    newgeno
}

# split a data frame with chr, pos (and marker names as rows)
#   into a list of chromosomes
split_map <-
function(map, chr_names=NULL)
{
    map <- reorder_map_table(map, chr_names=chr_names)

    pos <- as.numeric(map[,2])
    chr <- map[,1]
    uchr <- unique(chr)

    names(pos) <- rownames(map)

    lapply(split(pos, factor(chr, uchr)), sort)
}

# grab sex information
convert_sex <-
function(sex_control, covar, sep, comment.char, dir, quiet=TRUE)
{
    if(is.null(sex_control)) return(NULL)

    if("covar" %in% names(sex_control)) { # sex within the covariates
        if(!(sex_control$covar %in% colnames(covar)))
            stop('sex column "', sex_control$covar, '" not found in covar')
        sex <- covar[,sex_control$covar, drop=FALSE]
    }
    else if("file" %in% names(sex_control)) { # look for file
        if(!quiet) message(" - reading sex")
        file <- file.path(dir, sex_control$file)
        stop_if_no_file(file)
        sex <- read_csv(file, sep=sep, na.strings=NULL, comment.char=comment.char,
                        rownames_included=TRUE)
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
    codes <- convert_sexcodes(codes) # convert to 0 and 1
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

convert_sexcodes <-
function(codes)
{
    result <- match(substr(tolower(codes), 1, 1), c("f", "m")) - 1

    names(result) <- names(codes)
    result
}


# grab cross_info
convert_cross_info <-
function(cross_info_control, covar, sep, comment.char, dir, quiet=TRUE)
{
    if(is.null(cross_info_control)) return(NULL)

    if(!is.list(cross_info_control)) { # provided file name directly?
        if(file.exists(file.path(dir, cross_info_control)))
            cross_info_control <- list(file=cross_info_control)
    }

    if("file" %in% names(cross_info_control)) { # look for file
        if(!quiet) message(" - reading cross_info")

        file <- file.path(dir, cross_info_control$file)
        stop_if_no_file(file)
        cross_info <- read_csv(file, sep=sep, comment.char=comment.char, rownames_included=TRUE)

        if(any(is.na(cross_info)))
            stop(sum(is.na(cross_info)), " missing values in cross_info (cross_info can't be missing.")
    }
    else if("covar" %in% names(cross_info_control)) { # cross_info within the covariates
        if(!(cross_info_control$covar %in% colnames(covar)))
            stop('cross_info column "', cross_info_control$covar, '" not found in covar')
        cross_info <- covar[,cross_info_control$covar, drop=FALSE]
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

# is file on web (starts with http://, https://, or file://)
is_web_file <-
function(file)
{
    patterns <- c("^http://", "^https://", "^file://", "^ftp://")
    any(vapply(patterns, function(a,b) grepl(a, b), logical(1), file))
}

stop_if_no_file <-
function(filename)
{
    if(is_web_file(filename)) return(TRUE)
    if(!file.exists(filename))
        stop('file "', filename, '" does not exist.')
}


# read multiple CSV files and cbind results
read_mult_csv <-
    function(filenames, sep=",", na.strings=c("NA", "-"), comment.char="#", transpose=FALSE)
{
    result <- NULL
    for(file in filenames) {
        this <- read_csv(file, sep=sep, na.strings=na.strings,
                         comment.char=comment.char, transpose=transpose,
                         rownames_included=TRUE)

        if(is.null(result)) result <- this
        else result <- cbind_expand(result, this)
    }
    result
}

# read control file, as either YAML or JSON
read_control_file <-
function(filename)
{
    # ends in yaml?
    if(grepl("\\.yaml$", filename)) {
        control <- yaml::yaml.load_file(filename)
    }
    else if(grepl("\\.json$", filename)) {
        control <- jsonlite::fromJSON(readLines(filename))
    }
    else stop(paste('Control file', filename, 'should have extension ".yaml" or ".json"'))

    # default values for sep, na.strings, and comment.char
    if(is.null(control$sep)) control$sep <- ","
    if(is.null(control$na.strings)) control$na.strings <- "NA"
    if(is.null(control$comment.char)) control$comment.char <- "#"

    control
}

# read header lines
read_header <-
    function(filename, comment.char="#")
{
    con <- file(filename, "rt")
    on.exit(close(con))
    header <- NULL
    while(length(line <- readLines(con, 1))>0) {
        if(grepl(paste0("^", comment.char), line))
            header <- c(header, line)
        else break
    }
    header
}

# get expected dimensions from header
# (expecting header rows like "# nrow 5201" and "# ncol 59")
extract_dim_from_header <-
    function(header, comment.char="#")
{
    result <- c(NA, NA)
    nam <- c("nrow", "ncol")
    for(i in seq(along=nam)) {
        matchhead <- grepl(paste0("^", comment.char, "\\s*", nam[i], "\\s*\\d+"), header)
        if(sum(matchhead) > 1) # more than one matches
            stop('Multiple header lines with "', nam[i], '"')
        if(sum(matchhead) == 0) next
        spl <- strsplit(header[matchhead], "\\s+")[[1]]
        result[i] <- as.numeric(spl[length(spl)])
    }

    result
}

# drop markers that do not appear in all of geno, gmap, pmap, founder_geno
drop_incomplete_markers <-
    function( output )
{
    geno_mar <- colnames( output$geno )

    if(!is.null(output$gmap)) {
        pos <- suppressWarnings( as.numeric(output$gmap[,2]) )
        if(any(is.na(pos))) {
            n_na <- sum(is.na(pos))
            warning(n_na, ifelse(n_na==1, " marker", " markers"), " with missing genetic map position.")
        }
        gmap_mar <- rownames(output$gmap)[!is.na(pos)]
        output$gmap <- output$gmap[!is.na(pos),]
    }
    else gmap_mar <- geno_mar

    if(!is.null(output$pmap)) {
        pos <- suppressWarnings( as.numeric(output$pmap[,2]) )
        if(any(is.na(pos))) {
            n_na <- sum(is.na(pos))
            warning(n_na, ifelse(n_na==1, " marker", " markers"), " with missing physical map position.")
        }
        pmap_mar <- rownames(output$pmap)[!is.na(pos)]
        output$pmap <- output$pmap[!is.na(pos),]
    }
    else pmap_mar <- geno_mar

    if(!is.null(output$founder_geno))
        fg_mar <- colnames(output$founder_geno)
    else fg_mar <- geno_mar

    tot_mar <- unique(c(geno_mar, gmap_mar, pmap_mar, fg_mar))
    keep_mar <- find_common_ids(geno_mar, gmap_mar, pmap_mar, fg_mar)

    if(length(keep_mar) < length(tot_mar)) {
        n_drop <- length( tot_mar[!(tot_mar %in% keep_mar)] )
        warning("Omitting ", n_drop, " markers that are not in both genotypes and maps")

        output$geno <- output$geno[,geno_mar %in% keep_mar, drop=FALSE]
        if(!is.null(output$gmap))
            output$gmap <- output$gmap[gmap_mar %in% keep_mar, , drop=FALSE]
        if(!is.null(output$pmap))
            output$pmap <- output$pmap[pmap_mar %in% keep_mar, , drop=FALSE]
        if(!is.null(output$founder_geno))
            output$founder_geno <- output$founder_geno[, fg_mar %in% keep_mar, drop=FALSE]
    }

    output
}

reorder_map_table <-
    function(map_tab, chr_col=1, pos_col=2, chr_names=NULL)
{
    chr <- map_tab[,chr_col]

    if(is.null(chr_names)) chr_names <- unique(chr)

    chr <- factor(chr, levels=chr_names)
    pos <- suppressWarnings( as.numeric(map_tab[,pos_col]) )

    map_tab[ order(chr, pos, seq_len(nrow(map_tab))), , drop=FALSE]
}

# create unique temporary directory name
qtl2_temp_dir <-
    function(n_letters=15, initial_bit="qtl2")
{
    file.path(tempdir(),
              paste0(initial_bit, "_", paste(sample(c(letters, 0:9), n_letters, replace=TRUE), collapse="")) )

}

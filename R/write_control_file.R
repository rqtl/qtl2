# write_control_file
#' Write a control file for QTL data
#'
#' Write the control file (in [YAML](http://www.yaml.org) needed
#' by [read_cross2()] for a set of QTL data.
#'
#' @md
#'
#' @param output_file File name (with path) of the
#' [YAML](http://www.yaml.org) or [JSON](http://www.json.org) file to be created, as a character
#' string. If extension is `.json`, JSON format is used; otherwise, YAML is used.
#'
#' @param crosstype Character string with the cross type.
#' @param geno_file File name for genotype data.
#' @param founder_geno_file File name for the founder genotype data.
#' @param gmap_file File name for genetic map.
#' @param pmap_file File name for the physical map.
#' @param pheno_file File name for the phenotype data.
#' @param covar_file File name for the covariate data.
#' @param phenocovar_file File name for the phenotype covariate data
#' (i.e., metadata about the phenotypes).
#' @param sex_file File name for the individuals' sex. (Specify just
#' one of `sex_file` or `sex_covar`.)
#' @param sex_covar Column name in the covariate data that corresponds
#' to sex. (Specify just one of `sex_file` or `sex_covar`.)
#' @param sex_codes Named vector of character strings specifying the
#' encoding of sex. The names attribute should be the codes used in
#' the data files; the values within the vector should be
#' `"female"` and `"male"`.
#' @param crossinfo_file File name for the `cross_info` data. (Specify just
#' one of `crossinfo_file` or `crossinfo_covar`.)
#' @param crossinfo_covar Column name in the covariate data that
#' corresponds to the `cross_info` data. (Specify just one of
#' `crossinfo_file` or `crossinfo_covar`.)
#' @param crossinfo_codes In the case that `crossinfo_covar` is
#' provided (and not `crossinfo_file`; it would be left
#' untouched), a named vector of character strings specifying the
#' encoding of `cross_info`. The names attribute should be the
#' codes used in the covariate column; the values within the vector
#' should be the codes to which they will be converted (for example,
#' `0` and `1` for an intercross).
#' @param geno_codes Named vector specifying the encoding of
#' genotypes. The names attribute has the codes used within the
#' genotype and founder genotype data files; the values within the
#' vector should be the integers to which the genotypes will be
#' converted.
#' @param alleles Vector of single-character codes for the founder
#' alleles.
#' @param xchr Character string with the ID for the X chromosome.
#' @param sep Character string that separates columns in the data files.
#' @param na.strings Vector of character strings with codes to be
#' treated as missing values.
#' @param comment.char Character string that is used as initial
#' character in a set of leading comment lines in the data files.
#' @param geno_transposed If TRUE, genotype file is transposed (with markers as rows).
#' @param founder_geno_transposed If TRUE, founder genotype file is transposed (with markers as rows).
#' @param pheno_transposed If TRUE, phenotype file is transposed (with phenotypes as rows).
#' @param covar_transposed If TRUE, covariate file is transposed (with covariates as rows).
#' @param phenocovar_transposed If TRUE, phenotype covariate file is transposed (with phenotype covariates as rows).
#' @param description Optional character string describing the data.
#' @param comments Vector of character strings to be inserted as
#' comments at the top of the file (in the case of YAML), with each
#' string as a line. For JSON, the comments are instead included
#' within the control object.
#' @param overwrite If TRUE, overwrite file if it exists. If FALSE
#' (the default) and the file exists, stop with an error.
#'
#' @return (Invisibly) The data structure that was written.
#'
#' @details This function takes a set of parameters and creates the
#' control file (in [YAML](http://www.yaml.org) format) needed
#' for the new input data file format for
#' [R/qtl2](http://kbroman.org/qtl2).  See the
#' \href{http://kbroman.org/qtl2/pages/sampledata.html}{sample data
#' files} and the
#' \href{http://kbroman.org/qtl2/assets/vignettes/input_files.html}{vignette
#' describing the input file format}.
#'
#' @export
#' @keywords utilities
#' @seealso [read_cross2()], sample data files at
#' \url{http://kbroman.org/qtl2/pages/sampledata.html}
#' @examples
#' \dontrun{
#' # Control file for the sample dataset, grav2
#' write_control_file("~/my_data/grav2.yaml",
#'                    crosstype="riself",
#'                    geno_file="grav2_geno.csv",
#'                    gmap_file="grav2_gmap.csv",
#'                    pheno_file="grav2_pheno.csv",
#'                    phenocovar_file="grav2_phenocovar.csv",
#'                    geno_codes=c(L=1L, C=2L),
#'                    alleles=c("L", "C"),
#'                    na.strings=c("-", "NA"))
#'
#' # Control file for the sample dataset, iron
#' write_control_file("~/my_data/iron.yaml",
#'                    crosstype="f2",
#'                    geno_file="iron_geno.csv",
#'                    gmap_file="iron_gmap.csv",
#'                    pheno_file="iron_pheno.csv",
#'                    covar_file="iron_covar.csv",
#'                    phenocovar_file="iron_phenocovar.csv",
#'                    geno_codes=c(SS=1L, SB=2L, BB=3L),
#'                    sex_covar="sex",
#'                    sex_codes=c(f="female", m="male"),
#'                    crossinfo_covar="cross_direction",
#'                    crossinfo_codes=c("(SxB)x(SxB)"=0L, "(BxS)x(BxS)"=1L),
#'                    xchr="X",
#'                    alleles=c("S", "B"),
#'                    na.strings=c("-", "NA"))
#' }
write_control_file <-
function(output_file, crosstype=NULL, geno_file=NULL, founder_geno_file=NULL, gmap_file=NULL,
         pmap_file=NULL, pheno_file=NULL, covar_file=NULL, phenocovar_file=NULL,
         sex_file=NULL, sex_covar=NULL, sex_codes=NULL,
         crossinfo_file=NULL, crossinfo_covar=NULL, crossinfo_codes=NULL,
         geno_codes=NULL, alleles=NULL, xchr=NULL,
         sep=",", na.strings=c("-", "NA"), comment.char="#",
         geno_transposed=FALSE, founder_geno_transposed=FALSE,
         pheno_transposed=FALSE, covar_transposed=FALSE,
         phenocovar_transposed=FALSE,
         description=NULL, comments=NULL, overwrite=FALSE)
{
    output_file <- path.expand(output_file)
    if(!overwrite && file.exists(output_file))
        stop("The output file (", output_file, ") already exists. Remove it first (or use overwrite=TRUE).")

    result <- list(description="", # stub to be replaced or removed
                   comments="", # stub to be replaced or removed
                   crosstype=crosstype, sep=sep, na.strings=na.strings,
                   comment.char=comment.char)
    if(!is.null(description) && description!="") {
        paste(description, collapse="\n")
        result$description <- description
    } else {
        result$description <- NULL
    }

    if(!is.null(geno_file))
        result$geno <- geno_file
    if(!is.null(founder_geno_file))
        result$founder_geno <- founder_geno_file
    if(!is.null(gmap_file))
        result$gmap <- gmap_file
    if(!is.null(pmap_file))
        result$pmap <- pmap_file
    if(!is.null(pheno_file))
        result$pheno <- pheno_file
    if(!is.null(covar_file))
        result$covar <- covar_file
    if(!is.null(phenocovar_file))
        result$phenocovar <- phenocovar_file
    if(!is.null(alleles))
        result$alleles <- alleles
    if(!is.null(xchr))
        result$x_chr <- xchr
    if(!is.null(geno_codes)) {
        storage.mode(geno_codes) <- "integer"
        result$genotypes <- as.list(geno_codes)
    }

    # transposed file?
    if(!is.null(geno_transposed) && geno_transposed)
        result$geno_transposed <- geno_transposed
    if(!is.null(founder_geno_transposed) && founder_geno_transposed)
        result$founder_geno_transposed <- founder_geno_transposed
    if(!is.null(pheno_transposed) && pheno_transposed)
        result$pheno_transposed <- pheno_transposed
    if(!is.null(covar_transposed) && covar_transposed)
        result$covar_transposed <- covar_transposed
    if(!is.null(phenocovar_transposed) && phenocovar_transposed)
        result$phenocovar_transposed <- phenocovar_transposed

    # sex
    if(!is.null(sex_file)) {
        if(!is.null(sex_covar))
            stop("Specify just one of sex_file and sex_covar")
        if(is.null(sex_codes))
            stop("if sex_file is specified, sex_codes must be as well")

        result$sex <- c(list(file=sex_file),
                        as.list(sex_codes))
    }
    else if(!is.null(sex_covar)) {
        if(is.null(sex_codes))
            stop("if sex_covar is specified, sex_codes must be as well")
        result$sex <- c(list(covar=sex_covar),
                        as.list(sex_codes))
    }

    # cross_info
    if(!is.null(crossinfo_file)) {
        if(!is.null(crossinfo_covar))
            stop("Specify just one of crossinfo_file and crossinfo_covar")
        if(!is.null(crossinfo_codes))
            warning("if crossinfo_file is specified, crossinfo_codes is ignored")
        result$cross_info <- list(file=crossinfo_file)
    }
    else {
        if(!is.null(crossinfo_codes)) {
            if(is.null(crossinfo_covar))
                stop("if crossinfo_codes is provided, crossinfo_covar must also be provided")
            storage.mode(crossinfo_codes) <- "integer"
            result$cross_info <- c(list(covar=crossinfo_covar),
                                   as.list(crossinfo_codes))
        }
        else if(!is.null(crossinfo_covar)) {
            result$cross_info <- list(covar=crossinfo_covar)
        }
    }

    # JSON or YAML?
    if(grepl("\\.json$", output_file)) { # assume JSON
        if(is.null(comments))
            result$comments <- NULL
        else
            result$comments <- comments

        cat(jsonlite::toJSON(result, auto_unbox=TRUE, pretty=TRUE),
            file=output_file)
        cat("\n", file=output_file, append=TRUE) # add extra newline to avoid warning when reading
    }
    else {
        result$comments <- NULL # delete the placeholder

        # comments as a single string
        if(is.null(comments))
            comments <- ""
        else comments <- paste0("# ", comments, "\n", collapse="")

        # write data
        cat(comments, file=output_file)
        cat(yaml::as.yaml(result), file=output_file, append=TRUE)
    }
    invisible(result)
}

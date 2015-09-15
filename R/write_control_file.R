# write_control_file
#' Write a control file for QTL data
#'
#' Write the control file (in \href{http://www.yaml.org}{YAML} needed
#' by \code{\link{read_cross2}} for a set of QTL data.
#'
#' @param output_file File name (with path) of the
#' \href{http://www.yaml.org}{YAML} file to be created, as a character
#' string.
#' @param crosstype Character string with the cross type.
#' @param geno_file File name for genotype data.
#' @param foundergeno_file File name for the founder genotype data.
#' @param gmap_file File name for genetic map.
#' @param pmap_file File name for the physical map.
#' @param pheno_file File name for the phenotype data.
#' @param covar_file File name for the covariate data.
#' @param phenocovar_file File name for the phenotype covariate data
#' (i.e., metadata about the phenotypes).
#' @param sex_file File name for the individuals' sex. (Specify just
#' one of \code{sex_file} or \code{sex_covar}.)
#' @param sex_covar Column name in the covariate data that corresponds
#' to sex. (Specify just one of \code{sex_file} or \code{sex_covar}.)
#' @param sex_codes Named vector of character strings specifying the
#' encoding of sex. The names attribute should be the codes used in
#' the data files; the values within the vector should be
#' \code{"female"} and \code{"male"}.
#' @param crossinfo_file File name for the \code{cross_info} data. (Specify just
#' one of \code{crossinfo_file} or \code{crossinfo_covar}.)
#' @param crossinfo_covar Column name in the covariate data that
#' corresponds to the \code{cross_info} data. (Specify just one of
#' \code{crossinfo_file} or \code{crossinfo_covar}.)
#' @param crossinfo_codes In the case that \code{crossinfo_covar} is
#' provided (and not \code{crossinfo_file}; it would be left
#' untouched), a named vector of character strings specifying the
#' encoding of \code{cross_info}. The names attribute should be the
#' codes used in the covariate column; the values within the vector
#' should be the codes to which they will be converted (for example,
#' \code{0} and \code{1} for an intercross).
#' @param linemap_file File name for the \code{linemap} data
#' (indicating the mapping of individuals to lines). (Specify just one of
#' \code{linemap_file} or \code{linemap_covar}.)
#' @param linemap_covar  Column name in the covariate data that
#' corresponds to the \code{linemap} data. (Specify just one of
#' \code{linemap_file} or \code{linemap_covar}.)
#' @param geno_codes Named vector specifying the encoding of
#' genotypes. The names attribute has the codes used within the
#' genotype and founder genotype data files; the values within the
#' vector should be the integers to which the genotypes will be
#' converted.
#' @param alleles Vector of single-character codes for the founder
#' alleles.
#' @param xchr Character string with the ID for the X chromosome.
#' @param na.strings Vector of character strings with codes to be
#' treated as missing values.
#' @param geno_transposed If TRUE, genotype file is transposed (with markers as rows).
#' @param foundergeno_transposed If TRUE, founder genotype file is transposed (with markers as rows).
#' @param pheno_transposed If TRUE, phenotype file is transposed (with phenotypes as rows).
#' @param covar_transposed If TRUE, covariate file is transposed (with covariates as rows).
#' @param phenocovar_transposed If TRUE, phenotype covariate file is transposed (with phenotype covariates as rows).
#' @param comments Vector of character strings to be inserted as
#' comments at the top of the file, with each string as a line.
#'
#' @return (Invisibly) The data structure that was written.
#'
#' @details This function takes a set of parameters and creates the
#' control file (in \href{http://www.yaml.org}{YAML} format) needed
#' for the new input data file format for
#' \href{http://kbroman.org/qtl2}{R/qtl2}.  See the
#' \href{http://kbroman.org/qtl2/pages/sampledata.html}{sample data
#' files} and the
#' \href{http://kbroman.org/qtl2/assets/vignettes/input_files.html}{vignette
#' describing the input file format}.
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{read_cross2}}, sample data files at
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
function(output_file, crosstype, geno_file, foundergeno_file, gmap_file,
         pmap_file, pheno_file, covar_file, phenocovar_file,
         sex_file, sex_covar, sex_codes,
         crossinfo_file, crossinfo_covar, crossinfo_codes,
         linemap_file, linemap_covar, geno_codes, alleles, xchr,
         na.strings=c("-", "NA"),
         geno_transposed=FALSE, foundergeno_transposed=FALSE,
         pheno_transposed=FALSE, covar_transposed=FALSE,
         phenocovar_transposed=FALSE,
         comments)
{
    output_file <- path.expand(output_file)
    if(file.exists(output_file))
        stop("The output file (", output_file, ") already exists. Please remove it first.")

    result <- list(crosstype=crosstype, na.strings=na.strings)

    if(!missing(geno_file))
        result$geno <- geno_file
    if(!missing(foundergeno_file))
        result$foundergeno <- foundergeno_file
    if(!missing(gmap_file))
        result$gmap <- gmap_file
    if(!missing(pmap_file))
        result$pmap <- pmap_file
    if(!missing(pheno_file))
        result$pheno <- pheno_file
    if(!missing(covar_file))
        result$covar <- covar_file
    if(!missing(phenocovar_file))
        result$phenocovar <- phenocovar_file
    if(!missing(alleles))
        result$alleles <- alleles
    if(!missing(xchr))
        result$x_chr <- xchr
    if(!missing(geno_codes)) {
        storage.mode(geno_codes) <- "integer"
        result$genotypes <- as.list(geno_codes)
    }

    # transposed file?
    if(!missing(geno_transposed) && geno_transposed)
        result$geno_transposed <- geno_transposed
    if(!missing(foundergeno_transposed) && foundergeno_transposed)
        result$foundergeno_transposed <- foundergeno_transposed
    if(!missing(pheno_transposed) && pheno_transposed)
        result$pheno_transposed <- pheno_transposed
    if(!missing(covar_transposed) && covar_transposed)
        result$covar_transposed <- covar_transposed
    if(!missing(phenocovar_transposed) && phenocovar_transposed)
        result$phenocovar_transposed <- phenocovar_transposed

    # sex
    if(!missing(sex_file) && !is.null(sex_file)) {
        if(!missing(sex_covar) && !is.null(sex_covar))
            stop("Specify just one of sex_file and sex_covar")
        if(missing(sex_codes) || is.null(sex_codes))
            stop("if sex_file is specified, sex_codes must be as well")

        result$sex <- c(list(file=sex_file),
                        as.list(sex_codes))
    }
    else if(!missing(sex_covar) && !is.null(sex_covar)) {
        if(missing(sex_codes) || is.null(sex_codes))
            stop("if sex_covar is specified, sex_codes must be as well")
        result$sex <- c(list(covar=sex_covar),
                        as.list(sex_codes))
    }

    # cross_info
    if(!missing(crossinfo_file) && !is.null(crossinfo_file)) {
        if(!missing(crossinfo_covar) && !is.null(crossinfo_covar))
            stop("Specify just one of crossinfo_file and crossinfo_covar")
        if(!missing(crossinfo_codes) && is.null(crossinfo_codes))
            warning("if crossinfo_file is specified, crossinfo_codes is ignored")
        result$cross_info <- list(file=crossinfo_file)
    }
    else if(!missing(crossinfo_covar) && !is.null(crossinfo_covar)) {
        if(missing(crossinfo_codes) || is.null(crossinfo_codes))
            stop("if crossinfo_covar is specified, crossinfo_codes must be as well")
        storage.mode(crossinfo_codes) <- "integer"
        result$cross_info <- c(list(covar=crossinfo_covar),
                               as.list(crossinfo_codes))
    }

    # linemap
    if(!missing(linemap_file) && !is.null(linemap_file)) {
        if(!missing(linemap_covar) && !is.null(linemap_covar))
            stop("Specify just one of linemap_file and linemap_covar")
        result$linemap <- linemap_file
    }
    else if(!missing(linemap_covar) && !is.null(linemap_covar)) {
        result$linemap <- linemap_covar
    }

    # comments as a single string
    if(missing(comments) || is.null(comments))
        comments <- ""
    else comments <- paste0("# ", comments, "\n", collapse="")

    # write data
    cat(comments, file=output_file)
    cat(yaml::as.yaml(result), file=output_file, append=TRUE)

    invisible(result)
}

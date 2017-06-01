core <- c("qtl2geno", "qtl2scan", "qtl2plot")

qtl2_attach <- function() {

    for(pkg in core) {
        message("Loading ", pkg)
        library(pkg, character.only=TRUE)
    }

}

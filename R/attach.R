core <- c("qtl2geno", "qtl2scan", "qtl2plot")

qtl2_attach <- function() {

    lapply(core, library, character.only=TRUE)

}

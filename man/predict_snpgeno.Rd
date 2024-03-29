% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_snpgeno.R
\name{predict_snpgeno}
\alias{predict_snpgeno}
\title{Predict SNP genotypes}
\usage{
predict_snpgeno(cross, geno, cores = 1)
}
\arguments{
\item{cross}{Object of class \code{"cross2"}. For details, see the
\href{https://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.}

\item{geno}{Imputed genotypes, as a list of matrices, as from \code{\link[=maxmarg]{maxmarg()}}.}

\item{cores}{Number of CPU cores to use, for parallel calculations.
(If \code{0}, use \code{\link[parallel:detectCores]{parallel::detectCores()}}.)
Alternatively, this can be links to a set of cluster sockets, as
produced by \code{\link[parallel:makeCluster]{parallel::makeCluster()}}.}
}
\value{
A list of matrices with inferred SNP genotypes, coded 1/2/3.
}
\description{
Predict SNP genotypes in a multiparent population from inferred genotypes plus founder strains' SNP alleles.
}
\examples{
\dontrun{
# load example data and calculate genotype probabilities
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/main/DOex/DOex.zip")
DOex <- read_cross2(file)
probs <- calc_genoprob(DOex, error_prob=0.002)

# inferred genotypes
m <- maxmarg(probs, minprob=0.5)

# inferred SNP genotypes
inferg <- predict_snpgeno(DOex, m)
}

}
\seealso{
\code{\link[=maxmarg]{maxmarg()}}, \code{\link[=viterbi]{viterbi()}}, \code{\link[=calc_errorlod]{calc_errorlod()}}
}
\keyword{utilities}

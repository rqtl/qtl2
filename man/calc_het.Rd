% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_het.R
\name{calc_het}
\alias{calc_het}
\title{Calculate heterozygosities}
\usage{
calc_het(probs, by = c("individual", "marker"), omit_x = TRUE, cores = 1)
}
\arguments{
\item{probs}{List of arrays of genotype probabilities, as
calculated by \code{\link[=calc_genoprob]{calc_genoprob()}}.}

\item{by}{Whether to summarize by individual or marker}

\item{omit_x}{If TRUE, omit the X chromosome.}

\item{cores}{Number of CPU cores to use, for parallel calculations.
(If \code{0}, use \code{\link[parallel:detectCores]{parallel::detectCores()}}.)
Alternatively, this can be links to a set of cluster sockets, as
produced by \code{\link[parallel:makeCluster]{parallel::makeCluster()}}.}
}
\value{
The result is a vector of estimated heterozygosities
}
\description{
Calculate heterozygosites, by individual or by marker
}
\examples{
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
p <- calc_genoprob(iron, err=0.002)

# heterozygosities by individual
het_ind <- calc_het(p)

# heterozygosities by marker
het_mar <- calc_het(p, "marker")

}

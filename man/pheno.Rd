\name{pheno}
\alias{pheno}
\docType{data}
\title{
A dummy phenotype dataset that provides an example of the input phenotype data used by the package
}
\description{
A dummy dataset of 3 quantitative phenotypes measured in 150 individuals. The data has been randomly generated, for instructional purposes only, and does not yield a significant association with the genotype. The data is a matrix, with one phenotype for each column. 
}
\format{
  A matrix with 150 phenotype observations.
 \describe{
    \item{\code{testPhneo1}}{a numeric vector}
    \item{\code{testPhneo2}}{a numeric vector}
    \item{\code{testPhneo3}}{a numeric vector}
  }
}
\details{
Please note the following IMPORTANT issue: the 'pheno' matrix has both column names and row names! the column names MUST be the names of the phenotypes and the row names MUST be the codes representing each individual in the pheno matrix, one individual for each row. Both row names and column names are extracted by the main function and are therefore mandatory}
\examples{
data(pheno)
head(pheno)
dimnames(pheno)[[1]] # the row names
dimnames(pheno)[[2]] # the column names
}
\keyword{datasets}

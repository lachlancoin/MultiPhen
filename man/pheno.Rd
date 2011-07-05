\name{pheno}
\alias{pheno}
\docType{data}
\title{
A dummy dataset used as an example of the data used by the package
}
\description{
A dummy dataset of 3 phenotypes scored in 150 individuals. The data has been randomly generated, for instructional purposes only, and does not yield a significant association with the genotype. The data is a matrix, with one phenotype for each column. 
}
\usage{data(pheno)}
\format{
  The format is:
 num [1:150, 1:3] -0.4101 0.4434 -0.4763 -0.2449 -0.0915 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:150] "897211" "96511" "69611" "355811" ...
  ..$ : chr [1:3] "testPheno1" "testPheno2" "testPheno3"
}
\details{
Please note the following IMPORTANT issue: the pheno matrix has both column names and row names! the column names MUST be the names of the phenotypes and the row names MUST be the codes representing each individual in the pheno matrix, one individual for each row. Both row names and column names are extracted by the main function and are therefore mandatory}
\examples{
data(pheno)
head(pheno)
dimnames(pheno)[[1]] # the row names
dimnames(pheno)[[2]] # the column names
}
\keyword{datasets}

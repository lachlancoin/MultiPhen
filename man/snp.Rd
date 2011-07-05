\name{snp}
\alias{snp}
\docType{data}
\title{
Test genotype data
}
\description{
A toy genotype dataset of one single SNP, as a matrix of 1 column and 150 rows. The genotype is in the 0/1/2 format. The data has been randomly generated, for instructional purposes only. 
}
\usage{data(snp)}
\format{
  A data frame with 150 genotype observations.
  \describe{
    \item{\code{rsID}}{a numeric vector}
  }
}
\details{
The 150 genotypes in `snp' correspond to the 150 phenotypes in `pheno', i.e. one individual for each line. Please note the following important points: genotype data must be a matrix format, with one row for each individual and as many columns for each genotype. In case of one single genotype the data must still conform to this format, as a matrix of as many rows as individuals and one single column for the one genotype present. A second important point is that the column names must be the rsID of the SNP for genotypes in the 0/1/2 format. Further options concerning genotype format will be available in further releases but are yet not documented.
}
\examples{
data(snp)
dim(snp)
colnames(snp)
}
\keyword{datasets}

\name{MultiPhen-package}
\alias{MultiPhen-package}
\alias{MultiPhen}
\docType{package}
\title{
MultiPhen, a package for the genetic association testing of multiple phenotypes
}
\description{
MultiPhen is the package for the function mPhen, which performs association testing between genetic variants (SNPs/CNVs) and multiple phenotypes. The primary purpose is for modelling and testing multiple phenotypes jointly by performing an ordinal regression where SNPs are treated as the outcome and multiple phenotypes are predictors; this can have large increases in statistical power to detect genotype-phenotype associations over the univariate approach. However, mPhen can also be used to perform standard univariate linear regression (SNP as predictor) and univariate ordinal regression (SNP as outcome) on the phenotypes under study. mPhen can be applied to genotyped or imputed data.
}
\details{
\tabular{ll}{
Package: \tab MultiPhen\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2011-07-04\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Lachlan Coin, Paul O'reilly and Federico Calboli

Maintainer: Federico Calboli <f.calboli@imperial.ac.uk>}
\references{
O'Reilly et al. 2011. MultiPhen: Joint model of multiple phenotypes increases discovery in GWAS. Under review.
}
\keyword{ package }
\keyword{ GWAS}
\examples{
data(snp)
data(pheno)
mPhen(snp, pheno, phenotypes = c('testPheno1', 'testPheno2'), covariates = 'testPheno3') 
}

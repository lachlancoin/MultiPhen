\name{mPhen.plink}
\alias{mPhen.plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs mPhen on data in PLINK binary format
}
\description{
mPhen.plink opens a connection to a binary PLINK file and then streams in every marker, one at a time, running \code{mPhen} on it and writing the output to an external (to the R session) results file.  It is MANDATORY that in the directory containing the binary file also reside two accesory files, with the same name as the binary file but with extensions .fam and .bim, both produced by PLINK.
}
\usage{
mPhen.plink(root, results, phenoData, phenotypes = dimnames(phenoData)[[2]], 
            covariates = NULL, resids = NULL, strats = NULL, maf_thresh = 0.001, corThresh = 0, 
            inverseRegress = FALSE, JointModel = TRUE, multiGen = FALSE, 
            fillMissingPhens = FALSE, scoreTest = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root}{
filename of the dataset in PLINK binary format, WITHOUT the .bed extension.  
}
  \item{results}{
name of the file storing (as a text file) the results of \code{mPhen} on each locus of the ``.bed'' file.
}
  \item{phenoData}{
Data frame containing phenotype data, where each row corresponds to an individual and row.names are individual IDs. Each column contains data on a certain phenotype across the sample of individuals (can be quantitative, case/control or ordinal. Must be numeric); the column header provides the phenotype name. An example is provided by 'pheno'.}
  \item{phenotypes}{
Vector of phenotypes, from phenoData, to be tested (Default is all phenotypes).
}
  \item{covariates}{
Vector of phenotypes, from phenoData, to be considered as covariates to be controlled for in the regression (Default is no covariates).
}
  \item{resids}{
Vector of residuals, from phenoData, alternative way to adjust for covariates, which pre-calculates offset terms to use in the per SNP regression (Default is no residuals).
}
  \item{strats}{
Statification vector (i.e. cases/controls, exposed/not exposed, etc), from phenoData (Default is no stratification).
} 
  \item{maf_thresh}{
Only SNPs with minor allele frequency greater than this value (a float between 0 and 0.5) will be tested (Default is 0.01).
}
  \item{corThresh}{
Only phenotypes that have a pairwise Pearson's correlation coefficient greater than this value (a float between 0 and 1) will be tested (Default is 0).
}
  \item{inverseRegress}{
If TRUE then SNP is treated as the outcome in the (then Ordinal) regression, if FALSE then SNP is predictor and (one) phenotype is the outcome in the usual way. Default is FALSE.
}
  \item{JointModel}{
If TRUE then the phenotypes are modelled and tested jointly, with SNP as outcome and the phenotypes as predictors (overrides inverseREGRESS=FALSE); a likelihood ratio test of the full model vs the model without phenotypes is performed (F-test), providing a P-value to test the association between the genotypes and phenotypes at each marker (usual genome-wide significance threshold, 5e-8, is appropriate). Betas and P-values (from Wald test) relating to the association of each phenotype with the SNP from the model are also provided (NB. These are not the same as the corresponding univariate results since they are produced by the joint model). If FALSE the phenotypes are tested univariately (the order of the regression is then dependent on inverseREGRESS). Default is TRUE.
}
  \item{multiGen}{
Not yet implemented
}
  \item{fillMissingPhens}{
If TRUE missing phenotype data are imputed according to the mean of the phenotype, if FALSE then individuals with missing phenotype data are excluded (Default is FALSE).
}
  \item{scoreTest}{
If TRUE the score test is performed rather than the likelihood ratio test (Default is FALSE).
}
%  \item{imputed}{
%%     ~~Describe \code{imputed} here~~
%}
}
\details{
Please note that, if the binary file is listed a ``mydata.bed'', the filename is ``mydata'', and the extension is ``.bed''.  In this case ``mydata'' would be used as \code{root} value.  Also \code{results} must be the name of a file.  
}
\value{
Technically, \code{mPhen.plink} does not return anything to screen, but saves all the results to an external file, to avoid creating a gigantic results file that would eat up all the memory.
}
\references{
See the references in \code{\link{mPhen}} and \code{\link{read.plink}}
}
\author{
Federico Calboli \email{f.calboli@imperial.ac.uk}}

\seealso{
\code{\link{mPhen}}.
}
\examples{
## given a directory containing the following files:
## ~\$> ls
## chrom.bed  chrom.bim  chrom.dat  chrom.fam
## running mPhen on the genetic data in chrom.bed is achieved by:
data(pheno)

\dontrun{setwd("/path/to/the/binary/file/")}
\dontrun{mPhen.plink("chrom", "chrom-results.txt", pheno, 
        phenotypes = c('testPheno1', 'testPheno2'), 
        covariates = 'testPheno3',  JointModel = TRUE)}

## alternatively, without using setwd()

\dontrun{mPhen.plink("/path/to/the/binary/file/chrom", 
        "/path/to/the/binary/file/chrom-results.txt", 
        pheno, phenotypes = c('testPheno1', 'testPheno2'), 
        covariates = 'testPheno3',  JointModel = TRUE)}

## both methods would read the data in chrom.bed into the object called my.chrom
## under NO CIRCUMSTANCES you should run somethig like
\dontrun{results = mPhen.plink("chrom", "chrom-results.txt", pheno, 
        phenotypes = c('testPheno1', 'testPheno2'), 
         covariates = 'testPheno3',  JointModel = TRUE)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
\keyword{ file }
\keyword{ connection }

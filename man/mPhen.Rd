\name{mPhen}
\alias{mPhen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for the genetic association testing of multiple phenotypes
}
\description{
mPhen performs association testing between genetic variants (SNPs; CNVs to be added soon) and multiple phenotypes. The primary purpose is for modelling and testing multiple phenotypes jointly by performing an ordinal regression where SNPs are treated as the outcome and multiple phenotypes are predictors; this can have large increases in statistical power to detect genotype-phenotype associations over the univariate approach (method described in O'Reilly et al. 2011). However, mPhen can also be used to perform standard univariate linear regression (SNP as predictor) and univariate ordinal regression (SNP as outcome) on the phenotypes under study. mPhen can be applied to directly genotyped or imputed data.
}
\usage{
mPhen(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], 
      covariates = NULL, maf_thresh = 0.001, corThresh = 0, 
      inverseRegress = FALSE, JointModel = TRUE, multiGen = FALSE, 
      fillMissingPhens = FALSE, scoreTest = FALSE, imputed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoData}{
Data frame containing genotype data, where each row corresponds to an individual and row.names are individual IDs. When the genotype data are directly genotyped SNPs each column corresponds to a SNP, where the header is rsID and each data element is coded 0 (homozygote AA), 1 (heterozygote AB), or 2 (homozygote BB; A and B correspond, arbitrarily, to the two alleles) reflecting the genotype for each individual at the marker. If the genotype data contains imputed genotypes, then the data frame has three columns per SNP, with the header of the rsID suffixed by -0, -1, -2 corresponding to each of the three genotypes; the elements of the data frame are then probabilities (decimals/floats between 0 and 1) reflecting the probability estimates that the individual takes each genotype at the SNP according to imputation (for directly genotyped SNPs one genotype will have probability 1 and the other two 0). Examples are provided by 'snps' and 'snps.imputed'.
}
  \item{phenoData}{
Data frame containing phenotype data, where each row corresponds to an individual and row.names are individual IDs. Each column contains data on a certain phenotype across the sample of individuals (can be quantitative, case/control or ordinal. Must be numeric); the column header provides the phenotype name. An example is provided by 'pheno'. 
}
  \item{phenotypes}{
Vector of phenotypes, from phenoData, to be tested (Default is all phenotypes).
}
  \item{covariates}{
Vector of phenotypes, from phenoData, to be considered as covariates to be controlled for in the regression (Default is no covariates).
}
  \item{maf_thresh}{
Only SNPs with minor allele frequency greater than this value (a float between 0 and 0.5) will be tested (Default is 0.01).
}
  \item{corThresh}{
Only phenotypes that have a pairwise Pearson's correlation coefficient greater than this value (a float between 0 and 1) will be tested (Default is 0).
}
%  \item{rescale}{
%Integer scale factor used when expandData=TRUE (see below) by which to expand the phenotype and genotype data, and to convert genotype probabilities to integer weights (after rounding), in order to perform an ordinal regression that accounts for uncertainty in genotypes. Thus if rescale=10 and the genotype probabilities for a given individual at a certain SNP are 0.22,0.35,0.43 then the integer weights will be 2,4,4; this means that of the 10 replicates of the phenotype data (produced by expanding the data), 2 correspond to a genotype coded 0 at the SNP, 4 to a genotype coded 1, and 4 to a genotype coded 2. Default is 10.
%}
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
%  \item{exactTest}{
%Not yet implemented
%}
%  \item{exactMethod}{
%Defaults to the Wald test.
%}
%  \item{expandData}{
%Not yet implemented
%}
  \item{imputed}{
Set TRUE if the genotype data is imputed, FALSE otherwise
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A matrix where the first column contains the regression coefficients (betas) corresponding to the specified SNP/phenotype (if inverseREGRESS=FALSE then the beta corresponds to the SNP, otherwise the beta corresponds to the phenotype), the second column contains P-values for association between the SNP and phenotype. If MultiPhenTest=TRUE then there is also a row for each SNP containing the P-value for the multiple phenotype association test (the beta in this row is NA).
}
\references{
O'Reilly et al. 2011. MultiPhen: Joint model of multiple phenotypes increases discovery in GWAS. Under review.}
\author{
Lachlan Coin, Federico Calboli, Clive Hoggart, Paul O'Reilly, Yotsawat Pomyen.

Maintainer, Federico Calboli \email{f.calboli@imperial.ac.uk} 
}
\note{
The user should remember that the genotype data file is always a matrix of at least a column, hence if taking a subset of 1 SNP in the non-imputed genotype data matrix, the option drop = FALSE should be used (see the example below)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(snps); data(snps.imputed); data(pheno)

mPhen(snps, pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = 'testPheno3') 
# performs the standard regression (SNP as predictor) 
# separately on the first two phenotypes, with the third phenotype as a covariate.

mPhen(snps, pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = 'testPheno3',  JointModel = TRUE) 
# tests association of SNP genotypes and the first two phenotypes 
# in the joint phenotype model (SNP as outcome in an ordinal regression),
# with the third phenotype as a covariate.

mPhen(snps[,2, drop = FALSE], pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = 'testPheno3',  JointModel = TRUE) 
# please note the use use of drop = FALSE if analysing only one SNP

mPhen(snps.imputed, pheno, phenotypes = c('testPheno1', 'testPheno2'), 
      covariates = 'testPheno3',  JointModel = TRUE, imputed = TRUE) 
# for imputed data


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple phenotypes } 
\keyword{ association }
% __ONLY ONE__ keyword per line

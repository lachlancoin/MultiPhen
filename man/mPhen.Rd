\name{mPhen}
\alias{mPhen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mPhen, a function for  the genetic association testing of multiple phenotypes
}
\description{
mPhen performs association testing between genetic variants (SNPs/CNVs) and multiple phenotypes. The primary purpose is for modelling and testing multiple phenotypes jointly by performing an ordinal regression where SNPs are treated as the outcome and multiple phenotypes are predictors; this can have large increases in statistical power to detect genotype-phenotype associations over the univariate approach (method described in O'Reilly et al. 2011). However, mPhen can also be used to perform standard univariate linear regression (SNP as predictor) and univariate ordinal regression (SNP as outcome) on the phenotypes under study. mPhen can be applied to genotyped or imputed data.
}
\usage{
mPhen(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL, maf_thresh = 0.001, corThresh = 0, rescale = 10, inverseRegress = FALSE, multiPhenTest = FALSE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE, exactTest = FALSE, exactMethod = "wald", expandData = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoData}{
Data frame containing genotype data, where each row corresponds to an individual and row.names are individual IDs. When the genotype data are directly genotyped SNPs each column corresponds to a SNP, where the header is rsID and each data element is coded 0 (homozygote AA), 1 (heterozygote AB), or 2 (homozygote BB) reflecting the genotype for each individual at the marker. If the genotype data contains imputed genotypes, then the data frame has three columns per SNP, with the header of the rsID suffixed by -0, -1, -2 corresponding to each of the three genotypes; the elements of the data frame are then probabilities (decimals/floats between 0 and 1) reflecting the probability estimates that the individual takes each genotype at the SNP according to imputation (for directly genotyped numbers one genotype will have probability 1 and the other two 0). ??Raw geno data **ToDo**
}
  \item{phenoData}{
Data frame containing phenotype data, where each row corresponds to an individual and row.names are individual IDs. Each column contains data on a certain phenotype across the sample of individuals (can be quantitative, case/control or ordinal. Must be numeric); the column header provides the phenotype name. 
}
  \item{phenotypes}{
Vector of phenotypes to be tested from phenoData (Default is all phenotypes from phenoData).
}
  \item{covariates}{
Vector of phenotypes to be considered as covariates to be controlled for in the regression, chosen from phenoData (Default is no covariates).
}
  \item{maf_thresh}{
Only SNPs with minor allele frequency greater than this value (a float between 0 and 0.5) will be tested (default is 0.01).
}
  \item{corThresh}{
Only phenotypes that have a pairwise Pearson's correlation coefficient above this value (a float between 0 and 1) will be tested (default it 0).
}
  \item{rescale}{
just how much to multiply the weights by when analysing imputed data. Its set at 10 as a default, basically that means if the imputed genotype probability is at  0.91 , the weight will be set to 10.  The reason for this is that ordinal regression has to have integer weights.
}
  \item{inverseRegress}{
If TRUE then SNP is treated as the outcome in the (then Ordinal) regression, if FALSE SNP is predictor and (one) phenotype is the outcome in the usual way.
}
  \item{multiPhenTest}{
if TRUE the phenotypes are modelled and tested jointly, with SNP as outcome and the phenotypes as predictors (overrides inverseREGRESS=FALSE); a likelihood ratio test of the full model vs the model without phenotypes is performed (F-test), providing a P-value to test the association between the genotypes and phenotypes at each marker (usual genome-wide significance threshold, 5e-8, is appropriate). If FALSE the phenotypes are tested univariately (the order of the regression is then dependent on inverseREGRESS).
}
  \item{multiGen}{
Not yet implemented
}
  \item{fillMissingPhens}{
If TRUE missing phenotype data are imputed according to the mean of the phenotype, if FALSE then individuals with missing phenotype data are excluded (default is FALSE).
}
  \item{scoreTest}{
If TRUE the score test is performed rather than the likelihood ratio test (default is FALSE).
}
  \item{exactTest}{
Not yet implemented
}
  \item{exactMethod}{
Defaults to the Wald's test. Future changes are in the pipeline but not yet implemented
}
  \item{expandData}{
Not yet implemented
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix where the first column contains the regression coefficients (betas) corresponding to the specified SNP/phenotype (if inverseREGRESS=FALSE the beta corresponds to the SNP, otherwise the beta corresponds to the phenotype), the second column contains P-values for association between the SNP and phenotype. If MultiPhenTest=TRUE then there is a row for each SNP containing the P-value for the multiple phenotype association test (the beta in this row is NA).
}
\references{
O'Reilly et al. 2011. MultiPhen: Joint model of multiple phenotypes increases discovery in GWAS. Under review.}
\author{
Lachlan Coin, Federico Calboli, Clive Hoggart, Paul O'Reilly, Yotsawat Pomyen, 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(snp)
data(pheno)
mPhen(snp, pheno, phenotypes = c('testPheno1', 'testPheno2'), covariates = 'testPheno3') # test of association of SNP genotype with the first two phenotypes, using the third as a covariate

mPhen(snp, pheno, phenotypes = c('testPheno1', 'testPheno2'), covariates = 'testPheno3',  multiPhenTest = TRUE) # adding a test for the joint regression of the two phenotypes
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ GWAS }
\keyword{ regression }
\keyword{ multiple association } 
% __ONLY ONE__ keyword per line

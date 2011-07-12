mPhen <-
function(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL,maf_thresh = 0.001, corThresh = 0.0, inverseRegress = FALSE, multiPhenTest = FALSE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE, imputed = FALSE){
  rescale = 1
  exactTest = FALSE
  exactMethod = "wald"
  expandData = FALSE
  if(imputed){
    genoData = apply(matrix(1:dim(genoData)[2],nrow = 3),2, function(x){apply(genoData[,x],1,.mod)})
    if(inverseRegress) genoData = apply(genoData,c(1,2),round)
    dimnames(genoData) = list(NULL,rep('rsID', dim(genoData)[2]))
  }
  mPhenInternal = function(x){
    x = as.matrix(x)
    dimnames(x) = list(1:dim(genoData)[1], 'gtype')
    .mPhen(x, phenoData, phenotypes = phenotypes, covariates = covariates, maf_thresh = maf_thresh, corThresh = corThresh, inverseRegress = inverseRegress, multiPhenTest = multiPhenTest, multiGen = multiGen, fillMissingPhens = fillMissingPhens, scoreTest = scoreTest)
  }
  res = tapply(genoData, col(genoData), mPhenInternal)
  res
}


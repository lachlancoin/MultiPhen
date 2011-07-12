mPhen <-
function(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL, resids = NULL, maf_thresh = 0.001, corThresh = 0.0, inverseRegress = FALSE, JointModel = TRUE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE, imputed = FALSE){
  testgeno = !is.matrix(genoData)
  if(testgeno) {
    stop('ERROR! The genetic data is NOT in matrix format. Please reformat the genetic data as a matrix and rerun the test')
  }
  testpheno = !is.matrix(phenoData)
  if(testpheno){
    stop('ERROR! the phenotype data is NOT in matrix format. Please reformat. NOTE, if testing one single phenotype remember to add drop = F (see ?Extract)')
  }
  rescale = 1
  exactTest = FALSE
  exactMethod = "wald"
  expandData = FALSE
  if(multiGen){
    stop('ERROR! multiGen not yet implemented, resetting to FALSE')
    multiGen = FALSE
  }
  if(imputed){
    snps.names = colnames(genoData)
    snps.names = snps.names[seq(1, length(snps.names), by = 3)]
    if(expandData) weights = as.vector(t(apply(genoData,1,.mod1, rescale)))
    genoData = apply(matrix(1:dim(genoData)[2],nrow = 3),2, function(x){apply(genoData[,x],1,.mod)})
    if(inverseRegress) genoData = apply(genoData,c(1,2),round)
    dimnames(genoData)[[2]] = snps.names 
  }
  mPhenInternal = function(x){
    x = as.matrix(x)
    dimnames(x) = list(1:dim(genoData)[1], 'gtype')
    .mPhen(x, phenoData, phenotypes = phenotypes, covariates = covariates, resids = resids, maf_thresh = maf_thresh, corThresh = corThresh, inverseRegress = inverseRegress, JointModel = JointModel, multiGen = multiGen, fillMissingPhens = fillMissingPhens, scoreTest = scoreTest)
  }
  res = tapply(genoData, col(genoData), mPhenInternal)
  dimnames(res)[[1]] = colnames(genoData)
  res
}

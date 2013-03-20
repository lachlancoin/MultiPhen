mPhen.plink <-
function(root, results, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL, resids = NULL, strats = NULL, maf_thresh = 0.001, corThresh = 0.0, inverseRegress = FALSE, JointModel = TRUE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE){
  bed.file = paste(root, '.bed', sep = '')   
  bed.file.size = file.info(bed.file)$size # the bed file size in bytes
  sample.size = dim(read.table(paste(root, '.fam', sep = '')))[1] # number of individuals
  snp.size = ceiling(sample.size/4) # no. bytes storing the data for 1 SNP for all individuals; each byte stores 4 people
  n.snps = round((bed.file.size-3)/snp.size) # the total .bed file size is 3 plus size of ids x snps combo hence removing 3
  snp.names = as.factor(read.table(paste(root,'.bim', sep = ''), header = FALSE, colClasses = c('NULL', 'character', 'NULL', 'NULL', 'NULL', 'NULL'))[,1]) # the bim file has the snps names in the 2nd column
  snp.names = as.vector(snp.names)
  bin.connection = file(bed.file, 'rb') # opens a connection with .bed file
  test.bytes = readBin(bin.connection, what = "raw", n = 3) # bytes 1 and 2 store infor about file type, byte 3 about array type
  if(!identical(as.character(test.bytes), c('6c', '1b', '01'))) {
    stop('BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file')
  }
  genotype = matrix(ncol = 1, nrow = sample.size) # to store the genos (obviously)
  cat(paste(root, 'results', sep = ' '), sep = '\n', file = results) # inizialises the results file, all the results are appeneded
  for(i in 1:n.snps) {
    r.bin.snp = readBin(bin.connection, what = 'raw', n = snp.size)
    bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2, byrow = TRUE)[1:sample.size,]
    genotype[,1] = bin.snp[,1] + bin.snp[,2] - 10 * ((bin.snp[,1] == 1) & (bin.snp[,2] == 0))
    genotype[genotype == -9] = NA 
    g.res = try(mPhen(genotype, phenoData = phenoData, phenotypes = phenotypes, covariates = covariates, resids = resids, strats = strats, maf_thresh = maf_thresh, corThresh = corThresh, inverseRegress = inverseRegress, JointModel = JointModel, multiGen = multiGen, fillMissingPhens = fillMissingPhens, scoreTest = scoreTest, imputed = FALSE), silent = T)
    if(class(g.res) == 'try.error') next
    g.res1 = g.res[[1]]$Results
    g.res2 = g.res[[1]]$nobs
    r.names = dimnames(g.res1)[[1]]
    resmat = matrix(c(r.names, as.vector(g.res1)), nrow = 3, byrow =T)
    resmat = cbind(resmat, c('nobs','NA', g.res2 ))
    snp = snp.names[i]
    resmat = rbind(snp, resmat)
    cat(cbind(c('Marker', 'Phenotype', 'beta', 'p-value'), resmat), sep = c(', ', ', ', ', ', '\n'), append = T, file = results) # writes results
    cat(',,,', sep = '\n', append = T, file = results) # writes SNP name
  }
  close(bin.connection)
}

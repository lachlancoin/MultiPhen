mPhen <-
function(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL,maf_thresh = 0.001, corThresh = 0.0, rescale = 10, inverseRegress = FALSE, multiPhenTest = FALSE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE, exactTest = FALSE, exactMethod = "wald", expandData = FALSE){
  exactTest = FALSE
  exactMethod = "wald"
  expandData = FALSE
  multiPhen = multiPhenTest
  lmt1 = cbind(rep('pheno', length(phenotypes)), phenotypes)
  cvt1 = cbind(rep('covar', length(covariates)), covariates)
  limit = rbind(lmt1, cvt1)
  samples = as.matrix(dimnames(phenoData)[[1]]) 
  dimnames(samples) = list(NULL,"id")
  max_indiv = 100000
  ids = colnames(genoData)
  #ids = names[1,1:3] #c("Log","GType","Freq")
  baselevel=c(0,0,0,0)  ##this is used for maf calc.  For count all base (ie. normal) is 2
  baselevel[ids=="countAll"] = 2
  singleOutput=TRUE
  zipOutput=FALSE
  plate_id = "Plate"
  variance_id ="variance"
  varthresh =0.25
  step = 1
  #dimnames(samples)[[2]] = names[3,1:dim(samples)[2]]
  dataUnc = FALSE
  if(dataUnc){
    if(expandData) weights = as.vector(t(apply(genoData,1,.mod1, rescale)))
    genoData = as.matrix(apply(genoData,1,.mod))
    if(inverseRegress & !expandData) genoData = apply(genoData,c(1,2),round)
  }
  #dimnames(genoData) = list(1:length(samples[,1]),names[1,1:dim(genoData)[2]])
  # dimnames(genoData)[[2]] =names[1,1:dim(genoData)[2]]
  geno_header = dimnames(genoData)[[2]]
  geno = genoData
  genotmp = genoData[,1,drop=FALSE]
  if(max_indiv<dim(geno)[1]) geno = geno[1:max_indiv,,drop=FALSE]
  if(expandData & !dataUnc) weights = .convAvgToUncMat(geno,rescale)
  #if(dim(samples)[1]>max_indiv) samples = samples[1:max_indiv,,drop=FALSE]
  phenoData = .fixNonNumeric(phenoData)
  phenNames = dimnames(phenoData)[[2]]
  todo_ = limit
  todo = todo_[grep("^pheno",todo_[,1]),2]
  excl = todo_[grep("^exclpheno",todo_[,1]),2]
  covar = todo_[grep("^covar",todo_[,1]),2]
  stratify=todo_[grep("^strat",todo_[,1]),2]
  resid=todo_[grep("^resid",todo_[,1]),2]
  genocov=todo_[grep("^genocov",todo_[,1]),2]
  index = NULL
  index_cov = NULL
  index_resid = NULL
  index_strat = NULL
  exclindex = NULL
  for ( ik in todo) index =c(index,which(phenNames==ik))
  for ( ik in excl) exclindex =c(exclindex,grep(ik,phenNames))
  for ( ik in covar) index_cov = c(index_cov,which(phenNames==ik))
  for ( ik in resid) index_resid = c(index_resid,which(phenNames==ik))
  for ( ik in stratify) index_strat = c(index_strat,which(phenNames==ik))
  index = unique(index)
  exclindex = unique(exclindex)
  if(length(exclindex)>0)index = index[-match(exclindex,index)]
  index_cov = unique(index_cov)
  index_resid = unique(index_resid)
  index_strat = unique(index_strat)
  if(dim(todo_)[2]>2)phenoTrans = todo_[grep("^pheno",todo_[,1]),3] else phenoTrans = NULL
  if(length(genocov)>0) genoTrans = apply(todo_[grep("^genocov",todo_[,1]),,drop=FALSE],1,.findInd,geno_header,3) else genoTrans = NULL
  if(length(genocov)==0) genotmp=NULL
  if(length(index)==1) multiPhen = FALSE
  if(length(index_strat)>0 | length(index_cov) > 0) exactTest = FALSE
  if(exactTest){
    multiPhen = FALSE
    multiGen = FALSE
    inverseRegress = FALSE
  }
  funct = .getPvLrr
  pvFunct = .getPvLrrMult
  pvFunctMult = .getPvLrrAllGen
  if(inverseRegress) funct = .getPvRev
  if(multiPhen) pvFunct = .getPvRevPleio
  if(scoreTest & multiPhen) pvFunct = .ordTest1
  if(multiGen) pvFunctMult = .getPvLrrMultiGen
  if(exactTest) funct = .getPvTable
  toStrat = phenoData[,index_strat,drop=FALSE]
  stratNames = dimnames(phenoData)[[2]][index_strat]
  pheno_cov = phenoData[,index_cov,drop=FALSE]
  pheno_resid =phenoData[,index_resid,drop=FALSE]
  phenoData = phenoData[,index,drop=FALSE]
  dimnames(phenoData)[[2]] = phenNames[index]
  if(multiPhen){
    if(fillMissingPhens) phenoData = .fillMissing(phenoData)
    corel = (cor(phenoData,use="complete.obs"))^2
    index1 = corel[1,]>corThresh
    phenoData = phenoData[,index1,drop=FALSE]
    index = index[index1]
    dimnames(phenoData)[[2]] = phenNames[index]
  }
  var_id = which(dimnames(samples)[[2]]==variance_id)
  plate_ind = which(dimnames(samples)[[2]]==plate_id)
  plate = matrix(0, ncol=1, nrow = dim(samples)[1])
  if(length(plate_ind)>0){
    plate = .makeFactor(as.matrix(samples[,plate_ind]),c("plate"))
  }
  variance = samples[,var_id]
  if(length(var_id)>0) varna =  variance>varthresh else varna = rep(FALSE,dim(samples)[1]) 
  indiv = samples[,1]
  mat = match(samples[,1],dimnames(phenoData)[[1]])
  pheno1 = phenoData[mat,,drop=FALSE]
  dimnames(pheno1)[[2]] = dimnames(phenoData)[[2]]
  pheno_cov1 = pheno_cov[mat,,drop=FALSE]
  pheno_resid1 = pheno_resid[mat,,drop=FALSE]
  stratMatrix1 = .makeFactor(toStrat[mat,,drop=FALSE],stratNames)
  if(!is.null(phenoTrans)){
    #print(phenoTrans)
    qqfam =  grep("^quantile",phenoTrans)
    if(length(qqfam)>0) pheno1[,qqfam] = apply(pheno1[,qqfam,drop=FALSE],2,.makeQuantile)
    toptail = grep("^toptail",phenoTrans)
    if(length(toptail)>0) pheno1[,toptail] = apply(pheno1[,toptail,drop=FALSE],2,.makeTopTail,strsplit(phenoTrans[toptail[1]],"_")[[1]][2])
    thresh = grep("^thresh",phenoTrans)
    if(length(thresh)>0) pheno1[,thresh] = apply(pheno1[,thresh,drop=FALSE],2,.makeThresh,strsplit(phenoTrans[thresh[1]],"_")[[1]][1:2])
    dimnames(pheno1)[[2]] = dimnames(phenoData)[[2]]
  }
  #if(!is.null(varianceIsCov)){
   # if(varianceIsCov & length(var_id)>0) pheno_cov1 = cbind(pheno_cov1,as.numeric(variance))
   # if(varianceIsCov & length(plate_ind)>0) pheno_cov1 = cbind(pheno_cov1,plate)
   # if(!varianceIsCov & length(var_id)>0) pheno_resid1 = cbind(pheno_resid1,as.numeric(variance))
   # if(!varianceIsCov & length(plate_ind)>0) pheno_resid1 = cbind(pheno_resid1,plate)
  #}
  if(!is.null(genotmp)){
    #if(genoCovIsRes) pheno_resid1 = cbind(pheno_resid1,.fixNonNumeric(genotmp)) else 
    pheno_cov1 = cbind(pheno_cov1,.fixNonNumeric(genotmp))
  }
  if(dim(pheno_cov1)[2]==0) phenoCovna = rep(0,dim(samples)[1]) else phenoCovna = is.na(apply(pheno_cov1,1,sum))
  if(dim(pheno_resid1)[2]==0) phenoResna = rep(0,dim(samples)[1]) else phenoResna = is.na(apply(pheno_resid1,1,sum))
  families = rep(0, dim(pheno1)[2])
  inclM = as.matrix(apply(pheno1,2,.getIncl, phenoCovna | phenoResna | varna))
  for(i in 1:(dim(pheno1)[2])) families[i] = .getFamily(pheno1[,i],inclM[,i])
  binom = families=="binomial"
  caseInd = apply(pheno1,2,.iscase)
  for(i in 1:(dim(pheno1)[2])) if(binom[i]) pheno1[,i] = pheno1[,i]-min(pheno1[,i],na.rm=TRUE)
  if(dim(pheno_cov1)[2]>0)for(i in 1:(dim(pheno_cov1)[2])) pheno_cov1[,i] = .centralise(pheno_cov1[,i])
  inds = NULL
  for(i in 1:length(ids)) inds = c(inds,grep(ids[i],geno_header))
  ids = geno_header[inds]
  geno_header = geno_header[inds]
  alleleCols = grep("Allele[12]",geno_header)
  gCols = grep("GType",geno_header)
  if(length(alleleCols)>0){
    geno_header = c("geno",geno_header[!alleleCols])
  }
  if(!is.null(geno)) {
    geno = geno[,inds,drop=FALSE]
    if(length(alleleCols)>0) geno = .mergeAlleleCols(geno,alleleCols)
    if(length(gCols)>0) geno = .fixGenoCols(geno, gCols)
    maxg = apply(geno,2,max,na.rm=TRUE)
    ming = apply(geno,2,min,na.rm=TRUE)
  }
  nonint = maxg%%1 > 0
  nonint1 = ming%%1 > 0
  maxg[nonint] = floor(maxg[nonint])+1
  ming[nonint1] = floor(ming[nonint1])
  spl = seq(min(ming)+1,(max(maxg)),step) - 0.5
  typenames = c("pheno","incl","caseInd","offsets")
  phenN = phenNames[index]
  if(!is.null(phenoTrans)) phenN = apply(rbind(phenN,phenoTrans),2,paste,collapse=".")
  arraynames = list(phenN, typenames,indiv )
  phendata = .getArray(arraynames)
  phendata[,1,] = t(pheno1)
  phendata[,2,] = t(inclM)
  phendata[,3,] = t(caseInd)
  genoweights = .getArray(list(geno_header,c("geno","weights","include"),samples[,1]))
  if(!is.null(geno)) genoweights[,1,] = t(geno)
  genoweights[,2,] =  t( matrix(1,nrow = dim(samples)[1],ncol = length(geno_header)))
  datanme = list("pheno"=phendata,"pheno_cov" = pheno_cov1, "pheno_res" = pheno_resid1, "geno" = genoweights, "strat"=stratMatrix1, "alleleCols" = alleleCols,"gCols" = gCols, "ind" = inds)
  datanme2 = datanme
  if(length(geno_header)==1) multiGen = FALSE
  if(expandData){
    len = length(indiv)
    datanme2 = .expandDat(datanme,length(spl)+1)
    if(!is.null(weights))datanme2$geno[,2,] = t(weights[,inds])
    for(i in 1:(length(spl)+1)) datanme2$geno[,1,((i-1)*len +1):(i*len)] = rep((i-1),length(geno_header))
    datanme2$geno[,3,] = matrix(TRUE,ncol = length(geno_header), nrow = length(indiv)*(length(spl)+1))
  }
  if(length(index)>0) for(i in 1:(length(index))) datanme$pheno[i,4,] = .calcOffset(datanme$pheno[i,,],families[i], datanme$pheno_res)
  if(expandData){
    if(length(index)>0) for(i in 1:(length(index))) datanme2$pheno[i,4,] = .calcOffset(datanme2$pheno[i,,],families[i], datanme2$pheno_res)
  }
  if(!expandData) datanme2 = datanme
  fams = levels(as.factor(families))
  dimcounts = c(dim(datanme$strat)[2],length(index),length(geno_header),2)
  if(length(index_strat)>0) dimcounts[1] = dimcounts[1]+1
  if(multiPhen) dimcounts[2] = dimcounts[2]+1
  if(multiGen) dimcounts[3] = dimcounts[3]+1
  stratificNames = dimnames(datanme$strat)[[2]]
  if(length(stratificNames)==1) stratificNames = c("all")
  phenN = phenNames[index]
  if(multiPhen) phenN = c(phenN,"multiPhen")
  res = array(NA,dim=dimcounts,dimnames=list(stratificNames, phenN, geno_header,c("beta","pvalue")))
  geno1 = datanme$geno
  geno1[,3,] = !apply(as.matrix(geno1[,1,]),c(1,2),is.na)
  geno2 = geno1
  if(expandData) geno2 = datanme2$geno
  maf = rep(NA,dim(geno1)[1])
  for(k in 1:length(geno_header)) maf[k] = .getMaf(as.numeric(geno1[k,1,]),baselevel[k],2)
  for(j in 1:(dim(datanme$strat)[2])){
    inclu = datanme2$strat[,j]
    res[j,,,] =  pvFunctMult(geno2,datanme2$pheno,families,inclu,datanme2$pheno_cov,rescale,pvFunct,funct)
  }
  if(length(index_strat)>0) 
    res[dim(datanme$strat)[2]+1,,,] =apply(apply(res[1:(dim(datanme$strat)[2]-1),,,,drop=FALSE],c(2,3),.metaresInvVarianceFixed),c(3,1),t)
  dimcounts = c(dim(datanme$strat)[2],length(families),length(geno_header),(length(spl)+1))
  countsCase = array(NA,dim=dimcounts)
  countsControl = array(NA,dim=dimcounts)
  hwe_control = array(NA,dim=dimcounts[1:3])
  hwe_case = array(NA,dim=dimcounts[1:3])
  maf_control = array(NA,dim=dimcounts[1:3])
  maf_case = array(NA,dim=dimcounts[1:3])
  for(j in 1:(dimcounts[1])){
    stra = datanme$strat[,j]
    for(k in 1:length(families)){
       # print(paste(j,k))
      inclu = datanme$pheno[k,2,] & stra
      incluCont = inclu &  !datanme$pheno[k,3,]
      incluCase = inclu &  datanme$pheno[k,3,]
      if(length(which(incluCont))>0)  countsControl[j,k,,] = t(apply(geno1[,1,incluCont,drop=FALSE],1,.getHist,spl))
      if(length(which(incluCont))>0)  countsCase[j,k,,] = t(apply(geno1[,1,incluCase,drop=FALSE],1,.getHist,spl))
      hwe_control[j,k,] =  apply(countsControl[j,k,,,drop=FALSE],3,.hwepall)
      hwe_case[j,k,] =  apply(countsCase[j,k,,,drop=FALSE],3,.hwepall)
      maf_control[j,k,] =  apply(countsControl[j,k,,,drop=FALSE],3,.mafCount)
      maf_case[j,k,] =  apply(countsCase[j,k,,,drop=FALSE],3,.mafCount)
    }
  }
  res1 = res[1,,,]
  res1
  #print(res[dim(res)[1],,,])
}


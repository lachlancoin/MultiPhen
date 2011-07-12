.allFreq <-
function(y){
  vecm = c(2,1,0) 
  (vecm[1:length(y)] %*% y) / (2*(sum(y)))
}

.calcOffset <-
function(phenv,family,pheno_cov){
  offsets1 = rep(NA,dim(pheno_cov)[2])
  incl = phenv[2,]==1
  dimres = dim(pheno_cov)[2]
  if(dimres>0) 
    offsets1[incl] =  glm(phenv[1,incl]~pheno_cov[incl,],family=family)$lin else offsets1[incl] = 0
  offsets1
}

.centralise <-
function(vec){x = (vec-mean(vec,na.rm=TRUE))/sd(vec,na.rm=TRUE)}

.cnt <-
function(str) {
  if(length(grep("N",str))>0) 
    NA 
  else 
    nchar(gsub("A","",str))
}

.cntVec <-
function(str) apply(as.matrix(str,nrow = length(str)),1,.cnt)

.convAvgToUnc <-
function(vec,n,sc=1){
  res = rep(0,n)
  fl = floor(vec)
  if(fl==vec) 
    res[vec+1] = sc
  if(fl!=vec){
    ceil = fl+1
    val = round(sc*(ceil-vec))
    res[fl+1] = val
    res[ceil+1] = sc-val
  }
  res
}

.convAvgToUncM <-
function(vec,sc,maximum){
#print(vec)
 matr = apply(as.matrix(vec),1,.convAvgToUnc,maximum)
 as.vector(t(apply(as.matrix(vec),1,.convAvgToUnc,maximum,sc)))
}

.convAvgToUncMat <-
function(mat,sc){
 maximum = max(mat,na.rm=TRUE)+1
 if(maximum%%1>0) maximum = floor(maximum)+1
 apply(mat,2,.convAvgToUncM,sc,maximum)
}

.convertFactorToMatrix <-
function(f,nme){
   res= apply(as.matrix(levels(as.factor(f))),1,.getOnes,f)
   dimnames(res) = list(names(f),paste(nme,levels(as.factor(f)),sep="_"))
   res
}

.convertMatrixToFactor <-
function(mat){
 nme = as.matrix(data.frame(strsplit(dimnames(mat)[[2]],"_")))
 res = rep("",dim(mat)[1])
 for(i in 1:length(res)) 
   res[i] = paste(nme[2,mat[i,]==1],collase="")
 res
}

.estSigma <-
function(x){
  d = ncol(x)
  v = matrix(ncol=d,nrow=d)
  for( i in 1:d ){
    ptrI = !is.na(x[,i])
    v[i,i] = var(x[ptrI,i])
    if( i < d ){
      for( j in (i+1):d ){
        ptrJ = !is.na(x[,j])
        ptr = as.logical(ptrI*ptrJ)
        v[i,j] = cov( x[ptr,i], x[ptr,j] )
        v[j,i] = v[i,j]
      }
    }
  }
  return(v)
}

.expand <-
function(matr, repl){
 res = matr
 for(i in 2:repl) res = rbind(res,matr)
 res
}

.expand2 <-
function(matr, repl){
  res = matr
  for(i in 2:repl)
    res = abind(res,matr,along = 3)
  res
}

.expandCoeff <-
function(coeff,ind){
  diff =  max(ind) - dim(coeff)[1]
  if(diff>0) coeff = rbind(coeff,matrix(NA, nrow = diff,ncol = dim(coeff)[2]))
  coeff
}

.expandDat <-
function(data,repl){
  inds = sapply(data,.nullDim)
  data[inds==2] = lapply(data[inds==2],.expand,repl)
  data[inds==3] = lapply(data[inds==3],.expand2,repl)
  data
}

.expDist <-
function(y,p)   c(p^2, 2*p*(1-p), (1-p)^2)[1:length(y)]*sum(y)

.extractIds <-
function(vec,nme) {paste(nme[vec],collapse=",")[[1]]}

.extractIdsMat <-
function(mat,spl,caseInd) apply(mat,1,.extractIdsVec,spl,caseInd,dimnames(mat)[[3]])

.extractIdsVec <-
function(vec,spl,caseInd,nme) apply(.getin2(as.numeric(vec),spl,caseInd),2,.extractIds,nme)

.extractSig <-
function(res,ind = 2, totweight = 1){
  vec = summary(res)$coeff[ind,]
  if(length(vec)==3) vec = c(vec, 2*pt(abs(vec[3]),1,lower.tail=FALSE)) 
  vec[c(1,4)]
}

.extractSig1 <-
function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial"){
  ll1 = logLik(res)
  beta = summary(res)$coeff[ind,1]
  nullfam = is.null(family)
  if(hasCov & !nullfam) ll2 = logLik(update(res, ~covar,family=family)) else if(!nullfam) ll2 = logLik(update(res, ~1,family=family)) else if(hasCov) ll2 = logLik(update(res, ~covar)) else ll2 = logLik(update(res, ~1))
  df = attr(ll1,"df")[1]-attr(ll2,"df")[1]
  # print(paste("df",df))
  pv = pchisq((2*(ll1 - ll2))/(totweight),df,lower.tail=FALSE)
  c(beta,pv) 
}

.extractSig1P <-
function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial"){
  ll1 = logLik(res)
  coeff = summary(res)$coeff
  beta = coeff[ind,1]
  ve = abs(as.numeric(coeff[ind,3]))
  pvs = 2*pt(ve,1,lower.tail=FALSE)
  nullfam = is.null(family)
  if(hasCov & !nullfam) 
    ll2 = logLik(update(res, ~covar,family=family)) 
  else 
    if(!nullfam) 
    ll2 = logLik(update(res, ~1,family=family)) 
  else if(hasCov) 
    ll2 = logLik(update(res, ~covar)) 
  else 
    ll2 = logLik(update(res, ~1))
  pv = pchisq((2*(ll1 - ll2))/(totweight),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE)
  cbind(c(beta,NA),c(pvs,pv))
}

.extractSigExact <-
function(t){
  exactMethod = "wald"
  # print(t) note, not sure this is getting things in right order for good effect size estimate
  or1 = fisher.test(t)
  or = oddsratio(t,method = exactMethod)
  #  ind = grep(exactMethod,dimnames(or$p.value))[1]
  #  print(or)
  c(log(or$measure[2,1]),or1$p.value)
}

.extractSigP <-
function(res,ind,totweight = 1,hasCov=FALSE,family="binomial"){
  coeff = .expandCoeff(summary(res)$coeff,ind)
  vec  = coeff[ind,,drop=FALSE]
  ll1 = logLik(res)
  if(dim(vec)[2]==3) vec = cbind(vec, 2*pt(abs(vec[ind,3]),1,lower.tail=FALSE))
  if(hasCov)ll2 = logLik(update(res, ~covar,family=family)) else ll2 = logLik(update(res, ~1,family=family))
  pv = pchisq((2*(ll1 - ll2))/(totweight),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE)
  cbind(c(vec[,1],NA),c(vec[,4],pv))
}

.fillMissing <-
function( x){
  v = .estSigma(x)
  m = apply(x,2,mean,na.rm=TRUE)
  lambda = solve(v)
  xx = x
  ptr.ind = which( apply(is.na(x),1,sum)>0 )
  for( ii in 1:length(ptr.ind) ){
    i = ptr.ind[ii]
    ptr.miss = which(is.na(x[i,]))
    fill =  m[ptr.miss] - solve(lambda[ptr.miss,ptr.miss]) %*% lambda[ptr.miss,-ptr.miss] %*% as.matrix(x[i,-ptr.miss] - m[-ptr.miss])
    xx[i,ptr.miss] = fill
  }
  return(xx)
}

.findInd <-
function(vec,header,ind){x = grep(vec[ind],header)[1]; x}

.findmax <-
function(vec) max(0,max(as.numeric(vec),na.rm=TRUE))

.fisherPvalue <-
function(p){ pchisq(-2*sum(log(p)),df = 2*length(p),lower.tail=FALSE)}

.fixGenoCols <-
function(geno, gCols){
 geno1 = matrix(0,nrow=dim(geno)[1],ncol=dim(geno)[2])
 geno1[,gCols] = apply(geno[,gCols,drop=FALSE],2,.cntVec)
 geno1[,-gCols] = apply(geno[,-gCols,drop=FALSE],2,as.numeric)
 geno1
}

.fixNonNumeric <-
function(pheno){
 nonnumvec = .nonNumeric1(pheno)
 pheno1 = matrix(0,nrow=dim(pheno)[1],ncol=dim(pheno)[2])
 pheno1[,!nonnumvec] = apply(as.matrix(pheno[,!nonnumvec]),2,as.numeric)
 pheno1[,nonnumvec] = apply(as.matrix(pheno[,nonnumvec]),2,.makeNumeric)
 dimnames(pheno1) = dimnames(pheno)
 pheno = pheno1
}

.getArray <-
function(arraynames) array(NA,dim= .getDim(arraynames), dimnames = arraynames)

.getDim <-
function(arraynames) c(length(arraynames[[1]]), length(arraynames[[2]]), length(arraynames[[3]]))

.getFamily <-
function(phenv,incl){
 if(length(levels(as.factor(phenv[incl])))>2) "gaussian" else "binomial" 
}

.getGlm <-
function(phenvec,vars,covar, fam1,weights,offset){
  glm(phenvec~vars+covar,family = fam1,weights = weights,offset = offset)
}

.getGlmNoCovar <-
function(phenvec,vars, fam1,weights,offset){
  glm(phenvec~vars,weights = weights,offset = offset, family = fam1)
}

.getGlmRevBinom <-
function(vars,phenvec,covar,weights, family) glm(vars ~ phenvec + covar,family=family,weights=weights)

.getGlmRevBinomNoCov <-
function( vars, phenvec, weights, family){
  glm(vars ~ phenvec,family=family,weights=weights)
}

.getGlmRevOrd <-
function( vars, phenvec, covar,weights) {
 polr(vars ~ phenvec + covar,Hess=TRUE,method="logistic",weights = weights)
}

.getGlmRevOrdNoCov <-
function( vars, phenvec, weights){
  polr(vars ~ phenvec,Hess=TRUE,method="logistic",weights = weights)
}

.getHist <-
function(vec, spl){ 
  le = vec < spl[1]
  gt = vec > spl[length(spl)]
  #print(spl);print(max(vec));print(min(vec))
  c(length(which(le)),hist(as.numeric(vec[!le & !gt]),br=spl,plot=FALSE)$count,length(which(gt)))
}

.getIMatrix <-
function(d,phi,cats,num1,x){
  II = matrix(ncol=d,nrow=d)
  for(j in 1:d){
    for(k in j:d){
      num = 0
      for( i in 1:cats ) num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
      II[j,k] = sum(num1[[j]]*num1[[k]])
      II[j,k] = II[j,k] - sum(num)
      II[k,j] = II[j,k]
    }
  }
  II
}

.getin <-
function(stend,vec) vec>=stend[1] & vec<stend[2]   #NEW

.getin1 <-
function(genvec,spl) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,.getin, genvec) #NEW

.getin2 <-
function(genvec,spl,cc) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,.getin, genvec) & cc

.getIncl <-
function(phenv, na_cov){incl = !((is.na(phenv) | na_cov)); incl}

.getMaf <-
function(vec,base,noallele){
  .minOneMinus(sum(abs(vec[!is.na(vec)] -base))/(length(which(!is.na(vec)))*noallele))
}

.getOnes <-
function(vec,f) f==vec[1]

.getPvLrr <-
function(phend, gvar, family,inclu,covar,totweight){
  phenvec = phend[1,]
  vars = gvar[1,]
  weights = gvar[2,]
  offset = phend[4,]
  include = phend[2,] & inclu & weights > 0 & gvar[3,]>0 
  emptyCov = dim(covar)[2]==0
  todo = TRUE
  if(length(which(include))==0 | var(vars[include],na.rm=TRUE)==0) 
    todo=FALSE
  if(todo & emptyCov) 
    glm.res = .getGlmNoCovar(phenvec[include] ,as.numeric(vars[include]),family,weights[include],offset[include]) 
  else 
    if(todo) 
    glm.res = .getGlm(phenvec[include] , as.numeric(vars[include]), covar[include,],family,weights[include],offset[include])
  if(!todo) 
    res = rep(NA,2) 
  else 
    res = .extractSig(glm.res,ind=2,totweight = totweight)
  res
}

.getPvLrrAllGen <-
function(vars,phen,families,include,covar,totweight,functiAll,functi){
  x = t(apply(vars,1,functiAll,phen,families,include,covar,totweight,functi))
  x
}

.getPvLrrMult <-
function(vars,phen,families,include,covar,totweight, functi){
  res = matrix(NA, nrow = length(families),ncol = 2)
  fams = levels(as.factor(families))
  for(k in 1:length(fams))   
    if(length(which(families==fams[k]))>0) 
    res[families==fams[k],] = t(apply(phen[families==fams[k],,,drop=FALSE],1,functi,vars, fams[k],include,covar,totweight))
  res
}

.getPvLrrMultiGen <-
function(gvar,phend, family,inclu,covar,totweight,functiAll,functi){
  x = t(apply(phend,1,.getPvLrrMultiGen1,gvar,family,inclu,covar,totweight,functiAll,functi))
  x
}

.getPvLrrMultiGen1 <-
function(phend, gvar,family,inclu,covar,totweight,functiAll,functi){
  phenvec = phend[1,]
  vars = gvar[,1,]
  weights1 = gvar[1,2,]  
  nainds =  apply(as.matrix(gvar[,3,]),2,min)
  offset = phend[4,]
  include = phend[2,] & inclu & nainds >0 & weights1 > 0
  # & nainds # & weights > 0
  emptyCov = dim(covar)[2]==0
  ind = 1:(dim(vars)[1]) + 1
  todo = TRUE
  if(length(which(include))==0 ) 
    todo=FALSE
  if(todo & emptyCov) 
    glm.res = .getGlmNoCovar(phenvec[include],t(vars[,include]),family,weights[include],offset[include]) 
  else 
    if(todo) 
    glm.res = .getGlm(phenvec[include] , t(vars[,include]), covar[include,],family,weights[include],offset[include])
  if(!todo) rep(NA,2) 
  else 
    .extractSigP(glm.res,ind=ind,totweight=totweight,family=family,hasCov=!emptyCov)
}

.getPvRev <-
function(phend, gvar, family, inclu,covar,totweight, maxCatForOrd = 5){ 
  vars = gvar[1,]
  weights = as.numeric(gvar[2,])
  include = phend[2,]>0 & inclu & weights > 0 & gvar[3,]>0
  var1 = as.factor(vars[include])
  phenvec = phend[1,include] - phend[4,include]
  emptyCov = dim(covar)[2]==0
  binom =  length(levels(var1))<=2
  cont = length(levels(var1))>maxCatForOrd
  family = "binomial"
  ord = !binom & !cont
  if(cont) 
    var1 = vars[include]
  if(cont) 
    family = "gaussian"
  if(ord)
    family=NULL
  todo = TRUE
  # print(paste("running pv rev",binom,cont,family,ord))
  if(length(which(include))==0 | var(var1,na.rm=TRUE)==0) 
    todo=FALSE
  if(todo & !ord & !emptyCov)
    glm.res = .getGlmRevBinom(var1,phenvec,covar[include,],weights[include], family) 
  else 
    if(todo & !emptyCov) 
      glm.res = .getGlmRevOrd(var1,phenvec,covar[include,],weights[include]) 
    else 
      if(todo & !ord) 
        glm.res = .getGlmRevBinomNoCov(var1,phenvec,weights[include],family) 
      else 
        if(todo) 
          glm.res = .getGlmRevOrdNoCov(var1,phenvec,weights[include])
  if(!todo) rep(NA,2) 
  else 
    .extractSig1(glm.res,ind=1, totweight = totweight,hasCov=!emptyCov,family=family)
} #####FINISH

.getPvRevPleio <-
function(gvar,phen,families,inclu,covar,totweight,functi,maxCatForOrd = 5){
  phenvec = phen[,1,]-phen[,4,]
  vars = gvar[1,]
  weights = gvar[2,]
  nainds =  gvar[3,]
  phenincl = apply(as.matrix(phen[,2,]),2,min)
  include =  phenincl>0 & inclu & weights>0 & nainds >0
  var1 = as.factor(vars[include])
  binom = length(levels(var1))<=2
  cont = length(levels(var1))>maxCatForOrd
  family = "binomial"
  ord = !binom & !cont
  if(cont) family = "gaussian"
  if(ord)family=NULL
  if(cont) var1 = vars[include]
  emptyCov = dim(covar)[2]==0
  ind = 1:(dim(phenvec)[1])
  if(!ord) ind = ind+1
  todo = TRUE
  if(length(which(include))==0 | var(var1,na.rm=TRUE)==0) 
    todo=FALSE
  if(todo & !ord & emptyCov) 
    glm.res = .getGlmRevBinomNoCov(var1,t(phenvec[,include]),weights[include],family) 
  else 
    if(todo & emptyCov) 
    glm.res = .getGlmRevOrdNoCov(var1,t(phenvec[,include]),weights[include]) 
  else 
    if(todo & !ord) 
    glm.res = .getGlmRevBinom(var1,t(phenvec[,include]),covar[include,],weights[include],family) 
  else 
    if(todo)  
    glm.res = .getGlmRevOrd(var1,t(phenvec[,include]),covar[include,],weights[include])
  if(!todo) 
    rep(NA,2) 
  else 
    .extractSig1P(glm.res,ind=ind,totweight=totweight,family=family)
}

.getPvTable <-
function(phend, gvar, family,inclu,covar,totweight){
  phenvec = phend[3,]
  vars = as.factor(gvar[1,])
  weights = gvar[2,]
  offset = phend[4,]
  include = phend[2,] & inclu & weights > 0 & gvar[3,]>0
  emptyCov = dim(covar)[2]==0
  todo = TRUE
  if(length(which(include))==0 | var(vars[include],na.rm=TRUE)==0 | family!="binomial") 
    todo=FALSE
  if(!todo) 
    rep(NA,2) 
  else 
    .extractSigExact(.tabulateWithWeights(phenvec[include],vars[include],weights[include],totweight))
}

.getUNum1Matrix <-
function(d,y,x,phi,num1,cats){
  U = matrix(ncol=1,nrow=d)
  num1 = list(length=d)
  for( j in 1:d ){
    U[j] = sum(y*x[,j])    
    num1[[j]]=0 
    for( i in 1:cats ) 
      num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i]) 
    U[j,1] = U[j,1] - sum(num1[[j]])
  }
  list(U,num1)
}

.hwep <-
function(y) if(length(which(y>0))<=1) 1.0 else 1-pchisq(.hwestat(y,.expDist(y,.allFreq(y))),1)

.hwepall <-
function(y) min(.hwep(y[1:min(3,length(y))]), .hwep(y[3:min(5,length(y))]),na.rm=TRUE)

.hwestat <-
function(obs,exp) sum((obs-exp)^2/exp)

.iscase <-
function(vec) !is.na(vec) & vec>mean(vec,na.rm=TRUE)

.joinRes <-
function(mat,ncol,form){
 if(is.null(dim(mat))) 
   mat = matrix(mat,ncol=ncol)
 res1 =  apply(mat,2,.joinResVec,form)
 res1
}

.joinResVec <-
function(vec,form) paste(sprintf(form,vec),collapse=";")

.mafCount <-
function(y) max(.minOneMinus(.allFreq(y[1:min(3,length(y))])), .minOneMinus(.allFreq(y[3:min(5,length(y))])),na.rm=TRUE)

.makeFactor <-
function(mat1,nme){
  append = dim(mat1)[2]>0
  fact = factor(apply(mat1,1,paste,collapse=".")) 
  lev = levels(fact) 
  matrix = .convertFactorToMatrix(fact,paste(nme,sep="."))
  nmes = dimnames(matrix)[[2]]
  txt = dimnames(matrix)[[2]]
  index = setdiff(seq(length(txt)),grep("_NA",txt))
  res = matrix[,index,drop=FALSE]
  dimnames(res) = list(NULL,nmes[index])
  if(append)  
    res = cbind(res, rep(TRUE,dim(res)[1]))
  if(append)   
    dimnames(res) = list(NULL,c(nmes[index],"all")) 
  res 
}

.makeNumeric <-
function(vec) as.numeric(as.factor(vec))

.makeQuantile <-
function(vec)  qqnorm(vec,plot=FALSE)$x

.makeThresh <-
function(vec,thresh){
 vec1 = rep(NA,length(vec))
 vec1[vec<=thresh[1]] = 0
 vec1[vec>=thresh[2]] = 1
 vec1
}

.makeTopTail <-
function(vec,perc1){
 perc = as.numeric(perc1)/100
 vec1 = rep(NA,length(vec))
 quants = ecdf(sort(vec))(vec)
 vec1[quants<=perc] = 0
 vec1[quants>=1-perc] = 1
 vec1
}

.mergeAlleleCols <-
function(snp,cols){
  len = 1+ dim(snp)[2] - length(cols) 
  snp1 = matrix(0,nrow = dim(snp)[1],ncol =len) 
  snp1[,1] = .cntVec(apply(snp[,cols],1,paste,collapse=""))
  matr = snp[,-cols,drop=FALSE]
  nme1 = rep("",len)
  nme1[1] = "geno"
  if(len>1) 
    nme1[2:len] = dimnames(snp)[[2]][-cols]
  if(len>1) 
    for(i in 2:len) 
    snp1[,i] = matr[,i-1]
  dimnames(snp1) = list(dimnames(snp)[[1]],nme1)
  #print(snp1[1:5,])
  snp1
}

.metaresFisher <-
function(matr){
  inds = !is.na(matr[,2] )
  fisherp = .fisherPvalue(matr[inds,2])
  c(mean(matr[inds,1]),fisherp)
}

.metaresInvVarianceFixed <-
function(matr){
  inds = !is.na(matr[,2])
  metres = metagen(matr[inds,1,drop=FALSE],apply(matr[inds,,drop=FALSE],1,.seFromPBeta))
  pv = .pFromSEBeta(c(metres$TE.fixed,metres$seTE.fixed))
  ##c(metres$Q,pchisq(metres$Q,1,lower.tail=FALSE)) 
  c(metres$TE.fixed,pv)
}

.minOneMinus <-
function(x) min(x,1-x)

.mod <-
function(vec){
  #vec[vec==""] = 0; 
  #vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=TRUE); 
  (as.numeric(vec)%*% 0:2)#/1000
}

.mod1 <-
function(vec, rescale){
  #vec[vec==""] = 0
  #vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=TRUE);
  #vec[vec<200]  =0
  vec1 = floor(rescale * as.numeric(vec))#floor(rescale*(as.numeric(vec)/1000))
  ind = vec1==max(vec1)
  vec1[ind] = vec1[ind]+ rescale - sum(vec1)
  vec1
}

.mPhen <-
function(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL,maf_thresh = 0.001, corThresh = 0.0, inverseRegress = FALSE, multiPhenTest = FALSE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE, exactTest = FALSE, exactMethod = "wald", imputed = FALSE){
  rescale = 1
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
  if(imputed){
    if(expandData) weights = as.vector(t(apply(genoData,1,.mod1, rescale)))
    genoData = as.matrix(apply(genoData,1,.mod))
    if(inverseRegress & !expandData) genoData = apply(genoData,c(1,2),round)
    dimnames(genoData) = list(NULL,'rsID')
  } 
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
  #dimnames(genoData) = list(1:length(samples[,1]),names[1,1:dim(genoData)[2]])
  # dimnames(genoData)[[2]] =names[1,1:dim(genoData)[2]]
  geno_header = dimnames(genoData)[[2]]
  geno = genoData
  genotmp = genoData[,1,drop=FALSE]
  if(max_indiv<dim(geno)[1]) geno = geno[1:max_indiv,,drop=FALSE]
  if(expandData & !imputed) weights = .convAvgToUncMat(geno,rescale)
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
  #dimcounts = c(dim(datanme$strat)[2],length(families),length(geno_header),(length(spl)+1))
  #countsCase = array(NA,dim=dimcounts)
  #countsControl = array(NA,dim=dimcounts)
  #hwe_control = array(NA,dim=dimcounts[1:3])
  #hwe_case = array(NA,dim=dimcounts[1:3])
  #maf_control = array(NA,dim=dimcounts[1:3])
  #maf_case = array(NA,dim=dimcounts[1:3])
  #for(j in 1:(dimcounts[1])){
  #  stra = datanme$strat[,j]
  #  for(k in 1:length(families)){
  #     # print(paste(j,k))
  #    inclu = datanme$pheno[k,2,] & stra
  #    incluCont = inclu &  !datanme$pheno[k,3,]
  #    incluCase = inclu &  datanme$pheno[k,3,]
  #    if(length(which(incluCont))>0)  countsControl[j,k,,] = t(apply(geno1[,1,incluCont,drop=FALSE],1,.getHist,spl))
  #    if(length(which(incluCont))>0)  countsCase[j,k,,] = t(apply(geno1[,1,incluCase,drop=FALSE],1,.getHist,spl))
  #    hwe_control[j,k,] =  apply(countsControl[j,k,,,drop=FALSE],3,.hwepall)
  #    hwe_case[j,k,] =  apply(countsCase[j,k,,,drop=FALSE],3,.hwepall)
  #    maf_control[j,k,] =  apply(countsControl[j,k,,,drop=FALSE],3,.mafCount)
  #    maf_case[j,k,] =  apply(countsCase[j,k,,,drop=FALSE],3,.mafCount)
  #  }
  #}
  res1 = res[1,,,]
  res1
  #print(res[dim(res)[1],,,])
}

.nonNumeric <-
function(vec) sum(is.na(as.numeric(vec)))==length(vec)

.nonNumeric1 <-
function(mat) apply(mat,2,.nonNumeric)

.nullDim <-
function(matr){x = length(dim(matr)); x}

.numLev <-
function(vec) length(levels(as.factor(vec)))

.ordTes1 <-
function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  .ordTest(y,x)
}

.ordTest <-
function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  #print(is.numeric(x[,1]))
  cats = max(y) + 1
  x = x - apply(x,2,mean,na.rm=TRUE)
  d = ncol(x)
  gamma = as.vector(table(y)/length(y))
  phi = log(gamma)
  Unum1 = .getUNum1Matrix(d,y,x,phi,num1,cats)
  U = Unum1[[1]]
  num1 = Unum1[[2]]
  II = .getIMatrix(d,phi,cats,num1,x)
  zstats  = U/sqrt(diag(-II))
  pvs = 2*(1-pnorm(abs(U/sqrt(diag(-II)))))
  overallp = 	1 - pchisq(t(U) %*% solve(-II) %*% U, d)
  results = cbind(c(zstats, NA),c(pvs, overallp))
  results
}

.ordTest1 <-
function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  .ordTest2(y,x)
}

.ordTest2 <-
function(y,x){
  cats = max(y) + 1
  x = x - apply(x,2,mean,na.rm=TRUE)
  d = ncol(x)
  gamma = as.vector(table(y)/length(y))
  phi = log(gamma)
  U = matrix(ncol=1,nrow=d)
  num1 = list(length=d)
  for( j in 1:d ){
    U[j] = sum(y*x[,j])
    num1[[j]]=0
    for( i in 1:cats ){
      num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i])
    }
    U[j,1] = U[j,1] - sum(num1[[j]])
  }

  II = matrix(ncol=d,nrow=d)
  num2 = matrix(ncol=d,nrow=d)
  for( j in 1:d ){
    for( k in j:d ){
      num = 0
      for( i in 1:cats ){
        num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
      }
      II[j,k] = sum(num1[[j]]*num1[[k]])
      II[j,k] = II[j,k] - sum(num)
      II[k,j] = II[j,k]
    }
  }
  zstats  = U/sqrt(diag(-II))
  pvs = 2*(1-pnorm(abs(U/sqrt(diag(-II)))))
  overallp = 	1 - pchisq(t(U) %*% solve(-II) %*% U, d)
  results = cbind(c(zstats, NA),c(pvs, overallp))
  results
}

.pFromSEBeta <-
function(v){x =  pnorm(abs(v[1]),sd = v[2],lower.tail=FALSE)*2; x}

.proc <-
function(vec) (vec - 0.5)*vec*(1-vec)

.proc1 <-
function(vec) vec

.Random.seed <-
c(403L, 373L, 1550106100L, 838601621L, -1453538047L, -1218857010L, 
2091903414L, 1701334947L, -2071190038L, -1325172799L, -854370848L, 
-123597558L, 1009770699L, 969455859L, 740433579L, 67550253L, 
-1205714581L, 514627462L, 1353797431L, -848310372L, 1714975381L, 
1925093257L, -1106015188L, 247202605L, 1144584091L, -1020112686L, 
-716543581L, 133152719L, 1684125623L, -1076265275L, -1234986788L, 
-76629481L, -670330424L, -701724770L, 1347162808L, -1857942382L, 
586030930L, -165940417L, -625697718L, 937985191L, -293429048L, 
9831517L, 1826862792L, -478578024L, 2066725063L, 163297656L, 
-1863151304L, 306044288L, -666346730L, -2042926704L, -1255772666L, 
556348416L, 50058774L, -2093562466L, 1456067434L, 754442453L, 
-629114490L, 547598820L, -1941578998L, 150490190L, -1216853664L, 
1382779367L, 1673487256L, -2090953383L, 849035162L, 34772159L, 
-1932078408L, -335632525L, -1038233308L, -746967243L, 1135665566L, 
1490895439L, 339228423L, 1005578037L, 1054389384L, 1506275109L, 
2019926941L, -977290400L, 126799361L, 1559923915L, 729847863L, 
1696430802L, -1018143975L, 1554105059L, -944660222L, 1399493701L, 
-1745970034L, -495188578L, 1558069034L, 231344959L, 1976645281L, 
-940196548L, 1836250026L, -1701581766L, -1975018790L, 1372843048L, 
1726385151L, -292626507L, -1908932428L, 1108804204L, 1440018552L, 
-1682916752L, -264088150L, -350198328L, 511593804L, -1167867085L, 
261438471L, -2113638374L, 2007024934L, 891640134L, 568347596L, 
1898493480L, 1827884363L, -1093539901L, 1991923839L, 727476118L, 
2077156340L, -1094367422L, -1625392887L, -2141285655L, 1989605684L, 
-1479586090L, 51604506L, -914073462L, -257884920L, 1178632815L, 
647273840L, 818870120L, -30836433L, -1112676923L, 1713818208L, 
-137792171L, -1710102521L, -2068202152L, 16527693L, 1974980168L, 
1174511093L, 2133929714L, -1727591604L, 1362878030L, -1087283261L, 
-39189265L, 530922262L, 67468890L, 1621702634L, 1829410629L, 
-2037278381L, 233004634L, -1046507611L, 1034472363L, 1030963979L, 
1475101342L, 770309798L, -2049448647L, -133186321L, -1179211113L, 
-2081443L, -1174733217L, -1631412111L, 1613938712L, -562899475L, 
-1512935471L, -1843182048L, -2084847316L, 118915078L, 968512151L, 
2035276533L, -1997047764L, 2022700971L, 946758208L, 1804009662L, 
1847620496L, -1762856101L, -1296024150L, 2135215669L, 1544576245L, 
-765595663L, -410503643L, -1800458600L, -1975985274L, 1880775628L, 
-153488580L, -1969749179L, 444203072L, 504131503L, 352706357L, 
-1091225524L, 785235313L, -1004193413L, 234025727L, -57283920L, 
-803461240L, 1677239055L, 605229305L, 810956436L, -185634572L, 
-1928615981L, 971995325L, 70255738L, 980928928L, 2039178174L, 
-490288603L, 1645680650L, -388992414L, -262237464L, -1246229276L, 
-1340987279L, -1196968738L, -1751844676L, 723997992L, -354555659L, 
-383274852L, 550240783L, 1030053475L, -10503578L, 704889739L, 
-457527691L, -2004026515L, -1472450141L, 2006824471L, -791170005L, 
-815412056L, -885155338L, 649283124L, 1952722657L, -97296814L, 
-29952182L, -449541824L, -2126313265L, 1195820425L, -1980170402L, 
540131250L, -1129270288L, 117741799L, -850738149L, -1500752797L, 
6003877L, -1788283696L, 210392573L, 417154923L, 1066773880L, 
-1323372775L, 147683884L, -175418082L, 973072367L, 1867160023L, 
-292692127L, 67273051L, 759300819L, 1621833956L, 290905957L, 
-429310868L, -558501994L, 954892868L, -415670588L, 1547573890L, 
871646563L, -876073005L, 2062413604L, -1545008269L, -269118520L, 
1602700586L, 924480983L, 783927037L, 1895658725L, -1218652724L, 
1139103978L, 938946543L, 316879290L, -604446331L, 532770884L, 
2072721860L, -241410622L, -1285669116L, -465027045L, 1719497465L, 
2119509869L, 573687100L, 1197318344L, -658316666L, 1679209814L, 
-1810469252L, -153879405L, 1124587763L, 671539570L, 376144071L, 
1879698027L, 987554385L, -769396846L, -346962410L, 2007504627L, 
58267948L, 1129762787L, 924938983L, 151261963L, 999113930L, 287175682L, 
-2137172466L, -379010953L, 128844700L, 248338763L, -1851573103L, 
-856255517L, 715281585L, 54941629L, 2107239311L, 967881024L, 
885434066L, -1065933942L, 467422447L, 1713767251L, -88765030L, 
-384041140L, -1589755982L, -1343362784L, -1982270250L, -340305677L, 
1372394274L, 1463727022L, 59872232L, -1914195932L, 1272154810L, 
-1857541389L, 457658133L, 1745269658L, 1385285353L, -237557939L, 
265464586L, -694693793L, 1163302554L, -1312739001L, -1697806215L, 
1468096720L, 354106951L, 897209278L, 1957071404L, 2013317364L, 
-363726538L, -861355259L, -8556536L, -1101761621L, 1386923035L, 
1546516582L, -1966103392L, -1754769775L, -1999402816L, 1657898765L, 
-1794812079L, 959357537L, 158441694L, 1349436254L, -1270295700L, 
448847497L, 453202132L, 265284360L, 1341377640L, 1336943326L, 
1440214910L, -309302522L, -1988274396L, -290211525L, 1086333424L, 
1094500629L, -306574081L, -2023522344L, -147531937L, 387451779L, 
1872938427L, -716003157L, -1354116042L, -1322118926L, -779206342L, 
1481943318L, -336231506L, -749276453L, 623646443L, 424551442L, 
-1217766051L, -1172606325L, 294619234L, -1802494879L, 2110251567L, 
-1493374894L, 2110911274L, 307546360L, -1329495772L, 1041779858L, 
1330659343L, -712478049L, -912468139L, -951149257L, -1495860085L, 
-2025796336L, 25318215L, 888549828L, 1456295849L, 1351192537L, 
-2084771218L, 173339296L, 588868710L, 1503500359L, 338124489L, 
1292174979L, -1956819871L, 1512276386L, 1629882097L, 99781532L, 
706422146L, 1680574375L, 1592734487L, 312448783L, 433677495L, 
82767607L, 8296934L, 4012093L, 990159949L, -1320947075L, 1617991811L, 
1871516643L, 551948002L, 1205252114L, -175376780L, 1833324289L, 
931741400L, -2094756279L, -1650719881L, -1225948519L, 914313974L, 
-163457958L, -191448572L, 1571156189L, -631002367L, 443433088L, 
-437711196L, -1194324002L, -2070394526L, -683834907L, 1524956266L, 
-1879936991L, 1831692464L, -776887737L, -750117173L, 2072596860L, 
35468504L, 513710386L, -170448935L, -278201888L, -1322468224L, 
-715613251L, 962341447L, -1175928201L, -247507879L, -697057831L, 
-1797590324L, 1963683489L, -1297009733L, 104691873L, 1193723877L, 
-1775852176L, 574435067L, -1059225021L, 248368841L, -1461363093L, 
100629050L, -1322785911L, 201755008L, 1141582417L, 1819335659L, 
-37154701L, -2046065726L, 508559574L, -172263065L, -2057220016L, 
-284095054L, -310551804L, 224607573L, -223529307L, -257201737L, 
-648463021L, -76406836L, 264127073L, -1874733701L, 816073160L, 
1354926383L, -1997070151L, -100649966L, -1834819038L, 1092355979L, 
-980625154L, -1965331566L, -875215752L, 384014497L, -1173659838L, 
-2033780678L, -1628134434L, 1862365615L, -125651935L, 2122798495L, 
1461797692L, -1683376099L, 642928603L, -698927242L, 126111331L, 
-2128780632L, -1138064877L, -73265361L, 1623299920L, -1271671704L, 
2058675943L, -1583235948L, -1989732039L, 1456732496L, -1892536824L, 
1811712302L, 170925215L, 1522233118L, -1076111301L, -596284207L, 
198156487L, 1619307961L, 1844518593L, 617614110L, 2058429650L, 
-164692449L, 1850347398L, -1268645895L, -1441640255L, 571799119L, 
-384801109L, 683731115L, -1681201030L, 139437806L, 1600980004L, 
1646758649L, 881360486L, 1042109088L, -1983320444L, -759537848L, 
-12579803L, -1751631495L, 1524782044L, -1945007867L, 1294656484L, 
-1711070404L, 415298537L, -1056001045L, -1171920572L, -118775131L, 
294269164L, -190148136L, -1421687657L, -456267842L, -1220543646L, 
-658898740L, 317026859L, -445811892L, 966084703L, -200157255L, 
-124422270L, 1986819866L, 823118988L, -1308004996L, -457171659L, 
247804965L, 370292899L, 662441786L, -2010039331L, 1975707818L, 
-2145743588L, 1891845095L, -164574771L, 415747045L, -2106755805L, 
1937541943L, -774337647L, 1108012186L, -1283462801L, 1612920439L, 
1338054146L, -1213645723L, -1807206357L, -692046348L, 1317511062L, 
-1110510964L, -320165000L, 1468391L, 1428287965L, -1335001420L, 
-1635222328L, 1689925916L, -1885297598L, 1407684966L, 1655619786L, 
1289694378L, -1006846571L, 1885693037L, -1107439753L, 769201198L, 
1234826105L, -356754536L, 361192080L, 514115149L, 850228075L, 
1738289929L, -264537029L, 1925959733L, -1012083832L, -1586669457L, 
575436058L, 1170727247L, 457548198L, 103125867L, 2044260643L, 
1912546163L, -1364796340L, -580681423L, -435671296L, 189854745L, 
164092732L, -536126731L, -1491475882L, -1666004136L, -1925388631L, 
957672249L, -48987670L, 610708570L)
.seFromPBeta <-
function(v){x = abs(v[1])/qnorm(v[2]/2,lower.tail=FALSE); x}

.tabulateWithWeights <-
function(phen,vars,weights,totweight){
  tab = table(list(vars,phen))
  nme = dimnames(tab)
  nme1 = as.numeric(nme[[1]])
  nme2 = as.numeric(nme[[2]]) 
  tab1 = matrix(0,nrow = length(nme1), ncol  = length(nme2))
  for(i in 1:length(nme1)) 
    for(j in 1:length(nme2)) 
      tab1[i,j] = sum(weights[vars==nme1[i] & phen==nme2[j]])
  tab2 = tab1/totweight
  tab2
}


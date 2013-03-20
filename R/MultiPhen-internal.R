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
  if(sum(is.na(res)) == 2){
    res1 = res
  }
  else {
    ll1 = logLik(res)
    beta = summary(res)$coeff[ind,1]
    nullfam = is.null(family)
    if(hasCov & !nullfam) 
      ll2 = logLik(update(res, ~covar,family=family)) 
    else 
      if(!nullfam) 
        ll2 = logLik(update(res, ~1,family=family)) 
      else 
        if(hasCov) ll2 = logLik(update(res, ~covar)) 
        else 
          ll2 = logLik(update(res, ~1))
    df = attr(ll1,"df")[1]-attr(ll2,"df")[1]
    # print(paste("df",df))
    pv = pchisq((2*(ll1 - ll2))/(totweight),df,lower.tail=FALSE)
    res1 = c(beta,pv) 
  }
  res1
}
.extractSig1P <-
function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial"){
  if(sum(is.na(res)) == 2){
    res1 = res
  }
  else {
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
    res1 = cbind(c(beta,NA),c(pvs,pv))
  }
  res1
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
  if(hasCov) 
    ll2 = logLik(update(res, ~covar,family=family)) 
  else 
    ll2 = logLik(update(res, ~1,family=family))
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
 tryCatch(polr(vars ~ phenvec + covar,Hess=TRUE,method="logistic",weights = weights), error = function(e) print(c(NA, NA)))
}
.getGlmRevOrdNoCov <-
function( vars, phenvec, weights){
  tryCatch(polr(vars ~ phenvec,Hess=TRUE,method="logistic",weights = weights), error = function(e) print(c(NA, NA)))
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
function(stend,vec) vec>=stend[1] & vec<stend[2]
.getin1 <-
function(genvec,spl) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,.getin, genvec)
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
}
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
    .extractSig1P(glm.res, ind = ind, totweight = totweight, hasCov = !emptyCov, family = family)
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

.metaresHetMeasure<-function(matr){
  inds = !is.na(matr[,2])
  metres = metagen(matr[inds,1,drop=FALSE],apply(matr[inds,,drop=FALSE],1,.seFromPBeta))
  pv = .pFromSEBeta(c(metres$TE.fixed,metres$seTE.fixed))
  c(metres$Q,pchisq(metres$Q,1,lower.tail=FALSE)) 
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
function(genoData, phenoData, phenotypes = dimnames(phenoData)[[2]], covariates = NULL, resids = NULL, strats = NULL, maf_thresh = 0.001, corThresh = 0.0, inverseRegress = FALSE, JointModel = TRUE, multiGen = FALSE, fillMissingPhens = FALSE, scoreTest = FALSE, exactTest = FALSE, exactMethod = "wald", imputed = FALSE){
  rescale = 1
  exactTest = FALSE
  exactMethod = "wald"
  expandData = FALSE
  multiPhen = JointModel
  lmt1 = cbind(rep('pheno', length(phenotypes)), phenotypes)
  cvt1 = cbind(rep('covar', length(covariates)), covariates)
  rsd1 = cbind(rep('resid', length(resids)), resids)
  str1 = cbind(rep('strat', length(strats)), strats)
  limit = rbind(lmt1, cvt1, rsd1, str1)
  samples = as.matrix(dimnames(phenoData)[[1]]) 
  dimnames(samples) = list(NULL,"id")
  max_indiv = 100000
  #if(imputed){
  #  genoData = as.matrix(apply(genoData,1,.mod))
  #  if(inverseRegress & !expandData) genoData = apply(genoData,c(1,2),round)
  #  dimnames(genoData) = list(NULL,'rsID')
  #} 
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
  resid = todo_[grep("^resid",todo_[,1]),2]
  genocov = todo_[grep("^genocov",todo_[,1]),2]
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
  inds = match(ids, geno_header)
  #for(i in 1:length(ids)) inds = c(inds,grep(ids[i],geno_header))
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
  if(multiPhen) phenN = c(phenN,"JointModel")
  res = array(NA,dim=dimcounts,dimnames=list(stratificNames, phenN, geno_header,c("beta","pvalue")))
  geno1 = datanme$geno
  geno1[,3,] = !apply(as.matrix(geno1[,1,]),c(1,2),is.na)
  geno2 = geno1
  ssamples = NULL
  if(expandData) geno2 = datanme2$geno
  maf = rep(NA,dim(geno1)[1])
  for(k in 1:length(geno_header)) maf[k] = .getMaf(as.numeric(geno1[k,1,]),baselevel[k],2)
  for(j in 1:(dim(datanme$strat)[2])){
    inclu = datanme2$strat[,j]
    res[j,,,] =  pvFunctMult(geno2,datanme2$pheno,families,inclu,datanme2$pheno_cov,rescale,pvFunct,funct)
    ssamples = c(ssamples, sum(apply(cbind(as.numeric(is.na(geno2[1,1,])), as.vector(apply(is.na(datanme2$pheno[,1,]), 2, sum)), as.vector(apply(is.na(datanme2$pheno_cov), 1, sum))),1, sum) == 0))
  }
  if(length(index_strat)>0){ 
    res[dim(datanme$strat)[2]+1,,,] =apply(apply(res[1:(dim(datanme$strat)[2]-1),,,,drop=FALSE],c(2,3),.metaresInvVarianceFixed),c(3,1),t)
    ssamples = c(ssamples, sum(apply(cbind(as.numeric(is.na(geno2[1,1,])), as.vector(apply(is.na(datanme2$pheno[,1,]), 2, sum)), as.vector(apply(is.na(datanme2$pheno_cov), 1, sum))),1, sum) == 0))
}
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
  results = list(Results = res1, nobs = ssamples)
  return(results)
}
.nonNumeric <-
function(vec) sum(is.na(as.numeric(vec)))==length(vec)
.nonNumeric1 <-
function(mat) apply(mat,2,.nonNumeric)
.nullDim <-
function(matr){x = length(dim(matr)); x}
.numLev <-
function(vec) length(levels(as.factor(vec)))
.ordTest <-
function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  #print(is.numeric(x[,1]))
  cats = max(y) + 1
  m = apply(x,2,mean,na.rm=TRUE)
  x = t(t(x) - m)
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
  include =  apply(

    as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  .ordTest2(y,x)
}
.ordTest2 <-
function(y,x){
  cats = max(y) + 1
  m = apply(x,2,mean,na.rm=TRUE)
  x = t(t(x) - m)
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
c(403L, 10L, 1828790116L, -106126681L, 1781369853L, 1907359220L, 
1618369242L, 575036453L, -2100244161L, 551208646L, -1305195696L, 
1620063091L, 1467360657L, 1338402752L, -407563874L, 1163320169L, 
-255471717L, -158550358L, -1812528948L, 870089103L, -2106141243L, 
211777916L, 1557350866L, 765670733L, 772941831L, -102791186L, 
1511796424L, -922664501L, 923776041L, -1944870280L, 3379110L, 
748652993L, -993182125L, 896554242L, -1217702860L, -1306090057L, 
2099857837L, 407210692L, -1737355542L, 1201774933L, -538159377L, 
-553786122L, -1913441248L, -846060061L, -13789567L, -1530732304L, 
2066358158L, -1000476423L, -1517618549L, -808068422L, -434104900L, 
366411519L, -176328107L, 1197720940L, 541494722L, -1351586275L, 
-1199638889L, 1595621374L, -2062865480L, 883558299L, 48508409L, 
1971645768L, 1890658422L, -625486415L, 391830499L, 1401219250L, 
1693943044L, 627410503L, 826265885L, 836231060L, -1567440070L, 
1936527237L, 799186015L, 2111821542L, -798433936L, 214460051L, 
-619888847L, -1104476704L, 1268081278L, 975065353L, -227675333L, 
-1163343094L, 1183580908L, 742096751L, 752066725L, -1886844324L, 
-1424381390L, -1714749139L, -2044360729L, 1780256526L, 696557800L, 
-1060023573L, -1304265847L, 421308184L, 1772530246L, 1173047969L, 
436612403L, 907048098L, 510622612L, -99073513L, 588695053L, -21142748L, 
-186884342L, 1021375093L, 935653839L, -943932714L, -1855201664L, 
987873475L, -148496799L, -268937520L, -463310802L, -2094795303L, 
216467947L, 1178714L, -1310595556L, -220788961L, 256547701L, 
1313866252L, 58938338L, 1657077949L, 562966327L, 1956354014L, 
-2133269224L, -850507653L, -354292199L, -873658712L, 1861354582L, 
629260497L, -350462845L, -1155268846L, -1692808540L, -404089241L, 
-509179971L, 180214068L, -80487782L, -1605232411L, -1124676993L, 
548295686L, -604797040L, -1442723277L, -1314774575L, 271139712L, 
1914543710L, -1870977623L, -438590373L, -1825637270L, 1091986700L, 
1618650959L, -1953541115L, -2011927876L, 2141445010L, -718038899L, 
352493255L, -1450838098L, -98934264L, -1350973941L, -247853335L, 
-1708947016L, 1573512806L, -808721663L, 164767123L, -1373274558L, 
-979254796L, -579232777L, -236740627L, 873048196L, -135508310L, 
-608639979L, -478874449L, -392725834L, -1995065120L, -1945648093L, 
-1487013055L, -1015814608L, -2053000114L, -133243719L, 2117245131L, 
-1190533766L, 233822332L, -392090433L, -1893832299L, -859732L, 
731091842L, 1594415453L, 2098780247L, -309293250L, 213235576L, 
-390931365L, -2105439303L, 210172424L, -493133642L, 1119786993L, 
-2109586909L, 1614735986L, -2142539580L, -1736042873L, -1653567907L, 
-40243116L, -757099910L, 1442687429L, 962288415L, 220288678L, 
-2095527888L, 1185323731L, -1800736527L, -1510212576L, -881038658L, 
723377865L, 1452063867L, -1720675766L, 899458732L, 975348143L, 
-799651483L, -371927268L, -1924834446L, -1738934291L, -261392345L, 
-183861170L, -1935317848L, -1367930197L, 259866569L, -2025516072L, 
-946851194L, 702964577L, -115480845L, 915970402L, -1678104620L, 
-641232169L, 1432828109L, 1588173988L, 1355174584L, 378907482L, 
-244912752L, 1971642748L, -1715249900L, -1916810942L, -1587313920L, 
-1220977732L, -1079637168L, 1634935074L, -784724664L, -818554868L, 
-558778276L, -1417685806L, 1832355840L, -980579500L, -1375424024L, 
-263313414L, -1682839856L, 1755518252L, -82126604L, 1581143954L, 
1969654880L, 993758476L, 601548512L, 1731833186L, 117244584L, 
1848158796L, -956102788L, -1030862030L, -903672208L, -1188405948L, 
1808870360L, 2087575674L, -1421879760L, -1740449092L, -1466861740L, 
-257000894L, 1729813408L, -732276548L, 122778064L, -1493219998L, 
1270047720L, 22874060L, 1754171132L, -176476174L, 941980800L, 
-1432786700L, -1779190648L, -817253062L, -1883757744L, 2026843244L, 
862097556L, -1038128814L, -1826783520L, -1779121972L, -1426401120L, 
226531714L, -1983538968L, 271255788L, 1276282940L, 133980274L, 
-2112370576L, 1073144356L, -1525611976L, 452794906L, -762478384L, 
-123566532L, -615830380L, 1010746498L, -1723594816L, -57639620L, 
-1773605232L, -1395129310L, 602992648L, 394058892L, -1480273892L, 
-1999743854L, 1335471232L, 1998210068L, -1462066392L, 400266298L, 
844926160L, -601648020L, 770639540L, -280604654L, -1108058144L, 
1195748300L, 278093280L, -1909482462L, -530242264L, 1150146828L, 
-1894519236L, 646722802L, -1195464464L, 238741316L, 692441176L, 
1483227706L, 847804656L, -1144958276L, -2075887084L, 424504514L, 
9266784L, -1685228164L, -1098451760L, -134339550L, 347072680L, 
-857289332L, 2007257532L, 415012914L, -594705600L, -409137484L, 
1772183880L, 1236640442L, -582516976L, 1160601836L, -1628886956L, 
330368786L, -196009568L, -2121275444L, 1424031584L, -771990334L, 
425586216L, 965286316L, -794836164L, -2090478734L, 1202966320L, 
-1045405660L, 1558490680L, -194686246L, -1497729136L, -619645700L, 
1253571988L, -991743166L, -1176811008L, -290715588L, 1416933968L, 
-100274014L, -1520868536L, 2045348364L, 1119746780L, -1166135214L, 
503532672L, 63899476L, -1927794712L, -571833606L, 1668007760L, 
1823646764L, 1678170484L, 452791954L, 1791140448L, -4738932L, 
68628832L, -461760286L, 1481641128L, 1694287564L, 874420732L, 
1892891570L, 631243760L, -1714925116L, 1769721176L, -197320326L, 
-1398648400L, -1820137028L, 1889386068L, -1925604926L, 1583604512L, 
-312806468L, 815820880L, -229770270L, 872456040L, -834504884L, 
2018127612L, 440350450L, 520568064L, 972582004L, 1854803464L, 
-2003409478L, 1514501584L, 1426307052L, -910906476L, 2125056210L, 
2008925280L, 570740044L, 272438048L, -1205822846L, -1482980888L, 
-526509076L, -88637892L, -621902222L, -174847120L, -297124572L, 
1951472184L, -529813094L, 877096400L, -1678449220L, -2007175532L, 
1818272386L, -1837548096L, -2034690244L, -967765360L, -1418315998L, 
1563948552L, -1987861748L, -1613634660L, 2039700498L, 992816000L, 
675248916L, -820447960L, -823685062L, 84902864L, -1944880532L, 
1402227764L, -1540042606L, 954366304L, -1444373172L, 2021266912L, 
242326434L, 27900840L, -267903860L, -1647070276L, -2129304462L, 
-530438160L, -2044279740L, 1171236952L, 1201823802L, -517806928L, 
21230286L, -1107269637L, -2020904867L, -2019846262L, -461367384L, 
282926353L, -172383229L, 60722948L, -565042382L, -964116089L, 
-984471599L, 1511098710L, -1235720396L, 830425669L, 1796781711L, 
1142185992L, -2020883418L, 804019411L, -727269835L, -273582702L, 
-210400928L, -1091948343L, 287851515L, -2123291700L, 1447288410L, 
-1837317297L, 190667289L, -1850140626L, -143960516L, 2070916013L, 
-1911946153L, 1420489120L, 2105019134L, -132044341L, -934178867L, 
312902394L, -846268104L, 1263157985L, 634307987L, 717565364L, 
1748872898L, 1600356055L, 1723985697L, -1694742810L, 500908004L, 
626061077L, 912869183L, 1044533976L, -1902405898L, 149094659L, 
-1990458811L, -1699343774L, 1583550928L, 1142894329L, -1792086997L, 
1685564316L, -951172150L, -1268319041L, -540680695L, -1584359650L, 
740105932L, -1864180067L, -1232757977L, -62818992L, -1563362322L, 
-902871909L, -1572360451L, 1707505130L, -697727352L, -726151119L, 
-2016236893L, 927011492L, -1638071470L, 833915559L, 977008369L, 
-207858122L, -1080575276L, 868634213L, 1657633903L, -2107341144L, 
-1177123834L, -1523216077L, 1720618389L, 313238002L, 807374400L, 
264826793L, -747458021L, -732935188L, 923077626L, 2096850991L, 
-1278121927L, -1579550130L, -768779620L, -259637235L, -1779196809L, 
-264848256L, 12188382L, -1761859157L, 1233478317L, 1828206362L, 
-8591400L, -362781631L, 576762611L, 1649702164L, -1608671070L, 
1424555447L, 710046081L, 1055770502L, 1793762756L, -290458379L, 
-484868769L, 130265016L, 1929831062L, -1928737629L, -894728603L, 
-1228512382L, -1027522320L, -924211047L, -1085947765L, 1221929340L, 
1409980330L, 52873183L, 1444250985L, -1255106178L, -1627494036L, 
2118521277L, 437674631L, 227391088L, 1147608334L, 1772855867L, 
-2097583843L, 1388436426L, 1643812072L, 1497319377L, -1988175421L, 
125797060L, 1876795762L, 208817223L, -1378934767L, -436647530L, 
118074612L, 747079301L, 425430991L, -1115038008L, -1820453018L, 
95047059L, 1445646069L, -913502510L, 982776096L, -208234231L, 
90951227L, 1463313292L, 692838938L, 553525903L, 1327017433L, 
-1959169426L, 1006297596L, -1512598675L, 1061666071L, 1340897504L, 
771710654L, 10914955L, -1193645299L, 863370042L, 2144383864L, 
599515425L, 1357515219L, 603994356L, 1176567778L)
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

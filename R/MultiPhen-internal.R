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
    .getGlmRevBinom(var1,t(phenvec[,include]),covar[include],weights[include],family) 
  else 
    if(todo)  
    glm.res = .getGlmRevOrd(var1,t(phenvec[,include]),covar[include],weights[include])
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
function(vec){ vec[vec==""] = 0; vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=TRUE); (as.numeric(vec)%*% 0:2)/1000}

.mod1 <-
function(vec, rescale){
  vec[vec==""] = 0
  vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=TRUE);
  #vec[vec<200]  =0
  vec1 = floor(rescale*(as.numeric(vec)/1000))
  ind = vec1==max(vec1)
  vec1[ind] = vec1[ind]+ rescale - sum(vec1)
  vec1
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
c(403L, 469L, -768284368L, 688998009L, -909376839L, 159314433L, 
-369409609L, -2049069687L, -87854763L, 2061683073L, -1409622537L, 
1531233442L, -157842496L, 67001414L, 1626584440L, -1650516136L, 
1580500050L, -839455268L, -1105938706L, -550519929L, 1417373011L, 
-1200715273L, -1522848567L, -882407745L, 1757814119L, -1627337425L, 
-1637721627L, -891789472L, -1263805174L, -275488356L, -2001086174L, 
1368828350L, -928910972L, 667127978L, -225676196L, -1313521339L, 
52630341L, 669479301L, -155863085L, -1973030819L, -405271719L, 
-1565014667L, -1076203053L, 1062141054L, -1598396356L, 151134146L, 
-2122097244L, -1740896380L, -2012142002L, 496495792L, 275031418L, 
-1858090285L, -1879422401L, -542359829L, -688647835L, 1637988451L, 
57937675L, -1261253117L, -1787369151L, 864458044L, -1041205114L, 
-831034536L, -1894196786L, -92799590L, 1870784160L, -658659122L, 
1484572968L, 883812449L, 1834430049L, -1937366759L, -1038085201L, 
764159393L, -264718579L, 1957661609L, 1730687743L, 1800956538L, 
-1002133624L, -817513906L, 1434969216L, 1523275968L, -1344396326L, 
-1143181676L, 988919142L, 682537631L, -1663213205L, 2050197471L, 
399123633L, -1506174457L, 1589084447L, -50131785L, 689292813L, 
-1743397624L, 125184578L, 1621179908L, 1667897658L, -470686234L, 
1142658108L, 545765026L, 997407124L, -1733323699L, -1391009331L, 
-653316323L, 186692491L, 575771829L, -1470844431L, 1774191197L, 
-1689337701L, 1027295958L, 671088292L, 1464932074L, -1490301588L, 
31920972L, 1121460214L, 1654360296L, 202621650L, -1062239829L, 
2122631570L, -851084882L, 1327101005L, -1186035047L, 469257708L, 
-1344667660L, -1126942365L, -1285934448L, 99020479L, 835532033L, 
-226862966L, -1184869156L, 1778767747L, -1228120603L, -1311008714L, 
-2142447831L, -1837259644L, 1293814436L, 1754416735L, -945783245L, 
1926465186L, 435392342L, -1223306443L, 552512638L, 349255737L, 
394207855L, -695366772L, -1335249154L, 1064885009L, 1228930759L, 
-2055899648L, -85588721L, 1850946582L, 1754591562L, 177152841L, 
989223725L, 590286120L, -729774064L, -683935889L, -515540716L, 
-1716722101L, -1641497539L, -886452338L, -1278238176L, -259757537L, 
-1127598519L, 178963146L, 813096549L, -199709208L, -1685417040L, 
-2082071445L, -525930729L, -1536972354L, -29318414L, -1595731167L, 
1659117714L, 143786741L, -1384271045L, 658335904L, 2026717394L, 
-1468586563L, 1523457307L, 1349071652L, 1793818091L, -1995702454L, 
419760486L, 1041958037L, 1061189793L, 760937332L, -1589448884L, 
-375984965L, 1051469144L, -1420960553L, 1137776361L, -1667670862L, 
-652867980L, 1614396267L, 512446141L, 1014156462L, 925246001L, 
-14195748L, -1773549828L, -1626633145L, 325675387L, 259635962L, 
466178014L, -77921059L, -1759639130L, -738669807L, 1365748999L, 
-642321612L, -489828762L, -1100565959L, -119438545L, -707208984L, 
-1221292137L, -1416503506L, 110786L, -894724111L, 568553493L, 
-1476408048L, 1443948808L, -1767948377L, 537082460L, 1552453411L, 
274123749L, 1082236438L, 69921592L, -166350489L, -508668447L, 
-1190987742L, 540458044L, 2126688631L, -617693728L, 779008360L, 
-1658573773L, 859753446L, 1896858247L, 206441012L, 1219863266L, 
-388259171L, -1295609776L, 695816384L, -1740906221L, 166922106L, 
-1275573471L, 1231470196L, -1230441496L, 1619023191L, -1448262416L, 
-934291992L, 725131455L, -636706094L, 1585262223L, 956275712L, 
-122564446L, 481477481L, -1081903808L, -900554512L, -1503501193L, 
-1524735938L, 351073141L, -884983848L, 1124434388L, -1104305649L, 
-1921586608L, -364589032L, 870239603L, 1203885766L, 1710896063L, 
1048452268L, 662714778L, 1042576893L, -915512760L, 1056798904L, 
-44063365L, 1102589882L, 1989146905L, 1460788756L, 1464880040L, 
-1063275329L, 553816400L, -793284824L, -1759280409L, -844477158L, 
216767607L, -16759680L, -937223494L, 35136209L, 1088682645L, 
-64653015L, -1566759545L, -1667825153L, -1502841331L, -951787616L, 
-87983959L, 496258285L, -1177490337L, -891298233L, 1221659573L, 
-1628094103L, -978359083L, 229881530L, 1824419577L, -1163834067L, 
-1791566743L, -24767563L, -170434849L, -1420708177L, 1540860733L, 
1790221928L, 1997031361L, 1693600513L, -1708976245L, -579423045L, 
-1240669763L, 1895273097L, 740078489L, -1400081638L, -1793714003L, 
4149701L, -503730459L, 1265765065L, 122932991L, 950658463L, 328181461L, 
-181836112L, -378577255L, -971713139L, -725083977L, 976357343L, 
-1066825763L, -915250191L, 305807621L, 1579269338L, 864219689L, 
-252325643L, 20308121L, -2115624443L, 1648327919L, 1996312039L, 
-2139487571L, 1260096656L, -673052903L, -1914992831L, -1368233901L, 
760652931L, 715867928L, -1215876364L, 740506761L, -307700896L, 
-454706542L, 150498906L, 248482751L, 196964653L, 351948790L, 
2085698414L, -2026710671L, -2089623010L, -916144718L, 2144631230L, 
132401993L, 303441191L, -1926844984L, 1889820108L, 1911679421L, 
-159514288L, -1123123666L, -1657294470L, 582471739L, -2051578199L, 
1124753782L, 360729574L, -1486921607L, -1957911954L, -782882526L, 
-1945991126L, 1499743141L, 2064880883L, 1395154200L, 1605876388L, 
769960473L, 1711259080L, -1770345646L, 266768906L, 706839503L, 
1607937389L, -1390049090L, -1864003218L, -1293727143L, -236607250L, 
884431362L, 2065564734L, -1623798415L, -719483985L, -674231376L, 
-936740924L, -918357811L, -1085510064L, 615757406L, 701418130L, 
1109362219L, -793291367L, 1200821254L, -1546741634L, 679979625L, 
1464474918L, 557571946L, -67160502L, -1315425779L, -1100373893L, 
-591838536L, 354998500L, 612396521L, 849725328L, -1038687198L, 
977778122L, 1576068031L, -1416689763L, -710789626L, 485124702L, 
-485536607L, -1324089218L, -171135502L, -1084715250L, -1353012887L, 
-1546123417L, -1720791208L, -1045751268L, -1697182483L, 1469038608L, 
-1343853954L, -244377526L, 858166139L, 1494036889L, -1148324074L, 
301782742L, -1346183591L, -1054788418L, 650585234L, 1785422602L, 
32644629L, 207702371L, -113017432L, 881622308L, -491769047L, 
216562168L, 488546290L, -505806838L, 1912760447L, -471661651L, 
1900401518L, 1412271214L, 1701305929L, 1559828702L, 1307421874L, 
478196830L, 605926577L, 87331663L, -1064177312L, -1872402604L, 
-1309821747L, -881622592L, -841678402L, 220352962L, -1979346341L, 
1199287609L, 3783894L, -1036381426L, -277815623L, -524688874L, 
-1047983910L, -2121489126L, -1662280163L, -500276629L, -2050118536L, 
-750057260L, 1450013993L, 1147167520L, 1952654834L, 1859567674L, 
-891575393L, -180035635L, -450706314L, -87770130L, 1940529585L, 
-1283780386L, 909926258L, 1076459646L, 1199099977L, 1693986087L, 
1326525480L, -2045657332L, 572506781L, 722238864L, -255779410L, 
104800186L, -1843300613L, 1398439337L, -125031754L, -194025562L, 
-1946131815L, -1198467058L, -333342910L, 1868569194L, -1124341915L, 
152269747L, 1172666392L, 1087568100L, 820912345L, -2138341848L, 
615147442L, 923680810L, 2012808047L, -1703500883L, -959880930L, 
-874409266L, -1259943431L, -416675250L, 968904226L, 1455961086L, 
-896085199L, -178280465L, -1975618704L, 701121828L, -1543475475L, 
1234445904L, 1227846686L, 863541426L, 768879659L, 1658766361L, 
844835430L, 517394206L, -888806967L, -754651354L, 81521834L, 
-1731679126L, 2014749613L, 688632091L, -399278920L, 1765847172L, 
-1689080855L, -1098423152L, -1326848382L, -884095702L, -1455194689L, 
-8552067L, -1226808442L, 1106401310L, 865754529L, 66723038L, 
-1658205230L, -1803638130L, -311354551L, -106034617L, 1948783864L, 
1511898588L, 689559821L, -1994032720L, 692531518L, -1061945014L, 
802773659L, 714309593L, 342473270L, -110269674L, -1571858407L, 
417365214L, -1526358862L, 1269016234L, 45365845L, -172603229L, 
843765416L, 1608353988L, 1952472681L, 1736712413L, 195721263L, 
-1398362902L, -2031613811L, -1414087934L, -998516687L, -342186780L, 
27542045L, 324760223L, -449036205L, 668796338L, -1733726897L, 
-79824076L, -24660509L, -69181182L, 530056813L, -1297140563L, 
-1583964533L, -2139067238L, -217967759L, -128268674L, 1474153137L, 
-130240516L, -1149945771L, 184521823L, -311803837L, 638129022L, 
-1418718997L, -1613483624L, 754620307L, -1011536886L, 1687788273L, 
1299610069L, -1410772993L, -447314678L, 444998301L, -2060139694L, 
-2029234663L, -2028650252L, 1288415573L, -1689760049L, 1592762659L, 
1016629602L, -1895831897L, 1211190476L, -778804261L, 1638413322L, 
-1519823123L, -141189875L, 1624074315L, 2072709330L, -25262911L, 
698300606L, 986582694L, 290964512L, -1179836159L)
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


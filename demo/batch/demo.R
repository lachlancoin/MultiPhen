 ##Example running in batch mode

library(MultiPhen)

source("./init.txt")
mPhen.parseOpts()
print("check required params have been defined, or source from system")
mPhen.checkParams(c(".cohortDir",".analysisDir",".genoFile",".jointModel",".inverseRegress",".baseValue",".format",".phenoFiles",".excludeFile",".ending",".starting",".batch",".calcHWE",".writeQCFile",".plotQQ",".saveImage",".resultsName",".sigThresh",".fillMissingPhens",".thin",".quantileThresh",".writeBed"),     ls(all.names=T))


if(.calcHWE)library(HardyWeinberg)



.genoFile=paste(.cohortDir,.genoFile,sep="/")


print("NOW PREPARING PHENO FILE")
print(.phenoFiles)
print(.excludeFile)
pheno = mPhen.readPhenoFiles(.phenoFiles,.excludeFile,.naRowThresh = .naRowThresh, .naColThresh = .naColThresh,fillMissingPhens = .fillMissingPhens,.quantileThresh = .quantileThresh)
print("NOW READING LIMIT FILE")
limits = mPhen.readLimitFile("./limit.txt",dimnames(pheno)[[2]])
todo_ = limits$todo
suffix = limits$outn
if(.jointModel) suffix = paste(suffix,"J",sep="")




print("Now opening connection to geno file")
genoInput<-mPhen.openGenoConnection(.genoFile, indiv = rownames(pheno))
pheno_ = mPhen.preparePheno(pheno, genoInput$sampleids, limit=todo_,  fillMissingPhens =FALSE)
print("Now reading from geno file in first batch")
input<-mPhen.readGenoConnection(genoInput,.batch,.baseValue=.baseValue,.format=.format,.starting = .starting, .ending = .ending,.thin=.thin)

resultsDirParent=paste(.analysisDir,"results",sep="/")
dir.create(resultsDirParent)
resultsDir=paste(resultsDirParent,paste(.resultsName,suffix,sep="_"),sep="/")

###GET CHROM ###

resultsName=paste(input$chrom,input$firstPos,.resultsName, sep="_")
output<-mPhen.openOutputConnection(resultsDir,resultsName,limits,.writeQC=.writeQCFile,.writeBed = .writeBed)
sampleQC<-NULL




setwd("~/Documents/MRes/MyotonicDystrophyAffy2012")
#setwd("C:/JohnMcClure/People/DarrenMonckton/MyotonicDystrophyAffy2012")
#setwd("J:/CAMS/ResearchData/McClure/People/DarrenMonckton/MyotonicDystrophyAffy2012")
load("DM_MyoDys.RData")


#save.image("DM_MyoDys.RData") 

add.rownames <- function(x,nam="ID"){
  x1 <- cbind(rownames(x),x)
  colnames(x1) <- c(nam,colnames(x))
  x1
}
source("http://bioconductor.org/biocLite.R")
biocLite("beadarray")
a
biocLite("illuminaHumanv4.db")
a
biocLite("arrayQualityMetrics")
a
biocLite("lumi")
a
biocLite("affycoretools")
a
biocLite("GrapheR")
#not available for R 3.1.1
biocLite("ggplot2")
#not available for R 3.1.1
biocLite("RCurl")
#RCurl not available for R 3.1.1

library(beadarray)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
library(lumi)
library(limma)
library(affycoretools)
library(GrapheR)

#some error messages loading these libraries

DMdat <- read.table("Gene_level_expression.txt",sep="\t",header=T,as.is=T)
DMdat <- DMdat[,colnames(DMdat)!="X"]

#For duplicated entries:
DMnorm <- DMdat[!duplicated(DMdat[,1]),-1]
rownames(DMnorm) <- unique(DMdat[,1])
dups <- unique(DMdat[duplicated(DMdat[,1]),1])
for(i in 1:length(dups)){
  DMnorm[dups[i],] <- apply(DMdat[DMdat[,1]==dups[i],-1],2,median,na.rm=T)
}

sn <- colnames(DMnorm)
datIDs <- sub("X","",unlist(strsplit(unlist(strsplit(colnames(DMnorm),"HUEX1A11")),"HUEX1A21")))

colnames(DMnorm) <- sub("X","",colnames(DMnorm))
#this has been removed as it seems to cause problems by removing the X from the name

DMallinfo <- read.table("MASTERAnonymizedClinicalDataset23Feb2012.txt",sep="\t",header=T,skip=5,as.is=T)
rownames(DMallinfo) <- DMallinfo[,"DMBDI.PatID"]


#readDMb is a function for reading in data - read.table doesn't like data with empty columns at the end of excel files
readDMb <- function(rl){
  lrl <- length(rl)
  coln <- unlist(strsplit(rl[3],"\t"))
  nc <- length(coln)
  out <- matrix(NA,ncol=nc,nrow=lrl-3)
  for(i in 4:lrl){
    tmp <- unlist(strsplit(rl[i],"\t"))
    out[i-3,1:length(tmp)] <- tmp
  }
  colnames(out) <- coln
  out[out==""] <- NA

  out <- as.data.frame(out,stringsAsFactors=F)
  
  for(j in 1:nc){
    if(length(table(as.numeric(names(table(out[,j])))))>0) out[,j] <- as.numeric(out[,j])
  }
  out
}

rl <- readLines("DMBDI genotyping blood DNA Glasgow data with DMBDI ID for John.txt")

DMbloodInfo <- readDMb(rl)
# This didn't work - errors due to undefined columns. I tried to go through the backtrace and fix options but couldn't work it out


DMbloodInfo[,"Est.leng"] <- apply(DMbloodInfo[,c("Normal.allele.1","Normal.allele.2","Expanded.EstProgAllLeng")],1,FUN=function(x) max(c(mean(x[1:2],na.rm=T),x[3]),na.rm=T)) 
DMbloodInfo[!is.finite(DMbloodInfo[,"Est.leng"]),"Est.leng"] <- NA
DMbloodInfo[,"Modal.Alle"] <- apply(DMbloodInfo[,c("Normal.allele.1","Normal.allele.2","Expanded.ModalAllele")],1,max,na.rm=T) 
DMbloodInfo[!is.finite(DMbloodInfo[,"Modal.Alle"]),"Modal.Alle"] <- NA


alldat.inf <- DMbloodInfo[c(which(!is.na(DMbloodInfo[,"blood cel.file"])),which(!is.na(DMbloodInfo[,"muscle cel file"]))),]
rownames(alldat.inf) <- c(DMbloodInfo[which(!is.na(DMbloodInfo[,"blood cel.file"])),"blood cel.file"],
                          DMbloodInfo[which(!is.na(DMbloodInfo[,"muscle cel file"])),"muscle cel file"])
alldat.inf[,"sampletype"] <- c(rep("blood",length(which(!is.na(DMbloodInfo[,"blood cel.file"])))), 
                               rep("muscle",length(which(!is.na(DMbloodInfo[,"muscle cel file"])))))
alldat.inf <- alldat.inf[colnames(DMnorm),]

alldat.inf[,"Comb"] <- paste(alldat.inf[,"Patient no."],substr(alldat.inf[,"sampletype"],1,1),
                             alldat.inf[,"Est.leng"], substr(alldat.inf[,"Disease.Status"],1,1), alldat.inf[,"Modal.Alle"],sep="")

alldat.inf[,"group"] <- paste(substr(alldat.inf[,"sampletype"],1,1), substr(alldat.inf[,"Disease.Status"],1,1),sep="")


alldat.inf[,c("DOB","Date.of.Blood.draw")] <- NA
alldat.inf[alldat.inf[,"sampletype"]=="blood",c("DOB","Date.of.Blood.draw")] <- DMallinfo[as.character(alldat.inf[alldat.inf[,"sampletype"]=="blood","DMBDI ID"]),c("DOB","Date.of.Blood.draw")]

alldat.inf[,"Date.of.biopsy"] <- NA
alldat.inf[alldat.inf[,"sampletype"]=="muscle","Date.of.biopsy")] <- DMallinfo[as.character(alldat.inf[alldat.inf[,"sampletype"]=="blood","DMBDI ID"]),c("DOB","Date.of.Blood.draw")]


alldat.inf[,c("DOB","Date.of.Blood.draw","Date.of.biopsy")] <- NA
alldat.inf[,c("DOB","Date.of.Blood.draw","Date.of.biopsy")] <- DMallinfo[as.character(alldat.inf[,"DMBDI ID"]),c("DOB","Date.of.Blood.draw","Date.of.biopsy")]


dotm <- c(31,28,31,30,31,30,31,31,30,31,30,31)
pdotm <- c(0,cumsum(dotm[1:11]))

conv.date2num <- function(ddmmyyyy,sep="/"){
  sepdat <- strsplit(ddmmyyyy,sep)
  dats <- unlist(lapply(sepdat,FUN=function(x) as.numeric(x[3]) + (pdotm[as.numeric(x[2])] + as.numeric(x[1]))/365))
  dats
}

alldat.inf[,"Age.Blood.draw"] <- round(conv.date2num(alldat.inf[,"Date.of.Blood.draw"])-conv.date2num(alldat.inf[,"DOB"]),1)

alldat.inf[,"Age.biopsy"] <- round(conv.date2num(alldat.inf[,"Date.of.biopsy"])-conv.date2num(alldat.inf[,"DOB"]),1)


###############################
# Quality Checks
###############################

#Quality checking:

#min 0.71
boxplot(DMnorm,log="y")

#remember logshift is arbitrary chosen cutoff - 32
logshift <- 32
log2DMnorm <- log2(DMnorm + logshift)

pairs(log10(DMnorm[,1:6]),pch=".")

pairs(log10(DMnorm[,1:6]+logshift),pch=".")


log2DMnorm.aols <- log2DMnorm[apply(DMnorm,1,FUN=function(x) max(x)>=logshift),]



library(MASS)
#sammon
DMnorm.manh.dist <- dist(t(log2DMnorm.aols),method="manhattan")
DMnorm.sammon <- sammon(DMnorm.manh.dist)

#this chooses colours:
plot(c(rep(1:25,26),1:6),c(rep(1:26,each=25),rep(27,6)),pch=16,col=colors(),cex=2)

cols <- c("deepskyblue","hotpink","blue","red")
plot(1:13,col=cols,pch=16)
#chosen as two blues and two reds which are still easy to distinguish - i.e. two sets in two groups

xpand <- function(range,factor=1.2){
  width <- range[2]-range[1]
  mean <- mean(range)
  lo <- mean - factor*width/2 
  hi <- mean + factor*width/2     
  c(lo,hi)
}
#expand is a function to allow writing names under points so that points don't show up just as names - for clarity
reps <- rep(9,)

pdf("QC/SammonPlot_Group.pdf")
plot(DMnorm.sammon$points,type="n",xlim=xpand(range(DMnorm.sammon$points[,1]),factor=1.15))
text(DMnorm.sammon$points,alldat.inf[,"Comb"],col=cols[as.numeric(as.factor(alldat.inf[,"group"]))],cex=0.8)
dev.off()

pdf("QC/SammonPlot_Subject.pdf")
plot(DMnorm.sammon$points,type="n",xlim=xpand(range(DMnorm.sammon$points[,1]),factor=1.15))
text(DMnorm.sammon$points,as.character(alldat.inf[,"Patient no."]),col=alldat.inf[,"Patient no."],cex=0.8)
dev.off()

Bnorm <- DMnorm[,alldat.inf[,"sampletype"]=="blood"]
Mnorm <- DMnorm[,alldat.inf[,"sampletype"]=="muscle" & colnames(DMnorm)!="189821HUEX1A11"]

log2Bnorm <- log2DMnorm[,alldat.inf[,"sampletype"]=="blood"]
log2Mnorm <- log2DMnorm[,alldat.inf[,"sampletype"]=="muscle" & colnames(log2DMnorm)!="189821HUEX1A11"]
#189821HUEX1A11 is excluded as there is no allele length data


#####################################
# Useful functions
#####################################

# To calculate residuals:
  calcRes <- function(efit,data){
    print(dim(efit$design))
    print(dim(t(efit$coef)))
    pred <- t(efit$design %*% t(efit$coef))
    resid <- data - pred  
    resid
  }




#####################################
# Analysis
#####################################

############
# Blood data - Estimated Allele length + Age

ageBDpres <- rownames(alldat.inf[alldat.inf[,"sampletype"]=="blood" & !is.na(alldat.inf[,"Age.Blood.draw"]),])
designBlEstAge <- model.matrix(~alldat.inf[ageBDpres,"Est.leng"]+alldat.inf[ageBDpres,"Age.Blood.draw"],data=log2Bnorm[,ageBDpres])
colnames(designBlEstAge)[2:3] <- c("Est.leng","Age")
awBlEstAge <- arrayWeights(log2Bnorm[,ageBDpres],designBlEstAge)
fitBlEstAge <- lmFit(log2Bnorm[,ageBDpres], designBlEstAge, weights=awBlEstAge)

efitBlEstAge <- eBayes(fitBlEstAge)
topTable(efitBlEstAge,coef="Est.leng")
#ltoptabout(topTable(efitBlEstAge, coef="Est.leng",number=1000),nam="Blood.EstL.Age")
#ltoptableout again
BlEstAge.Limma <- topTable(efitBlEstAge, coef="Est.leng",number=nrow(log2Bnorm))

write.table(topTable(efitBlEstAge,coef="Est.leng",n=100),file="Limma/LimmaBlood_EstLeng&Age_EstLeng_Top100.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(topTable(efitBlEstAge,coef="Age",n=100),file="Limma/LimmaBlood_EstLeng&Age_Age_Top100.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(topTableF(efitBlEstAge,n=100),file="Limma/LimmaBlood_EstLeng&Age_Ftest_Top100.csv",sep=",",row.names=F,col.names=T,quote=F)

############
# Blood data - Estimate Allele length + Age with GSE

##### below here we are trying to create estlength/age/gse analysis
residBlEstAge <- calcRes(efitBlEstAge,log2Bnorm[,ageBDpres])
BlEstAge.gse <- apply(residBlEstAge,2,mean)

log2BnormEstAgeGSE <- sweep(log2Bnorm[,ageBDpres],2,BlEstAge.gse)


awBlEstAgeGSE <- arrayWeights(log2BnormEstAgeGSE[,ageBDpres],designBlEstAge)

fitBlEstAgeGSE <- lmFit(log2BnormEstAgeGSE[,ageBDpres], designBlEstAge, weights=awBlEstAgeGSE)

efitBlEstAgeGSE <- eBayes(fitBlEstAgeGSE)
topTable(efitBlEstAgeGSE,coef="Est.leng")
#ltoptabout(topTable(efitBlEstAge, coef="Est.leng",number=1000),nam="Blood.EstL.Age")
#ltoptableout again
BlEstAgeGSE.Limma <- topTable(efitBlEstAgeGSE, coef="Est.leng",number=nrow(log2Bnorm))

################
# 11 March 2015


# Muscle data - Age + Est Length

# Need to sort this  bit out
ageEstMBpres <- rownames(alldat.inf[alldat.inf[,"sampletype"]=="muscle" & !is.na(alldat.inf[,"Age.biopsy"]) & !is.na(alldat.inf[,"Est.leng"]),])
designMuEstAge <- model.matrix(~alldat.inf[ageEstMBpres,"Est.leng"]+alldat.inf[ageEstMBpres,"Age.biopsy"],data=log2Mnorm[,ageEstMBpres])
colnames(designMuEstAge)[2:3] <- c("Est.leng","Age")
awMuEstAge <- arrayWeights(log2Mnorm[,ageEstMBpres],designMuEstAge)
fitMuEstAge <- lmFit(log2Mnorm[,ageEstMBpres], designMuEstAge, weights=awMuEstAge)

efitMuEstAge <- eBayes(fitMuEstAge)
topTable(efitMuEstAge,coef="Est.leng")

MuEstAge.Limma <- topTable(efitMuEstAge, coef="Est.leng",number=nrow(log2Mnorm[,ageEstMBpres]))

write.table(topTable(efitMuEstAge,coef="Est.leng",n=100),file="Limma/LimmaMuscle_EstLeng&Age_EstLeng_Top100.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(topTable(efitMuEstAge,coef="Age",n=100),file="Limma/LimmaMuscle_EstLeng&Age_Age_Top100.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(topTableF(efitMuEstAge,n=100),file="Limma/LimmaMuscle_EstLeng&Age_Ftest_Top100.csv",sep=",",row.names=F,col.names=T,quote=F)


# Muscle data - Age + Est Length with GSE

residMuEstAge <- calcRes(efitMuEstAge,log2Mnorm[,ageEstMBpres])
MuEstAge.gse <- apply(residMuEstAge,2,mean)

log2MnormEstAgeGSE <- sweep(log2Mnorm[,ageEstMBpres],2,MuEstAge.gse)

awMuEstAgeGSE <- arrayWeights(log2MnormEstAgeGSE[,ageEstMBpres],designMuEstAge)

fitMuEstAgeGSE <- lmFit(log2MnormEstAgeGSE[,ageEstMBpres], designMuEstAge, weights=awMuEstAgeGSE)

efitMuEstAgeGSE <- eBayes(fitMuEstAgeGSE)
topTable(efitMuEstAgeGSE,coef="Est.leng")
MuEstAgeGSE.Limma <- topTable(efitMuEstAgeGSE, coef="Est.leng",number=nrow(log2Mnorm))



#####################################
# New bit 1:

# Blood data - Estimate Allele length + Age with DSE

# For a defined cut-off:
#    Top 100 most significant probes:
topB100IDs <- rownames(BlEstAge.Limma[1:100,]) # get their IDs
BlEstAge.dse100 <- apply(residBlEstAge[topB100IDs,],2,mean)

log2BnormEstAgeDSE100 <- sweep(log2Bnorm[,ageBDpres],2,BlEstAge.dse100)


awBlEstAgeDSE100 <- arrayWeights(log2BnormEstAgeDSE100[,ageBDpres],designBlEstAge)

fitBlEstAgeDSE100 <- lmFit(log2BnormEstAgeDSE100[,ageBDpres], designBlEstAge, weights=awBlEstAgeDSE100)

efitBlEstAgeDSE100 <- eBayes(fitBlEstAgeDSE100)
topTable(efitBlEstAgeDSE100,coef="Est.leng")
#ltoptabout(topTable(efitBlEstAge, coef="Est.leng",number=1000),nam="Blood.EstL.Age")
#ltoptableout again
BlEstAgeDSE100.Limma <- topTable(efitBlEstAgeDSE100, coef="Est.leng",number=nrow(log2Bnorm))


# For a defined cut-off:
#    Top 800 most significant probes (where Nplot lines from DSE100 cross line GSE):
##top800IDs <- BlEstAge.Limma[1:800,"ID"] # get their IDs [didn't work]
top800BIDs <- rownames(BlEstAge.Limma[1:800,]) # get their IDs
BlEstAge.dse800 <- apply(residBlEstAge[top800BIDs,],2,mean)

log2BnormEstAgeDSE800 <- sweep(log2Bnorm[,ageBDpres],2,BlEstAge.dse800)


awBlEstAgeDSE800 <- arrayWeights(log2BnormEstAgeDSE800[,ageBDpres],designBlEstAge)

fitBlEstAgeDSE800 <- lmFit(log2BnormEstAgeDSE800[,ageBDpres], designBlEstAge, weights=awBlEstAgeDSE800)

efitBlEstAgeDSE800 <- eBayes(fitBlEstAgeDSE800)
topTable(efitBlEstAgeDSE800,coef="Est.leng")
#ltoptabout(topTable(efitBlEstAge, coef="Est.leng",number=8000),nam="Blood.EstL.Age")
#ltoptableout again
BlEstAgeDSE800.Limma <- topTable(efitBlEstAgeDSE800, coef="Est.leng",number=nrow(log2Bnorm))

# Plots looking at correlation between GSE, DSE100 and DSE800
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(cbind(BlEstAge.gse,BlEstAge.dse100,BlEstAge.dse800), upper.panel=panel.cor)

#pairs(cbind(BlEstAge.gse,BlEstAge.dse100,BlEstAge.dse800))
#plot(BlEstAge.gse,BlEstAge.dse100)
#plot(BlEstAge.dse100,BlEstAge.dse800)

# End of new bit 1
#####################################

################
# 12 March 2015

# Muscle data - Estimate Allele length + Age with DSE

# For a defined cut-off:
#    Top 100 most significant probes:
top100MIDs <- rownames(MuEstAge.Limma[1:100,]) # get their IDs
MuEstAge.dse100 <- apply(residMuEstAge[top100MIDs,],2,mean)

log2MnormEstAgeDSE100 <- sweep(log2Mnorm[,ageEstMBpres],2,MuEstAge.dse100)



awMuEstAgeDSE100 <- arrayWeights(log2MnormEstAgeDSE100[,ageEstMBpres],designMuEstAge)

fitMuEstAgeDSE100 <- lmFit(log2MnormEstAgeDSE100[,ageEstMBpres], designMuEstAge, weights=awMuEstAgeDSE100)

efitMuEstAgeDSE100 <- eBayes(fitMuEstAgeDSE100)
topTable(efitMuEstAgeDSE100,coef="Est.leng")
MuEstAgeDSE100.Limma <- topTable(efitMuEstAgeDSE100, coef="Est.leng",number=nrow(log2Mnorm))


# For a defined cut-off:
#    Top 800 most significant probes:
top800MIDs <- rownames(MuEstAge.Limma[1:800,]) # get their IDs
MuEstAge.dse800 <- apply(residMuEstAge[top800MIDs,],2,mean)

log2MnormEstAgeDSE800 <- sweep(log2Mnorm[,ageEstMBpres],2,MuEstAge.dse800)


awMuEstAgeDSE800 <- arrayWeights(log2MnormEstAgeDSE800[,ageEstMBpres],designMuEstAge)

fitMuEstAgeDSE800 <- lmFit(log2MnormEstAgeDSE800[,ageEstMBpres], designMuEstAge, weights=awMuEstAgeDSE800)

efitMuEstAgeDSE800 <- eBayes(fitMuEstAgeDSE800)
topTable(efitMuEstAgeDSE800,coef="Est.leng")
MuEstAgeDSE800.Limma <- topTable(efitMuEstAgeDSE800, coef="Est.leng",number=nrow(log2Mnorm))


#####################################
# 4 March 2015

# For a defined cut-off:
#    Top 800 most significant probes (where Nplot lines from DSE100 cross line GSE):
##top800IDs <- BlEstAge.Limma[1:800,"ID"] # get their IDs [didn't work]
top101to800BIDs <- rownames(BlEstAge.Limma[101:800,]) # get their IDs
BlEstAge.dse101800 <- apply(residBlEstAge[top101to800BIDs,],2,mean)

log2BnormEstAgeDSE101800 <- sweep(log2Bnorm[,ageBDpres],2,BlEstAge.dse101800)


awBlEstAgeDSE101800 <- arrayWeights(log2BnormEstAgeDSE101800[,ageBDpres],designBlEstAge)

fitBlEstAgeDSE101800 <- lmFit(log2BnormEstAgeDSE101800[,ageBDpres], designBlEstAge, weights=awBlEstAgeDSE101800)

efitBlEstAgeDSE101800 <- eBayes(fitBlEstAgeDSE101800)
topTable(efitBlEstAgeDSE101800,coef="Est.leng")
#ltoptabout(topTable(efitBlEstAge, coef="Est.leng",number=8000),nam="Blood.EstL.Age")
#ltoptableout again
BlEstAgeDSE101800.Limma <- topTable(efitBlEstAgeDSE101800, coef="Est.leng",number=nrow(log2Bnorm))

pairs(cbind(BlEstAge.gse,BlEstAge.dse100,BlEstAge.dse800,BlEstAge.dse101800), upper.panel=panel.cor)

# to get an idea of how correlated top100 is with each next 100 up to 701to800

top101to200BIDs <- rownames(BlEstAge.Limma[101:200,]) # get their IDs
BlEstAge.dse101200 <- apply(residBlEstAge[top101to200BIDs,],2,mean)
top201to300BIDs <- rownames(BlEstAge.Limma[201:300,]) # get their IDs
BlEstAge.dse201300 <- apply(residBlEstAge[top201to300BIDs,],2,mean)
top301to400BIDs <- rownames(BlEstAge.Limma[301:400,]) # get their IDs
BlEstAge.dse301400 <- apply(residBlEstAge[top301to400BIDs,],2,mean)
top401to500BIDs <- rownames(BlEstAge.Limma[401:500,]) # get their IDs
BlEstAge.dse401500 <- apply(residBlEstAge[top401to500BIDs,],2,mean)
top501to600BIDs <- rownames(BlEstAge.Limma[501:600,]) # get their IDs
BlEstAge.dse501600 <- apply(residBlEstAge[top501to600BIDs,],2,mean)
top601to700BIDs <- rownames(BlEstAge.Limma[601:700,]) # get their IDs
BlEstAge.dse601700 <- apply(residBlEstAge[top601to700BIDs,],2,mean)
top701to800BIDs <- rownames(BlEstAge.Limma[701:800,]) # get their IDs
BlEstAge.dse701800 <- apply(residBlEstAge[top701to800BIDs,],2,mean)


pairs(cbind(BlEstAge.dse100,BlEstAge.dse101200,BlEstAge.dse201300,BlEstAge.dse301400,
            BlEstAge.dse401500,BlEstAge.dse501600,BlEstAge.dse601700,BlEstAge.dse701800), upper.panel=panel.cor)


# Two-stage Blood GSE then DSE800:
#    Top 800 most significant probes (where Nplot lines from DSE100 cross line GSE):
##top800BIDs <- BlEstAge.Limma[1:800,"ID"] # get their IDs [didn't work]
GSEtop800BIDs <- rownames(BlEstAgeGSE.Limma[1:800,]) # get their IDs
BlEstAge.gsedse800 <- apply(residBlEstAge[GSEtop800BIDs,],2,mean)

log2BnormEstAgeGSEDSE800 <- sweep(log2Bnorm[,ageBDpres],2,BlEstAge.gsedse800)

awBlEstAgeGSEDSE800 <- arrayWeights(log2BnormEstAgeGSEDSE800[,ageBDpres],designBlEstAge)
fitBlEstAgeGSEDSE800 <- lmFit(log2BnormEstAgeGSEDSE800[,ageBDpres], designBlEstAge, weights=awBlEstAgeGSEDSE800)
efitBlEstAgeGSEDSE800 <- eBayes(fitBlEstAgeGSEDSE800)
topTable(efitBlEstAgeGSEDSE800,coef="Est.leng")
BlEstAgeGSEDSE800.Limma <- topTable(efitBlEstAgeGSEDSE800, coef="Est.leng",number=nrow(log2Bnorm))

pairs(cbind(BlEstAge.gse,BlEstAge.dse100,BlEstAge.dse800,BlEstAge.dse101800,BlEstAge.gsedse800), upper.panel=panel.cor)


################
# 12 March 2015

# Two-stage Musle GSE then DSE800:
#    Top 800 most significant probes (where Nplot lines from DSE100 cross line GSE):
GSEtop800MIDs <- rownames(MuEstAgeGSE.Limma[1:800,]) # get their IDs
MuEstAge.gsedse800 <- apply(residMuEstAge[GSEtop800MIDs,],2,mean)

log2MnormEstAgeGSEDSE800 <- sweep(log2Mnorm[,ageEstMBpres],2,MuEstAge.gsedse800)

awMuEstAgeGSEDSE800 <- arrayWeights(log2MnormEstAgeGSEDSE800[,ageEstMBpres],designMuEstAge)
fitMuEstAgeGSEDSE800 <- lmFit(log2MnormEstAgeGSEDSE800[,ageEstMBpres], designMuEstAge, weights=awMuEstAgeGSEDSE800)
efitMuEstAgeGSEDSE800 <- eBayes(fitMuEstAgeGSEDSE800)
topTable(efitMuEstAgeGSEDSE800,coef="Est.leng")
MuEstAgeGSEDSE800.Limma <- topTable(efitMuEstAgeGSEDSE800, coef="Est.leng",number=nrow(log2Mnorm))

pairs(cbind(MuEstAge.gse,MuEstAge.dse100,MuEstAge.dse800,MuEstAge.gsedse800), upper.panel=panel.cor)


# End 4 March 2015
#####################################

# Blood data - Estimate Allele length + Age with GSE
# N plot:

plot(c(1,1000),c(1,1000),type="l",lty=3,ylab="N (GSE)",xlab="N (noGSE)")
lines(ncomp[1:2000,1],ncomp[1:2000,2],lty=1)

Nplot.pre <- function(pv1,pv2){
    print(date())
  pvcomp <- cbind(sort(pv1),sort(pv2)) # this will be two sets of ordered p-values in the real data
  pvcomb <- sort(unique(as.vector(pvcomp)))
  ncomp <- matrix(NA,nrow=length(pvcomb),ncol=3)
  for(i in 1:length(pvcomb)){
    ncomp[i,] <- c(sum(pvcomp[,1] <= pvcomb[i]) , sum(pvcomp[,2] <= pvcomb[i]) , pvcomb[i])
  }
    print(date())
  ncomp
}
ncomp <- Nplot.pre(BlEstAge.Limma[,"P.Value"],BlEstAgeGSE.Limma[,"P.Value"])


Nplot <- function(ncomp,N=nrow(ncomp)/2,add=F,colour="red"){
    if(add==F)  plot(c(1,N),c(1,N),type="l",lty=3,ylab="N (ISE)",xlab="N (noISE)")
    lines(ncomp[1:(2*N),1],ncomp[1:(2*N),2],lty=1,col=colour)
}
Nplot(ncomp,colour="grey70") 

Nplot(ncomp) 

Nplot(ncomp,1000) 

#different numbers here for different numbers of genes, 1000-5000 
# show most obvious, up to 10000 shows flattening out seen when whole gene set is run
#for write up maybe have 4000 and 10000 to show difference

summary(ncomp[,1]-ncomp[,2])

summary(ncomp[1:1000,1]-ncomp[1:1000,2]) # try different values e.g. 10000



# To run simulation of randomly permuted GSE values:

GSE <- BlEstAge.gse
NOGSE.Limma.pv <- BlEstAge.Limma[,"P.Value"]
design <- designBlEstAge
dataset <- log2Bnorm[,ageBDpres]
coeff="Est.leng"
N=10000

NplotSim <- function(NOGSE.Limma.pv, GSE, dataset, design, coeff, nsim, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    ncompl <- vector("list",length=nsim)
    for(i in 1:nsim){
        GSEdata <- sweep(dataset,2,sample(GSE))
        GSEaw <- arrayWeights(GSEdata,design)
        GSEfit <- lmFit(GSEdata, design, weights=GSEaw)
        GSEefit <- eBayes(GSEfit)
        GSE.Limma <- topTable(GSEefit, coef=coeff, number=nrow(dataset))
        ncompi <- Nplot.pre(NOGSE.Limma.pv,GSE.Limma[,"P.Value"])[,1:2] # gets rid of p-value column to save space
        ncompl[[i]] <- ncompi
    }
    ncompl
}
# this takes about 4 mins on John's machine: (2 mins per iteration):
Sim2Nplots <- NplotSim(NOGSE.Limma.pv, GSE, dataset, design, coeff, nsim=2, seed=NULL)


Sim2Nplots

#NplotSim(NOGSE.Limma.pv, GSE, dataset, design, coeff, N=10000, nsim=2, bgdcol="grey80", seed=NULL)


Nplotmult <- function(ncompl,N=nrow(ncomp)/2,add=F, colour="grey80"){
    plot(c(1,N),c(1,N),type="l",lty=3,ylab="N (GSE)",xlab="N (noGSE)")
    for(i in 1:length(ncompl)) lines(ncompl[[i]][1:(2*N),1],ncompl[[i]][1:(2*N),2],lty=1,col=colour)
}
Nplotmult(Sim2Nplots,N=10000,add=F, colour="grey80")
Nplot(ncomp,N=10000,add=T)

#
Sim100N <- NplotSim(NOGSE.Limma.pv, GSE, dataset, design, coeff, nsim=100, seed=NULL)
Nplotmult(Sim100N,N=10000,add=F, colour="grey80")
Nplot(ncomp,N=10000,add=T)


nulldist <- unlist(lapply(Sim2Nplots,FUN=function(x) max(x[,2] - x[,1])))
nulldist <- unlist(lapply(Sim100N,FUN=function(x) max(x[,2] - x[,1])))
pv = min(2*sum(nulldist>=max(ncomp[,2]-ncomp[,1]))/length(nulldist), 1)


 #       if(i==1) Nplot(ncompi,N) else Nplot(ncompi,N,add=T,colour=bgdcol)

test <- rbind(c(0,0),Sim2Nplots[[1]][1:20,])

#####################################
# New bit 2:

# Nplots with DSE100

ncompDSE100 <- Nplot.pre(BlEstAge.Limma[,"P.Value"],BlEstAgeDSE100.Limma[,"P.Value"])

Nplot(ncomp,colour="blue") 
Nplot(ncompDSE100,colour="red",add=T) 
title("N plots for GSE (blue) and DSE top 100 (red)")
# much worse for DSE100 over most genes

Nplot(ncomp,colour="blue",N=2000) 
Nplot(ncompDSE100,colour="red",N=2000,add=T) 
title("N plots for GSE (blue) and DSE top 100 (red)")
#

Nplot(ncomp,colour="blue",N=800) 
Nplot(ncompDSE100,colour="red",N=800,add=T) 
title("N plots for GSE (blue) and DSE top 100 (red)")

# Nplots with DSE800

ncompDSE800 <- Nplot.pre(BlEstAge.Limma[,"P.Value"],BlEstAgeDSE800.Limma[,"P.Value"])

Nplot(ncomp,colour="blue") 
Nplot(ncompDSE800,colour="red",add=T) 
title("N plots for GSE (blue) and DSE top 800 (red)")

Nplot(ncomp,colour="blue",N=1000) 
Nplot(ncompDSE800,colour="red",N=1000,add=T) 
title("N plots for GSE (blue) and DSE top 800 (red)")


Nplot(ncomp,colour="blue") 
Nplot(ncompDSE100,colour="red",add=T) 
Nplot(ncompDSE800,colour="green",add=T) 
title("N plots for GSE (blue) DSE top 100 (red) and DSE top 800 (green)")


Nplot(ncomp,colour="blue",N=1000) 
Nplot(ncompDSE100,colour="red",N=1000,add=T) 
Nplot(ncompDSE800,colour="green",N=1000,add=T) 
title("N plots for GSE (blue) DSE top 100 (red) and DSE top 800 (green)")

########################
# 4 March 2015 (part II)

# Nplots with DSE800

ncompDSE101800 <- Nplot.pre(BlEstAge.Limma[,"P.Value"],BlEstAgeDSE101800.Limma[,"P.Value"])

ncompGSEDSE800 <- Nplot.pre(BlEstAge.Limma[,"P.Value"],BlEstAgeGSEDSE800.Limma[,"P.Value"])

Nplot(ncomp,colour="blue") 
Nplot(ncompDSE100,colour="red",add=T) 
Nplot(ncompDSE800,colour="green",add=T) 
Nplot(ncompDSE101800,colour="purple",add=T) 
title("N plots for GSE (blue) DSE top 100 (red), DSE top 800 (green)\nand DSE top 101to800 (purple)")


Nplot(ncomp,colour="blue",N=1000) 
Nplot(ncompDSE100,colour="red",N=1000,add=T) 
Nplot(ncompDSE800,colour="green",N=1000,add=T) 
Nplot(ncompDSE101800,colour="purple",N=1000,add=T) 
title("N plots for GSE (blue) DSE top 100 (red), DSE top 800 (green)\nand DSE top 101to800 (purple)")


Nplot(ncomp,colour="blue") 
Nplot(ncompDSE100,colour="red",add=T) 
Nplot(ncompDSE800,colour="green",add=T) 
Nplot(ncompDSE101800,colour="purple",add=T) 
Nplot(ncompGSEDSE800,colour="orange",add=T)
title("N plots for GSE (blue) DSE top 100 (red), DSE top 800 (green)\nDSE top 101to800 (purple) and GSEDSE800 (orange)")

Nplot(ncomp,colour="blue",N=1000) 
Nplot(ncompDSE100,colour="red",N=1000,add=T) 
Nplot(ncompDSE800,colour="green",N=1000,add=T) 
Nplot(ncompDSE101800,colour="purple",N=1000,add=T) 
Nplot(ncompGSEDSE800,colour="orange",N=1000,add=T)
title("N plots for GSE (blue) DSE top 100 (red), DSE top 800 (green)\nDSE top 101to800 (purple) and GSEDSE800 (orange)")

Nplot(ncomp,colour="blue",N=3000) 
Nplot(ncompDSE100,colour="red",N=3000,add=T) 
Nplot(ncompDSE800,colour="green",N=3000,add=T) 
Nplot(ncompDSE101800,colour="purple",N=3000,add=T) 
Nplot(ncompGSEDSE800,colour="orange",N=3000,add=T)
title("N plots for GSE (blue) DSE top 100 (red), DSE top 800 (green)\nDSE top 101to800 (purple) and GSEDSE800 (orange)")

# To run simulation of randomly permuted DSE800 values:

NplotSim <- function(NOGSE.Limma.pv, GSE, dataset, design, coeff, nsim, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    ncompl <- vector("list",length=nsim)
    for(i in 1:nsim){
        GSEdata <- sweep(dataset,2,sample(GSE))
        GSEaw <- arrayWeights(GSEdata,design)
        GSEfit <- lmFit(GSEdata, design, weights=GSEaw)
        GSEefit <- eBayes(GSEfit)
        GSE.Limma <- topTable(GSEefit, coef=coeff, number=nrow(dataset))
        ncompi <- Nplot.pre(NOGSE.Limma.pv,GSE.Limma[,"P.Value"])[,1:2] # gets rid of p-value column to save space
        ncompl[[i]] <- ncompi
    }
    ncompl
}
# this takes about 4 mins on John's machine: (2 mins per iteration):
Sim2N800 <- NplotSim(NOGSE.Limma.pv=BlEstAge.Limma[,"P.Value"], GSE=BlEstAge.dse800, dataset=log2Bnorm[,ageBDpres], design=designBlEstAge, coeff="Est.leng", nsim=2, seed=NULL)

Sim2N800

#NplotSim(NOGSE.Limma.pv, GSE, dataset, design, coeff, N=10000, nsim=2, bgdcol="grey80", seed=NULL)


Nplotmult <- function(ncompl,N=nrow(ncomp)/2,add=F, colour="grey80"){
    plot(c(1,N),c(1,N),type="l",lty=3,ylab="N (GSE)",xlab="N (noGSE)")
    for(i in 1:length(ncompl)) lines(ncompl[[i]][1:(2*N),1],ncompl[[i]][1:(2*N),2],lty=1,col=colour)
}
Nplotmult(Sim2N800,N=10000,add=F, colour="grey80")
Nplot(ncomp,colour="blue",N=10000,add=T) 
Nplot(ncompDSE100,colour="red",N=10000,add=T) 
Nplot(ncompDSE800,colour="green",N=10000,add=T) 
Nplot(ncompDSE101800,colour="purple",N=10000,add=T) 

Nplotmult(Sim2N800,N=1000,add=F, colour="grey80")
Nplot(ncomp,colour="blue",N=1000,add=T) 
Nplot(ncompDSE100,colour="red",N=1000,add=T) 
Nplot(ncompDSE800,colour="green",N=1000,add=T) 
Nplot(ncompDSE101800,colour="purple",N=1000,add=T) 



# For Wilf to do:
Sim100N800 <- NplotSim(NOGSE.Limma.pv=BlEstAge.Limma[,"P.Value"], GSE=BlEstAge.dse800, dataset=log2Bnorm[,ageBDpres], design=designBlEstAge, coeff="Est.leng", nsim=100, seed=NULL)
Nplotmult(Sim100N800,N=10000,add=F, colour="grey80")
Nplot(ncomp,N=10000,add=T)



# Line plots and Venns
BlEstAgeTop100 <- topTable(efitBlEstAge,coef="Est.leng",n=100,sort.by="p")[1:100,]
pvcut100 <- topTable(efitBlEstAge,coef="Est.leng",n=100,sort.by="p")[100,4]
n100equivGSE <- max(which(topTable(efitBlEstAgeGSE,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
BlEstAgeGSETop100equiv <- topTable(efitBlEstAgeGSE,coef="Est.leng",n=n100equivGSE,sort.by="p")
n100equivDSE100 <- max(which(topTable(efitBlEstAgeDSE100,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
BlEstAgeDSE100Top100equiv <- topTable(efitBlEstAgeDSE100,coef="Est.leng",n=n100equivDSE100,sort.by="p")
n100equivDSE800 <- max(which(topTable(efitBlEstAgeDSE800,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
BlEstAgeDSE800Top100equiv <- topTable(efitBlEstAgeDSE800,coef="Est.leng",n=n100equivDSE800,sort.by="p")


source("venn234way.R")
v4wBlEstAgeGSEDSE100800 <- venn4way(rownames(BlEstAgeTop100),rownames(BlEstAgeGSETop100equiv),rownames(BlEstAgeDSE100Top100equiv),rownames(BlEstAgeDSE800Top100equiv),names=c("no GSE","GSE","DSE100","DSE800"),
         main="4-way Venn Diagram of genes with pv < noGSE Top100 pv\n from Limma analysisof Blood samples",
         sub=paste("pv cutpoint is ",signif(pvcut100,2),".  Model is Est Length + Age, coefficient of interest is Est Length",sep=""))


BlEstAgeTop100d100 <- topTable(efitBlEstAgeDSE100,coef="Est.leng",n=100,sort.by="p")[1:100,]
BlEstAgeTop100d100n <- rownames(BlEstAgeTop100d100)
pvcut100DSE100 <- topTable(efitBlEstAgeDSE100,coef="Est.leng",n=100,sort.by="p")[100,4]
n100D100equivDSE800 <- max(which(topTable(efitBlEstAgeDSE800,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100DSE100))
BlEstAgeDSE800Top100equiv <- topTable(efitBlEstAgeDSE800,coef="Est.leng",n=n100D100equivDSE800,sort.by="p")

D800equivposD100 <- charmatch(BlEstAgeTop100d100n,rownames(BlEstAgeDSE800Top100equiv))

v2wBlEstAgeDSE100800 <- venn2way(rownames(BlEstAgeTop100d100),rownames(BlEstAgeDSE800Top100equiv),names=c("DSE100","DSE800"),
         main="2-way Venn Diagram of genes with pv < DS100 Top100 pv\n from Limma analysis of Blood samples",
         sub=paste("pv cutpoint is ",signif(pvcut100DSE100,2),".  Model is Est Length + Age, coefficient of interest is Est Length",sep=""))
plot(rep(0:1,c(100,100)),c(1:100,1:100),type="n")
segments(rep(0,100),1:100,rep(1,100),D800equivposD100)
points(rep(0,length(v2wBlEstAgeDSE100800$A.only)),
       charmatch(v2wBlEstAgeDSE100800$A.only,BlEstAgeTop100d100n),col="red",pch=16)
points(rep(1,length(v2wBlEstAgeDSE100800$B.only)),
       charmatch(v2wBlEstAgeDSE100800$B.only,rownames(BlEstAgeDSE800Top100equiv)),col="blue",pch=16)


# 11 March 2015:

lrng <- log10(range(c(BlEstAgeTop100d100[,"P.Value"], BlEstAgeDSE800Top100equiv[,"P.Value"])))

rng <- range(c(BlEstAgeTop100d100[,"P.Value"], BlEstAgeDSE800Top100equiv[,"P.Value"]))

plot(0:1,rng,type="n",log="y",ylab="log10 p-value",xlab="",xaxt="n")
segments(rep(0,100),BlEstAgeTop100d100[v2wBlEstAgeDSE100800$AB,"P.Value"],rep(1,100),BlEstAgeDSE800Top100equiv[v2wBlEstAgeDSE100800$AB,"P.Value"])
points(rep(0,length(v2wBlEstAgeDSE100800$A.only)),
       BlEstAgeTop100d100[v2wBlEstAgeDSE100800$A.only,"P.Value"],col="red",pch=16)
points(rep(1,length(v2wBlEstAgeDSE100800$B.only)),
       BlEstAgeDSE800Top100equiv[v2wBlEstAgeDSE100800$B.only,"P.Value"],col="blue",pch=16)
mtext(c("DSE100","DSE800"),1,line=1,at=0:1)


BlEstAgeTop100d100
BlEstAgeDSE800Top100equiv

v2wBlEstAgeDSE100800$AB





plot(rep(0:1,c(100,100)),c(1:100,1:100),type="n")
segments(rep(0,100),1:100,rep(1,100),D800equivposD100)
points(rep(0,length(v2wBlEstAgeDSE100800$A.only)),
       charmatch(v2wBlEstAgeDSE100800$A.only,BlEstAgeTop100d100n),col="red",pch=16)
points(rep(1,length(v2wBlEstAgeDSE100800$B.only)),
       charmatch(v2wBlEstAgeDSE100800$B.only,rownames(BlEstAgeDSE800Top100equiv)),col="blue",pch=16)



# End 4 March 2015 (part II)
########################


# Plots of most significant p-values:
plot(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[,"P.Value"]-BlEstAgeGSE.Limma[,"P.Value"],type="l",ylim=c(-0.06,0.01),col="blue")
lines(c(0,44000),c(0,0),lty=2,col="grey")
lines(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[,"P.Value"]-BlEstAgeDSE100.Limma[,"P.Value"],lty=2,col="red")
lines(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[,"P.Value"]-BlEstAgeDSE800.Limma[,"P.Value"],lty=3,col="green")
title("Differences between noISE ordered pvales and\nGSE (blue), DSE100 (red) and DSE800 (green)")

plot(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[,"P.Value"]-BlEstAgeGSE.Limma[,"P.Value"],type="l",ylim=c(-0.003,0.002),xlim=c(0,1500),col="blue")
lines(c(0,44000),c(0,0),lty=2,col="grey")
lines(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[,"P.Value"]-BlEstAgeDSE100.Limma[,"P.Value"],lty=2,col="red")
lines(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[,"P.Value"]-BlEstAgeDSE800.Limma[,"P.Value"],lty=3,col="green")
title("Differences between noISE ordered pvales and\nGSE (blue), DSE100 (red) and DSE800 (green)")


plot(1:length(BlEstAge.Limma[,"P.Value"]),log2(BlEstAge.Limma[,"P.Value"]/BlEstAgeGSE.Limma[,"P.Value"]),type="l",ylim=c(-1.5,5),col="blue")
lines(c(0,44000),c(0,0),lty=2,col="grey")
lines(1:length(BlEstAge.Limma[,"P.Value"]),log2(BlEstAge.Limma[,"P.Value"]/BlEstAgeDSE100.Limma[,"P.Value"]),lty=2,col="red")
lines(1:length(BlEstAge.Limma[,"P.Value"]),log2(BlEstAge.Limma[,"P.Value"]/BlEstAgeDSE800.Limma[,"P.Value"]),lty=3,col="green")
title("Differences between noISE ordered pvales and\nGSE (blue), DSE100 (red) and DSE800 (green)")

plot(1:length(BlEstAge.Limma[,"P.Value"]),log2(BlEstAge.Limma[,"P.Value"]/BlEstAgeGSE.Limma[,"P.Value"]),type="l",ylim=c(-0.5,5),xlim=c(0,1500),col="blue")
lines(c(0,44000),c(0,0),lty=2,col="grey")
lines(1:length(BlEstAge.Limma[,"P.Value"]),log2(BlEstAge.Limma[,"P.Value"]/BlEstAgeDSE100.Limma[,"P.Value"]),lty=2,col="red")
lines(1:length(BlEstAge.Limma[,"P.Value"]),log2(BlEstAge.Limma[,"P.Value"]/BlEstAgeDSE800.Limma[,"P.Value"]),lty=3,col="green")
title("Differences between noISE ordered pvales and\nGSE (blue), DSE100 (red) and DSE800 (green)")

# Tried to look at pairing the p-values using the probe ID, but it didn't really work
#plot(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[rownames(BlEstAgeGSE.Limma),"P.Value"]-BlEstAgeGSE.Limma[rownames(BlEstAgeGSE.Limma),"P.Value"],type="l")
#plot(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[rownames(BlEstAgeGSE.Limma),"P.Value"]-BlEstAgeDSE100.Limma[rownames(BlEstAgeGSE.Limma),"P.Value"],type="l")
#plot(1:length(BlEstAge.Limma[,"P.Value"]),BlEstAge.Limma[rownames(BlEstAgeGSE.Limma),"P.Value"]-BlEstAgeDSE800.Limma[rownames(BlEstAgeGSE.Limma),"P.Value"],type="l")

# End of new bit 2
#####################################


################
# 12 March 2015

# Nplots for Muscle
 

ncompMu <- Nplot.pre(MuEstAge.Limma[,"P.Value"],MuEstAgeGSE.Limma[,"P.Value"])

Nplot(ncompMu,colour="grey70") 




# Nplots with DSE100

ncompMuDSE100 <- Nplot.pre(MuEstAge.Limma[,"P.Value"],MuEstAgeDSE100.Limma[,"P.Value"])

Nplot(ncompMu,colour="blue") 
Nplot(ncompMuDSE100,colour="red",add=T) 
title("Muscle ~ Est + Age; coef=Est\nN plots for GSE (blue) and DSE top 100 (red)")
#  worse for DSE100 over most genes, but better at the beginning

Nplot(ncompMu,colour="blue",N=2000) 
Nplot(ncompMuDSE100,colour="red",N=2000,add=T) 
title("Muscle ~ Est + Age; coef=Est\nN plots for GSE (blue) and DSE top 100 (red)")
#

Nplot(ncompMu,colour="blue",N=5000) 
Nplot(ncompMuDSE100,colour="red",N=5000,add=T) 
title("Muscle ~ Est + Age; coef=Est\nN plots for GSE (blue) and DSE top 100 (red)")
#

Nplot(ncompMu,colour="blue",N=800) 
Nplot(ncompMuDSE100,colour="red",N=800,add=T) 
title("Muscle ~ Est + Age; coef=Est\nN plots for GSE (blue) and DSE top 100 (red)")

# Nplots with DSE800

ncompMuDSE800 <- Nplot.pre(MuEstAge.Limma[,"P.Value"],MuEstAgeDSE800.Limma[,"P.Value"])

# Nplots with GSEDSE800

ncompMuGSEDSE800 <- Nplot.pre(MuEstAge.Limma[,"P.Value"],MuEstAgeGSEDSE800.Limma[,"P.Value"])

Nplot(ncompMu,colour="blue") 
Nplot(ncompMuDSE100,colour="red",add=T) 
Nplot(ncompMuDSE800,colour="green",add=T) 
Nplot(ncompMuGSEDSE800,colour="orange",add=T)
title("Muscle ~ Est + Age; coef=Est N plots for\nGSE (blue) DSE top 100 (red), DSE top 800 (green) and GSEDSE800 (orange)")

Nplot(ncompMu,colour="blue",N=1000) 
Nplot(ncompMuDSE100,colour="red",N=1000,add=T) 
Nplot(ncompMuDSE800,colour="green",N=1000,add=T)  
Nplot(ncompMuGSEDSE800,colour="orange",N=1000,add=T)
title("Muscle ~ Est + Age; coef=Est N plots for\nGSE (blue) DSE top 100 (red), DSE top 800 (green) and GSEDSE800 (orange)")

Nplot(ncompMu,colour="blue",N=3000) 
Nplot(ncompMuDSE100,colour="red",N=3000,add=T) 
Nplot(ncompMuDSE800,colour="green",N=3000,add=T) 
Nplot(ncompMuGSEDSE800,colour="orange",N=3000,add=T)
title("Muscle ~ Est + Age; coef=Est N plots for\nGSE (blue) DSE top 100 (red), DSE top 800 (green) and GSEDSE800 (orange)")

#
# DSE blood vs DSE muscle

BlMuIC <- intersect(alldat.inf[ageBDpres,"DMBDI ID"], alldat.inf[ageEstMBpres,"DMBDI ID"])

plot(BlEstAge.gse[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.gse[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] ,
     xlab="Blood GSE",ylab="Muscle GSE")
     
cor.test(BlEstAge.gse[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.gse[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] )


plot(BlEstAge.dse100[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.dse100[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] ,
     xlab="Blood DSE 100",ylab="Muscle DSE 100")
     
cor.test(BlEstAge.dse100[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.dse100[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] )


plot(BlEstAge.dse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.dse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] ,
     xlab="Blood DSE 800",ylab="Muscle DSE 800")
     
cor.test(BlEstAge.dse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.dse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] )


plot(BlEstAge.gsedse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.gsedse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] ,
     xlab="Blood GSEDSE 800",ylab="Muscle GSEDSE 800")
     
cor.test(BlEstAge.gsedse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="blood",])],
     MuEstAge.gsedse800[rownames(alldat.inf[is.element(alldat.inf[,"DMBDI ID"],BlMuIC) & alldat.inf[,"sampletype"]=="muscle",])] )

#venn2way(alldat.inf[ageBDpres,"DMBDI ID"], alldat.inf[ageEstMBpres,"DMBDI ID"])

ageEstMBpres

MuEstAge.dse100

BlEstAge.dse100

# Est leng vs GSE regression

plot(designMuEstAge[,"Est.leng"],MuEstAge.dse100)
lm1 <-  lm(MuEstAge.dse100~designMuEstAge[,"Est.leng"])
lines(designMuEstAge[,"Est.leng"],fitted(lm1))

###########################
# 18 March 2015

  calcRes <- function(efit,data){
    print(dim(efit$design))
    print(dim(t(efit$coef)))
    pred <- t(efit$design %*% t(efit$coef))
    resid <- data - pred  
    resid
  }
  

residBlEstAge <- calcRes(efitBlEstAge,log2Bnorm)
residMuEstAge <- calcRes(efitMuEstAge,log2Mnorm)


ttsum <- function(x){
  tt <- t.test(x,mu=0)
  output <- c(tt$p.value,tt$est,tt$conf)
  output
}
ttsum(residBlEstAge[,1])

ttresBlEstAge <- t(apply(residBlEstAge,2,FUN=function(x) ttsum(x)))
colnames(ttresBlEstAge) <- c("p.value","mean","CI.low","CI.high")
write.table(add.rownames(ttresBlEstAge),file="Residuals/LimmaResiduals_BloodAgeEstResults.xls",sep="\t",row.names=F,col.names=T,quote=F)

ttresMuEstAge <- t(apply(residMuEstAge,2,FUN=function(x) ttsum(x)))
colnames(ttresMuEstAge) <- c("p.value","mean","CI.low","CI.high")
write.table(add.rownames(ttresMuEstAge),file="Residuals/LimmaResiduals_MuscleAgeEstResults.xls",sep="\t",row.names=F,col.names=T,quote=F)

cumav <- function(x) cumsum(x)/(1:length(x))
# plot(log10(BlEst.Limma[,"P.Value"]),residBlEst[BlEst.Limma[,1],1],pch=".",col="red")
# lines(log10(BlEst.Limma[,"P.Value"]),cumav(residBlEst[BlEst.Limma[,1],1]))
#these two lines won't work - ask John
lines(par("usr")[1:2],c(0,0),lty=2)

for(j in 1:4){
  npl <- if(j<4) 9 else 8
  pdf(paste("Residuals/LimmaResults_Pval.vs.Resid_BloodEstAge_Est_",j,".pdf",sep=""))
  par(mfrow=c(3,3))
  for(i in 1:npl){
    plot(log10(BlEstAge.Limma[,"P.Value"]),residBlEstAge[rownames(BlEstAge.Limma),(j-1)*9+i],xlab="log10 p-value",ylab="Limma residual",
         main=colnames(residBlEstAge)[(j-1)*9+i],ylim=c(-0.4,0.4),pch=".",col="red")
    lines(log10(BlEstAge.Limma[,"P.Value"]),cumav(residBlEstAge[rownames(BlEstAge.Limma),(j-1)*9+i]))
    lines(c(-100,0),c(0,0),lty=2)
  }
  dev.off()
}
#error in the for loop. All the errors just above are to do with
#"only 0s may be mixed with negative subscripts"

for(j in 1:4){
  npl <- if(j<4) 9 else 2
  pdf(paste("Residuals/LimmaResults_Pval.vs.Resid_MuscleEstAge_Est_",j,".pdf",sep=""))
  par(mfrow=c(3,3))
  for(i in 1:npl){
    plot(log10(MuEstAge.Limma[,"P.Value"]),residMuEstAge[rownames(MuEstAge.Limma),(j-1)*9+i],xlab="log10 p-value",ylab="Limma residual",
         main=colnames(residMuEstAge)[(j-1)*9+i],ylim=c(-0.4,0.4),pch=".",col="red")
    lines(log10(MuEstAge.Limma[,"P.Value"]),cumav(residMuEstAge[rownames(MuEstAge.Limma),(j-1)*9+i]))
    lines(c(-100,0),c(0,0),lty=2)
  }
  dev.off()
}

# line plots NoISE vs GSE

#BlEstAgeTop100 <- topTable(efitBlEstAge,coef="Est.leng",n=100,sort.by="p")[1:100,]
#pvcut100 <- topTable(efitBlEstAge,coef="Est.leng",n=100,sort.by="p")[100,4]
#n100equivGSE <- max(which(topTable(efitBlEstAgeGSE,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
#BlEstAgeGSETop100equiv <- topTable(efitBlEstAgeGSE,coef="Est.leng",n=n100equivGSE,sort.by="p")
#n100equivDSE100 <- max(which(topTable(efitBlEstAgeDSE100,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
#BlEstAgeDSE100Top100equiv <- topTable(efitBlEstAgeDSE100,coef="Est.leng",n=n100equivDSE100,sort.by="p")
#n100equivDSE800 <- max(which(topTable(efitBlEstAgeDSE800,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
#BlEstAgeDSE800Top100equiv <- topTable(efitBlEstAgeDSE800,coef="Est.leng",n=n100equivDSE800,sort.by="p")


v2wBlEstAgeNoISEvGSE <- venn2way(rownames(BlEstAgeTop100),rownames(BlEstAgeGSETop100equiv),names=c("NoISE","GSE"),
         main="2-way Venn Diagram of genes with pv < NoISE Top100 pv\n from Limma analysis of Blood samples",
         sub=paste("pv cutpoint is ",signif(pvcut100,2),".  Model is Est Length + Age, coefficient of interest is Est Length",sep=""))

rng <- range(c(BlEstAgeTop100[,"P.Value"], BlEstAgeGSETop100equiv[,"P.Value"]))

pdf("LinePlots/BloodEstAge_Est_NoISEvsGSE_lineplot.pdf")
plot(0:1,rng,type="n",log="y",ylab="log10 p-value",xlab="",xaxt="n", main="Blood Est Age, Est coef top100 genes' p-values from NoISE and GSE")
segments(rep(0,100),BlEstAgeTop100[v2wBlEstAgeNoISEvGSE$AB,"P.Value"],rep(1,100),BlEstAgeGSETop100equiv[v2wBlEstAgeNoISEvGSE$AB,"P.Value"])
points(rep(0,length(v2wBlEstAgeNoISEvGSE$A.only)),
       BlEstAgeTop100[v2wBlEstAgeNoISEvGSE$A.only,"P.Value"],col="red",pch=16)
points(rep(1,length(v2wBlEstAgeNoISEvGSE$B.only)),
       BlEstAgeGSETop100equiv[v2wBlEstAgeNoISEvGSE$B.only,"P.Value"],col="blue",pch=16)
mtext(c("NoISE","GSE"),1,line=1,at=0:1)
dev.off()


v2wBlEstAgeNoISEvDSE800 <- venn2way(rownames(BlEstAgeTop100),rownames(BlEstAgeDSE800Top100equiv),names=c("NoISE","DSE800"),
         main="2-way Venn Diagram of genes with pv < NoISE Top100 pv\n from Limma analysis of Blood samples",
         sub=paste("pv cutpoint is ",signif(pvcut100,2),".  Model is Est Length + Age, coefficient of interest is Est Length",sep=""))

rngDSE800 <- range(c(BlEstAgeTop100[,"P.Value"], BlEstAgeDSE800Top100equiv[,"P.Value"]))

pdf("LinePlots/BloodEstAge_Est_NoISEvsDSE800_lineplot.pdf")
plot(0:1,rngDSE800,type="n",log="y",ylab="log10 p-value",xlab="",xaxt="n", main="Blood Est Age, Est coef top100 genes' p-values\nfrom NoISE and DSE800")
segments(rep(0,100),BlEstAgeTop100[v2wBlEstAgeNoISEvDSE800$AB,"P.Value"],rep(1,100),BlEstAgeDSE800Top100equiv[v2wBlEstAgeNoISEvDSE800$AB,"P.Value"])
points(rep(0,length(v2wBlEstAgeNoISEvDSE800$A.only)),
       BlEstAgeTop100[v2wBlEstAgeNoISEvDSE800$A.only,"P.Value"],col="red",pch=16)
points(rep(1,length(v2wBlEstAgeNoISEvDSE800$B.only)),
       BlEstAgeDSE800Top100equiv[v2wBlEstAgeNoISEvDSE800$B.only,"P.Value"],col="blue",pch=16)
mtext(c("NoISE","DSE800"),1,line=1,at=0:1)
dev.off()



n100equivGSEDSE800 <- max(which(topTable(efitBlEstAgeGSEDSE800,coef="Est.leng",n=1000,sort.by="p")[1:1000,4] <= pvcut100))
BlEstAgeGSEDSE800Top100equiv <- topTable(efitBlEstAgeGSEDSE800,coef="Est.leng",n=n100equivGSEDSE800,sort.by="p")

v2wBlEstAgeNoISEvGSEDSE800 <- venn2way(rownames(BlEstAgeTop100),rownames(BlEstAgeGSEDSE800Top100equiv),names=c("NoISE","GSEDSE800"),
         main="2-way Venn Diagram of genes with pv < NoISE Top100 pv\n from Limma analysis of Blood samples",
         sub=paste("pv cutpoint is ",signif(pvcut100,2),".  Model is Est Length + Age, coefficient of interest is Est Length",sep=""))

rngGSEDSE800 <- range(c(BlEstAgeTop100[,"P.Value"], BlEstAgeGSEDSE800Top100equiv[,"P.Value"]))

pdf("LinePlots/BloodEstAge_Est_NoISEvsGSEDSE800_lineplot.pdf")
plot(0:1,rngGSEDSE800,type="n",log="y",ylab="log10 p-value",xlab="",xaxt="n", main="Blood Est Age, Est coef top100 genes' p-values\nfrom NoISE and GSEDSE800")
segments(rep(0,100),BlEstAgeTop100[v2wBlEstAgeNoISEvGSEDSE800$AB,"P.Value"],rep(1,100),BlEstAgeGSEDSE800Top100equiv[v2wBlEstAgeNoISEvGSEDSE800$AB,"P.Value"])
points(rep(0,length(v2wBlEstAgeNoISEvGSEDSE800$A.only)),
       BlEstAgeTop100[v2wBlEstAgeNoISEvGSEDSE800$A.only,"P.Value"],col="red",pch=16)
points(rep(1,length(v2wBlEstAgeNoISEvGSEDSE800$B.only)),
       BlEstAgeGSEDSE800Top100equiv[v2wBlEstAgeNoISEvGSEDSE800$B.only,"P.Value"],col="blue",pch=16)
mtext(c("NoISE","GSEDSE800"),1,line=1,at=0:1)
dev.off()



#######################
# 19 March 2015

topTable(efitBlEst,coef="Est.leng")

designBlEst <- model.matrix(~alldat.inf[alldat.inf[,"sampletype"]=="blood","Est.leng"],data=log2Bnorm)
colnames(designBlEst)[2] <- "Est.leng"
awBlEst <- arrayWeights(log2Bnorm,designBlEst)
fitBlEst <- lmFit(log2Bnorm, designBlEst, weights=awBlEst)

efitBlEst <- eBayes(fitBlEst)
topTable(efitBlEst,coef="Est.leng")
#ltoptabout(topTable(efitBlEst, coef="Est.leng",number=1000),nam="Blood")
# not working because of ltoptabout isn't run
BlEst.Limma <- topTable(efitBlEst, coef="Est.leng",number=nrow(log2Bnorm))


#plot(alldat.inf[alldat.inf[,"sampletype"]=="blood","Est.leng"] , log2Bnorm["239362_at",],type="n")
#text(alldat.inf[alldat.inf[,"sampletype"]=="blood","Est.leng"] , log2Bnorm["239362_at",],1:ncol(log2Bnorm))
# 27th row/patient is good to choose
pdf("Residuals/BloodEst_239362_at_NAPA_RegressionPlotWithResidual.pdf")
plot(alldat.inf[alldat.inf[,"sampletype"]=="blood","Est.leng"] , log2Bnorm["239362_at",],pch=rep(c(1,16,1),c(26,1,8)),ylab="log2 Expression for NAPA gene",xlab="Estimated Allele Length")
cfs <- fitBlEst["239362_at",]$coef
lines(par("usr")[1:2], cfs[1] + par("usr")[1:2] * cfs[2])
xpp <- alldat.inf[alldat.inf[,"sampletype"]=="blood","Est.leng"][27]
lines(rep(xpp,2), c(log2Bnorm["239362_at",27], cfs[1] + xpp * cfs[2]),lty=2)
dev.off()

# NAPA   N-ethylmaleimide-sensitive factor attachment protein, alpha

rngtt <- range(ttresBlEstAge[,3:4])
pdf("Residuals/BloodEstAge_IndividualsMeanResidualsWithCI.pdf")
plot(ttresBlEstAge[,"mean"], 1:nrow(ttresBlEstAge),pch=16,xlim=c(rngtt[1],rngtt[2]),ylab="Patient",xlab="Average Residual",yaxt="n")
segments(par("usr")[1],1:nrow(ttresBlEstAge) ,par("usr")[2], 1:nrow(ttresBlEstAge),lty=3,col="grey")
segments(ttresBlEstAge[,"CI.low"],1:nrow(ttresBlEstAge),ttresBlEstAge[,"CI.high"],1:nrow(ttresBlEstAge),lwd=2)
mtext(1:nrow(ttresBlEstAge),2,0.5,at=1:nrow(ttresBlEstAge),cex=0.7,las=1)
lines(c(0,0),par("usr")[3:4],lty=3)
dev.off()





###########################
# leftovers

awBlEstAgeGSE <- arrayWeights(log2BnormEstAgeGSE[,ageBDpres],designBlEstAge)

fitBlEstAgeGSE <- lmFit(log2BnormEstAgeGSE[,ageBDpres], designBlEstAge, weights=awBlEstAgeGSE)

efitBlEstAgeGSE <- eBayes(fitBlEstAgeGSE)
topTable(efitBlEstAgeGSE,coef="Est.leng")
#ltoptabout(topTable(efitBlEstAge, coef="Est.leng",number=1000),nam="Blood.EstL.Age")
#ltoptableout again
BlEstAgeGSE.Limma <- topTable(efitBlEstAgeGSE, coef="Est.leng",number=nrow(log2Bnorm))

library(affy)
library(limma)

targetsFile <- "estrogen/estrogen.txt"
targetsFile

pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)
pData(pd)

raw <-ReadAffy(celfile.path = "estrogen", filenames=rownames(pData(pd)),phenoData = pd)
raw


boxplot(raw,col="red",las=2)

par(mfrow=c(2,1))
hist(log2(pm(raw[,1])),breaks=100,col="steelblue",main="PM",xlim=c(4,14))
hist(log2(mm(raw[,1])),breaks=100,col="steelblue",main="MM",xlim=c(4,14))

mva.pairs(pm(raw)[,1:4], plot.method="smoothScatter")
mva.pairs(pm(raw)[,5:8], plot.method="smoothScatter")

#Normalised Unscaled Standard Error (NUSE)
library(affyPLM)
plmset <- fitPLM(raw)
NUSE(plmset,las=2)
RLE(plmset,las=2)

#QC plots
library(arrayQualityMetrics)
arrayQualityMetrics(eset)

#Design matrix
ER <- pData(pd)$estrogen
Time <- pData(pd)$time.h
design <- model.matrix(~ER + Time)
design2 <- model.matrix(~ER*Time)

#Read and normalize estrogen data
eset <- rma(raw)
head(exprs(eset))
#Fit the model
fit1 <- lmFit(eset, design)
fit1 <- eBayes(fit1)
topTable(fit1, coef = 2)

fit2 <- lmFit(eset, design2)
fit2 <- eBayes(fit2)
topTable(fit2, coef = 2)
head(fit2$coefficients)

#Decide DE 
decideTests(fit1)
table(decideTests(fit1))
sum(abs(decideTests(fit1))==1)
volcanoplot(fit1, highlight = 10, names = rownames(fit1$coefficients))

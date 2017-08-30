#set the working directory
setwd("C:/Users/cherukuri/Desktop/affymetrics/")
#ftp link to download R-packages
source("http://bioconductor.org/biocLite.R")
#download GEOquery package to retrive datasets from geo database
biocLite("GEOquery")
#Install required packages
biocLite("affy")
biocLite("limma")
#call the installed packages
library(affy)
library(GEOquery)
#download the completedataset of the study from GEO database.
getGEOSuppFiles("GSE47516")
###Initialize the path, where your CEL files are located### 
##make sures the CEL files are unzipped##
setwd("C:/Users/cherukuri/Desktop/affymetrics/GSE47516/")
##Read the CEL files##
affy.data = ReadAffy()
##normalize the expression values ## 
eset.mas5 = mas5(affy.data)
##get the expression values##
exprSet.nologs = exprs(eset.mas5)
##list the column names ##
exprSet_colnames= colnames(exprSet.nologs)
##log2 transformation ## ## or can use RMA for normalization which outputs log2 transformation values ##
exprSet = log(exprSet.nologs, 2)
## save the expression matrix to an outputfile ##
write.table(exprSet,file = "C:/Users/cherukuri/Desktop/affymetrics/expression_matrix.txt",sep = " ", quote =F, col.names=T)
write.table(exprSet_colnames,file = "C:/Users/cherukuri/Desktop/affymetrics/pheno_data.txt",sep = " ", quote =F, col.names=T)        
########################## DIFFERENTIAL EXPRESSIONS ##########################################
setwd("C:/Users/cherukuri/Desktop/affymetrics/")
##load package ##
library(limma)
## read phenotypefile ## ## defines the categories in the data ##
pheno =read.table("pheno_matrix.txt",row.names=1)
## load the expression file ##
raw.data <- read.delim("expression_data.txt")
## prepare a design matrix ##
Group<-factor(pheno$V2,levels=levels(pheno$V2))
data = raw.data[,2:22]
rownames(data) <- raw.data[, 1]
design<-model.matrix(~0+Group)
colnames(design)<-c("mutant", "wild")
## fit the linear model ##
fit = lmFit(data,design) 
## perform the statistical Analysis ##
fit = eBayes(fit)
results <- decideTests(fit, p.value=0.05, method="global",
adjust.method="BH")
## report the top 20 Differentially Expressed genes ##
x<-topTable(fit,number=20,coef =2)
## write the Differentially Expressed genes to the output file ## 
write.table(x,file="DE_result_top20.txt",sep="\t")
write.table(results,file="DE_result_pval_0.05_cutoff.txt",sep="\t")



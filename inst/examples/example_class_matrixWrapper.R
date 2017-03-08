#####
# Fake AUC, normally this object is returned by calcAUC()
aucMatrix <- matrix(sample(1/1:30, 120, replace=TRUE), nrow=8)
 rownames(aucMatrix) <- paste("geneSet",1:nrow(aucMatrix), sep="")
 colnames(aucMatrix) <- paste("cell",1:ncol(aucMatrix), sep="")
fakeAUC <- matrixWrapper(matrix=aucMatrix, rowType="gene-set", colType="cell", matrixType="AUC")
#####

fakeAUC

# To subset & access the AUC or ranking slot:
fakeAUC[1:2,]
fakeAUC[,3:4]
getAuc(fakeAUC)[,1:5] # fakeAUC@matrix[,1:5]

# To subset the AUC or ranking, but keeping the class:
subset(fakeAUC, 1:2) # Default subsets on rows
subset(fakeAUC, 3:4, "cell")

# These methods are also available for this class:
dim(fakeAUC)
nrow(fakeAUC)
ncol(fakeAUC)
colnames(fakeAUC)
rownames(fakeAUC)


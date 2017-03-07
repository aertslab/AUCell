
# write.gmt.single (one gene set)
aGeneSet <- c("gene1", "gene2", "gene3")
write.gmt.single(aGeneSet, setName="geneSet1", fileName="aMiniGeneSet.gmt")


# write.gmt (multiple gene sets)
anotherGeneSet <- paste("gene", sample(1:100, 10), sep="")
geneSets <- list(geneSet1=aGeneSet,
                 geneSet2=anotherGeneSet)
geneSets
write.gmt(geneSets, fileName="bothGeneSets.gmt", gmtDir=".")


# read.gmt
read.gmt("bothGeneSets.gmt")


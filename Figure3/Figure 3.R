#R script for plots in Figure 3
#load libraries
library(WGCNA)
library(ggplot2)
library(reshape)

#DIA protein abundance data - see Figure1 directory for file
egg.dat<-read.csv('dia protein abundances revised.csv', header=T, row.names=1)

#log transform data
df = egg.dat
df.log=log2(df)

dat.t<-t(df.log)
rownames(dat.t) <- gsub(".Protein.Abundance","", rownames(dat.t))

sampleTree = hclust(dist(dat.t), method='average')
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree)

SampleData<-read.csv('coral.eggs.metadata2.csv')

AllSamples<-rownames(dat.t)
SampleRows <- match(AllSamples, SampleData$file.name)
datSamples <- SampleData[SampleRows, -1]
rownames(datSamples)<-SampleData[SampleRows, 1]

collectGarbage()
allowWGCNAThreads()
#for RStudio: allwWGCNAThreads; for R enableWGCNAThreads

powers<-c(c(1:10), seq(from=12,to=20, by=2))
sft<-pickSoftThreshold(dat.t, powerVector=powers, verbose=5)

sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab='Soft Threshold (power)', ylab='Scale Free Topology Model Fit, signed R^2', type='n', main=paste('Scale Independence'))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col='red')
abline(h=0.9, col='red')
#soft power of 16
plot(sft$fitIndices[,1],sft$fitIndices[,5], xlab='Soft Threshold (power)', ylab='Mean Connectivity', type='n', main=paste('Mean Connectivity'))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers, cex=cex1, col='red')

bwnet = blockwiseModules(dat.t, power=16, TOMType='signed', minModuleSize=30, reassignThreshold=0, mergeCutHeight=0.25, numericLabels=T, saveTOMs=T, saveTOMFileBase='Mcap-eggs-TOM', verbose=3)

softPower=16
adjacency<-adjacency(dat.t, power=softPower, type='signed')
TOM = TOMsimilarity(adjacency, TOMType="signed")
save(TOM, file="TOM.MCAPEggs.RData", compress = F)

moduleLabels = bwnet$colors
bwLabels = matchLabels(bwnet$colors, moduleLabels)
bwModuleColors = labels2colors(bwLabels)

#relating modules to clinical traits
nGenes = ncol(dat.t)
nSamples=nrow(dat.t)

MEs0 = moduleEigengenes(dat.t, bwModuleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datSamples, use='p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sept="")
dim(textMatrix) = dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
labeledHeatmap(Matrix=moduleTraitCor, xLabels=names(datSamples), yLabels=names(MEs), ySymbols=names(MEs), colorLabels=F, colors=blueWhiteRed(50), textMatrix=textMatrix, setStdMargins=F, cex.text=0.5, zlim=c(-1,1), main=paste("Module-trait relationships"))

####Look more closely at variables with significant correlation to modules
#define GS (gene significance) and MM (module membership) for Genet 13
Gen13 = as.data.frame(datSamples$Genet13)
names(Gen13) = 'Gen13'

modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(dat.t, MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="s")

geneTraitSignificance = as.data.frame(cor(dat.t, Gen13, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Gen13), sep="")
names(GSPvalue) = paste("p.GS.", names(Gen13), sep="")

#Genet26
Gen26 = as.data.frame(datSamples$Genet26)
names(Gen26) = 'Gen26'

geneTraitSignificance2 = as.data.frame(cor(dat.t, Gen26, use="p"))
GSPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance2), nSamples))

names(geneTraitSignificance2) = paste("GS.", names(Gen26), sep="")
names(GSPvalue2) = paste("p.GS.", names(Gen26), sep="")

#define GS and MM for bleaching
bleach = as.data.frame(datSamples$Bleaching)
names(bleach) = 'bleach'

geneTraitSignificance3 = as.data.frame(cor(dat.t, bleach, use='p'))
GSPvalue3 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance3), nSamples))

names(geneTraitSignificance3) = paste("GS.", names(bleach), sep="")
names(GSPvalue3) = paste("p.GS.", names(bleach), sep="")

#create export files
proteins = colnames(dat.t)
protInfo0 = data.frame(clusters = proteins, moduleColor = bwModuleColors, geneTraitSignificance, geneTraitSignificance2, geneTraitSignificance3, GSPvalue, GSPvalue2, GSPvalue3)

modOrder = order(-abs(cor(MEs, Gen13, use="p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(protInfo0)
  protInfo0 = data.frame(protInfo0, geneModuleMembership[,modOrder[mod]], MMPvalue[, modOrder[mod]]);
  names(protInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

protOrder = order (protInfo0$moduleColor, -abs(protInfo0$GS.Gen13));
protInfo = protInfo0[protOrder, ]

write.csv(protInfo, file='protInfoGen13.csv')

modOrder2 = order(-abs(cor(MEs, Gen26, use="p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(protInfo0)
  protInfo0.2 = data.frame(protInfo0, geneModuleMembership[,modOrder2[mod]], MMPvalue[, modOrder2[mod]]);
  names(protInfo0.2) = c(oldNames, paste("MM.", modNames[modOrder2[mod]], sep=""), paste("p.MM.", modNames[modOrder2[mod]], sep=""))
}

protOrder2 = order (protInfo0.2$moduleColor, -abs(protInfo0.2$GS.Gen26));
protInfo2 = protInfo0.2[protOrder2, ]

write.csv(protInfo2, file='protInfoGen26.csv')

modOrder3 = order(-abs(cor(MEs, bleach, use="p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(protInfo0)
  protInfo0.3 = data.frame(protInfo0, geneModuleMembership[,modOrder3[mod]], MMPvalue[, modOrder3[mod]]);
  names(protInfo0.3) = c(oldNames, paste("MM.", modNames[modOrder3[mod]], sep=""), paste("p.MM.", modNames[modOrder3[mod]], sep=""))
}

protOrder3 = order (protInfo0.3$moduleColor, -abs(protInfo0.3$GS.bleach));
protInfo3 = protInfo0.3[protOrder3, ]

write.csv(protInfo3, file='protInfobleach.csv')

#investigating gene significance and module membership for highly significant modules
###blue
module = 'blue'
column = match(module, modNames)
moduleProteins = bwModuleColors == module

sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleProteins, column]),
                   abs(geneTraitSignificance2[moduleProteins,1]),
                   xlab=paste("Module Membership in", module, "module"),
                   ylab="Gene Significance for Bleaching",
                   main = paste("Module membership vs. genen significant\n"),
                   cex.main = 1.2, cex.lab=1.2, cex.axis=1.2, col=module)


###Turquoise
module = 'turquoise'
column = match(module, modNames)
moduleProteins = bwModuleColors == module
sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleProteins, column]),
                   abs(geneTraitSignificance2[moduleProteins,1]),
                   xlab=paste("Module Membership in", module, "module"),
                   ylab="Gene Significance for Bleaching",
                   main = paste("Module membership vs. genen significant\n"),
                   cex.main = 1.2, cex.lab=1.2, cex.axis=1.2, col=module)

###Grey
module = 'grey'
column = match(module, modNames)
moduleProteins = bwModuleColors == module
sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleProteins, column]),
                   abs(geneTraitSignificance[moduleProteins,1]),
                   xlab=paste("Module Membership in", module, "module"),
                   ylab="Gene Significance for Spawn Date",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab=1.2, cex.axis=1.2, col=module)

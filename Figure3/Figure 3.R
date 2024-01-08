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

SampleData<-read.csv('coral.eggs.metadata.csv')

AllSamples<-rownames(dat.t)
SampleRows <- match(AllSamples, SampleData$file.name)
datSamples <- SampleData[SampleRows, -1]
rownames(datSamples)<-SampleData[SampleRows, 1]

collectGarbage()
allowWGCNAThreads()
#for RStudio: allwWGCNAThreads; for R enableWGCNAThreads

SampleData$Bleaching<-as.numeric(SampleData$Bleaching)
sampleTree2<-hclust(dist(dat.t), method='average')
traitColors = numbers2colors(SampleData[,2:4], signed=F)
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(SampleData[,2:4]))

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

#define GS (gene significance) and MM (module membership) for spawn date
spawn = as.data.frame(datSamples$Spawn.Date)
names(spawn) = 'spawn'

modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(dat.t, MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="s")

geneTraitSignificance = as.data.frame(cor(dat.t, spawn, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(spawn), sep="")
names(GSPvalue) = paste("p.GS.", names(spawn), sep="")

#define GS and MM for bleaching
bleach = as.data.frame(datSamples$Bleaching)
names(bleach) = 'bleach'

geneTraitSignificance2 = as.data.frame(cor(dat.t, bleach, use='p'))
GSPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance2), nSamples))

names(geneTraitSignificance2) = paste("GS.", names(bleach), sep="")
names(GSPvalue2) = paste("p.GS.", names(bleach), sep="")

#create export files
proteins = colnames(dat.t)
protInfo0 = data.frame(clusters = proteins, moduleColor = bwModuleColors, geneTraitSignificance, geneTraitSignificance2, GSPvalue, GSPvalue2)

modOrder = order(-abs(cor(MEs, spawn, use="p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(protInfo0)
  protInfo0 = data.frame(protInfo0, geneModuleMembership[,modOrder[mod]], MMPvalue[, modOrder[mod]]);
  names(protInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

protOrder = order (protInfo0$moduleColor, -abs(protInfo0$GS.spawn));
protInfo = protInfo0[protOrder, ]

write.csv(protInfo, file='protInfoSpawn.csv')

modOrder2 = order(-abs(cor(MEs, bleach, use="p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(protInfo0)
  protInfo0.2 = data.frame(protInfo0, geneModuleMembership[,modOrder2[mod]], MMPvalue[, modOrder2[mod]]);
  names(protInfo0.2) = c(oldNames, paste("MM.", modNames[modOrder2[mod]], sep=""), paste("p.MM.", modNames[modOrder2[mod]], sep=""))
}

protOrder2 = order (protInfo0.2$moduleColor, -abs(protInfo0.2$GS.bleach));
protInfo2 = protInfo0.2[protOrder2, ]

write.csv(protInfo2, file='protInfoBleaching.csv')

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

#Boxplots of module proteins associated with significantly enriched GO terms for each module
mod.mem<-read.csv('module membership coral proteins only.csv')

#Significant by bleaching: blue and turquoise modules
#Blue Module
colnames(egg.dat)
bleached<-as.data.frame(cbind(egg.dat$X10B_32.Protein.Abundance, egg.dat$X13B_46.Protein.Abundance, egg.dat$X21B_43.Protein.Abundance, egg.dat$X26B_28.Protein.Abundance, egg.dat$X66B_44.Protein.Abundance, egg.dat$X74B_31.Protein.Abundance))
nonbleached<-as.data.frame(cbind(egg.dat$X10NB_47.Protein.Abundance, egg.dat$X13NB_33.Protein.Abundance, egg.dat$X21NB_30.Protein.Abundance, egg.dat$X26NB_45.Protein.Abundance, egg.dat$X66NB_29.Protein.Abundance, egg.dat$X74NB_48.Protein.Abundance))
bleached$avg.bleached<-rowMeans(data.frame(bleached))
nonbleached$avg.nonbleached<-rowMeans(data.frame(nonbleached))

rownames(bleached)<-rownames(egg.dat)
colnames(bleached)<-c('B.10_32', 'B.13_46', 'B.21_43', 'B.26_28', 'B.66_44', 'B.74_31', 'B.avg')

rownames(nonbleached)<-rownames(egg.dat)
colnames(nonbleached)<-c('NB10_47', 'NB13_33', 'NB21_30', 'NB26_45', 'NB66_29', 'NB74_48', 'NB.avg')

blue.all<-subset(mod.mem, moduleColor=='blue', select=c(Proteins))
merge.blue<-merge(x=blue.all, y=nonbleached, by.x='Proteins', by.y='row.names')
merge.blue2<-merge(x=merge.blue, y=bleached, by.x='Proteins', by.y='row.names')

blue.noavg<-merge.blue2[-c(8, 15)]

blue.melt<-melt(blue.noavg, id.vars=c('Proteins'))
blue.melt$Colony<-blue.melt$variable
blue.melt$Colony<-gsub("NB", "", blue.melt$Colony)
blue.melt$Colony<-gsub("B", "", blue.melt$Colony)
blue.melt$Colony<-gsub("\\.", "", blue.melt$Colony)

blue.melt$variable<-gsub("_.*", "", blue.melt$variable)
blue.melt$variable<-gsub("\\.", "", blue.melt$variable)
blue.melt$variable<-gsub("10", "", blue.melt$variable)
blue.melt$variable<-gsub("13", "", blue.melt$variable)
blue.melt$variable<-gsub("21", "", blue.melt$variable)
blue.melt$variable<-gsub("26", "", blue.melt$variable)
blue.melt$variable<-gsub("66", "", blue.melt$variable)
blue.melt$variable<-gsub("74", "", blue.melt$variable)

ggplot(blue.melt) +
  geom_line(aes(x=variable, y=log(value), group=Proteins), alpha=0.3, color='blue') +
  geom_boxplot(aes(x=variable, y=log(value)), color='blue') +
  theme_bw() +
  scale_x_discrete('Colony Bleaching Status', limits=c('NB', 'B')) +
  theme(legend.position = "none") +
  labs(x="Parent Colony Bleaching History", y='Log(Avg. Area Under the Curve)')

##Turquoise
turq.all<-subset(mod.mem, moduleColor=='turquoise', select=c(Proteins))
merge.turq<-merge(x=turq.all, y=nonbleached, by.x='Proteins', by.y='row.names')
merge.turq2<-merge(x=merge.turq, y=bleached, by.x='Proteins', by.y='row.names')

turq.noavg<-merge.turq2[-c(8, 15)]

turq.melt<-melt(turq.noavg, id.vars=c('Proteins'))
turq.melt$Colony<-turq.melt$variable
turq.melt$Colony<-gsub("NB", "", turq.melt$Colony)
turq.melt$Colony<-gsub("B", "", turq.melt$Colony)
turq.melt$Colony<-gsub("\\.", "", turq.melt$Colony)

turq.melt$variable<-gsub("_.*", "", turq.melt$variable)
turq.melt$variable<-gsub("\\.", "", turq.melt$variable)
turq.melt$variable<-gsub("10", "", turq.melt$variable)
turq.melt$variable<-gsub("13", "", turq.melt$variable)
turq.melt$variable<-gsub("21", "", turq.melt$variable)
turq.melt$variable<-gsub("26", "", turq.melt$variable)
turq.melt$variable<-gsub("66", "", turq.melt$variable)
turq.melt$variable<-gsub("74", "", turq.melt$variable)

ggplot(turq.melt) +
  geom_line(aes(x=variable, y=log(value), group=Proteins), alpha=0.3, color='turquoise') +
  geom_boxplot(aes(x=variable, y=log(value)), color='turquoise') +
  theme_bw() +
  scale_x_discrete('Colony Bleaching Status', limits=c('NB', 'B')) +
  theme(legend.position = "none") +
  labs(x="Parent Colony Bleaching History", y='Log(Avg. Area Under the Curve)')

#Significant by spawn date: grey module
##Grey
grey.all<-subset(mod.mem, moduleColor=='grey', select=c(Proteins))

June13<-as.data.frame(cbind(egg.dat$X10NB_47.Protein.Abundance, egg.dat$X26NB_45.Protein.Abundance))
June14<-as.data.frame(cbind(egg.dat$X13NB_33.Protein.Abundance, egg.dat$X13B_46.Protein.Abundance))
June16<-as.data.frame(cbind(egg.dat$X74NB_48.Protein.Abundance))
July10<-as.data.frame(cbind(egg.dat$X66NB_29.Protein.Abundance))
July12<-as.data.frame(cbind(egg.dat$X21NB_30.Protein.Abundance, egg.dat$X10B_32.Protein.Abundance, egg.dat$X26B_28.Protein.Abundance))
July13<-as.data.frame(cbind(egg.dat$X21B_43.Protein.Abundance, egg.dat$X74B_31.Protein.Abundance))
July14<-as.data.frame(cbind(egg.dat$X66B_44.Protein.Abundance))

rownames(June13)<-rownames(egg.dat)
colnames(June13)<-c('June13NB10_47', 'June13NB26_45')

rownames(June14)<-rownames(egg.dat)
colnames(June14)<-c('June14NB13_33', 'June14B13_46')

rownames(June16)<-rownames(egg.dat)
colnames(June16)<-c('June16NB74_48')

rownames(July10)<-rownames(egg.dat)
colnames(July10)<-c('July10NB66_29')

rownames(July12)<-rownames(egg.dat)
colnames(July12)<-c('July12NB21_30', 'July12B10_32', 'July12B26_28')

rownames(July13)<-rownames(egg.dat)
colnames(July13)<-c('July13B21_43', 'July13B74_31')

rownames(July14)<-rownames(egg.dat)
colnames(July14)<-c('July14B66_44')


grey.merge<-merge(x=grey.all, y=June13, by.x='Proteins', by.y='row.names')
grey.merge2<-merge(x=grey.merge, y=June14, by.x='Proteins', by.y='row.names', all.x=T)
grey.merge3<-merge(x=grey.merge2, y=June16, by.x='Proteins', by.y='row.names', all.x=T)
grey.merge4<-merge(x=grey.merge3, y=July10, by.x='Proteins', by.y='row.names', all.x=T)
grey.merge5<-merge(x=grey.merge4, y=July12, by.x='Proteins', by.y='row.names', all.x=T)
grey.merge6<-merge(x=grey.merge5, y=July13, by.x='Proteins', by.y='row.names', all.x=T)
grey.merge7<-merge(x=grey.merge6, y=July14, by.x='Proteins', by.y='row.names', all.x=T)


grey.melt<-melt(grey.merge7, id.vars=c('Proteins'))

grey.melt$Colony<-grey.melt$variable
grey.melt$Colony<-gsub(".*B", "", grey.melt$Colony)
grey.melt$Colony<-gsub("_.*", "", grey.melt$Colony)

grey.melt$variable<-gsub("B.*", "", grey.melt$variable)
grey.melt$variable<-gsub("N", "", grey.melt$variable)

ggplot(grey.melt) +
  geom_line(aes(x=variable, y=log(value), group=Proteins), alpha=0.3, color='grey') +
  geom_boxplot(aes(x=variable, y=log(value)), color='grey') +
  theme_bw() +
  scale_x_discrete('Spawn Date', limits=c('June13', 'June14', 'June16', 'July10', 'July12', 'July13', 'July14')) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x="Spawn Date", y='Log(Avg. Area Under the Curve)')


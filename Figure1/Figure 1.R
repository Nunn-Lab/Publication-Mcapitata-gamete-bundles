#R script for comparison of workflows
#load libraries
library(reshape)
library(ggplot2)

#1A - DDA AUC vs. DDA spectral counts
#read in the spectral count data
dda.spn<-read.csv('dda spectral counts.csv')
#remove unnecessary column
dda.spn<-dda.spn[,2:14]

#read in the AUC data
dda.auc<-read.csv('dda protein abundances pivot.csv')

#rename to shorten column names
colnames(dda.spn)<-c('Protein', '10B', '10NB', '13B', '13NB', '21B', '21NB', '26B', '26NB', '66B','66NB', '74B', '74NB')
colnames(dda.auc)<-c('Protein', '10B', '10NB', '13B', '13NB', '21B', '21NB', '26B', '26NB', '66B','66NB', '74B', '74NB')

#reformat the data
spn.melt<-melt(dda.spn, id.vars='Protein')
auc.melt<-melt(dda.auc, id.vars='Protein')

#merge DDa AUC and spectral counts
spn.auc.merge<-merge(x=spn.melt, y=auc.melt, by=c('Protein', 'variable'))
spn.auc.merge$value.y<-as.numeric(spn.auc.merge$value.y)

ggplot(spn.auc.merge) +
  geom_point(aes(x=log(value.x), y=log(value.y)), alpha=0.5) +
  geom_smooth(aes(x=log(value.x), y=log(value.y)), method='lm', se=FALSE, color='dodgerblue', formula=y~x)+
  labs(x="Log(DDA Spectral Counts)", y="Log(DDA Area Under the Curve)") +
  theme_bw()

#1B - DIA AUC vs. DDA AUC
#read in the DIA dataset and rename columns
dia.dat<-read.csv('dia protein abundances revised.csv')
colnames(dia.dat)<-c('Protein', '10B', '10NB', '13B', '13NB', '21B', '21NB', '26B', '26NB', '66B', '66NB', '74B', '74NB')

#reformat
melt.dia<-melt(dia.dat, id.vars="Protein")

#mrege with DDA AUC
dda.dia<-merge(x=melt.dia, y=auc.melt, by=c('Protein', 'variable'))
dda.dia$value.y<-as.numeric(dda.dia$value.y)

ggplot(dda.dia) +
  geom_point(aes(x=log(value.y), y=log(value.x)), alpha=0.5) +
  geom_smooth(aes(x=log(value.y), y=log(value.x)), method='lm', se=FALSE, color='dodgerblue', formula=y~x)+
  labs(x="Log(DDA Area Under the Curve)", y="Log(DIA Area Under the Curve)") +
  theme_bw()

#1C - DIA AUC vs. DDA spectral counts
#merge datasets
dda.dia2<-merge(x=melt.dia, y=spn.melt, by=c('Protein', 'variable'))

#reformat
dda.dia2.refmt<-melt(dda.dia2, id.vars=c('Protein', 'variable'))
colnames(dda.dia2.refmt)<-c('Protein', 'Colony', 'Acquisition', 'Abundance')
dda.dia2.refmt$Abundance<-as.numeric(dda.dia2.refmt$Abundance)

ggplot(dda.dia2) +
  geom_point(aes(x=value.y, y=log(value.x)), alpha=0.5) +
  geom_smooth(aes(x=value.y, y=log(value.x)), method='lm', se=FALSE, color='dodgerblue', formula=y~x)+
  labs(x="DDA Spectral Counts", y="Log(DIA Area Under the Curve)") +
  theme_bw()

#1D - Protein counts Venn diagram
#subset proteins from colony 10, bleached, from each dataset
DDA.spn.10B<-subset(spn.melt, variable=='10B')
DDA.spn.10B<-subset(DDA.spn.10B, value>0)
DDA.auc.10B<-subset(auc.melt, variable=='10B')
DDA.auc.10B<-subset(DDA.auc.10B, value>0)
DIA.10B<-subset(melt.dia, variable=='10B')
DIA.10B<-subset(DIA.10B, value>0)

#create a list for each data type
list.B<-list('DDA Spec.Counts'=DDA.spn.10B$Protein, 'DDA AUC'=DDA.auc.10B$Protein, 'DIA AUC'=DIA.10B$Protein)

venn.pl<-ggVennDiagram(list.B, label_alpha = 0, ) +
  scale_fill_distiller(palette = "BuGn", direction = 1)

#1E - Density plot of DIA and DDA
#reformat data
melt.both<-melt(dda.dia, id.vars=c('Protein', 'variable'))
colnames(melt.both)<-c('Protein', 'Colony', 'Acquisition', 'Abundance')
melt.both$Abundance<-as.numeric(melt.both$Abundance)

ggplot(melt.both, aes(x=Abundance, fill=Acquisition)) +
  geom_density(alpha=0.7) +
  theme_bw() +
  xlim(0,5e7) +
  scale_fill_manual("Acquisition Type", values=c('red', 'blue'), labels=c('DIA', 'DDA'))

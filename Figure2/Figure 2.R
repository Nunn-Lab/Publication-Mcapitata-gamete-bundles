#R script for plots in Figure 2
#load libraries
source("Biostats.r")
library(dplyr)
library(vegan)

# NMDS : DDA - Spectral Library
#read in data
mcap.dat<-read.csv('dda protein abundances pivot.csv', header=T, row.names=1)
mcap.dat[, 1] <- sapply(mcap.dat[, 1], as.numeric)
mcap.dat <- na.omit(mcap.dat)
#transpose data
mcap.t<-as.data.frame(t(mcap.dat))

#log-transform data
mcap.tra<-data.trans((mcap.t+1), method='log', plot=F)

#run NMDS on a bray-curtis dissimilarity matrix
nmds.mcap<-metaMDS(mcap.tra, distance='bray', trymax=10, autotransform=F)

#make object grps that is a list of group assignments for each of sample
grps<-c("Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached")

#NMDS plot
fig.mcap<-ordiplot(nmds.mcap, choices=c(1,2), type='none', display='sites')
points(fig.mcap, 'sites', pch=c(rep(21,9)) , bg=c("gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue") , col='grey50', cex=1.5)
ordihull(fig.mcap, groups=grps, draw='lines',col='grey75', label=T)

#NMDS : DIA - Spectral Library

#read in data
mcap.dat<-read.csv('dia protein abundances revised.csv', header=T, row.names=1)
mcap.dat[, 1] <- sapply(mcap.dat[, 1], as.numeric)
mcap.dat <- na.omit(mcap.dat)
#transpose data
mcap.t<-as.data.frame(t(mcap.dat))

#log-transform data
mcap.tra<-data.trans((mcap.t+1), method='log', plot=F)

#run NMDS on a bray-curtis dissimilarity matrix
nmds.mcap<-metaMDS(mcap.tra, distance='bray', trymax=10, autotransform=F)

# groups are the same as before
grps<-c("Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached")

#make and export NMDS plot
fig.mcap<-ordiplot(nmds.mcap, choices=c(1,2), type='none', display='sites')
points(fig.mcap, 'sites', pch=c(rep(21,9)) , bg=c("gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue") , col='grey50', cex=1.5)
ordihull(fig.mcap, groups=grps, draw='lines',col='grey75', label=T)

# NMDS: DDA - Database Search

#read in data
mcap.dat<-read.csv('combined db file adjnsaf only.csv', header=T, row.names=1)
mcap.dat[, 1] <- sapply(mcap.dat[, 1], as.numeric)
mcap.dat <- na.omit(mcap.dat)
#transpose data
mcap.t<-as.data.frame(t(mcap.dat))

#log-transform data
mcap.tra<-data.trans((mcap.t+1), method='log', plot=F)

#run NMDS on a bray-curtis dissimilarity matrix
nmds.mcap<-metaMDS(mcap.tra, distance='bray', trymax=10, autotransform=F)

# once again same groups
grps<-c("Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached", "Bleached", "Not Bleached")

#make and export NMDS plot
fig.mcap<-ordiplot(nmds.mcap, choices=c(1,2), type='none', display='sites')
points(fig.mcap, 'sites', pch=c(rep(21,9)) , bg=c("gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue", "gold", "blue") , col='grey50', cex=1.5)
ordihull(fig.mcap, groups=grps, draw='lines',col='grey75', label=T)

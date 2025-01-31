#R script for volcano plot (adult tissue vs. gamete bundle)
#load libraries
library(ggplot2)

#read in data
q.dat<-read.csv('diff abund with annotations.csv')

#assigned broad category colors to GO terms
GO.colors<-c('antioxidant' = 'red', 'chromatin condensation' = 'orange', 'chromatin structure' = 'orange', 'cilia & flagella cytoskeleton' = 'yellow', 'cilia cytoskeleton' = 'yellow', 'cilium' = 'yellow', 'cleaves C-O bonds' = 'green', 'cytoskeleton' = 'yellow', 'ER membrane' = 'blue', 'iron storage' = 'pink', 'post-translational modification, O-glycosylation' = 'purple', 'protein cross-linking' = 'turquoise', 'unknown' ='brown')

ggplot(dat=q.dat) +
  geom_point(aes(x=LogFoldChange, y=abs(Zstatistic), color=Function, size=abs(LogFoldChange)), alpha=0.5) +
  theme_bw() +
  geom_text(aes(x=LogFoldChange, y=abs(Zstatistic), label=Name), check_overlap=T, hjust='inward', vjust='inward')+
  scale_color_manual(values=GO.colors) +
  labs(x='Log Fold Change', y='|Z-statistic|') +
  geom_hline(yintercept=2, linetype='dashed') +
  geom_vline(xintercept=0.5, linetype='dashed') +
  geom_vline(xintercept=-0.5, linetype='dashed') +
  theme(legend.position="none")

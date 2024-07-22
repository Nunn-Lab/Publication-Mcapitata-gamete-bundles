### Calculating C to D symbiont ratios from qPCR Cq data and generating figures for Timmins-Schiffman et al. 
### author: Jenna Dilworth, CEE Lab, University of Southern California

# references: #
#Cunning et al. 2012 - Excess algal symbionts increase the susceptibility of reef corals to bleaching
#Cunning et al. 2017 - Patterns of bleaching and recovery of Montipora capitata in Kaneohe Bay, Hawaii, USA
#Innis et al. 2018- Coral color and depth drive symbiosis ecology of Montipora capitata in K ̄ ane‘ohe Bay, O‘ahu, Hawai‘i

#load required packages 
library(tidyverse)
library(readxl)
library(gtools)
library(lme4)
library(car)

options(scipen=999) # To avoid scientific notation

#read in and format raw qPCR data ####
rawdata<-read_excel("UW_allsamples.xlsx")%>%
  filter(well_type!= "NTC")

#set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
fluo.C<-0.55                                             #probe correction for C
fluo.D<-0                                                

data<-rawdata%>%  
  select(well_name, Target, Cq)%>%
  mutate_at("Cq", as.numeric)%>%
  na.replace(41)%>%
  group_by(well_name, Target)%>%                         ###group technical replicates
  summarise(meanCq = mean(Cq, na.rm = TRUE),             ###take mean of technical replicates
            stdev = sd(Cq, na.rm = TRUE), 
            .groups = "keep") 

C<-filter(data,Target=="VIC")
D<-filter(data,Target=="FAM")

#final formatting of data
final_data<-left_join(C,D,by="well_name")%>%
  select(-Target.x,-Target.y)%>%                                              
  rename(d_mean=meanCq.y)%>%                                           
  rename(d_sd=stdev.y)%>%
  rename(c_mean=meanCq.x)%>%
  rename(c_sd=stdev.x)%>%
  mutate(c_mean=c_mean-fluo.C)%>% 
  mutate_at("c_mean", round, 0)%>%
  mutate_at("c_mean", round, 0)%>%  
  mutate(presence=case_when(d_mean<36.9 && c_mean<36.9 ~"CD",   ###set presence/absence in new column based on ct values - Ct for a single molecule is 37
                            d_mean>36.9~"C", 
                            c_mean>36.9~"D"))%>%  
  mutate(cd_ratio=(2^(c_mean-d_mean))/copy.n.ratioCD)%>%        ###calculate ratio, log of 2^(difference in ct values) from Innis et al. 2018                                                             #from Cunning et al. 2012
  mutate(cd_ratio=case_when(d_mean>36.9 && c_mean>36.9 ~"0", 
                            d_mean>36.9 && c_mean<36.9~"Inf", 
                            c_mean>36.9 && d_mean<36.9~"-Inf",
                            TRUE ~ as.character(cd_ratio)))%>%   ###from Innis et al. 2018
  mutate(cd_ratio=as.numeric(cd_ratio))%>%
  mutate(prop_d=cd_ratio/(cd_ratio+1))%>%
  mutate(prop_d=case_when(cd_ratio == Inf ~"0", 
                          cd_ratio == -Inf ~ "1", 
                          cd_ratio == 0 ~"NA", 
                          TRUE ~ as.character(prop_d)))%>%       ###from Innis et al. 2018
  mutate(prop_d=as.numeric(prop_d))%>%
  mutate(prop_c=1-prop_d)%>%                                     ###set proportion C off proportion D
  mutate(prop_c=case_when(cd_ratio == -Inf ~"0", 
                          cd_ratio == Inf ~ "1", 
                          cd_ratio == 0 ~"NA",
                          TRUE ~ as.character(prop_c)))%>%       ###from Innis et al. 2018
  mutate(prop_c=as.numeric(prop_c))%>%
  mutate(abundance=case_when(prop_d>prop_c~"D>C",
                             prop_c>prop_d~"C>D"))%>%            ###set abundance based on which proportion is dominant, create character factor
  mutate(sd_warning=case_when((c_sd>1&d_sd>1)~"cd*",
                              c_sd>1~"c*",
                              d_sd>1~"d*"))%>%                   ###set warnings where standard deviation of tech replicates is >1
  mutate(ct_warning=case_when((c_mean>36.9&d_mean>36.9)~"cd*",
                              c_mean>36.9~"c*",
                              d_mean>36.9~"d*"))%>%              ###set warnings where ct values is later than 37
  mutate(rep=case_when((is.na(d_sd)&d_mean>0)|(is.na(c_sd)&c_mean>0)~"1replicate"))%>%     ###set warning if only one replicate amplified (no standard deviation but mean present)
  mutate(C=as.integer(prop_c*100))%>%
  mutate(D=100-C)

write.csv(final_data, file = "UW_final_datav3.csv")  

## pull out subset of colonies and their eggs, add bleaching/lifestage info
combined_data <- final_data%>%
  mutate(bleaching = case_when(grepl("NB", well_name) ~ "NB",
                               grepl("B", well_name) ~"B"))%>%
  mutate(bleaching = ifelse(is.na(bleaching), "NB", bleaching))%>%
  mutate(type = case_when(grepl("egg", well_name) ~ "egg",
                          TRUE ~ "parent"))%>%
  mutate_at("well_name", stringr::str_remove_all, pattern = 'egg')

# add additional info to parent data
# get subset of symbiont data that was used for determining community type - T6 when possible, but substitute with T2 when not available
idx <- c(5,6,10,11,14,15,30,31,40,43:54)

regexp <- "[[:digit:]]+"
parent_data <- combined_data%>%
  filter(type == "parent")%>%
  separate(well_name, into = c("Timepoint", "Colony"), sep = " ")%>%
  mutate(Colony= str_extract(Colony, regexp))%>%
  mutate_at("Timepoint", as.factor)%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("bleaching", as.factor)%>%
  select(Timepoint, Colony, bleaching, type, prop_d)%>%
  na.omit()
  filter(Timepoint!="T1")%>%
  select(-Timepoint)

paired_parents <- parent_data[idx,]

#same info for egg bundles
egg_data <- combined_data%>%
  filter(type == "egg")%>%
  rename(Colony=well_name)%>%
  mutate(Colony= str_extract(Colony, regexp))%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("bleaching", as.factor)%>%
  select(Colony, bleaching, type, prop_d)%>%
  na.omit()

#join the two datasets back together
paired_data <- paired_parents%>%
  full_join(egg_data)

#check out the data
library(ggpubr)
bxp <- ggboxplot(
  paired_data, x = "bleaching", y = "prop_d", color = "type",
  facet.by = "Colony", short.panel.labs = FALSE
)
bxp

paired_data %>%
  group_by(Colony, type) %>%
  get_summary_stats(prop_d, type = "mean_sd")

#going into analysis - are there differences in eggs vs parents?
#get rid of combinations we don't have paired egg-parent samples for
parent.analysis.data <- paired_data %>%
  filter(Colony!=15)%>%
  filter(Colony !=32)%>%
  filter(Colony !=41)%>%
  filter(Colony !=69)%>%
  filter(Colony !=73)%>%
  tibble::rownames_to_column("ID")%>%
  filter(ID !=1)%>%
  filter(ID !=3)%>%
  filter(ID !=5)

#check out the data
bxp.parent <- ggboxplot(
  parent.analysis.data, x = "bleaching", y = "prop_d", color = "type",
  facet.by = "Colony", short.panel.labs = FALSE
)
bxp.parent

#paired model accounting for repeated measures of the same individuals (as parents, and as eggs)
parent.model =  glmer(prop_d ~ bleaching*type + (1 | Colony), family=binomial, data = parent.analysis.data)
summary(parent.model)
Anova(parent.model, type ="III")
# significant effect of type: 0.021337 * 

#Analysis of Deviance Table (Type III Wald chisquare tests)

#Response: prop_d
#                 Chisq Df Pr(>Chisq)   
#(Intercept)    10.0487  1   0.001525 **
#  bleaching       0.0195  1   0.889082   
#type            5.2990  1   0.021337 * 
#  bleaching:type  0.2290  1   0.632240  

#are there differences in B vs NB eggs?
#get rid of combinations we don't have paired B-NB egg samples for
egg.analysis.data <- paired_data %>%
  filter(Colony!=20)%>%
  filter(Colony !=25)%>%
  filter(Colony !=28)%>%
  filter(Colony !=5)%>%
  filter(Colony !=55)%>%
  filter(Colony !=65)%>%
  filter(Colony !=32)%>%
  filter(Colony !=73)

#check out the data
bxp.egg <- ggboxplot(
  egg.analysis.data, x = "bleaching", y = "prop_d",
  facet.by = "Colony", short.panel.labs = FALSE
)
bxp.egg

#paired model accounting for repeated measures of the same individuals (B and NB)
egg.model =  glmer(prop_d ~ bleaching + (1 | Colony), family=binomial, data = egg.analysis.data)
summary(egg.model)
Anova(egg.model, type ="III")

#Analysis of Deviance Table (Type III Wald chisquare tests)

#Response: prop_d
#             Chisq Df Pr(>Chisq)  
#(Intercept)  3.9194  1    0.04773 *
 # bleaching   0.5800  1    0.44632  

#create a figure
parent.fig.data <- parent.analysis.data%>%
  mutate(prop_c = 1-prop_d)%>%
  pivot_longer(cols = c(prop_d, prop_c), names_to = "symbiont")

#reorder colonies from least to highest proportion D
#parents vs eggs
parent.fig.data$Colony <- factor(parent.fig.data$Colony, levels=c("28", "13", "21", "5", "20","10","25","26","55","65","66","74"))

parent_plot <-ggplot(parent.fig.data, aes(fill=symbiont, x=type, y=value))+
  geom_bar(position="fill", stat="identity")+
  ylab("Proportion") +
  xlab("Life Stage") +
  scale_fill_brewer(palette="Blues")+
  theme_bw()+
  guides(fill = guide_legend(title = "Symbiont Genus"))+
  facet_wrap(~Colony, nrow = 1)
parent_plot

#eggs - B vs NB
egg.fig.data <- egg.analysis.data%>%
  mutate(prop_c = 1-prop_d)%>%
  pivot_longer(cols = c(prop_d, prop_c), names_to = "symbiont")

#reorder colonies from least to highest proportion D
egg.fig.data$Colony <- factor(egg.fig.data$Colony, levels=c("13", "21", "69", "10","15","26","41","66","74"))

egg_plot <-ggplot(egg.fig.data, aes(fill=symbiont, x=bleaching, y=value))+
  geom_bar(position="fill", stat="identity")+
  ylab("Proportion") +
  xlab("Bleaching Status") +
  scale_fill_brewer(palette="Blues")+
  theme_bw()+
  guides(fill = guide_legend(title = "Symbiont Genus"))+
  facet_wrap(~Colony, nrow = 1)
egg_plot

#create composite figure
library(cowplot)

composite <- plot_grid(parent_plot, egg_plot, labels = c('A', 'B'), label_size = 12, nrow = 2)
composite

## calculation of labeling efficiency of Laco.M

setwd("E:/DDA-PASEF/20230717_Microbiology_culture_proteome/Murinus.L/Output1/Murinus_Heavy")

data <- read_excel("psm_for_Labeling_efficiency.xlsx", sheet = "Sheet5")

f1 <- ggplot(data) + geom_histogram(aes(x=Percentage, y=..count..), bins = 50, fill="lightblue", colour="black")
f1+theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) +  theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Incorporation efficiency")+ ylab("PSM count")


## tier 2 ID summary

Tier2_ID_Summary <- read_excel("Tier2_ID_Summary.xlsx", sheet = "Sheet2")
library(reshape2)
data <- melt(Tier2_ID_Summary, id.vars = "Sample", variable.name = "Group", value.name = "IDs")

library(tidyr)
data1 <- data %>% separate(Group, c("Species", "Level"), "_")

data1$Level <- factor(data1$Level,levels=c("Precursor", "Peptide", "ProteinGroup"))

data1$Sample <- factor(data1$Sample,levels=c("BacPep_1pg" ,  "BacPep_2pg",   "BacPep_5pg" ,  "BacPep_10pg" , "BacPep_20pg" 
                                             ,"BacPep_50pg" , "BacPep_100pg", "BacPep_150pg" ,"BacPep_250pg" ,"BacPep_500pg"))

ggplot(data1) + geom_bar(aes(x=Sample, y=IDs, fill=Species),width=0.9,alpha=0.8,stat = "identity",position = "dodge", color="black") + facet_wrap(~Level, ncol = 1) + scale_fill_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust =0.5, angle = 90),axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"))


## adding error bar

library(reshape2)
data2 <- melt(STD, id.vars = "Sample", variable.name = "Group", value.name = "IDs")

library(tidyr)
data3 <- data2 %>% separate(Group, c("Species", "Level"), "_")

data3$Level <- factor(data3$Level,levels=c("Precursor", "Peptide", "ProteinGroup"))

data3$Sample <- factor(data3$Sample,levels=c("BacPep_1pg" ,  "BacPep_2pg",   "BacPep_5pg" ,  "BacPep_10pg" , "BacPep_20pg" 
                                             ,"BacPep_50pg" , "BacPep_100pg", "BacPep_150pg" ,"BacPep_250pg" ,"BacPep_500pg"))


data4 <- cbind(data1, data3)

data4 <- data4[, c(1:4,8)]

ggplot(data4, aes(x=Sample, y=IDs, fill=Species)) + 
  geom_bar(stat="identity", color=NA, width=0.9,alpha=0.8,
           position=position_dodge()) +
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.9)) + facet_wrap(~Level, ncol = 1, scales = "free_y")+ scale_fill_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust =0.5, angle = 90),axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"))


ggplot(data4[21:60,], aes(x=Sample, y=IDs, color=Species), group=1) + 
  geom_point(size=2,alpha=0.8, position = position_dodge(0.5)) + geom_line(group=1)+
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.5)) + facet_wrap(~Level, ncol = 1, scales = "free_y")+ scale_color_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust =0.5, angle = 90),axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"))


## Tier3 ID summary
library(readxl)
ID_Summary_Tier3 <- read_excel("ID_Summary_Tier3.xlsx", sheet = "Average")
library(reshape2)
data <- melt(ID_Summary_Tier3, id.vars = "Sample", variable.name = "Group", value.name = "IDs")

STD <- read_excel("ID_summary_Tier3.xlsx", sheet = "STD")

data1 <- data %>% separate(Group, c("Species", "Level"), "_")
library(reshape2)
data2 <- melt(STD, id.vars = "Sample", variable.name = "Group", value.name = "IDs")
data3 <- data2 %>% separate(Group, c("Species", "Level"), "_")
library(tidyr)
data4 <- cbind(data1, data3)
data4 <- data4[, c(1:4,8)]
fix(data4)

row.names(data4) <- data4$Sample
data4$Level <- factor(data4$Level,levels=c("Precursor", "Peptide", "ProteinGroup"))

data4$Sample <- factor(data4$Sample,levels=c(      
                                             "BacCell_1B", "BacCell_100M", "BacCell_10M","BacCell_5M", "BacCell_1M",  "BacCell_0.5M","BacCell_0.1M",  "BacCell_50000","BacCell_10000"  ))


library(ggplot2)
library(ggbreak)
f1 <- ggplot(data4[1:18,], aes(x=Sample, y=IDs, fill=Species)) + 
  geom_bar(stat="identity", color=NA, width=0.8,alpha=0.8,
           position=position_dodge()) +
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.8)) + scale_fill_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ xlab(NULL)+
  scale_y_cut(breaks=c(200,1500), which = c(1,2,3),scales=c(0.5,0.3,1), expand = F) + ylim(0,10000)

f2 <- ggplot(data4[19:36,], aes(x=Sample, y=IDs, fill=Species)) + 
  geom_bar(stat="identity", color=NA, width=0.8,alpha=0.8,
           position=position_dodge()) +
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.8)) + scale_fill_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ xlab(NULL)+
  scale_y_cut(breaks=c(200,1500), which = c(1,2,3),scales=c(0.5,0.3,1), expand = F) + ylim(0,9000)

f3 <- ggplot(data4[37:54,], aes(x=Sample, y=IDs, fill=Species)) + 
  geom_bar(stat="identity", color=NA, width=0.8,alpha=0.8,
           position=position_dodge()) +
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.8)) + scale_fill_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust =0.5, angle = 90),panel.border = element_blank(),axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ xlab(NULL)+
  scale_y_cut(breaks=c(250), which = c(1,2),scales=c(0.5,1), expand = F) + ylim(0,1500)



## Tier2 peptide average intensity correlation


library(corrplot)
row.names(data) <- data$ModifiedSequence
data <- as.matrix(data)
data1<- data[,-1]
data2 <- log2(data1)
M <- cor(data2, use = "pairwise.complete.obs", method = "spearman")

corrplot(M, method = "ellipse", is.corr = F, tl.col = 'black', tl.cex = 0.6)


corrplot.mixed(M, is.corr=F, lower = 'number', upper = 'square',col.lim=c(0.80, 1), tl.cex = 0.75)



corrplot(M, is.corr=F,method = 'square',type = "lower",col.lim=c(0.3, 1),addCoef.col = 'black', tl.pos = 'd',tl.cex = 1,
         cl.pos = 'b', col = COL1("Greens"), tl.col = 'black')

corrplot(M, order = 'hclust', type="lower", is.corr = F)


## Tier3 Taxa annotation and quantification 1B cells

##improt the extracted SNI data, here we have pre and 14D
library(readr)
BuiltIn_pepTaxa <- read_csv("BuiltIn.pepTaxa.csv")

data <- BuiltIn_pepTaxa[,c(1:12,18)]

##keep only rows without NA at specific taxa level, e.g. species, family...(since some peptides are only annotated at higher taxa level)
## for species level
data1 <- data[!(is.na(data$Species) | data$Species==""),]

## now we can count how many peptides are annotated at specific taxa level
data2 <- as.data.frame(table(data1$Species))


##keep species with unique peptides >= 3 
library(dplyr)
data3 <- data2 %>% 
  filter(if_all(starts_with('Freq'), ~ . > 2)) ## after this filter, all species have at least 3 peptides annotated

##make a list of remaining species
species_unique <- as.list(unique(data3$Var1))

## substract the orignal data (with peptide info) based on the unique species
data4 <- data1[data1$Species %in% unlist(species_unique),]

## extract the minimum info for the following analysis, peptides, species and then the intensity columns
data4 <- data4[,c(2, 12:13)]
data5 <- data4

## change the data into long format
library(reshape2)
data5 <- melt(data4, id.var=c("Sequence", "Species"), variable.name="Experiment", value.name = "Intensity")


##Sum the peptides (all) intensity for each species
library(dplyr)
data6 <- data5 %>% group_by(Species, Experiment) %>% summarise(SumIntensity = sum(Intensity))

library(tidyr)
## change the summarized data into wide format
data7 <- spread(data6,key=Experiment, value=SumIntensity)
## here you can export the raw intensity to the folder
write.csv(data7, "Species_Minimum_3Peptides_SumIntensity_1B_Cells.csv", sep = ",", eol = "\n")


## Tier3 Taxa annotation and quantification 100M and NegCtrl cells

##improt the extracted SNI data, here we have pre and 14D
library(readr)
BuiltIn_pepTaxa <- read_csv("BuiltIn.pepTaxa.csv")

data <- BuiltIn_pepTaxa[,c(1:12,16)]

##keep only rows without NA at specific taxa level, e.g. species, family...(since some peptides are only annotated at higher taxa level)
## for species level
data1 <- data[!(is.na(data$Species) | data$Species==""),]

## now we can count how many peptides are annotated at specific taxa level
data2 <- as.data.frame(table(data1$Species))


##keep species with unique peptides >= 3 
library(dplyr)
data3 <- data2 %>% 
  filter(if_all(starts_with('Freq'), ~ . > 2)) ## after this filter, all species have at least 3 peptides annotated

##make a list of remaining species
species_unique <- as.list(unique(data3$Var1))

## substract the orignal data (with peptide info) based on the unique species
data4 <- data1[data1$Species %in% unlist(species_unique),]

## extract the minimum info for the following analysis, peptides, species and then the intensity columns
data4 <- data4[,c(2, 12:13)]
data5 <- data4

## change the data into long format
library(reshape2)
data5 <- melt(data4, id.var=c("Sequence", "Species"), variable.name="Experiment", value.name = "Intensity")


##Sum the peptides (all) intensity for each species
library(dplyr)
data6 <- data5 %>% group_by(Species, Experiment) %>% summarise(SumIntensity = sum(Intensity))

library(tidyr)
## change the summarized data into wide format
data7 <- spread(data6,key=Experiment, value=SumIntensity)
## here you can export the raw intensity to the folder
write.csv(data7, "Species_Minimum_3Peptides_SumIntensity_100M_Cells.csv", sep = ",", eol = "\n")

## Tier3 Taxa annotation and quantification all

##improt the extracted SNI data, here we have pre and 14D
library(readr)
BuiltIn_pepTaxa <- read_csv("BuiltIn.pepTaxa.csv")

data <- BuiltIn_pepTaxa

##keep only rows without NA at specific taxa level, e.g. species, family...(since some peptides are only annotated at higher taxa level)
## for species level
data1 <- data[!(is.na(data$Species) | data$Species==""),]

## now we can count how many peptides are annotated at specific taxa level
data2 <- as.data.frame(table(data1$Species))


##keep species with unique peptides >= 3 
library(dplyr)
data3 <- data2 %>% 
  filter(if_all(starts_with('Freq'), ~ . > 2)) ## after this filter, all species have at least 3 peptides annotated

##make a list of remaining species
species_unique <- as.list(unique(data3$Var1))

## substract the orignal data (with peptide info) based on the unique species
data4 <- data1[data1$Species %in% unlist(species_unique),]

## extract the minimum info for the following analysis, peptides, species and then the intensity columns
data4 <- data4[,c(2, 12:22)]
data5 <- data4

## change the data into long format
library(reshape2)
data5 <- melt(data4, id.var=c("Sequence", "Species"), variable.name="Experiment", value.name = "Intensity")


##Sum the peptides (all) intensity for each species
library(dplyr)
data6 <- data5 %>% group_by(Species, Experiment) %>% summarise(SumIntensity = sum(Intensity))

library(tidyr)
## change the summarized data into wide format
data7 <- spread(data6,key=Experiment, value=SumIntensity)
## here you can export the raw intensity to the folder
write.csv(data7, "Species_Minimum_3Peptides_SumIntensity_All.csv", sep = ",", eol = "\n")


## 100M species abundance

Species_Minimum_3Peptides_SumIntensity_100M_Cells <- read_excel("Species_Minimum_3Peptides_SumIntensity_100M_Cells.xlsx", sheet = "Sheet1")

data <- Species_Minimum_3Peptides_SumIntensity_100M_Cells[,c(1,5)]
library(treemapify)
library(ggplot2)
fix(data)
library(dplyr)

data$Species <- factor(data$Species, levels = data$Species)
ggplot(data[1:11,], aes(area = Percentage, fill = Species, label = paste(Percentage, "%"))) +
  geom_treemap(alpha=0.7) +
  geom_treemap_text(colour = "black",
                    place = "centre",
                    size = 12)+ scale_fill_brewer(palette = "RdYlBu")
f1



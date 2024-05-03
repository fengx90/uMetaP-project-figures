## Tier1 ID summary

setwd("E:/DIA-PASEF/20240112_LOD_Test_Ultra_Level1/Dilution_1")
library(readxl)
ID_Summary_Tier1 <- read_excel("ID_Summary.xlsx", sheet = "Average")
library(reshape2)
data <- melt(ID_Summary_Tier1, id.vars = "Sample", variable.name = "Group", value.name = "IDs")

STD <- read_excel("ID_Summary.xlsx", sheet = "STD")
library(dplyr)
data1 <- data %>% separate(Group, c("Species", "Level"), "_")
library(reshape2)
data2 <- melt(STD, id.vars = "Sample", variable.name = "Group", value.name = "IDs")
data3 <- data2 %>% separate(Group, c("Species", "Level"), "_")
library(tidyr)
data4 <- cbind(data1, data3)
data4 <- data4[, c(1:4,8)]
fix(data4)
data4$Level <- factor(data4$Level,levels=c("Precursor", "Peptide", "ProteinGroup"))

data4$Sample <- factor(data4$Sample,levels=c(      
  "500fg", "1pg", "2.5pg","5pg", "10pg",  "25pg","50pg",  "75pg","125pg","250pg"))
library(ggplot2)
library(ggbreak)
ggplot(data4[c(11:20,41:50),], aes(x=Sample, y=IDs, fill=Species)) +
  geom_linerange(aes(x = Sample, ymin = 0, ymax = IDs, color = Species), position = position_dodge(.6),lwd = 2) +
  geom_point(position = position_dodge2(width=0.6, preserve="single"),size = 5, aes(shape=Species), stat = "identity")+ scale_shape_manual(values = c(21,23))+
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.6)) + facet_wrap(~Level, ncol = 1, scales = "free_y")+ scale_fill_manual(values = c("lightblue","salmon"))+scale_color_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ xlab(NULL)+scale_y_cut(breaks=c(900,5000), which = c(1,2,3),scales=c(0.5,0.4,1), expand = F)+ylim(0,12000)

ggplot(data4[c(21:30,51:60),], aes(x=Sample, y=IDs, fill=Species)) +
  geom_linerange(aes(x = Sample, ymin = 0, ymax = IDs, color = Species), position = position_dodge(.6),lwd = 2) +
  geom_point(position = position_dodge2(width=0.6, preserve="single"),size = 5, aes(shape=Species), stat = "identity")+ scale_shape_manual(values = c(21,23))+
  geom_errorbar(aes(ymin=IDs-Sd, ymax=IDs+Sd), width=.5,
                position=position_dodge(.6)) + facet_wrap(~Level, ncol = 1, scales = "free_y")+ scale_fill_manual(values = c("lightblue","salmon"))+scale_color_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust =0.5, angle = 90),panel.border = element_blank(),axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("IDs") + theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ xlab(NULL)+ scale_y_cut(breaks=c(450), which = c(1,2),scales=c(0.5,1), expand = F)+ ylim(0,1800)

library(ggpubr)
ggarrange(f1,f2,nrow = 2)
## quantitative accuracy peptide

data <- read_excel("Peptide_Intensity_quantitative_accuracy.xlsx", sheet = "Average")
library(reshape2)
data1 <- melt(data, id.vars = c("ModifiedPeptide","Species"), variable.name = "Injection", value.name = "Intensity")
fix(data1)
data1$Injection <- factor(data1$Injection,levels=c(      
  "500fg", "1pg", "2.5pg","5pg", "10pg",  "25pg","50pg",  "75pg","125pg","250pg"))
ggplot(data1)+geom_boxplot(aes(x=Injection, y=log2(Intensity), fill=Species)) + facet_wrap(~Species, ncol = 2)



##L.Murinus Protein IDs venn 
library(eulerr)
fit <- euler(c("250pg&1pg&500fg" = 124, "250pg&1pg" = 133,"250pg&500fg" = 8, 
               "500fg&1pg" = 0, "250pg"=924, "1pg"=4
), shape = "ellipse", control = list(extraopt = T)
)

plot(fit, quantities = list(font=7),edges = list(lty = 1),
     fills = list(fill = c("lightblue", 
                           "gray89", "white"
                           
     )),legend = list(side = "top"), hjust=1)

##L.Murinus Protein IDs venn 
library(eulerr)
fit <- euler(c("250pg&1pg&500fg" = 128, "250pg&1pg" = 130,"250pg&500fg" = 5, 
               "500fg&1pg" = 2, "250pg"=1474, "1pg"=1, "500fg"=1
), shape = "ellipse",control = list(extraopt = T)
)

plot(fit, quantities = list(font=7),edges = list(lty = 1),
     fills = list(fill = c("salmon", 
                           "gray89", "white"
                           
     )),legend = list(side = "top"), hjust=1)


## quantitative accuracy 500fg-250pg

data <- read_excel("Peptide_Intensity_quantitative_accuracy.xlsx", sheet = "Median")
library(reshape2)
data1 <- melt(data, id.vars = "Sample", variable.name = "Species", value.name = "Median")
library(ggplot2)
f1 <- ggplot(data1,aes(x=Sample, y=Median, fill=Species)) + 
  geom_point(position = position_dodge2(width=0.6, preserve="single"),size = 5, aes(shape=Species), stat = "identity")+ scale_shape_manual(values = c(21,23)) + scale_fill_manual(values = c("lightblue","salmon"))+ geom_line(stat = "identity", size=1)+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.ticks.x = element_blank(),panel.border = element_rect(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  ylab("Median-Exp.ratio/Expected ratio") + theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+
  facet_wrap(~Species,ncol = 1)+ylim(0,3)+ geom_abline(intercept = 1, slope = 0, color="orange",linetype="dashed", size=1)+ xlab("Loaded peptide (pg)")

## Cv distribution of common peptides
library(readxl)
CV <- read_excel("Peptide_Intensity_quantitative_accuracy.xlsx", sheet = "Technical_CV")
library(reshape2)
data2 <- melt(CV, id.vars = c("ModifiedPeptide", "Species"), variable.name = "Peptide", value.name = "CV")
library(ggplot2)
f2<- ggplot(data2)+ geom_violin(aes(x=Peptide, y=CV, color=Species), size=1,fill="white",trim=T)+ geom_boxplot(aes(x=Peptide, y=CV),width=0.08,outlier.size = 1.2, outlier.color = "black")+scale_color_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ theme(panel.border = element_rect())+
  facet_wrap(~Species,ncol = 1)+ geom_abline(intercept = 0.2, slope = 0, color="red",linetype="dashed", size=1)+ xlab("Loaded peptide (pg)")

library(ggpubr)
ggarrange(f1,f2, nrow = 1, widths = c(0.7,0.5))


CV <- read_excel("Peptide_Intensity_quantitative_accuracy.xlsx", sheet = "Sheet6")
data2 <- melt(CV, id.vars = c("ModifiedPeptide", "Species"), variable.name = "Peptide", value.name = "CV")

f2<- ggplot(data2)+ geom_violin(aes(x=Peptide, y=CV, color=Species), size=1,fill="white",trim=T)+ geom_boxplot(aes(x=Peptide, y=CV),width=0.08,outlier.size = 1.2, outlier.color = "black")+scale_color_manual(values = c("lightblue","salmon"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank()) + 
  theme(strip.text = element_text(size=12, face = "plain"), axis.line = element_line())+ theme(panel.border = element_rect())+
  facet_wrap(~Species,ncol = 1)+ geom_abline(intercept = 0.2, slope = 0, color="red",linetype="dashed", size=1)+ xlab("Loaded peptide (pg)")





##improt the extracted SNI data, here we have pre and 14D
library(readr)
data <- read_csv("BuiltIn.pepTaxa.csv")

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
write.csv(data7, "Species_Minimum_3Peptides_SumIntensity_25ng.csv", sep = ",", eol = "\n")


## 25ng species abundance

data <- read_excel("Species_Minimum_3Peptides_SumIntensity_25ng.xlsx", sheet = "Sheet2")

data <- data[,c(1,3)]
library(treemapify)
library(ggplot2)


f1 <- ggplot(data[1:11,], aes(area = Percentage, fill = Species, label = paste(Percentage, "%"))) +
  geom_treemap(alpha=0.7) +
  geom_treemap_text(colour = "black",
                    place = "centre",
                    size = 12)+ scale_fill_brewer(palette = "RdYlBu")
f1

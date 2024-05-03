setwd("D:/Data_Bruker/DIA/DIA_0.01-100ng_PD48970_search")

## dynamic range peptide 50ng-30min
data <- read_excel("Peptide_intensity_50ng-30min_for_Dynamic_range.xlsx")
data1 <- data[,c(1,9)]
library(ggplot2)
f1 <- ggplot(data, aes(x=PeptideRank, y=Log10Intensity)) + geom_line(colour="#CDC673",size=2)+ 
  geom_abline(slope = 0, intercept = 5.2344, color="black", linetype=2, size=1)+
  geom_abline(slope = 0, intercept = 3.8373, color="black", linetype=2, size=1)

f1<- f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Peptide rank")+ ylab("Log10(Peptide intensity)")+
  theme(panel.grid = element_blank())
f1
library(tidyr)
data$Cato <- factor(data$Cato,levels=c("High", "Medium", "Low"))
f2 <- ggplot(data, aes(x=CV, y=..scaled..)) + geom_density(aes(fill=Cato), alpha=0.8)+facet_wrap(~Cato, nrow = 3) + geom_vline(xintercept = 0.2, color="red", linetype=2,size=1)+
  scale_fill_manual(values = c("khaki", "darkkhaki", "khaki4")) + xlim(0,1)
f2 <- f2+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain"), legend.position = "none") + xlab("Coefficient of varation")+ ylab("Density")+
  theme(panel.grid = element_blank()) + theme(strip.text = element_text(size=11))
f2
library(ggpubr)
f3 <- ggarrange(f1,f2, ncol = 2, widths = c(1,0.7))
f3


## dynamic range peptide 50ng-66min
data <- read_excel("Peptide_intensity_50ng-66min_for_Dynamic_range.xlsx")
data1 <- data[,c(1,9)]
library(ggplot2)
f1 <- ggplot(data, aes(x=PeptideRank, y=Log10Intensity)) + geom_line(colour="#CDC673",size=2)+ 
  geom_abline(slope = 0, intercept = 5.4094, color="black", linetype=2, size=1)+
  geom_abline(slope = 0, intercept = 3.9076, color="black", linetype=2, size=1)

f1<- f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Peptide rank")+ ylab("Log10(Peptide intensity)")+
  theme(panel.grid = element_blank())
f1
library(tidyr)
data$Cato <- factor(data$Cato,levels=c("High", "Medium", "Low"))
f2 <- ggplot(data, aes(x=CV, y=..scaled..)) + geom_density(aes(fill=Cato), alpha=0.8)+facet_wrap(~Cato, nrow = 3) + geom_vline(xintercept = 0.2, color="red", linetype=2,size=1)+
  scale_fill_manual(values = c("khaki", "darkkhaki", "khaki4")) + xlim(0,1)
f2 <- f2+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain"), legend.position = "none") + xlab("Coefficient of varation")+ ylab("Density")+
  theme(panel.grid = element_blank()) + theme(strip.text = element_text(size=11))
f2
library(ggpubr)
f3 <- ggarrange(f1,f2, ncol = 2, widths = c(1,0.7))
f3


## peptide abundance fold change
setwd("D:/Data_Bruker/DIA/DIA_0.01-100ng_PD48970_search")
data <- read_excel("StrippedPeptide_PepAbundanceChange.xlsx", sheet = "Sheet3")
data1 <- data[,7:12]
data1 <- data[,c(1,8:13)]
library(reshape2)
data2 <- melt(data1, id.vars = "StrippedSequence", variable.name = "Peptide", value.name = "FC")
data2$Peptide <- factor(data2$Peptide,levels=c("1ng", "2.5ng", "5ng","10ng","25ng", "50ng"))
library(ggplot2)
f1 <- ggplot(data2, aes(x=Peptide, y=FC)) + geom_violin(aes(fill=Peptide), colour="black") + geom_boxplot(width=.1, fill="white") + scale_fill_brewer(palette = 1)
f1 + ylim(0,50)


## protein quantification CV cross 3 replicates
library(readxl)
data <- read_csv("ProteinAbundance_forCV.csv")
library(reshape2)
data1 <- melt(data, id.vars = "ProteinGroup", variable.name = "Peptide", value.name = "CV")
data1$Peptide <- factor(data1$Peptide,levels=c("10pg","25pg","50pg","100pg","250pg","500pg","1ng", "2.5ng", "5ng","10ng","25ng", "50ng"))
library(ggplot2)
library(dplyr)
data2 <- na.omit(data1)
dummy <- data2 %>%
  group_by(Peptide) %>%
  summarize(median = median(CV/100))
f1 <- ggplot(data1,aes(x=CV/100)) +geom_density(fill="#71b7be",aes(y=..scaled..,),position = "identity",stat = "density",colour="gray", alpha=0.7,size=0.7)+ 
  scale_fill_brewer(palette = "Spectral")+ facet_wrap(~Peptide, nrow = 12) + scale_x_continuous(position = "top")


# Add median lines
f1 <- f1+geom_vline(data=dummy, aes(xintercept=median),size=0.7,color="red",
                linetype="dashed") + geom_text(data=dummy,aes(label=median, x=median+0.08, y=0.6)) + xlim(0,1)

f1 + theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain"), legend.position = "none") +ylab("Density")+ xlab("Coefficient of varation")+
  theme(panel.grid = element_blank()) + theme(strip.text = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

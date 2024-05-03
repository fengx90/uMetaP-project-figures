data <- read_excel("Protein_Intensity_All_For_coreelation.xlsx", sheet = 2)
 row.names(data) <- data$ProteinGroup
 data1 <- as.matrix(data)
 data1 <- data1[,-1]
 data1 <- as.data.frame(data1)
 fix(data1)
 data2 <- log2(data1)
 data3 <- cor(data2, use = "pairwise.complete.obs", method = "pearson")
 write.csv(data3, "Protein_quantification_sample_correlations.csv", eol = "\n")
 
 ##dynamic range
 data <- read_excel("Peptide_intensity_25ng_for_Dynamic_range.xlsx")
 data1 <- data[,c(1,9)]
 library(ggplot2)
f1 <- ggplot(data1, aes(x=PeptideRank, y=Log10Intensity)) + geom_point(size=3, alpha=0.5, color="lightblue")+ 
  geom_abline(slope = 0, intercept = 5.001, color="black", linetype=2, size=1)+
  geom_abline(slope = 0, intercept = 3.454, color="black", linetype=2, size=1)

f1<- f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Peptide rank")+
  theme(panel.grid = element_blank())

library(tidyr)
data$Cato <- factor(data$Cato,levels=c("High", "Medium", "Low"))
f2 <- ggplot(data, aes(x=CV, y=..scaled..)) + geom_density(aes(fill=Cato), alpha=0.7)+facet_wrap(~Cato, nrow = 3) + geom_vline(xintercept = 0.2, color="red", linetype=2,size=1)+
  scale_fill_brewer(palette = 1) 
f2 <- f2+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain"), legend.position = "none") + xlab("Coefficient of varation")+ ylab("Density")+
  theme(panel.grid = element_blank()) + theme(strip.text = element_text(size=11))
f2
library(ggpubr)
f3 <- ggarrange(f1,f2, ncol = 2, widths = c(1,0.5))
f3


data <- read_excel("Peptide_intensity_10ng_for_Dynamic_range.xlsx")
data1 <- data[,c(1,9)]
library(ggplot2)
f1 <- ggplot(data1, aes(x=PeptideRank, y=Log10Intensity)) + geom_line(colour="#CDC673",size=2)+ 
  geom_abline(slope = 0, intercept = 4.948, color="black", linetype=2, size=1)+
  geom_abline(slope = 0, intercept = 3.417, color="black", linetype=2, size=1)

f1<- f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Peptide rank")+ ylab("Log10(Intensity)")+
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
f3 <- ggarrange(f1,f2, ncol = 2, widths = c(1,0.5))
f3

##protein boxplot
data <- read_excel("Protein_averaged_intensity_Boxplot.xlsx", sheet = 1)

library(reshape2)

data<- data[,-1]
data1 <- melt(data[,1:13], id.vars = c("ProteinGroup"), variable.name = "Peptide", value.name = "Intensity")


f1 <- ggplot(data=data1, aes(x=Peptide,y=log2(Intensity)))+geom_boxplot(width=0.3, outlier.alpha = 0.1) + scale_fill_brewer(palette = "paired")
f1

##correlation
data <- read_excel("Protein_10ng_Correlation_30min-66min.xlsx")
library(ggplot2)
f1 <- ggplot(data, aes(x=data$`10ng-30min`, y=data$`10ng-60min`, color=Org))+ geom_point(alpha=0.5, size=2)+scale_color_manual(values = c("lightsteelblue4","darkslategray4"))+xlab("Log2(Intensity) 10ng_30min")+ylab("Log2(Intensity) 10ng_66min")+geom_smooth(method=lm, se=FALSE, linetype="dashed",
                                                                                                                                                                                                          color="darkred", size=1.2)
f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain"), legend.position = "top")




write.table(data1, "HumanGutMicrobiome_Database_1.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE, fileEncoding = "utf8")



## plot simulation data
data <- read_excel("Taxa_Summary_simulation.xlsx")
View(data)                                                                                                                                 
library(reshape2)
data1 <- melt(data, id.vars = "NO.Simulation", variable.name = "Sample", value.name = "Number")
fix(data1)
library(ggplot2)
data1$Taxa <- factor(data1$Taxa,levels=c("Phylum","Class","Order","Family","Genus","Species"))
f1 <- ggplot(data=data1, aes(x=Taxa,y=Number))+geom_boxplot(aes(fill=Taxa), alpha=0.6,width=0.8, outlier.color = "gray") + scale_fill_brewer(palette = 1)+ facet_wrap(~Taxa, nrow = 1, scales = "free")+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "plain", vjust = 1, hjust = 1, angle = 45)) + 
  theme(title= element_text(size=12, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5)) + 
  theme(panel.grid.major = element_line(color = "black", linetype = 2)) + geom_jitter(position=position_jitter(), pch=22,color="black", size=2.0, alpha=0.9)
f1

f1 <- f1 + theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain.italic", colour="black", vjust = 0.5, hjust = 0.5)) + theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain.italic"),legend.title = element_text(size = 12, face = "plain")) + theme(axis.text = element_text(size = 10, face = "plain", color = "black"), axis.title = element_text(size = 10, face = "plain"), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black", size = 1), axis.ticks.x = element_blank())

A_new <- c(38)
B_new <- c(61)
C_new <- c(96)
D_new <- c(141)
E_new <- c(172)
F_new <- c(237)
new_vals <- c(A_new, B_new, C_new, D_new, E_new, F_new)
new_data <- data.frame(name=c("Phylum","Class","Order","Family","Genus","Species"), value=new_vals)
new_data$Taxa <- factor(new_data$Taxa,levels=c("Phylum","Class","Order","Family","Genus","Species"))
f2 <- ggplot()+ geom_point(data=new_data, aes(x=Taxa, y=Number), fill="orange", color="blue", size=2, shape=22)+ facet_wrap(~Taxa, nrow = 1, scales = "free")+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "plain", vjust = 1, hjust = 1, angle = 45)) + 
  theme(title= element_text(size=12, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5)) + 
  theme(panel.grid.major = element_line(color = "black", linetype = 2))
f2

f2 <- f2 + theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain.italic", colour="black", vjust = 0.5, hjust = 0.5)) + theme(legend.position = "none", legend.text = element_text(size = 10, face = "plain.italic"),legend.title = element_text(size = 12, face = "plain")) + theme(axis.text = element_text(size = 10, face = "plain", color = "black"), axis.title = element_text(size = 10, face = "plain"), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black", size = 1), axis.ticks.x = element_blank())



## Ultra_DDA protein dynamic range

setwd("D:/Data_Bruker/Ultra20ng_Pro100ng_DDA_comaprison")
library(readxl)
data <- read_excel("DynamicRange_Ultra_proteins.xlsx", sheet = "Sheet2")

library(ggplot2)
f1 <- ggplot(data, aes(x=Rank, y=Log2Intensity)) + geom_point(shape=21,aes(colour=Cat),size=3, fill=NA)+scale_colour_manual(values = c("#83699b","mediumorchid3"))

f1<- f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Protein rank")+ ylab("Log2Intensity)")+
  theme(axis.line = element_line(colour = "black"),panel.grid = element_blank(), legend.position = "none", panel.border = element_blank())
f1

f2 <-ggplot(data, aes(y=Log2Intensity, x=..count..)) + geom_histogram(aes(fill=Cat),colour="black",alpha=0.7, bins = 30, position = "identity") + scale_fill_manual(values = c("#83699b","mediumorchid3"))

f2<- f2+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Protein count")+
  theme(axis.line.x = element_line(colour = "black"),panel.grid = element_blank(), legend.position = "none") + theme(panel.border = element_blank(),axis.title.y = element_blank(),axis.text.y= element_blank(), axis.ticks.y = element_blank())

library(ggpubr)
ggarrange(f1,f2,ncol = 2, widths = c(1,0.3))



## PD48970 in silico peptides taxa annotation
setwd("D:/Data_Bruker/DIA/Taxa_analysis_InSilico_Peptides")
data <- read_excel("Summary_taxa.xlsx")
library(ggplot2)
data$Rank <- factor(data$Rank,levels=c("Phylum","Class","Order","Family","Genus","Species"))
f1 <- ggplot(data) + geom_bar(stat = "identity",position = "dodge", aes(x=Rank, y= Number, fill=Rank))+ scale_fill_brewer(palette = "OrRd")
f1+theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain"), legend.position = "none")+
  theme(axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(), legend.position = "none") + theme(panel.border = element_blank())


f1

data1 <- read_excel("Summary_taxa.xlsx", sheet = "Sheet2")
data1 <- data1[c(1:11), c(1,3)]
data2 <- data1 %>% 
  arrange(desc(Species)) %>%
  mutate(prop = Percentage / sum(data1$Percentage) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

f2 <- ggplot(data2, aes(x="", y=Percentage, fill=Species)) +
  geom_bar(stat="identity", width=1, alpha=0.7) +
  coord_polar("y", start=0)+theme_void() + scale_fill_brewer(palette = "RdYlBu")+
  geom_text(aes(y = ypos, label = Percentage), color = "black", size=5)
f2


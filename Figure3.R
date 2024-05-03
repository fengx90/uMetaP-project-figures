## comparison of taxa annotation between 25ng and previous 250ng 
library(readxl)
data <- read_excel("25ng_taxa_annotation_extracted.xlsx", sheet = "Sheet8")
library(reshape2)
library(ggplot2)
library(tidyverse)
data1 <- melt(data, id.vars = "Taxa", variable.name = "Sample", value.name = "Number")

data1 %>%
  mutate(across(Taxa, factor, levels=c("Phylum","Class", "Order", "Family", "Genus","Species"))) %>%
  ggplot() + geom_bar(aes(x=Sample, y=Number, fill=Sample),width=0.9,alpha=0.8,stat = "identity",position = "dodge", color="black") + facet_wrap(~Taxa, ncol = 6) + scale_fill_manual(values = c("#CDC673","lavender"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = c(0,50, 100,200)) + ylab("Annotated taxa")+ xlab(NULL) + theme(strip.text = element_text(size=12, face = "plain"))


## biomass analysis of 25ng

##improt the extracted SNI data, here we have pre and 14D
library(readxl)
data <- read_excel("25ng_taxa_annotation_extracted.xlsx")

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


## 25ng COG analysis

data <- read_excel("D:/Data_Bruker/DIA/DIA_0.01-100ng_PD48970_search/DIA_0.01-100ng_PD48970_Taxa_Function_analysis/functional_annotation/25ng_Function_annotation_extracted.xlsx", sheet = "Sheet3")
library(pheatmap)

data1 <- data[,c(2,4)]
row.names(data1) <- data1$Description
data2 <- data1[,-1]
row.names(data2) <- data1$Description

data2 <- data1 %>% 
  arrange(desc(Species)) %>%
  mutate(prop = Percentage / sum(data1$Percentage) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

f2 <- ggplot(data2, aes(x="", y=Percentage, fill=Species)) +
  geom_bar(stat="identity", width=1, alpha=0.7) +
  coord_polar("y", start=0)+theme_void() + scale_fill_brewer(palette = "RdYlBu")+
  geom_text(aes(y = ypos, label = Percentage), color = "black", size=5)
f2





# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)

# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
library(packcircles)
library(ggplot2)
library(viridis)
# Create data
data <- read_excel("D:/Data_Bruker/DIA/DIA_0.01-100ng_PD48970_search/DIA_0.01-100ng_PD48970_Taxa_Function_analysis/functional_annotation/25ng_Function_annotation_extracted.xlsx", sheet = "Sheet3")
data1 <- data[,c(1,2,4)]

# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(data1$Number, sizetype='area')

# We can add these packing information to the initial data frame
data2 <- cbind(data1, packing)

# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
# Make the plot
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
  
  # Add text in the center of each bubble + control its size
  geom_text(data = data2, aes(x, y, size=Number, label = Cluster)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()



edges <- read_excel("D:/Data_Bruker/DIA/DIA_0.01-100ng_PD48970_search/DIA_0.01-100ng_PD48970_Taxa_Function_analysis/functional_annotation/25ng_Function_annotation_extracted.xlsx", sheet = "Sheet5")
vertices <- read_excel("D:/Data_Bruker/DIA/DIA_0.01-100ng_PD48970_search/DIA_0.01-100ng_PD48970_Taxa_Function_analysis/functional_annotation/25ng_Function_annotation_extracted.xlsx", sheet = "Sheet6")


p <- ggraph(mygraph, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = as.factor(depth), color = as.factor(depth) , alpha=0.5)) +
  scale_fill_manual(values=c("0" = "white", "1" = "#440d54", "2" = "yellow", "3" = viridis(4)[3], "4"=viridis(4)[4])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black") ) +
  theme_void() + 
  theme(legend.position="FALSE")
library(data.tree)
tree <- FromDataFrameNetwork(edges)
mylevels <- data.frame( name=tree$Get('name'), level=tree$Get("level") )

vertices <- vertices %>% 
  left_join(., mylevels, by=c("name"="name"))

vertices <- vertices %>% 
  mutate(new_label=ifelse(!level==1, name, NA))

mygraph <- graph_from_data_frame( edges, vertices=vertices )

ggraph(mygraph, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = as.factor(depth), color = as.factor(depth) , alpha=0.5)) +
  scale_fill_manual(values=c("0" = "white", "1" = "#440d54", "2" = "yellow", "3" = viridis(4)[3], "4"=viridis(4)[4])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black") ) +
  geom_node_label(aes(label=new_label), size=4) +
  theme_void() + 
  theme(legend.position="FALSE", plot.margin = unit(rep(0,5), "cm"))


## simulation data taxa
data <- read_excel("Taxa_Summary_simulation.xlsx", sheet = "Sheet2")

library(ggplot2)
data$Rank <- factor(data$Rank,levels=c("Phylum","Class","Order","Family","Genus","Species"))
data1 <- melt(data, id.vars = "Rank", variable.name = "Group", value.name = "Number")
library(tidyverse)
library(dplyr)

fix(data1) ## change each workflow into numbers 1-4

data1$Group <- factor(data1$Group,levels=c("Simulation_Max","Simulation_Min","DIA-PASEF"))
fix(data1)
data2 <- filter(data1, Group %in% c("Simulation_Max", "Simulation_Min")) %>% group_by(data2$Rank) %>% arrange(Group) %>% mutate(pos2 = sum(Number) - Number / 2)               

data3 <- filter(data1, Group %in% c("DIA-PASEF")) %>% group_by(data3$Rank) %>% arrange(Group) %>% mutate(pos2 = sum(Number) - Number / 2)
barwidth = 0.4
library(ggplot2)
f1 <- ggplot() + 
  geom_bar(data = data2, mapping=aes(y = Rank, x = Number, fill = Group),  colour="black",stat="identity", position='stack', width = barwidth) +
  geom_bar(data = data3, aes(y = Rank-barwidth, x = Number,fill = Group),colour="black",stat="identity", position='dodge', width = barwidth) + 
  theme(axis.text.y = element_blank())+ scale_fill_manual(values = c("#999999","#666666","#CDC673"))
f1 <- f1 + theme_bw() + theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.text.x = element_text(angle=0, vjust = 0, hjust = 0.5),axis.title = element_text(size = 12, face = "plain", color = "black"),legend.position = "top") +xlab("Annotated taxa") +
  scale_y_discrete(limits=c("Phylum","Class","Order","Family","Genus","Species"))
f1


##Theoretical, simulation, DIA-PASEF comparison, species


library(eulerr)
fit <- euler(c("Protein_database" = 219, "DIA-PASEF" = 0, "Simulation" = 0,
               "Protein_database&DIA-PASEF" = 2, "Protein_database&Simulation&DIA-PASEF" = 218,
              "Protein_database&Simulation" = 373, "Simulation&DIA-PASEF"=0), shape = "ellipse", control = list(extraopt = T)
)
plot(fit, quantities = TRUE
)

plot(fit, quantities = list(font=2),edges = list(lty = 1),
     fills = list(fill = c( "#bc4a64","#bfaf7f", "#549e39"
     )),legend = list(side = "top"), hjust=1)


##histogram of detected and undeteced species
data <- read_excel("Taxa_Summary_simulation.xlsx", sheet = "Sheet5")

f1 <- ggplot(data, aes(x=Rank, y=log10(Percentage))) + geom_point(shape=21,aes(colour=Group),size=5, fill=NA)+scale_colour_manual(values = c("#bfaf7f","#549e39"))

f1<- f1+ theme_bw()  + theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), axis.title = element_text(size = 12,face = "plain")) + xlab("Species abundance rank")+ ylab("Log10 (Abundance %)")+
  theme(axis.line = element_line(colour = "black"),panel.grid = element_blank(), legend.position = "top", panel.border = element_blank())+
  geom_abline(slope = 0, intercept = -1.301, color="black", linetype=2, size=1)+geom_abline(slope = 0, intercept = -2, color="black", linetype=2, size=1)+ylim(-2.5,1.5)
f1


## proteins with specific functions

data <- read_excel("PUFs_25ng_30min_Database.xlsx", sheet = "Sheet1")
data1 <- melt(data, id.vars = "Category", variable.name = "Group", value.name = "Number")
data1 %>%
  mutate(across(Category, factor, levels=c("PUFs","Small proteins", "AMPs"))) %>%
  ggplot() + geom_bar(aes(x=Category, y=Number, fill=Group),width=0.9,alpha=0.8,stat = "identity",position = "dodge", color="black") + scale_fill_manual(values = c("#CDC673","#bc4a64"))+
  theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 10, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.ticks.x = element_blank(),axis.text = element_text(size = 14, face = "plain", color = "black"), axis.title = element_text(size = 12, face = "plain"), panel.grid.minor = element_blank())+
  ylab("Number") + theme(strip.text = element_text(size=12, face = "plain"))


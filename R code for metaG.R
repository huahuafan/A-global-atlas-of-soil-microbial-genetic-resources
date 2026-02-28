# 1. calculating the Jaccard distance
#for the 200 sampling plot
library(vegan)
Ko <- read.table("1.Jaccard distance.zip/KEGG_200_2025.txt",header = TRUE)
df.abundances=data.frame(Ko)
# dataframe with presence/absence
df.binary = df.abundances
#df.binary[df.binary > 0] = 1
dis=vegdist(t(df.binary), method = "jaccard", binary = T)
dis=as.matrix(dis)
write.table(dis,"dis_KEGG_200_2025.txt",sep="\t")
#for the 1609 sampling plot
Ko <- read.table("1.Jaccard distance.zip/KEGG_1609_2025.txt",header = TRUE)
df.abundances=data.frame(Ko)
# dataframe with presence/absence
df.binary = df.abundances
#df.binary[df.binary > 0] = 1
dis=vegdist(t(df.binary), method = "jaccard", binary = T)
dis=as.matrix(dis)
write.table(dis,"dis_KEGG_1609_2025.txt",sep="\t")

#2.Stacked Bar Chart for protected areas across different climate regions and vegetation types
library (reshape2)
library(ggplot2) # CRAN v3.4.2
library(ggprism) # CRAN v1.0.4
library(plyr)    # CRAN v1.8.8
library(ggthemes)
data <- read.table("2.Stacked Bar Chart for protected areas.zip/protected_areas_traits.txt", header = T, sep = '\t', check.names = FALSE)
DJ <- melt(data,id.vars='x_axis')
names(DJ)[1:2] <- c("Taxonomy","sample")  
DJ$sample <- factor(DJ$sample, levels = unique(DJ$sample))       
DJ$Taxonomy <- factor(DJ$Taxonomy, levels = unique(DJ$Taxonomy)) 

p1<-ggplot(DJ, aes( x = sample,y = value, fill = Taxonomy))+  
  geom_col(position = 'stack', width = 0.5, show.legend = F)+  #width = 0.5 bar width，show.legend = F remove legend
  scale_y_continuous(expand = c(0,0))+  
  labs(x="",y="Proportion of areas(%)", 
       fill="")+   
  coord_flip()+    
  theme_few()+     
  scale_fill_manual(values=c("#FFFEE0","#C0C000"))  
p1
saveRDS(p1, paste0("/2.Stacked Bar Chart for protected areas", "/traits_enzyme.rds"))
ggsave(plot = p1, filename = paste0("/2.Stacked Bar Chart for protected areas", "/traits_enzyme.png"),width = 5, height = 7, dpi = "retina")

data <- read.table("2.Stacked Bar Chart for protected areas.zip/protected_areas_gene_richness_Fig1_1.txt", header = T, sep = '\t', check.names = FALSE)
DJ <- melt(data,id.vars='x_axis')
names(DJ)[1:2] <- c("Taxonomy","sample")  
DJ$sample <- factor(DJ$sample, levels = unique(DJ$sample))       
DJ$Taxonomy <- factor(DJ$Taxonomy, levels = unique(DJ$Taxonomy)) 

p2<-ggplot(DJ, aes( x = sample,y = value, fill = Taxonomy))+  
  geom_col(position = 'stack', width = 0.5, show.legend = F)+  #width = 0.5 bar width，show.legend = F remove legend
  scale_y_continuous(expand = c(0,0))+  
  labs(x="",y="Proportion of areas(%)", 
       fill="")+   
  coord_flip()+    
  theme_few()+     
  scale_fill_manual(values=c("#FFFEE0","#C0C000"))  
p2
saveRDS(p2, paste0("/2.Stacked Bar Chart for protected areas", "/protected_areas_gene_richness_Fig1_1.rds"))
ggsave(plot = p2, filename = paste0("/2.Stacked Bar Chart for protected areas", "/protected_areas_gene_richness_Fig1_1.png"),width = 4, height = 6.4, dpi = "retina")

data <- read.table("2.Stacked Bar Chart for protected areas.zip/protected_areas_gene_dissimilarity_Fig1_2.txt", header = T, sep = '\t', check.names = FALSE)
DJ <- melt(data,id.vars='x_axis')
names(DJ)[1:2] <- c("Taxonomy","sample")  
DJ$sample <- factor(DJ$sample, levels = unique(DJ$sample))       
DJ$Taxonomy <- factor(DJ$Taxonomy, levels = unique(DJ$Taxonomy)) 

p3<-ggplot(DJ, aes( x = sample,y = value, fill = Taxonomy))+  
  geom_col(position = 'stack', width = 0.5, show.legend = F)+  #width = 0.5 bar width，show.legend = F remove legend
  scale_y_continuous(expand = c(0,0))+  
  labs(x="",y="Proportion of areas(%)", 
       fill="")+   
  coord_flip()+    
  theme_few()+     
  scale_fill_manual(values=c("#FFFEE0","#C0C000"))  
p3
saveRDS(p3, paste0("/2.Stacked Bar Chart for protected areas", "/protected_areas_gene_dissimilarity_Fig1_2.rds"))
ggsave(plot = p3, filename = paste0("/2.Stacked Bar Chart for protected areas", "/protected_areas_gene_dissimilarity_Fig1_2.png"),width = 4, height = 6.4, dpi = "retina")

data <- read.table("2.Stacked Bar Chart for protected areas.zip/protected_areas_gene_richness_dissimilarity_Fig1_3.txt", header = T, sep = '\t', check.names = FALSE)
DJ <- melt(data,id.vars='x_axis')
names(DJ)[1:2] <- c("Taxonomy","sample")  
DJ$sample <- factor(DJ$sample, levels = unique(DJ$sample))       
DJ$Taxonomy <- factor(DJ$Taxonomy, levels = unique(DJ$Taxonomy)) 

p4<-ggplot(DJ, aes( x = sample,y = value, fill = Taxonomy))+  
  geom_col(position = 'stack', width = 0.5, show.legend = F)+  #width = 0.5 bar width，show.legend = F remove legend
  scale_y_continuous(expand = c(0,0))+  
  labs(x="",y="Proportion of areas(%)", 
       fill="")+   
  coord_flip()+    
  theme_few()+     
  scale_fill_manual(values=c("#FFFEE0","#C0C000"))  
p4
saveRDS(p4, paste0("/2.Stacked Bar Chart for protected areas", "/protected_areas_gene_richness_dissimilarity_Fig1_3.rds"))
ggsave(plot = p4, filename = paste0("/2.Stacked Bar Chart for protected areas", "/protected_areas_gene_richness_dissimilarity_Fig1_3.png"),width = 4, height = 6.4, dpi = "retina")

#3.boxplot for the relative abundance of microbial traits

library(ggplot2)
data=read.table("3.Boxplot for genes.zip/box_climate_richness_dissimilarity.txt",header=T)
meandf=read.table("3.Boxplot for genes.zip/meandf_climate_richness_dissimilarity.txt",header=T)
pbar_richnness <- ggplot(data,aes(Response,Richness))+
  geom_boxplot(aes(fill=Response),outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.25) +
  geom_point(data = meandf, aes(x = Response, y = Richness), size = 5, color = "cyan", pch = 18) +
  labs(x = "Response",y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))
pbar_richnness
saveRDS(pbar_richnness, paste0("/3.Boxplot for genes", "/gene_richness.rds"))
ggsave(plot = pbar_richnness, filename = paste0("/3.Boxplot for genes", "/gene_richness.png"),width = 10, height = 6.4, dpi = "retina")

data=read.table("3.Boxplot for genes.zip/box_climate_richness_dissimilarity.txt",header=T)
meandf=read.table("3.Boxplot for genes.zip/meandf_climate_richness_dissimilarity.txt",header=T)
pbar_Dissimilarity <- ggplot(data,aes(Response,Dissimilarity1))+
  geom_boxplot(aes(fill=Response),outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.25) +
  geom_point(data = meandf, aes(x = Response, y = Dissimilarity), size = 5, color = "cyan", pch = 18) +
  labs(x = "Response",y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))
pbar_Dissimilarity
saveRDS(pbar_Dissimilarity, paste0("/3.Boxplot for genes", "/gene_Dissimilarity.rds"))
ggsave(plot = pbar_Dissimilarity, filename = paste0("/3.Boxplot for genes", "/gene_Dissimilarity.png"),width = 10, height = 6.4, dpi = "retina")

#drawing the relative abundance of microbial traits across biomes
data=read.table("3.Boxplot for genes.zip/box_plot.txt",header=T)
meandf=read.table("3.Boxplot for genes.zip/meandf.txt",header=T)
pbar_microbial_traits <- ggplot(data,aes(Response,multigenes))+
  geom_boxplot(aes(fill=Response),outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.25) +
  geom_point(data = meandf, aes(x = Response, y = multigenes), size = 5, color = "cyan", pch = 18) +
  labs(x = "Response",y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))
pbar_microbial_traits
saveRDS(pbar_microbial_traits, paste0("/3.Boxplot for genes", "/relative_abundance_of_multigenes.rds"))
ggsave(plot = pbar_microbial_traits, filename = paste0("/3.Boxplot for genes", "/relative_abundance_of_multigenes.png"),width = 10, height = 6.4, dpi = "retina")

#drawing the relative abundance of microbial enzyme profiles across biomes
data=read.table("3.Boxplot for genes.zip/box_plot.txt",header=T)
meandf=read.table("3.Boxplot for genes.zip/meandf.txt",header=T)
pbar_enzyme <- ggplot(data,aes(Response,multienzyme))+
  geom_boxplot(aes(fill=Response),outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.25) +
  geom_point(data = meandf, aes(x = Response, y = multienzyme), size = 5, color = "cyan", pch = 18) +
  labs(x = "Response",y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))
pbar_enzyme
saveRDS(pbar_enzyme, paste0("/3.Boxplot for genes", "/relative_abundance_of_multienzyme.rds"))
ggsave(plot = pbar_enzyme, filename = paste0("/3.Boxplot for genes", "/relative_abundance_of_multienzyme.png"),width = 10, height = 6.4, dpi = "retina")

##drawing the relative abundance of microbial traits across different climate
data=read.table("3.Boxplot for genes.zip/box_climate.txt",header=T)
meandf=read.table("3.Boxplot for genes.zip/meandf_climate.txt",header=T)
pbar1 <- ggplot(data,aes(Response,multigenes))+
  geom_boxplot(aes(fill=Response),outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.25) +
  geom_point(data = meandf, aes(x = Response, y = multigenes), size = 5, color = "cyan", pch = 18) +
  labs(x = "Response",y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))
pbar1
saveRDS(pbar1, paste0("/3.Boxplot for genes", "/relative_abundance_climate_traits.rds"))
ggsave(plot = pbar1, filename = paste0("/3.Boxplot for genes", "/relative_abundance_climate_traits.png"),width = 10, height = 6.4, dpi = "retina")

##drawing the relative abundance of microbial enzyme profiles across different climate
data=read.table("3.Boxplot for genes.zip/box_climate.txt",header=T)
meandf=read.table("3.Boxplot for genes.zip/meandf_climate.txt",header=T)
pbar2 <- ggplot(data,aes(Response,multienzyme))+
  geom_boxplot(aes(fill=Response),outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.25) +
  geom_point(data = meandf, aes(x = Response, y = multienzyme), size = 5, color = "cyan", pch = 18) +
  labs(x = "Response",y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))
pbar2
saveRDS(pbar2, paste0("/3.Boxplot for genes", "/relative_abundance_climate_enzyme.rds"))
ggsave(plot = pbar2, filename = paste0("/3.Boxplot for genes", "/relative_abundance_climate_enzyme.png"),width = 10, height = 6.4, dpi = "retina")

#4.Spearman_correlations and drawing the heatmaps and bubble plot
##Spearman correlations
library(tidyverse)##
library(psych)
library(corrplot)
all <- read.table("4.Spearman_correlations_Bubble_plot.zip/Spearman_env_filter_VPA_to_genes.txt", header = TRUE,row.names = 1,sep = "\t")
abc <- corr.test(all,method="spearman",adjust="BH")##
write.table(t(abc$r),"spearmanr_spearman_correlations_1609.txt",sep="\t")
write.table(t(abc$p),"spearmanp_spearman_correlations_1609.txt",sep="\t")

#From these results, we could find the spearman correlations between env with the gene richnness/dissimilarity, the relative abundance of microbial traits and enzyme profiles
#Drawing the bubble plot
library(ggplot2)
data <- read.table("4.Spearman_correlations_Bubble_plot.zip/bubble_abun_new.txt",header = TRUE)
bubble=ggplot(data, aes(x = y, y = x, size = size, color = category)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(0.1, 9)) + # Adjust the range of bubble sizes
  scale_color_manual(values = c("C" = "pink", "B" = "Grey90", "A" = "cyan"))+
  theme_minimal() +
  labs(
    title = "Bubble Chart Example",
    x = "X Axis Label",
    y = "Y Axis Label",
    size = "Bubble Size",
    color = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    legend.position = "right"
  )
bubble
saveRDS(bubble, paste0("/4.Spearman_correlations_Bubble_plot", "/env_traits.rds"))
ggsave(plot = bubble, filename = paste0("/4.Spearman_correlations_Bubble_plot", "/env_traits.png"),width = 10, height = 6.4, dpi = "retina")

#5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables
library(vegan)
library(MASS)
library(MuMIn)
library(olsrr)
library(performance)
library(car)
library(GGally)
library(tidyverse)
library(Hmisc)

rm(list=ls())

#Co-linearity analysis
varechem<-read.table("5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables.zip/Spearman_env_filter_VPA_to_genes.txt", header = TRUE,row.names = 1,sep = "\t")
varechem=varechem[c(1:19)]
env<-as.matrix(varechem)
co_linear <- varclus(env, similarity="spear") # spearman is the default
co_linear
p1=plot(co_linear)#一般来说当spearman ρ2 > 0.7时则表明分支上的变量可能存在严重的共线性，此时应该根据研究的目的去除其中一个变量。
saveRDS(p1, paste0("/5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables/", "/all.rds"))
ggsave(plot = p1, filename = paste0("/5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables/", "/all.png"),width = 5, height = 7, dpi = "retina")

varechem=varechem[c(2:4,6:11,13,15:16,18)]
env<-as.matrix(varechem)
co_linear <- varclus(env, similarity="spear") # spearman is the default
co_linear
p2=plot(co_linear)
saveRDS(p2, paste0("/5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables/", "/selected_variables.rds"))
ggsave(plot = p2, filename = paste0("/5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables/", "/selected_variables.png"),width = 5, height = 7, dpi = "retina")

#6.VPA_analysis and plot
#We show one example from our datasets
library("vegan")
a=read.table("6.VPA_analysis_plot.zip/Spearman_env_filter_VPA_to_genes.txt",header = TRUE,row.names = 1,sep = "\t")
a=a[1:1529,]
b=a[1:1529,42]
#c=a[c(2:4,6:11,13,15:16,18)]
mod1 <- varpart(b, ~ Elevation+Slope,
                   ~ MAT+PSEA+MDR+TSEA+Aridity,
                   ~ Plant_Cover, 
                   ~ Bulk_Density+Soil_C_g_kg+Soil_CN+Fine_texture+Soil_pH,data=a)
plot(mod1)
mod1

topography <- model.matrix(~ Elevation+Slope, data=a)
aFrac <- rda(b, topography, data=a)
anova(aFrac, step=9999, perm.max=9999)

climate <- model.matrix(~ MAT+PSEA+MDR+TSEA+Aridity, data=a)
aFrac <- rda(b, climate, data=a)
anova(aFrac, step=9999, perm.max=9999)

plant <- model.matrix(~ Plant_Cover, data=a)
aFrac <- rda(b, plant, data=a)
anova(aFrac, step=9999, perm.max=9999)

soil <- model.matrix(~ Bulk_Density+Soil_C_g_kg+Soil_CN+Fine_texture+Soil_pH,data=a)
aFrac <- rda(b, soil, data=a)
anova(aFrac, step=9999, perm.max=9999)

# Drawing the VPA plot
library(RColorBrewer)
library(reshape2)
library(ggplot2)
data<-read.csv("6.VPA_analysis_plot.zip/duiji_new.csv",sep=",",na.strings="NA",stringsAsFactors=FALSE)
data<-melt(data,id.vars='taxa')
sorder = factor(data$variable,levels=unique(data$variable),order=TRUE)
porder = factor(data$taxa,levels=unique(data$taxa),order=TRUE)
windowsFonts(A=windowsFont("Times New Roman"))
P <- ggplot(data=data,aes(x=sorder,y=value,fill=porder)) +
  geom_bar(stat="identity",position="stack") +
  labs(x="Sample Names",y="Relative Abundance (%)",fill="Class",title="funguild")
P+ theme_classic (base_size = 18)+theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,family="A",face="bold",size=3))+
  scale_fill_manual(values = (c( '#F0E442','#0072B2', '#B3DE69','#D55E00', '#CCEBC5', 'grey')))

saveRDS(P, paste0("~/6.VPA_analysis_plot/", "/VPA.rds"))
ggsave(plot = P, filename = paste0("~/6.VPA_analysis_plot/", "/VPA.png"),width = 5, height = 7, dpi = "retina")

#R code for maps
### Part Random forest map of Global distribution of soil microbial genetic resources


# This part of code performs a spatial prediction analysis using a # random forest model and then exports the results to a Gtiff file. #

# It first loads several packages, including the raster, randomForest, # sf and rgdal libraries. The raster data files (in GeoTIFF format) are then loaded using the "raster" function, and their pixel sizes # are changed to 25 km using the "aggregate" function.

# A raster stack is created with the "stack" function and the # coordinates of a data frame are loaded and used to extract values # from the raster stack with the "extract" function.

# Then, the resistance variable "datapoint" is selected from # the data frame "Database " and a database is # created using this variable and the extracted raster values. #

# A random forest model is created using the function "randomForest" # with the variable "datapoint" as response variable and the extracted # raster values as predictor variables. The model is trained using the # database, with 999 trees and repeated 100 times.

# Finally, the "predict" function is used to generate a map of predicted values using the trained random forest model, and the resulting map is written to a GeoTIFF file using the "writeRaster" # function, with the "options" argument specifying that a .tfw file should also be created. #

### Part Random forest map 

#Loading libraries
library(raster)
library(randomForest)
library(sf)
library(rgdal)
#
#Loading available raster data AI <- raster ("path/to/AI.tif")
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\ELEVATION25.tif")->ELEVATION25
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\SLOPE25.tif")->SLOPE25

raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\MAT25.tif")->MAT25
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\MAP25.tif")->MAP25
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\PSEA25.tif")->PSEA25
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\MDR25.tif")->MDR25

raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\PH25.tif")->PH25
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\FINE_TEXTURE25.tif")->FINE_TEXTURE25
raster("\\Users\\usuario\\Desktop\\CAPAS_EMILIO\\CAPAS_25km\\SOC25.tif")->SOC25

#group list
rasters <- list(SLOPE25, ELEVATION25, MAT25,
                MAP25, PSEA25, MDR25, PH25, FINE_TEXTURE25, SOC25)
sng <- stack(rasters)

#Load coordinates and extract values
read.csv(file="\\Users\\usuario\\Desktop\\Kunkun_metagenomic.csv", 
         header=T, sep=";", dec=".")->MEGA_ENV
coordinates <-cbind(MEGA_ENV[,4:3])
d <- extract(sng, MEGA_ENV[,4:3])

#Select resistance of datapoints
datapoint <- MEGA_ENV[,16]
database <- cbind(datapoint=datapoint, d)

#Random Forest 
rfng <-randomForest(datapoint~ . , data=database, importance=TRUE, ntree=999, na.action = na.omit, nrep = 50)
p <- predict(sng, rfng)
rfng
rfng$importance

#Export the map
writeRaster(p,'\\Users\\usuario\\Desktop\\Sporulation.tif',options=c('TFW=YES'))

### Part Map of outliers by Mahalanobis distance.

# This part of the code performs an analysis to identify outliers # in a data set based on Mahalanobis distances, and then exports # the results to a file.

# This code imports the dismo library and selects a subset of the data frame "Data_for_Model_stressors1" containing the coordinates # in the second and first columns. Next, the "m" object is created # using the "mahal" function of the dismo package, which calculates
# the Mahalanobis distances between the "sng" object (which is a set # of pixels in space) and the "points" object (which is the set of  reference points).

# The resulting Mahalanobis distances are stored in the object # "mahatotal" using the function "predict" with object "m" as argument. Finally, the "mahatotal" object is written to a GeoTIFF # file named "outliers_stressors1.tif" using the "writeRaster" function, with the "options" argument specifying that a # .tfw file should also be created.


### Part Map of outliers by Mahalanobis distance

#Load library
library(dismo)

#Load coordinates and calculate Mahalanobis distances
points <- subset(MEGA_ENV[,4:3])
m<- mahal(sng, points)

#Store distances
mahatotal <- predict(m, sng)

#Export the map
writeRaster(mahatotal,'\\Users\\usuario\\Desktop\\Mahalanobis_kunkun.tif',options=c('TFW=YES'))

beep(1)

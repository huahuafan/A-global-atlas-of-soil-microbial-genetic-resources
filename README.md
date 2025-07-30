# A-global-atlas-of-soil-microbial-genetic-resources
This R code includes the R packages (including versions) and code used in this manuscript.
We also include dataset to demo the code, and also provide instructions to run the code on my data.
In detail, the R code incude the following codes:
# 1. calculating the Jaccard distance
# examples for the 200 sampling plot
# 2.Stacked Bar Chart for protected areas across different climate regions and vegetation types
# 3.boxplot for the relative abundance of microbial traits
#drawing the relative abundance of microbial traits across biomes
#drawing the relative abundance of microbial enzyme profiles across biomes
##drawing the relative abundance of microbial traits across different climate
##drawing the relative abundance of microbial enzyme profiles across different climate
# 4.Spearman_correlations and drawing the heatmaps and bubble plot
#From these results, we could find the spearman correlations between env with the gene richnness/dissimilarity, the relative abundance of microbial traits and enzyme profiles
# 5.Hierarchical clustering to evaluate the inter-dependencies between pairwise variables
# 6.VPA_analysis and plot
#We show one example from our datasets
#R code for maps
# 7. Part Random forest map of Global distribution of soil microbial genetic resources
# This part of code performs a spatial prediction analysis using a # random forest model and then exports the results to a Gtiff file. #

# It first loads several packages, including the raster, randomForest, # sf and rgdal libraries. The raster data files (in GeoTIFF format) are then loaded using the "raster" function, and their pixel sizes # are changed to 25 km using the "aggregate" function.

# A raster stack is created with the "stack" function and the # coordinates of a data frame are loaded and used to extract values # from the raster stack with the "extract" function.

# Then, the resistance variable "datapoint" is selected from # the data frame "Database " and a database is # created using this variable and the extracted raster values. #

# A random forest model is created using the function "randomForest" # with the variable "datapoint" as response variable and the extracted # raster values as predictor variables. The model is trained using the # database, with 999 trees and repeated 100 times.

# Finally, the "predict" function is used to generate a map of predicted values using the trained random forest model, and the resulting map is written to a GeoTIFF file using the "writeRaster" # function, with the "options" argument specifying that a .tfw file should also be created. #

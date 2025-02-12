##### Title: Microbial Community Composition Analysis in Inflammatory Bowel Disease (IBD) #####

----------------------------------------------------------------------------------------------
# Data Overview
# -------------
# Data Source: 
#   - The data is used to study the role of gut bacteria in Inflammatory Bowel Disease (IBD).
#
# Files Provided:
#   1. `ibd_taxa`:
#      - Contains bacterial abundances (taxon x sample matrix).
#      - Each row represents a bacterial taxon, and each column represents a sample.
#   2. `ibd_lineages`:
#      - Contains taxonomic lineages of the bacteria (lineage x taxon matrix).
#      - Provides hierarchical classification (e.g., phylum, class, order, family, genus, species).
#   3. `ibd_metadata`:
#      - Contains sample properties (metadata x sample matrix).
#      - Includes clinical and demographic information (e.g., diagnosis, age, antibiotic use).
#
# - File Format:
#   - Files are provided in CSV format.
#   - Browsers like Chrome or Edge may change the extension to `.xls`. If this happens, adjust the script to read `.xls` instead of `.csv`.
#
# ----------------------------------------------------------------------------------------------

## 1: Loading Required Libraries and Data

# Install and load required packages
install.packages('ggplot2')  # For advanced plotting
install.packages('vegan')    # For ecological diversity analysis
install.packages('reshape2') # For data reshaping
install.packages('devtools') # For installing GitHub packages
library(ggplot2)
library(vegan)
library(reshape2)
library(devtools)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomformat")  # For handling BIOM files
BiocManager::install("phyloseq")    # For microbiome analysis
library(biomformat)
library(phyloseq)

# Install GitHub packages
install_github("SpiecEasi")  # For microbial network analysis
library(SpiecEasi)

# Set working directory
setwd("~/Microbiome Analysis/Bowel Disease (IBD)")

# Create Output Directory
output_dir <- "Analysis_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load data
ibd_taxa <- read.csv('ibd_taxa.xls', row.names = 1)      # Taxon x sample matrix
ibd_lineages <- read.csv('ibd_lineages.xls', row.names = 1) # Lineage x taxon matrix
ibd_metadata <- read.csv('ibd_metadata.xls', row.names = 1) # Metadata x sample matrix

# Save Loaded Data to Output Directory
write.csv(ibd_taxa, file.path(output_dir, "ibd_taxa_loaded.csv"))
write.csv(ibd_lineages, file.path(output_dir, "ibd_lineages_loaded.csv"))
write.csv(ibd_metadata, file.path(output_dir, "ibd_metadata_loaded.csv"))

## 2: Microbial Community Composition Visualization

# Generate colors for phyla
colours <- rainbow(length(unique(ibd_lineages$Phylum)))
color.indices <- match(ibd_lineages$Phylum, unique(ibd_lineages$Phylum))
colvec <- colours[color.indices]

# Stacked bar plot
barplot(as.matrix(ibd_taxa), col = colvec, xlab = "Samples", ylab = "Relative Abundance")
legend('topright', fill = colours, legend = unique(ibd_lineages$Phylum), cex = 0.8)
dev.off()

## 3: Alpha Diversity Analysis

# Define Shannon Diversity function
shannon_div <- function(x) {
  x <- x[x > 0]  # Remove zeros
  p <- x / sum(x)  # Calculate proportions
  H <- -sum(p * log(p))  # Shannon index
  return(H)
}

# Calculate Shannon diversity for each sample
ibd_metadata$shannon <- apply(ibd_taxa, 2, shannon_div)

# Boxplot of Shannon diversity by diagnosis
ggplot(ibd_metadata, aes(x = Diagnosis, y = shannon, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity by Diagnosis", y = "Shannon Diversity Index", x = "Diagnosis") +
  theme_minimal()
dev.off()

# Save Shannon Diversity Results
write.csv(ibd_metadata, file.path(output_dir, "shannon_diversity_results.csv"))

## 4: Beta Diversity Analysis

# Perform PCoA
pcoa.res.ibd <- capscale(t(ibd_taxa) ~ 1, distance = 'bray', na.action = na.omit)
eigenvalues <- eigenvals(pcoa.res.ibd)

# Plot PCoA
colvec <- rainbow(length(unique(ibd_metadata$Diagnosis)))[as.numeric(as.factor(ibd_metadata$Diagnosis))]
plot(pcoa.res.ibd$CA$u[, c(1, 2)], col = colvec,
     xlab = paste('PCoA1 (', round(eigenvalues[1], 2), ')', sep = ''),
     ylab = paste('PCoA2 (', round(eigenvalues[2], 2), ')', sep = ''),
     main = "PCoA of Microbial Communities")

legend("topleft", legend = unique(ibd_metadata$Diagnosis), col = unique(colvec), pch = 1)
dev.off()

# Save PCoA Results
write.csv(pcoa.res.ibd$CA$u, file.path(output_dir, "pcoa_results.csv"))

## 5: Differential Abundance Analysis

# Find index of Faecalibacterium prausnitzii
fp.index <- which(rownames(ibd_taxa) == "Faecalibacterium_prausnitzii")

# Perform Wilcoxon test
wilcox_result <- wilcox.test(
  t(ibd_taxa[fp.index, ibd_metadata$Diagnosis == "Control"]),
  t(ibd_taxa[fp.index, ibd_metadata$Diagnosis %in% c("CD", "UC")])
)
print(wilcox_result)

# Save Wilcoxon Test Results
sink(file.path(output_dir, "wilcoxon_test_results.txt"))
print(wilcox_result)
sink()

# Visualize distributions
hist(t(ibd_taxa[fp.index, ibd_metadata$Diagnosis == "Control"]), 
     breaks = "FD", col = rgb(0, 1, 0, 0.5), xlab = "Abundance", main = "Faecalibacterium prausnitzii")
hist(t(ibd_taxa[fp.index, ibd_metadata$Diagnosis %in% c("CD", "UC")]), 
     breaks = "FD", col = rgb(1, 0, 0, 0.5), add = TRUE)
legend("topright", legend = c("Control", "IBD"), fill = c("green", "red"))
dev.off()

## 6: Correlation Analysis

# Correlation between Shannon diversity and Age
cor_age <- cor(ibd_metadata$shannon, ibd_metadata$Age, method = "spearman")
print(paste("Correlation between Shannon diversity and Age:", cor_age))

# Save Correlation Results
sink(file.path(output_dir, "correlation_results.txt"))
print(paste("Correlation between Shannon diversity and Age:", cor_age))
sink()

# Scatter Plot
png(file.path(output_dir, "shannon_vs_age_scatterplot.png"), width = 600, height = 400)
plot(shannon ~ Age, data = ibd_metadata, ylab = "Shannon Diversity", xlab = "Age", main = "Shannon Diversity vs Age")
dev.off()

## 7: Machine Learning for Disease Prediction

# Install and load randomForest package
install.packages("randomForest")
library(randomForest)

# Prepare data
X <- t(ibd_taxa)  # Transpose taxa matrix
Y <- as.factor(ibd_metadata$Diagnosis)  # Diagnosis as labels

# Train random forest model
rf_model <- randomForest(X, Y, ntree = 500, importance = TRUE)

# Save Random Forest Model
saveRDS(rf_model, file.path(output_dir, "random_forest_model.rds"))

# Evaluate model
sink(file.path(output_dir, "random_forest_results.txt"))
print(rf_model) # Plot variable importance
sink()

## 8: Network Analysis

# Use SpiecEasi for network analysis
spiec_net <- spiec.easi(t(ibd_taxa), method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20)

# Load the igraph package for network visualization
install.packages("igraph")  # Install if not already installed
library(igraph)
library(SpiecEasi)

# Assuming `spiec_net` is your SpiecEasi network object
opt_beta <- as.matrix(getRefit(spiec_net))  # Use the correct function

# Extract the adjacency matrix from the SpiecEasi output
adj_matrix <- as.matrix(getOptBeta(spiec_net))  # Get the optimal beta matrix
adj_matrix <- as.matrix(adj_matrix)  # Ensure it's a matrix
adj_matrix[adj_matrix != 0] <- 1  # Binarize the matrix (optional, depending on your needs)

# Create an igraph object from the adjacency matrix
network <- graph.adjacency(adj_matrix, mode = "undirected", weighted = TRUE)

# Simplify the network (remove loops and multiple edges)
network <- simplify(network)

# Plot the network
plot(network, 
     vertex.size = 3,          # Size of nodes
     vertex.label.cex = 0.7,   # Size of labels
     vertex.color = "lightblue",  # Node color
     edge.color = "gray",      # Edge color
     layout = layout_with_fr(network))  # Layout algorithm (Fruchterman-Reingold)
dev.off()

# Save Adjacency Matrix
write.csv(adj_matrix, file.path(output_dir, "adjacency_matrix.csv"))

##### Save the Entire Script #####
script_name <- "Microbial_Community_Analysis_Script.R"
file.copy("Microbial_Community_Analysis_Script.R", file.path(output_dir, script_name))

# Print Completion Message
print(paste("Analysis complete! All outputs saved to:", output_dir))

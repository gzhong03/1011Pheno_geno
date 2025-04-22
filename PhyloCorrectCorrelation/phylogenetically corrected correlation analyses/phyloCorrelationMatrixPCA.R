library(ape)
library(phytools)
library(corrplot)
library(ggplot2)

# Load tree and data
mytree <- read.tree("S.c.tree")
mytree <- multi2di(mytree)
mydata <- read.csv("1000GenomeS.cerevisiaeMulticellularMasterForTree.csv", header = TRUE)
rownames(mydata) <- mydata$Standardized_name

# Select only the traits of interest
mydat2 <- mydata[, c("Inv", "PSH", "LD", "Let", "Flor", "Adh")]

# Drop tree tips not in data
mytree <- drop.tip(mytree, setdiff(mytree$tip.label, rownames(mydat2)))

# Match data order to tree
mydat2 <- mydat2[match(mytree$tip.label, rownames(mydat2)), ]

# Apply Felsenstein’s Phylogenetically Independent Contrasts (PIC) to each trait
pic_matrix <- apply(mydat2, 2, function(trait) pic(trait, mytree))

# Compute phylogenetically corrected correlation matrix
phylo_corr_matrix <- cor(pic_matrix, use = "pairwise.complete.obs")
print(round(phylo_corr_matrix, 3))

#visualize plot
corrplot(phylo_corr_matrix,
         method = "color",       # Use color shading
         type = "upper",          # Show full matrix, not just upper
         addCoef.col = "maroon",  # Display correlation values
         number.cex = 0.7,       # Size of the numbers
         tl.cex = 0.8,           # Size of trait labels
         tl.col = "black",       # Label color
         diag = TRUE)            # Show diagonal (can set FALSE if you want it blank)

#PCA
pic_matrix_clean <- na.omit(pic_matrix)

pca_pic <- prcomp(pic_matrix_clean, scale. = TRUE)
summary(pca_pic)

# Convert scores to data frame
pca_df <- as.data.frame(pca_pic$x)

# Add strain names if needed
pca_df$Strain <- rownames(pca_df)

# Plot first two PCs
# Scores (strain positions in PCA space)
pca_df <- as.data.frame(pca_pic$x)
pca_df$Strain <- rownames(pca_df)

# Loadings (how much each phenotype contributes to each PC)
loadings <- as.data.frame(pca_pic$rotation)
loadings$Phenotype <- rownames(loadings)

# Scale the loadings so arrows fit nicely in plot space
arrow_scale <- 5  # You can tweak this value if arrows are too long/short
loadings_scaled <- loadings
loadings_scaled[,1:2] <- loadings[,1:2] * arrow_scale

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "firebrick") +
  geom_text(data = loadings_scaled,
            aes(x = PC1, y = PC2, label = Phenotype),
            color = "firebrick", vjust = -0.5, size = 4) +
  theme_minimal() +
  labs(title = "PCA of PIC-Corrected Traits",
       x = paste0("PC1 (", round(summary(pca_pic)$importance[2, 1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_pic)$importance[2, 2]*100, 1), "%)")) +
xlim(-3, 7) +      # <- adjust these values as needed
ylim(-3, 4)

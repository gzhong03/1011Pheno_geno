setwd("/Users/helldawg/Helen/Work/Projects/1000genomes/Tree")

#This rscript is based on: https://www.zoology.ubc.ca/~schluter/R/Phylogenetic.html#:~:text=in%20trait%20%22x%22-,Independent%20contrasts,data%20points%20resulting%20from%20phylogeny.

#Get the libraries
library(ape)
library(phytools)
library(nlme)
library(visreg)

#Get the tree
mytree <- read.tree("S.c.tree")

#Make it dichotomous (necessary for ape to calculate things); not necessary for other packages
mytree=multi2di(mytree)

#If we want any info on the tree
#mytree               # prints basic tree information
#mytree$tip.label     # species names
#mytree$edge.length   # lengths of all tree branches
#mytree$edge          # identity of branches (from node, to node)

mytree$edge.length        # lengths of all tree branches
range(mytree$edge.length) # shortest and longest branch


#Look at the tree
plot(mytree, cex=0.7)

mydata <-read.csv("1000GenomeS.cerevisiaeMulticellularMasterForTree.csv", header=TRUE)

#Set the row names to be the strain names
row.names(mydata) <- mydata$Standardized_name

#If we want to prune the tree to only have the strains for which there is data
mytree<-drop.tip(mytree, mydata$Missing)

#Re-order the dataframe of phenotypes to match the strain order in the tree
mydata <- mydata[match(mytree$tip.label,rownames(mydata)),]

#Make a dataframe that contains only the phenotypes of interest
mydat2 <-mydata[, c("Inv", "PSH", "LD", "Let", "Flor", "Adh")]

#Test for a phylogenetic signal using Pagel's Lambda (ranges from 0-1, with 1 being maximum phylogenetic signal)
phylosig(mytree, as.matrix(mydat2)[,c("Flor")], test = TRUE)


#FELSENSTEIN PHYLOGENETICALLY INDEPENDENT CONTRASTS
#This is a method that estimates and tests a regression or correlation between two variables 
#while correcting for the nonindependence of data points resulting from phylogeny. 
#It assumes that trait evolution mimics a Brownian motion process with unchanging rates 
#through time and along all branches of the tree. This is a difficult assumption to verify.

#First, use the pic() function to transform the two variables x and y separately into 
#phylogenetically independent contrasts (PICs), x1 and y1. 
# Keep only complete cases for both traits
# Step 1: Find complete cases
complete_cases <- complete.cases(mydat2$LD, mydat2$PSH)

# Step 2: Subset phenotype data
LD_clean <- mydat2$LD[complete_cases]
PSH_clean <- mydat2$PSH[complete_cases]

# Step 3: Subset the tree to match the filtered strains
filtered_strains <- rownames(mydat2)[complete_cases]
tree_clean <- drop.tip(mytree, setdiff(mytree$tip.label, filtered_strains))

# Step 4: Make sure the data order matches the tree
LD_clean <- LD_clean[match(tree_clean$tip.label, filtered_strains)]
PSH_clean <- PSH_clean[match(tree_clean$tip.label, filtered_strains)]

# Step 5: Compute PICs
LD_pic <- pic(LD_clean, tree_clean)
PSH_pic <- pic(PSH_clean, tree_clean)

# Step 6: Combine and clean PICs
pic_df <- data.frame(LD = LD_pic, PSH = PSH_pic)
pic_df_clean <- na.omit(pic_df)

# Step 7: Compute phylogenetically corrected correlation
r <- cor(pic_df_clean$LD, pic_df_clean$PSH)
print(r)


#GLM LINEAR REGRESSION

#Get the correlation matrix, which is based on the evolutionary distances
v_matrix <-vcv(mytree, corr = TRUE)

# Ensure phenotypes have no NA for traits used
complete_cases <- complete.cases(mydat2$LD, mydat2$Let)

# Filter tree and data
filtered_strains <- rownames(mydat2)[complete_cases]
tree_clean <- drop.tip(mytree, setdiff(mytree$tip.label, filtered_strains))
mydata_clean <- mydata[filtered_strains, ]  # Use full data frame (not just mydat2)

# Reorder rows to match tip labels
mydata_clean <- mydata_clean[match(tree_clean$tip.label, rownames(mydata_clean)), ]

# Get the covariance matrix again
v_matrix <- vcv(tree_clean, corr = TRUE)

#Model the two variables. By calling the corBrownian correlation structure
z <- gls(PSH ~ LD, data = mydata_clean, correlation = corBrownian(1, tree_clean))
z <- gls(Inv ~ PSH, data = mydata_clean, correlation = corBrownian(1, tree_clean))
summary(z)
confint(z)

#-----------------------------------------------------------------------------
z <- lm.gls(x=mydat2$LD, y=mydat2$Let, v=v_matrix)

#Test the function here
z <- lm (mydat2$Let ~ mydat2$LD)
X <-model.matrix(z)
colnames(X)[1] <- "intercept"
colnames(X)[2] <- "LD"
ve<- eigen (v_matrix)
#Had to add this to make one very small negative value a positive value
ve[["values"]]=abs(ve[["values"]])
pv <-ve$vectors %*% sqrt(diag(ve$values)) %*% t(ve$vectors)
X <- as.data.frame(solve(pv) %*% X)
ynew <-solve(pv) %*% mydat2$Let

#colnames(X, c("intercept", "LD"))
f <-as.formula(paste("ynew ~", paste(names(X), collapse = " + "), "-1"))
lm.fit(ynew ~ intercept + mydat2$LD, data = X)

lm.fit(ynew ~ X$intercept + X$LD -1)

lm.gls <- function(x, y, v){
  # Function to implement general least squares
  # Assumes no missing values
  
  # x is a univariate vector or a matrix of variables
  # y is a vector, the response variable
  # v is the phylogenetic covariance matrix
  
  # There is a transformation of the original x and y variables
  # that would fulfill the assumptions of uncorrelated residuals having
  # equal variances (see Draper and Smith p 108 for this logic in the
  # case of weighted least squares for uncorrelated data; p156 for
  # the case with serially correlated data). The transformation (P)
  # is obtained by solving P'P=V as follows:
  
  # Generate the design matrix X, including column of 1's for intercept
  z <- lm(y ~ x)
  X <- model.matrix(z)
  colnames(X)[1] <- "intercept"
  
  
  # Calculate pv, the "square root" of v
  ve <- eigen(v)
  pv <- ve$vectors %*% sqrt(diag(ve$values)) %*% t(ve$vectors)
  
  # Transform y and X by pv
  X <- as.data.frame(solve(pv) %*% X)
  ynew <- solve(pv) %*% y
  
  # Fit the gls model, corrected for phylogeny
  f <- as.formula(paste("ynew ~", paste(names(X), collapse = " + "), "-1"))
  lm.fit <- lm(f, data = X)
  
  # Collect results
  yhat <- predict(lm.fit)
  resid <- residuals(lm.fit)  
  result <- data.frame(X, ynew, yhat, resid)
  
  return(list(lm.fit = lm.fit, result = result))
}

library(ggtree)
ggtree(tr=mytree, layout="circular", options(ignore.negative.edge=TRUE))

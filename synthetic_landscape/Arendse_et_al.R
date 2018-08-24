# Matthew Arendse
# 3440197

# BCB743: Quantitative Ecology practical (96 Marks)
# A fun practical conducted to explore whether meaningful patterns between species dispersal ability
# and environmental variables could be visulaised and quantified.

# Load Packages:

library(tidyverse)
library(ggpubr)
library(vegan)
library(readr)
library(gclus)
library(FactoMineR)

# Load additionnal functions
source("Numerical_Ecology_with_R/All Data/evplot.R")
source("Numerical_Ecology_with_R/All Data/hcoplot.R")
source("Numerical_Ecology_with_R/All Data/triplot.rda.R")
source("Numerical_Ecology_with_R/All Data/plot.lda.R")
source("Numerical_Ecology_with_R/All Data/cleanplot.pca.R")
source("Numerical_Ecology_with_R/All Data/PCA.newr.R")
source("Numerical_Ecology_with_R/All Data/CA.newr.R")
source("Numerical_Ecology_with_R/All Data/coldiss.R")

# Load Datasets

spe <- read_csv("spe.csv")
spe_mass <- read_csv("spe_mass.csv")
spa <- read_csv("spa.csv")
env <- read_csv("env.csv")

# Question 1 (2 Marks) - Dimensions of the datasets------------------

# Give dimensions

dim(env)
# This data set has 28 rows and 3 columns 
dim(spe)
# This data set has 28 rows and 17 columns 
dim(spe_mass)
# This data set has 28 rows and 17 columns 
dim(spa)
# This data set has 28 rows and 2 columns 

# Question 2 (6 Marks) - Spatial Plot ---------------------------------------------------------

# (i) Spatial co-ords of sites and transects

spa1 <- spa %>% 
  group_by(x = 50)
spa2 <- spa %>% 
  group_by(x = 100)
spa3 <- spa %>% 
  group_by(x = 150)
spa4 <- spa %>% 
  group_by(x = 200)

dev.new(title="Site Locations")
plot(spa, type="n", main="Site Locations", xlab="x 	coordinate (cm)", ylab="y coordinate (cm)")
lines(spa1, col="light blue")
lines(spa2, col="light blue")
lines(spa3, col="light blue")
lines(spa4, col="light blue")
text(spa, row.names(spa), cex=0.8, col="red")
text(50, 25, "Transect 1", cex=0.6, col="red")
text(100, 25, "Transect 2", cex=0.6, col="red")
text(150, 25, "Transect 3", cex=0.6, col="red")
text(200, 25, "Transect 4", cex=0.6, col="red")

# These lines of code provide a plot visualising the way in which our experimental landscape was set out.
# We see that our landscape was 2m by 3.6m. To sample our landscape, 4 transect lines were set out at 50cm intervals
# along the bredth of the landscape. Each transect ran the length of the landscape (360cm) and consisted of 7 sites.
# Sites were spaced 60cm from eachother starting at 0 and ending at 360cm.
# A 30x30 cm quadrat was placed at each site, and was used to sample the number of individuals present at each site.


# To compute the number of species at each site
# sum by rows
sit.pres <- apply(spe > 0, 1, sum)
sit.pres
# This gives us an indication of the total number of species present at each site.

# Sort the results in increasing order
sort(sit.pres)


# (i a) Species Richness Map

dev.new(title="Species Richness", width=10, height=5)
plot(spa, main="Map of Species Richness", pch=21, col="white", 
     bg="red", cex=5*sit.pres/max(sit.pres), xlab="x coordinate (cm)", 
     ylab="y coordinate (cm)")
lines(spa1, col="light blue")
lines(spa2, col="light blue")
lines(spa3, col="light blue")
lines(spa4, col="light blue")
text(spa, row.names(spa), cex=0.8,pos = 2, col="black")
text(50, 25, "Transect 1", cex=0.6, col="black")
text(100, 25, "Transect 2", cex=0.6, col="black")
text(150, 25, "Transect 3", cex=0.6, col="black")
text(200, 25, "Transect 4", cex=0.6, col="black")

# These few lines of code provide a plot visualising our experimental landscape.
# Here bubbles have been plotted at each site to provide a representation of species richness at each site.
# Larger bubbles indicate sites with a greater species richness. 
# Here we see that site 5 had the most species (11) followed by site 19 (8) and then site 11 (7).


# (i b) Shannon-Weaver Map

H <- diversity(spe)            # Shannon entropy
N1 <- exp(H) # Shannon diversity (number of abundant species)
N1
# This line of code provides an indication diversity of species at each site.


plot(spa, main="Map of Shannon diversity", pch=21, col="white", 
     bg="red", cex=5*N1/max(N1), xlab="x coordinate (cm)", 
     ylab="y coordinate (cm)")
lines(spa1, col="light blue")
lines(spa2, col="light blue")
lines(spa3, col="light blue")
lines(spa4, col="light blue")
text(spa, row.names(spa), cex=0.8,pos = 2, col="black")
text(50, 25, "Transect 1", cex=0.6, col="black")
text(100, 25, "Transect 2", cex=0.6, col="black")
text(150, 25, "Transect 3", cex=0.6, col="black")
text(200, 25, "Transect 4", cex=0.6, col="black")

# Here we are provided with a plot of the Shannon diversity at each site.
# Once again, site 5 appears to have the highest diversity, followed by site 19 and site 11.

# (i c) Simpsons Map
?diversity
N2 <- diversity(spe, index = "simpson")  # Simpson diversity (number of dominant species)
N2

dev.new(title="Simpsons Diversity", width=10, height=5)
plot(spa, main="Map of Simpson diversity", pch=21, col="white", 
     bg="red", cex=(N2*4), xlab="x coordinate (cm)", 
     ylab="y coordinate (cm)")
lines(spa1, col="light blue")
lines(spa2, col="light blue")
lines(spa3, col="light blue")
lines(spa4, col="light blue")
text(spa, row.names(spa), cex=0.8,pos = 2, col="black")
text(50, 25, "Transect 1", cex=0.6, col="black")
text(100, 25, "Transect 2", cex=0.6, col="black")
text(150, 25, "Transect 3", cex=0.6, col="black")
text(200, 25, "Transect 4", cex=0.6, col="black")

# Here we a provided with a plot that displays the Simpsons diversity (1-D) per site.
# A value closer to 1 is interpreted as having a relativley high diversity. Once again, our sites which showed a relativley high
# species richness (site 5, 11 and 19) also show a relativley high diversity.

# Question 3 (8 Marks) - Association Matrix--------------------------------------------------------

# (i) Association Matrix for Species Data
# Transpose the dataset
spe.pa <- decostand(spe, method="pa")


spe.t <- t(spe.pa)
spe.t[1:17, 1:28]

spe.jac.AM <- vegdist(spe.t, "jac", binary=TRUE)
spe.jac.AM
round(as.matrix(spe.jac.AM)[1:16, 1:15], 3)

# (ii) Plot the Matrix
# Magenta = dissimilarity close to 0 (max similarity)
# Cyan = dissimilarity close to 1 (min similarity)
dev.new(title="Jaccard", width=10, height=5)
coldiss(spe.jac.AM, byrank=FALSE, diag=TRUE)

# (iii) Briefly describe the obvious positive and negative species associations.

# The plotted heat map shows that sites containing big crumpled (b_c) 
# squares and b_c triangles & b_c half (half circles) and b_c squares
# show strong similarity.
# Sites such as those containing gums and pink_el show high dissimilarity.

# These results suggest that b_c triangles and b_c squares show similar dispersal
# ability across the landscape and as such end up occupying very similar sites.

# Similarity Matrix (Non-Binary)
spe.t2 <- t(spe)
spe.t2[1:17, 1:28]

spe.bray.AM <- vegdist(spe.t, "bray", binary=FALSE)
spe.bray.AM
round(as.matrix(spe.jac.AM)[1:16, 1:15], 3)

# (ii) Plot the Matrix
# Magenta = dissimilarity close to 0 (max similarity)
# Cyan = dissimilarity close to 1 (min similarity)
dev.new(title="Bray-Curtis", width=10, height=5)
coldiss(spe.bray.AM, byrank=FALSE, diag=TRUE)

# The Jaccard Presence/Absence association matrix produces a simmilar visual output to the Bray-Curtis output,
# but produces a more discriminatory result with fewer species similarities being present. However, the major similarites persist
# across both methods.













# Question 4 (30 Marks) - PCA----------------------------------------------

# PCA on species data

# Hellinger pre-transformation of the species data
spe.h <- decostand(spe.pa, "hellinger")
(spe.h.pca <- rda(spe.h))
summary(spe.h.pca)
# the total iniertia is the sum of the diagonal 
# of the correlation matrix that feeds into the PCA

# The eigenvalues indicate the importance of each of the PC axes, 
# with the first PC axis always having the highest eigenvalue. 
# Values diminish until the smallest one is found for the last PC axis.
# The sum of all the eigenvalues is equal to the total inertia.

# Here we can see that the first 4 PC axes explain 61% of the variation

# PCA biplots
source("Numerical_Ecology_with_R/All Data/cleanplot.pca.R")
spe.rda.sc1 <- scores(spe, display="species", scaling=1)
spe.rda.sc2 <- scores(spe, display="species", scaling=2)

par(mfrow = c(2, 2))
biplot(spe.h.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.h.pca, main = "PCA - scaling 2 ('species')") # default scaling 2
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04) # species scaling


#PC 1 V 2
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.01,ax1 = 1,ax2 = 2) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 1,ax2 = 2) 

#PC 1 V 3
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.01,ax1 = 1,ax2 = 3) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 1,ax2 = 3) 

#PC 1 V 4
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 1,ax2 = 4) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 1,ax2 = 4) 

#PC 1 V 5
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 1,ax2 = 5) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 1,ax2 = 5) 

#PC 1 V 6
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 1,ax2 = 6) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 1,ax2 = 6) 

#PC 2 V 3
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 2,ax2 = 3) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 2,ax2 = 3) 

#PC 2 V 4
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 2,ax2 = 4) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 2,ax2 = 4) 

#PC 2 V 5
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 2,ax2 = 5) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 2,ax2 = 5) 

#PC 2 V 6
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 2,ax2 = 6) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 2,ax2 = 6) 

#PC 3 V 4
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 3,ax2 = 4) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 3,ax2 = 4) 

#PC 3 V 5
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 3,ax2 = 5) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 3,ax2 = 5) 

#PC 3 V 6
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 3,ax2 = 6) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 3,ax2 = 6) 

#PC 4 V 5
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 4,ax2 = 5) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 4,ax2 = 5) 

#PC 4 V 6
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 4,ax2 = 6) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 4,ax2 = 6) 

#PC 5 V 6
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.001,ax1 = 5,ax2 = 6) # sites scaling
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.04,ax1 = 5,ax2 = 6) 


# Patterns evident in biplots:

# 3 main patterns are shown by the biplots (PC1 v PC2):
# (1) site 8 and 26 are similar based on the presence of the beans species.
# (2) sites 3, 18, 24 etc are similar based on the presence of the small_rectangle
# species at those sites.
# (3) sites 19, 20, 21 etc are similar in that these sites share many similar species.
# this pattern suggests that those species share similar dispersal ability, and were able
# to occupy similar sites, whereas beans were only able to disperse effectivley to sites 8 and 26
# and small_rectangles to sites 3, 18, 24 etc.

# How many axes must be retained, and why?

# Given that the first 3 PC axes only explain 51% of the total variation,
# it would be advisable to retain at least the first 4 PC axes (61%)
# but ideally you'd want to retain the first 6 PC axes as these explain >75% of the variation.
# By retaining more of the PC axes and plotting them against one another (PC1 v PC3, PC1 v PC4 etc)
# we get a better understanding of how our species are distributed among the various sites.


# Hellinger pre-transformation of the species mass data
spe.m.h <- decostand(spe_mass, "hellinger")
spe.m.pca <- rda(spe.m.h)
summary(spe.m.pca)

# PCA biplots
source("Numerical_Ecology_with_R/All Data/cleanplot.pca.R")
spe.rda.sc1 <- scores(spe, display="species", scaling=1)
spe.rda.sc2 <- scores(spe, display="species", scaling=2)

par(mfrow = c(1, 2))
biplot(spe.m.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.m.pca, main = "PCA - scaling 2 ('species')") # default scaling 2
par(mfrow = c(1, 2))
cleanplot.pca(spe.m.pca, scaling = 1, mar.percent = 0.001) # sites scaling
cleanplot.pca(spe.m.pca, scaling = 2, mar.percent = 0.04) # species scaling

#Here the biplots show the association between species mass and the sites that species were most likely associated with as a product of their mass.
# Here we observe strong associations by the pink_ele , beans, and s_rectangle species with their respective sites.
# This would suggest that the mass of species played a key role in determining which sites they were able to disperse to.

# Comparing Scaling 1 for spe data with Scaling 1 for spe_mass data
dev.new(title="Comparing species data (left) with species mass data (right)", width=10, height=5)
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.01,ax1 = 1,ax2 = 2) # Species data
cleanplot.pca(spe.m.pca, scaling = 1, mar.percent = 0.001) # Species mass data

# Question 5 (12 Marks) - NMDS-------------------------------------------------

# NMDS applied to the QE_Prac species - percentage difference
# dissimilarity matrix

spe.nm <- spe[-16,]
spe.nm <- spe.nm[-26:-27,]

spe.nmds <- metaMDS(spe.nm, distance="bray")
spe.nmds

spe.nmds.j <- metaMDS(spe.nm, distance="jaccard")
spe.nmds.j

#Plot (bray)
spe.nmds$stress
dev.new(title="NMDS on QE_Prac species - Percentage difference")
plot(spe.nmds, type="t", main=paste("NMDS/Percentage difference - Stress =", round(spe.nmds$stress,3)))
# Shepard plot and goodness of fit (Bray)
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe.nmds, main="Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type="t", main="Goodness of fit")
points(spe.nmds, display="sites", cex=gof*300)


#Plot (jacc)
spe.nmds.j$stress
dev.new(title="NMDS on QE_Prac species - Percentage difference")
plot(spe.nmds.j, type="t", main=paste("NMDS/Percentage difference - Stress =", round(spe.nmds.j$stress,3)))
# Shepard plot and goodness of fit (Jaccard)
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe.nmds.j, main="Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds.j, type="t", main="Goodness of fit")
points(spe.nmds.j, display="sites", cex=gof*300)

# Here, provided with the results from the Bray Curtis NMDS analysis, there is a clear
# pattern present. In both the NMDS/Percentage difference and Goodness of Fit plots, we see that
# both the beans and the big crumpled (b_c) rectangle species are stongly associated with
# sites 9, 10, and 25 and site 7 respectivley. The rest of our species have seemingly been able to distribute
# themselves fairly evenly across the remaining sites. This pattern is reflected in the goodness
# of fit plot, given poorly fitted sites have larger bubbles, we see that the beans species for exapmle has been well fitted to site 9.
# This type of analysis plends itself as another tool to be used by ecologists when exploring the relationship between species and sites.
# From these analysis we can begin to ask questions as: What variables have contributed to the beans species being so strongly associated with site 9,
# whilst the rest of the species have distributed relativley evenly across sites.


# Question 6 (40 Marks) - Constrained RDA analysis---------------------------------

# Hellinger-transform the species dataset
spe.hel <- decostand(spe, "hellinger")
spe.hel
# Redundancy analysis (RDA)
# *************************

# RDA of the Hellinger-transformed QE_Prac species data, constrained
# by all the environmental variables contained in env
(spe.rda <- rda(spe.hel ~ ., env)) # Observe the shortcut
# formula
summary(spe.rda)	# Scaling 2 (default)


#(i) - Summary(spe.rda):

# the summary output of our constrained RDA analysis shows the following:
# (1) our constrained vaariables only account for 7% of the total inertia whilst the unconstrained variables
# account for 92% of the total Inertia.
# (2) we see that the cumulative contribution of the constrained axes (RDA1 & RDA2) again only accounts for 
# 7% of the total variance.
# (3) Once again 6 principle component (PC) axes would need to be retained for the analysis to explain at least 75% of the variation.
# (4) When only considering the constrained axes (RDA1 & RDA2), RDA1 accounts for 79% of the total variance in terms of 
# the accumulated constrained eigenvalues.


# Canonical coefficients from the rda object
coef(spe.rda)
# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(spe.rda)$r.squared)
# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

# Scaling 1: distance biplot
dev.new(title="RDA scaling 1 + wa")
par(mfrow=c(1,2))
plot(spe.rda, scaling=1, 
     main="Biplot RDA spe.hel ~ env - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")
# Scaling 2 (default): correlation biplot
plot(spe.rda, main="Biplot RDA spe.hel ~ env - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")
## Biplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
plot(spe.rda, scaling=1,
     main="Biplot RDA spe.hel ~ env - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")
# Scaling 2
plot(spe.rda, 
     main="Biplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length=0, lty=1, col="red")

#(ii) Biplot Interpretation:

# RDA spe~env wa and lc scores Scaling 1 plots:
# in the above mentioned plots we see that of the three environmental variables (carpet, cardboard, and wood)
# carpet facilitated dispersal far  more than that of cardboard and wood. The large majority of our species show
# strong associations with the carpet variable, whilst only our "pink elephant", "half circle", and "gums" species show some association
# with the cardboard variable.
# Based on the Scaling 2 plots, we see that site (row) 18 is strongly associated with cardboard and is a reflection of the sampling that took place 
# on top of our incline. The rest of the sites seem to be more strongly associated with the carpet variable, with sites 17 most likely being associated with 
# the wooden incline itself that lead up to the cardboard at site 18.
# Overall the species associations with sites and specific environmental variables seems to reflect the distribution patterns 
# observed throughout our practical experiment.

# Global test of the RDA result
anova(spe.rda, permutations=how(nperm=999))
# with a significance score of f = 0.9 we accept H0 and state that our environmental variables
# had no significant influence on the dispersal ability of our organisms across our landscape.

# Tests of all canonical axes
anova(spe.rda, by="axis", permutations=how(nperm=999))
# with a significance score of f = 1.5 and f = 0.3 we can confirm that neither of the 
# RDA axes had a significant influence on the variation explaining the distribution of our species across our landscape.

anova(spe.rda, by="terms", permutations=how(nperm=999))
# with a significance score of f = 1.5 and f = 0.3 we can confirm that neither of the environmental
#variables (Carpet/Cardboard) had a significant influence on the dispersal ability of our organisms across the 
#landscape.

# Variance inflation factors (VIF)
vif.cca(spe.rda)
# Based on the VIF scores of our environmental variables ( x < 10 ) we would discard both 
# variables as neither of them contribute to the variance explaining the dispersal ability of 
# our species across the landscape.



# Redundancy analysis (RDA) - spe_mass
# *************************
spe.m.hel <- decostand(spe_mass, "hellinger")
spe.m.hel


# RDA of the Hellinger-transformed QE_Prac species data, constrained
# by all the environmental variables contained in env
(spe.m.rda <- rda(spe.m.hel ~ ., env)) # Observe the shortcut
# formula
summary(spe.m.rda)	# Scaling 2 (default)


# Scaling 1: distance biplot (mass)
dev.new(title="RDA scaling 1 + wa")
par(mfrow=c(2,2))
plot(spe.m.rda, scaling=1, 
     main="Biplot RDA spe.m.hel ~ env - scaling 1 - wa scores")
spe.sc1 <- scores(spe.m.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")
# Scaling 2 (default): correlation biplot (mass)
plot(spe.m.rda, main="Biplot RDA spe.hel ~ env - scaling 2 - wa scores")
spe.sc2 <- scores(spe.m.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")
# Scaling 1: species data
plot(spe.rda, scaling=1, 
     main="Biplot RDA spe.hel ~ env - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")
# Scaling 2 (default): correlation biplot
plot(spe.rda, main="Biplot RDA spe.hel ~ env - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")


#(ii) Biplot Interpretation:

# RDA spe_mass~env vs spe~env Scaling 1 plots:
# when comparing the scaling 1 biplot for our  RDA of species mass by environmental variables and RDA of just species by environmental variables.
# In both biplots, we see that beans have been strongly associated with those sites that were on the carpet.
# However, when comparing the those sites and species that have been associated with the cardboard environment some differences become apparent.
# When taking species mass into account we see that in the spe_mass scaling 2 plot that its only site 18 and the gums species that are associated with the cardboard environment.
# This is contrasted with the fact that in the spe scaling 2 site 18 and 3 other sites had been associated with the carpet environment, as well as gums and 3 other species.
# These results suggest that species mass may have an influence on which sites species were able to disperse to.

# Question 7 (10 Marks)---------------------------------------------------

# Using the Doubs data as an example, the above set of analysis produced plots that helped researchers
# in establishing transition zones along the river, each with a unique indicator species and environmental variables
# that contributed to the occupation of each zone by a given species. The Doubs data set provides an excellent example of 
# how useful these types of analysis can be in visualising  how environmental variables can influence the occurence and success 
# of various species at various sites along a river. These types of analysis are incredibly useful to ecologists in that it 
# allows them to narrow in on specific patterns seen in the analysis that lend themselves well to further
# quantitative analysis whereby the influence of environmental variables can prove to be statistically significant in influencing
# how species are distributed.

# Our experiment was an attempt to assess how the dispersal ability, across a landscape, 
# of various species was affected by their shape (flat shapes, crumpled shapes, sweets and foam peanuts)
# and the envrionmental variables (carpet, wood, cardboard) of the landscape. 
# Although the experiment worked, the data collected did not
# lend itslef well to this type of analysis and proved to be unsuccessful in highlighting any meaningful patterns.
# In the Doubs data, clear associations could be drawn between site locations, levels of pollution at those sites, and tolerant/hardy
# species occuring at those sites - Logical ecological conclusions could be drawn with regard to the environmenta variables and species present at each site. 
# Here however, the affinity of most of our species to be associated
# with carpet is more likely a result of the fact that our landscape was made up of at least 85% carpet rather than carpet facilitating/hindering the dispersal ability of any given species.
# Had we placed a fan at a set entry point to our landscape, had set introduction periods for each species, and accounted for wind speed, we would have most likely
# observed more of a gradient in terms of how species had been dispersed. 
# Heavier species would not have made it as far from the point of entry compared to lighter species,
# and through the above analysis more meaningful patterns of dispersal ability would have arisen as a product of wind speed. 
# We would then have been able to quatify these patterns as a product of species mass and wind speed at the time of entry into the landscape,
# and as such  would have been able to draw 
# meaningful associations between mass, wind speed, and dispersal ability of various species.

# Based on the above analysis, my understanding of the ecology of our habitat is that it is simply a random assortment of species dispersed across a landscape.
# Conducting the experiment over many trials would not result in a more significant result being produced.
# In order to improve the experiment, stricter environmental variables and conditions need to be set so as to allow for more dispersal patterns to present themselves.
# By "allowing" certain species to play to their strengths, rather than giving each a fair chance of dispersal, we are more likely to see dispersal patterns emerging.
# These emergent patterns can then be explored in light of various environmental (wind, substrate, etc) or species specific (mass, size,etc) variables and the
# magnitude and significance of their influence quantified.

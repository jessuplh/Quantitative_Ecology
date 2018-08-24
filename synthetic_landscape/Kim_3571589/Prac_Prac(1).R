# Quantitative Ecology Practical Assigment
# 16 August 2018
# Kim Scholtz

##load neccessary libraries
library(tidyverse)
library(ggpubr)
library(vegan)

## load data
env <- read_csv("env.csv")
spe <- read_csv("spe.csv")
spa <- read_csv("spa.csv")
spe_mass <- read_csv("spe_mass.csv")

#####

# Question 1 (2 marks)
# Specify the dimensions of the three data sets.

# The dimensions are given as rows x columns.
dim(env)
dim(spe)
dim(spa)

# Question 2 (6 marks)
# i.	Provide a plot on Cartesian (spatial) coordinates of a) the species richness,
# the b) Shannon-Weaver, and c) the Simpson’s Index.

# Site location maps 
mutate_spa <- spa %>% 
  mutate(N = 1:28)

### Plotting the site location
Figure1 <- ggplot(mutate_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_text(aes(label = N), colour = "blue4") +
  labs(title = "Site Location", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
Figure1

# Compute the number of species at each site
# To sum by rows, the second argument of apply(), MARGIN, is set to 1
sit.pres <- apply(spe > 0, 1, sum)

site_df <- data.frame(spa = sit.pres)

# Results to be sorted in increasing order

sit.pres_df <- data.frame(species_richness = sit.pres)

species_richness <- sit.pres_df %>% 
  mutate(N = 1:28)

# Species richness on spatial data
Figure2 <-  ggplot(mutate_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=species_richness$species_richness), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Map showing Species Richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
Figure2

# Shannon  index 
Shanny <- as.tibble(diversity(spe, index = "shannon"))

Figure3 <- ggplot(mutate_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=Shanny$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Shannon-Weaver Diversity Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
Figure3

# Simpson index
Simpson <- as.tibble(diversity(spe, index = "simpson"))

Figure4 <- ggplot(mutate_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=Simpson$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Simpson Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
Figure4 

### Plotting the univariate measures of diversity
plot <- ggarrange(Figure1, Figure2, Figure3, Figure4)
###

# ii.	Explain each figure using your knowledge of how the ‘landscape was created’.

# Man introduced a range of different species into the environment (made up of a 
# carpet floor, wood and cardboard structures as well as a stainless steel bowl). 
# Individual species were introduced via wind dispersal using a man-made fan board.
# When looking at the species richness graph most of the sites have fewer than 3 
# species present. The highest number of species  were present at site 3. 
# At sites 27 and 28 no species were present (this being the point furthest  from 
# the point of dispersal). Wind being the environmental factor has played a role in the aggregate of species in the middle region of the study site.
# Gum and jelly species were able to disperse further from their point of release as they are heavier
# than the paper species and the effects of the wind has strongly influenced there distribution 
# and richness in the various sites under study. The lowest species richness was seen
# at sites 27 and 28. This could be as a result of the light weight of the paper as well as
# the intensity of the wind and the point of dispersal.In the Shanon-weaver index
# the eveness of species diversity is higher in the middle of the study area. This 
# could be as a result of the mass of the species as well as point of dispersal and
# wind intensity.  Sites 1, 2, 22 and 29 were of the sites that had a lower evenness.  
# In the Simpsons index species diversity is greater towards the middle and upper 
# regions of the environment than in the lower regions (site 1 and 2). Site 6 and 7 
# has the higher species diversity. Site 7 is dominated by paper species whereas site
# 6 is dominated by the gum and jelly species. They are able to dispersed further away from their point of release as a result of having a greater mass.
# The general trend shown by the figures are that diversity increases as species richness increase.

# Question 3 (8 marks)

# Load the required packages
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(FD)

# Load additionnal functions
# (files must be in the working directory)
# source("coldiss.R")
# source("panelutils.R")

# i.Calculate an association matrix for the species data. (1)
spe.t <- t(spe)
spe.t
spe.t.S7 <- vegdist(spe.t)
round(as.matrix(spe.t.S7), 2)

# ii.Produce a plot that is suitable for showing the kinds of information contained 
# in the association matrix. (3)

## Percentage difference (Bray-Curtis) dissimilarity matrix on raw species abundance data
dev.new(title="Percentage difference (Bray-Curtis), raw data", width=10, height=5)
coldiss(spe.t.S7, byrank=FALSE, diag=TRUE)
# Same but on log-transformed data

dev.new(title="Percentage difference (Bray-Curtis), ln(y+1) data", width=10, height=5)
coldiss(spe.t.S7, byrank=FALSE, diag=TRUE)

# iii.Briefly describe a few of the more ‘obvious’ positive and negative species
# associations. If possible, add an ‘ecological’ explanation for why this is so. (4)

# Positive associations between two species can occur when both select the same
# habitat or have the same environmental requirements. Conversely, negative
# associations can occur if the species have differing ecological requirements.
# Similarites are seen between flattened shapes and crumbled shapes irrespective of there
# weight or shape. One of the weakest association is seen between the gum species and any of the other species.
# Gum species are one of the heaviest species (in mass) and thus its weight could have had
# an effect on its abilty to disperse. The major diver of dispersal was wind.
# A strong correlation is evident in both the squares, triangles and crumbled squares 
# and crumbled triangles. In addition, a correlation could be seen between small triangles 
# and big triangles. Big pink elephants are correlated with half circle. 
# A 82% dissimilarity is seen between small rectangles and big rectangles. 

# Question 4 (20)
# i.	Undertake a PCA on the species data. Comprehensively explain the various (all)
# components of the summary() of the PCA object. (20)

# Load the required packages
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(FactoMineR)

# Load additionnal functions

source("evplot.R")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

spe.pca <- rda(spe, scale=TRUE) 
summary(spe.pca) # Default scaling 2

# The principal components are the linear combinations of the original variables that
# account for the variance in the data. The maximum number of components extracted 
# always equals the number of variables. The eigenvectors, which are comprised of 
# coefficients corresponding to each variable, are used to calculate the principal 
# component scores. The coefficients indicate the relative weight of each variable
# in the component.In these results, the first three principal components ( PC1-3.791, PC2-2.4752, PC3-2.124) 
# have eigenvalues greater than 1. These three components explain 49.3% of the variation
# in the data. If 49.3% is an adequate amount of variation explained in the data
# then you should use the first three principal components. The proportion is used to
# determine which principal components explain most of the variability in the data.
# The higher the proportion, the more variability that the principal component
# explains. The size of the proportion can help you decide whether the principal 
# component is important enough to retain. A principal component (PC7) with a proportion 
# of 0.062 explains 62.15% of the variability in the data. Therefore, this component 
# important to include. Another component (PC1) has a proportion of 0.223, and thus 
# explains only 0.2% of the variability in the data. This component may not be important enough to include.
# The cumulative proportion can help you determine the number of principal components
# to use. Retain the principal components that explain an acceptable level 
# of variance. For example, you may only need 80% of the variance explained by the principal
# components if you are only using them for descriptive purposes. However, if you want to 
# perform other analyses on the data, you may want to have at least 90% of the 
# variance explained by the principal components. Scores are linear combinations of the data that are determined by the coefficients 
# for each principal component. Species scores explains 12.7% (small triangle) of PC1,
# 75% explained for big triangles whilst for gumsthe species score explains -7%.
# Most of the site scores have negative values for most of the principal components,
# with the exception of certain PCA's (PC6-2.46% at site8 , PC4-31.46 at site 9% and
# PC1-14,95% at site 6).

# ii.	Provide plots of the PCA. How many axes must be retained, and why?
par(mfrow = c(2, 2))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") # default scaling 2
cleanplot.pca(spe.pca, scaling = 1, mar.percent = 0.08) # sites scaling
cleanplot.pca(spe.pca, scaling = 2, mar.percent = 0.04) # species scaling


# ii. What patterns become visible from the biplot? Explain your findings in detail. (10)

#  PC1-PC6 axes must be retained as their cummalative proportion explained adds up 
# to 75.96%.PCA Axis 1 to be on the X-axis and PCA Axis 2 to be on the Y-axis to 
# get the diagram above. Positively correlated variables are grouped together.
# Negatively correlated variables are positioned on opposite sides of the plot 
# origin (opposed quadrants).
# The angle between the arrows of the environmental variables can be used to
# approximate the correlations among those variables in the scaling focused on
# species correlations. Angles between the species arrows in a diagram from a linear
# ordination method. The angle between the species suggests that they are almost
# non-correlated.Species are present only at certain sites. Many species have significant abundances at 
# sites where arrows extend further than the circle. Sites scaling (scaling 1) is used in the
# figure to the left. Sites are scaled by eigenvalues. The distances in multidimensional space are more accurately reflected
# on the graph plane so it shows the relationships between sites better. The figure on
# the left show how closely together sites aggregate in the left side of the graph, while 
# fewer sites are seen in the right part of the graph. Sites are more scattered in the
# right of the figure than the left. The figure shows how sites sharing similar species plot closer together and the ones
# that are further apart on the graph are also more dissimilar in their actual multidimensional space. 
# The circle of equilibrium on the bottom left panel shows that species (gums, big half circles
# and pink elephant) influence the spread of species across sites. The presence of species are generally associated with the spread of sites, although
# their influence is not perpectly alligned with the direction of spread in most sites aside from
# at site 4 and 18. These patterns of association between species at sites could indicate some degree of
# collinearity. In the figure to the bottom left (scaling 2) species present at the various sites 
# in the bottom left is low and no correlation appears to be evident between the species. 


# Question 5 (10)
# i.	Undertake an NMDS on the species data, and provide plots of the NMDS. (2)
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(FactoMineR)

# Load additionnal functions
# (files must be in the working directory)
source("evplot.R")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

spe <- read_csv("~/SMIT HONOURS/PRACTICAL/spe.csv")
spe_1 <- spe[-16,]
r_species <- spe_1[-26:-27,]

spe.nmds <- metaMDS(r_species, distance="bray")
spe.nmds
spe.nmds$stress
dev.new(title="NMDS on the different species - Percentage difference")
plot(spe.nmds, type="t", main=paste("NMDS/Percentage difference - Stress =", round(spe.nmds$stress,3)))

# Shepard plot and Goodness of fit
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe.nmds, main="Shepard plot")
G <- goodness(spe.nmds)
plot(spe.nmds, type="t", main="Goodness of fit")
points(spe.nmds, display="sites", cex=gof*300)

# ii.	What patterns become visible from the ordination plot? Explain your findings in detail. (10)

# Shepard diagram shows the relationship between NMDS ordination distance and
# original observed distance.Shepard diagram shows that the observed dissimilarities and the 
# ordination distances were up to 98% correlated (P,0.0001,Figure to the left). The 
# Shepard stress plot shows the relationship between the actual dissimilarities
# between species (from the original dissimilarity matrix) and the ordination distances (i.e. the distances on the final plot). Species that are 
# ordinated closer to one another are likely to be more similar than those further 
# apart.In the diagram there is a large amount of scatter indicating a poor linear 
# relationship. Therefore, the ordination is not representative of the original distances. 
# Occasionally, specific objects may be ordinated poorly despite the overall solution 
# being acceptable.
# In the Goodness fit figure to the right, sites that are poorly fitted have larger
# bubbles. As seen in the ordination plot it is clear that most of the sites if not 
# all are poorly fitted. The site with the smallest circle is site 10, which may be
# a better fit than the rest of the sites,  but it is still not a very strong fit.

# Question 6 (40)
# i.Do a constrained analysis on the data, and comprehensively explain all 
# the findings as per the of the summary() of the NMDS object. (20)
spe.hel <- decostand(spe, "hellinger")
(spe.rda <- rda(spe.hel ~ ., env)) # Observe the shortcut
                                    # formula
summary(spe.rda)	# Scaling 2 (default)
# RDA computes axes that are linear combinations of the explanatory variables.
# In RDA the axes explain the variation of the response matrix.
# The overall variance is partitioned into constrained and unconstrained
# fractions. The constrained fraction (0.051) is the amount of variance of the Y matrix.
# The unconstrained fraction yielded a value of 0.6551.
# This analysis yielded 2 canonical axes (with eigenvalues labelled RDA1 to RDA2) 
# and 16 additional,unconstrained axes for the residuals (with eigenvalues labelled
# PC1 to PC16). The results give the eigenvalues themselves, as well as the cumulative proportion
# of variance explained (for the RDA axes) or represented (for the residual axes).
# The last cumulative value is therefore 1. The cumulative contribution to
# the variance obtained by the 2 canonical axes is the proportion of the total 
# variance of the response data explained by the RDA. It is the same value as the
# “Proportion constrained” presented above (0.0726). The canonical eigenvalues RDA1 to RDA2 increase in value. The first
# residual eigenvalue (PC1), on the contrary, is larger than the last canonical
# eigenvalue and is  larger than most RDA eigenvalues. This means that the first residual structure (axis) of the data has more variance
# than some of the structures that can be explained by the explanatory variables in X.
# The eigenvalue for RDA1 = 0.04097 and RDA2 = 0.01037.  The results give the eigenvalues themselves, as well as the cumulative proportion 
# of variance. 

# ii.	Run the necessary permutation tests. Explain the outcomes. (6)

# Global test of the RDA result
anova(spe.rda, permutations=how(nperm=999))
# Tests of all canonical axes
anova(spe.rda, by="term", permutations=how(nperm=999))
# Variance inflation factors (VIF)
vif.cca(spe.rda)

# Hypothesis 0: The response (species richness) is independent on the site.
# Hypothesis 1: The response (species richness) is dependent on the site.
# For each of the tests there was no significant difference as seen by the
# p value(p >0.05). Tests of significance showed that there was no significant difference
# between matrices X and Y after 999 permutations. The  are all 
# P>0.05 for all of the canonical axes and shows a weak species environment 
# correlation.The canonical axes explain 10% and 40% of the variance. 
# Important variables were missing from the data (such as taken into account wind intensity)
# This couldve influenced the outcome of the tests.
# Global test of the RDA results indicate that the model had a P-value > 0.05,df = 2 and 
# variance = 0.051. P value indicates that there is no significant difference.
# While the residual had a df = 25 with a variance of 0.6551. 
# Tests of all canonical axes showed that the carpet had a df = 1, 
# variance = 0.0408, P value = 0.080.
# Cardboard had a df = 1, variance = 0.0104 and p value = 0.956.
# Residual had a df = 25 and variance = 0.6551.
# The variance inflation factors (VIF)for both carpet and cardboard structures 
# yielded a variance of 1.9285, whilst wood had no variance value (N/A).

# iii. Produce the necessary biplots. Explain the findings. (14)
# Scaling 1: distance triplot
dev.new(title="RDA scaling 1 + wa") 
par(mfrow = c(2, 2))
plot(spe.rda, scaling=1, 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2 (default): correlation triplot
# dev.new(title="RDA scaling 2 + wa")
plot(spe.rda, main="Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
# dev.new(title="RDA scaling 1 + lc")
plot(spe.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2
# dev.new(title="RDA scaling 2 + lc")
plot(spe.rda, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length=0, lty=1, col="red")

# On the bottom left the figure illustrates how most species are present on the carpet.
# A small proportion of species are found on the cardboard region but none of the species ubder
# study are evident on the wood area. The graph clearly shows the presence of only two 
# environmental variables rather than three. Thus we can deduce that the variable that is not
# present (wood) does not have a significant influence on the variability within the study area.
# The position at which species are found  relative to the environmental variables 
# shows how strong the association between a species and a particular environmental
# variable is. The environmental variables that were able to significantly explain variability 
# in the environmental structure were fitted to the ordinations . Hence, only 
# two (carpet and cardboard) of the three environmental variables were present. 
# We can therefore deduce that wood has no significant influence on the variability
# within the study area. The different species present at a specific site 
# (n = 28) is represented by the points on the biplot.The length of the arrow can be
# related to the rate of change.  The figure on the left (Scaling 1) shows how the rate of change of the environmental variable
#(carpet) is higher than the rate of change of the cardboard. 
# For scaling 1 species that were found closer together were expected to have similar
# variable values. The biplot for scaling 1 shows how more species are correlated with 
# the carpet and fewer species are correlated with the cardboard. When looking at the figure to the bottom right (scaling 2) the obtuse angle between 
# carpet and cardboard suggests that there is a strong negative correlation between
# the two variables.


# Question 7 (10)
# Write down your understanding of the ‘ecology’ of the landscape that was sampled 
# and analysed the various steps, above.

# Ecology underpins the study of living organisms (plants and animals)in the environment
# and their physical surroundings. Lanscape ecology focuses explicitly on the ecological
# understanding of spatial heterogeneity.The landscape in this study consisted of cardboard,
# wood and carpet. The landscape that was sampled was flat and spacious to sustain a wide range of species.
# Species were randomly distributed amoung sites and between carpet , wood and cardboard 
# environments. The carpet was an ideal environment in which to do sampling as 
# there were no other species present on landscape prior to the introduction of the species under study. 
# The only steep area was the wood and carboard, which was arranged in the form of a mountain.
# The wood, cardboard and bowl (environmental variables) influences the distribution
# of the different species as it acted as a barrier  seperating species from eachother and isolating
# some species from the other. The area beneath the object created a microenvironment for some species.
# Species were dispersed by wind onto the land. There was also a stainless steel bowl which prevented 
# species from inhabiting that area. This object also acted as a barrier and species werent able to occupy that space.
# Certain species were able to disperse further than others and this dispersal rate was influenced by the strength of the wind as well as the wood and cardboard acting 
# as a barrier, and the weight of the species. Multivariate analysis is based on the
# statistical principle of multivariate statistics, which involves observation and 
# analysis of more than one statistical outcome variable at a time.

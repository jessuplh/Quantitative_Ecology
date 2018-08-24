# Quantitative Ecology
# Practical R assignment
# Shaundré Hofstander
# Due 24 August 2018
###########################################################################


# Let's begin -------------------------------------------------------------

# load libraries
library(tidyverse)
library(ggpubr)
library(vegan)
library(gclus)
library(cluster)


# Load additionnal functions
# (files must be in the working directory)
source("coldiss.R")
source("panelutils.R")
source("cleanplot.pca.R")
source("evplot.R")
source("PCA.newr.R")
source("CA.newr.R")

# load data

env <- read_csv("env.csv")
spa <- read_csv("spa.csv")
spe <- read_csv("spe.csv")
spe_mass <- read_csv("spe_mass.csv")



# Question 1 --------------------------------------------------------------


dim(env)
dim(spa)
dim(spe)
dim(spe_mass)

# Question 2  -------------------------------------------------------------

# i) 

# Number of species at each site
# create site location maps 
mut_spa <- spa %>% 
  mutate(N = 1:28)

fig1 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_text(aes(label = N), colour = "blue4") +
  labs(title = "Site Location", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()

# Compute the number of species at each site
# To sum by rows, the second argument of apply(), MARGIN, is set to 1
sit.pres <- apply(spe > 0, 1, sum)

site_df <- data.frame(spa = sit.pres)

# Sort the results in increasing order
# sort_sit.pres <- sort(sit.pres$apply.spe...0..1..sum.)
sit.pres_df <- data.frame(spe_rich = sit.pres)
# sit.pres_df$spe_rich <- as.numeric(sit.pres_df$spe_rich)

spe_richness <- sit.pres_df %>% 
  mutate(N = 1:28)

fig2 <- ggplot(spe_richness, aes(x = N, y = spe_rich)) +
  geom_step() +
  geom_text(aes(label = N), colour = "red") +
  labs( x = "Positions of sites", y = "Species Richness", 
        title = "Species Richness") +
  theme_minimal()

# species richness on spatial
fig3 <-  ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=spe_richness$spe_rich), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Map of Species Richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()



# species mass data
spe_mass_mutate <- spe_mass %>% 
  mutate(N = 1:28)

try1 <-  ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=spe_mass_mutate$spe_rich), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Map of Species Richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()


# shannon-weaver 
H <- as.tibble(diversity(spe, index = "shannon"))

fig4 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=H$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Shannon-Weaver Diversity Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()


# simpson index
simp <- as.tibble(diversity(spe, index = "simpson"))

fig5 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=simp$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Simpson Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()

ggarrange(fig2, fig3, fig4, fig5)


# ii)
# The species richness and map of species richness plot are used to show the
# abundance of each species. There were two main dispersal points, in which wind was the vector. Sites 4 & 5, 11 & 12, and 18 & 19, were the sites found in or close to the centre and had the highest species richness. These sites were 1 and 7. The species richness of each species differed, because of their random dispersal points as well as their physiological differences. These sites were the dispersal points for the gummy and paper species. This serves as a possible explanantion for the high species richness in the centre of the study rgion, with the wind acting as a force of dispersal. The physiological traits of the gummy and bean species, in this case specifically heavier mass  un comparison to the paper species, allowed them to disperse independently of the wind. 
# There are less than three different species present within most of the
# sites.  The landscape of sites 17 and 18 comprised of wood and cardboard, which formed an incline
# the species richness here are of an average standard. Another site of dispersal was site 5, and it had the greatest species richness of all the sites. The sites 27 and 28 found on the outskirts of the study region, had the lowest species richness as it did not comprise of any species.
# This is possibly due to these sites being a greater distance from the dispersal points, in comparison with the other sites and with wind as the vector and the physiological constraints of some species, it made it difficult for species to disperse to these sites. 
# The species from the dispersal points, generally dispersed to the centre of the study region, this can is an explanation as to why species did not disperse to these sites. 

# Shannon - Weaver - The Shannon index is used to show the species evenness in an environment. #The sites found in the centre of the study region had the greatest species richness, this is due to the great species diversity and species richness in these sites. This is possibly due to the points of dispersal 
# On the other hand, sites located the furthest from the centre of the study region, generally found on the 
# outskirts, had lower species richness and species diversity, resulting in a lower evenness. These sites were 1 & 2, 27 & 29, and 22 & 23.

# Simpson's index - This index focuses on species richness and the 
# abundance of relative species. The species diversity and species richness were closely correlated. 
# The species diversity was highest in the middle of the study region and had
# species is greatest in the centre of the study area, the areas possessing the 
# highest species richness. Sites 1 and 2 had the lowest diversity of species,
# the gummy and jelly species (species with the greater mass) and shared site 1 as a dispersal point. #However, sites 6 and 7 have a higher species diversity, where site 7 was a dispersal point for the 
# paper species that had little to no mass. This allowed the wind to facilitate their dispersal further from their dispersal point. The general trend shows that the species diversity increases with the increase in species diversity and species evenness


# Question 3 --------------------------------------------------------------
# i)
spe.t <- t(spe)
spe.t[1:17, 1:17]
spe.t.S7 <- vegdist(spe.t, binary = FALSE) 
round(as.matrix(spe.t.S7), 2)

# or
cor(spe)

# ii)
coldiss(spe.t.S7, byrank=FALSE, diag=TRUE)

# Those with dissimilarities of 0 represents maximum similarity, and are individual of the same species. The large and small paper individuals of the same species, had  
# dissimilarities that varied. This was seen with the small and large triangle that had a 
# 47 % dissimilarity. There was a 82% dissimilarity for the small and large rectangle. 
# The individuals of the crumbled and uncrumbled species of the same shape, had greater
# dissimilarities. This was evident in the 100 % dissimilarity of the small crumbled half circle and a small half circle, as well as an 87% dissimilarity presented by the small crumbled triangle and a small triangle.
# There was great dissimilarity between the crumbled shapes and
# the gums , with dissimilarities of 1. This was seen with the dissimilarity of 1 between the
# small crumbled triangle and the gums, and the big crumbled square and the gums.
# In terms of ecology, individuals of different species have a greater dissimilarity value than individuals of similar species.


# Question 4 --------------------------------------------------------------

# i) 
spe.pca <- rda(spe, scale=TRUE)
spe.pca 
summary(spe.pca) 
summary(spe.pca, scaling=1)

# Total inertia -  is the sum of the diagonal of the correlation matrix that gets fed into the PCA. PC1 #explains about 22% of the  variance. PC6, including the other principal components, explains #approximately 75.96%  of the existing variance. A great inertia total shows that the variance is high. #The inertia is not related to the species abundance variance, but rather to species scores variance.
# Species scores – represented by the position of the arrow heads on the PCA biplot. Gives an #indication of the both the magnitude and direction of the species in association with a gradient linear #to the landscape.
# Sites scores - most influential variables that form part of the PC1 and show the
# coordinates of the sites. Here PC1 shows that sites 4, 5, 11, 12, and 19 share the most similarities. 
# Eigenvalues - these are measures of the variance of
# the axes. It can be expressed as Proportions Explained, or proportions of
# variation, when divided by the total inertia.
# Eigenvalues shows the level of importance of the PC axes, first PC axis has the highest #eigenvalue, and descends to the last PC axis, which has the lowest eigenvalue. The total inertia is #equivalent to the sum of all the eigenvalues.
# Eigenvectors -  represent the 'loadings' of the original, untransformed environmental variables and #corresponds to a particular eigenvalue that is found on the PC axis. Furthermore, indicates the #original variable’s importance along the new reduced dimension. 


# ii
par(mfrow = c(2, 2))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") 
cleanplot.pca(spe.pca, scaling = 1, mar.percent = 0.08) 
cleanplot.pca(spe.pca, scaling = 2, mar.percent = 0.04)

# ii)
# In scaling 1, scales the sites via eigenvalues. The distances in this scaling reflect multidimensional 
# space more accurately. This is best used to show the relationships between sites more accurately. 
# Sites that are plotted closer together are more similar in terms of environmental characteristic, 
# than those spaced further apart. Scaling 1 focuses on the at sites while, scaling 2 focuses on species 
# biplots: 
# The circle of equilibrium on the plot indicates that certain species, such as gums, big half circles and pink elephants affect the spread of sites. 
# The presence of species are generally associated with the spread of sites
# The movement to the left of the vertical zero line, usually associates the distribution of sites with the species that are present. This brings to light that a significant variable for the distribution of sites is missing, however the species present can be explained by the distribution of sites.
# The axes that must be retained are PC1-PC6, as the PC1-PC6 cumulative proportion explains approximately 75.96%. 

# Scaling 1: from left to right, the first group consisting of sites 6, 14, 17, 24, and 25 
# with greatest abundance of small half rectangle and gums. Moving away from these sites, 
# results in the decrease of these species and the lowest abundances for all the other species. 
# The sites in the top left quadrant are not associated with any of the species. 
# A group consisting of sites 11, 12 and 19 have a high abundance for many species independent 
# of the size and type of species. These sites also have the lowest abundance in the gummy species. 
# The big crumbled half circle small crumbled, half circle and big square species are most abundant in sites 4 and 5. 

# Scaling 2: there are no specie present in the group found in the top left quadrant. 
# The bottom left has a very strong and positive correlation with the small half circles and gum species. 
# However, they have a very strong and negative correlation with the group in top right of the plot. 
# The top right quadrant contains all the species that are well correlated with one another. 
# The top right well correlated with the bottom right quadrant. The orthogonal arrows of the 
# distribution of the species in the right quadrants of the plot show a correlation close to 0. 
# This can be seen with the big half circles and big crumbled squares. 
# Small half circles, big square and small triangle have short arrows in comparison with others, 
# showing it is less significant in terms of the sites ordination.
# PCA is one of the many ways to analyse a correlation matrix. 
# The first principal axis maximizes the variance that is reflected by its eigenvalue, 
# and the second one is orthogonal to the first principal axis, and continues to maximize the 
# remaining variance. PCA is based on eigenvectors. The eigenvalues gives an indication of the 
# amount of the original variance that is captured with regards to the corresponding eigenvector. 
# Furthermore giving an indication of the importance of each eigenvector. 
# The first axis (PC1) has the greatest explanatory power, as it has the greatest eigenvalue of 3.791.


# Question 5 --------------------------------------------------------------

# Exclude rows without data
spe <- read_csv("spe.csv")

spe.1 <- spe[-16,]

remove_spe <- spe.1[-26:-27,]

spe_prac.nmds <- metaMDS(remove_spe, distance="bray")
spe_prac.nmds
spe_prac.nmds$stress
dev.new(title="NMDS on species - Percentage difference")
plot(spe_prac.nmds, type = "t", main = paste("NMDS/Percentage difference", round(spe_prac.nmds$stress, 3)))

# Shepard plot and goodness of fit 
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe_prac.nmds, main="Shepard plot")
gof <- goodness(spe_prac.nmds)
plot(spe_prac.nmds, type="t", main="Goodness of fit")
points(spe_prac.nmds, display="sites", cex=gof*300)


# ii)
# In the stress plot, only the site scores are present,as the dissimilarities did not have the necessary 
# information about the species. It also shows two correlation like statistics of goodness of fit.

# Shepard diagram - shows the relationship between NMDS ordination distance and the actual
# distance that was observed. It plots the distances of the ordination against the dissimilarities of the 
# community. The red monotone step line represents the fit. This diagram shows that about 98% of the 
# observed dissimilarities and the ordination distances was correlated (P = 0.0001). The 
# the species are more similar to one another, when they are ordinated closer together.
# The level of scattering in this diagram is high, therefore a poor linear relationship exists.
# Thus, the ordination is not a good representative of the actual observed distances. 

# The goodness of fit of the transformation function can be characterized by a 
# residual sum of squares (RSS), just as in regression.
# In the goodness plot, sites with larger bubbles are poorly fitted than those with smaller bubbles.
# This shows that majority of the sites if not are poorly fitted. Although not it is not a very strong fit, 
# site 10 has a better fit out of all of the sites as it has the smallest circle.
# The crumbled species mostly occur close to or in sites 4, 5, 12 and 20. 
# None of the species occur around site 3, 13, 19, 21 and 22. 
# The gum species are separated from the other species. Small 
# triangles are strongly and positively associated with sites 15 and 23. 
# Although no species are found in sites 7 and 8, they are most likely the 
# most positively fitted sites in comparison with their bubble sizes to other sites. 


# Question 6 --------------------------------------------------------------

#i)

spe.hel <- decostand(spe, "hellinger")
(spe.rda <- rda(spe.hel ~ ., env))
summary(spe.rda)
sum(spe.rda$CCA$eig) 

# RDA computes axes that are linear combinations of the explanatory variables.
# In RDA the axes explain the variation of the response matrix.
# The overall variance is partitioned into constrained and unconstrained
# fractions. The constrained fraction which = 0.051, is the amount of variance of the Y matrix.
# The unconstrained fraction yielded a value of 0.6551.
# This analysis yielded 2 canonical axes (with eigenvalues labelled RDA1 to RDA2) 
# and 16 additional,unconstrained axes for the residuals (with eigenvalues labelled
# PC1 to PC16).
# The results give the eigenvalues themselves, as well as the cumulative proportion
# of variance explained (for the RDA axes) or represented (for the residual axes).
# The last cumulative value is therefore 1. The cumulative contribution to
# the variance obtained by the 2 canonical axes is the proportion of the total 
# variance of the response data explained by the RDA. It is the same value as the
# “Proportion constrained” presented above.
# The canonical eigenvalues RDA1 to RDA2 increase in value. The first
# residual eigenvalue (PC1), on the contrary, is larger than the last canonical
# eigenvalue and is  larger than most RDA eigenvalues.
# This means that the first residual structure (axis) of the data has more variance
# than some of the structures that can be explained by the explanatory variables in X.

# The overall variance is divided into constrained and unconstrained sections. 
# The eigenvalues are categorised into RDA1 to RDA2 and 16 additional, unconstrained axes for the 
# residual components. 
# The value 0.07267 describes the unconstrained proportion. The amounts of variance explained by 
# the RDA model is measured by the canonical eigenvalues. The amounts of variance represented by the 
# residual axes is measure by the residual, but the RDA model l does not explain it. 
# The accumulated constrained eigenvalues are where the proportions of the total explained variance 
# are an expression of the cumulative amounts of variance.
# Species scores are the vectors’ coordinates tips that represent the response variables in the bi- or 
# triplots. 
# Site scores are the sites coordinates as expressed in the space of the response variables Y. 
# Site constraints are the site coordinates found in the space of the explanatory variables X. 
# Biplot scores for constraining variables are the vectors’ coordinates tips that represents the explanatory variables. 
# Centroids for factor constraints are the centroids’ tips of coordinates of levels of the different factor variables.


# ii)
# Global test of the RDA result
anova(spe.rda, permutations=how(nperm=999)) 
# Tests of all canonical axes
anova(spe.rda, by="term", permutations=how(nperm=999)) 

anova(spe.rda, by="axis", permutations=how(nperm=999)) 
# Variance inflation factors (VIF)
vif.cca(spe.rda)

# Global test of the RDA showed a  p value of 0.597, df (degrees of freedom) = 2 and 
# variance = 0.051. the p value showed no significant difference.
# The residual had a df  of 25 and the variance was 0.655. 
# The canonical axes showed that carpet had a df of 1, with a variance of 0.0408, a p value of 0.081.
# Cardboard region had a df of 1, with a variance of 0.0104 and a p value of 0.968.
# For both carpet and cardboard the VIF or variance inflation factors had
# a variance of 1.928, whereas the wood had did not have a variance at all.

# The dataset’s total variance is divided into an constrained and unconstrained analysis. 
# Most of the variation in the response data may be explained by the explanatory variables, 
# when the unconstrained variance is less than the constrained variance. Therefore, the majority 
# of the variation in this analysis can be explained by the explanatory variables. The unconstrained 
# variance is higher than the constrained variance. The explanatory variance explains approximately 
# 0.05134 of the  amount of variance of the y matrix. There are two canonical axes in this analysis. 
# The analysis includes the eigenvalues labelled RDA1 and RDA2 , with 14 additional, unconstrained axes 
# for the residuals. RDA1 has an eigenvalue of 0.04097  in this analysis, whereas RDA2 = 0.01037. 
# The sum of the eigenvalues is equivalent to the total variance, for both constrained and unconstrained. 
# The results calculates the eigenvalues and the cumulative proportion of variance explained. 
# The last cumulative value = 1. The “Proportion constrained” is 0.07267 and is equivalent to the 
# proportion of the total variance of the response data, which the RDA explains.

# Permutation: There are two main tests that are done. Done relative to a permutated matrix, 
# the overall model fit of the dataset is examined by the first test. Next the partial 
# effects of each variable is examined in the second test. An overall test to show the 
# significance indicated that canonical between matrices X and Y become greatly insignificant 
# after 999 permutations have been reached. The canonical axes have a p value of greater than 0.05, 
# showing a poor correlation between species and its environment.

# iii)
coef(spe.rda)

# Scaling 1: distance triplot
dev.new(title="RDA scaling 1 + wa") 
par(mfrow = c(2, 2))
plot(spe.rda, scaling=1, 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2 (default): correlation triplot
# dev.new(title="RDA scaling 2 + wa")
plot(spep.rda, main="Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")

# Triplots of the rda results (lc scores)
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

# The plots display the environmental variables, as well as their loading on the RDA axes. 
# Cardboard explains a great section of the variance on axis 1, due to its heavy loading on 
# RDA1 and the opposite is true for carpet for RDA2 on axis 1. The closer or further the species 
# is located to or from an environmental variable, indicates how strongly a species and the environmental 
# variable is associated. The biplot shows only two of the three environmental variables present, 
# wood is excluded. The exclusion of wood shows that it is insignificant in explaining the community’s variance. 
# The different species are represented by the different points at specific sites (n = 28). 
# The length of the rate of change is in proportion with the length of the arrow. 
# The triplot1, Scaling 1 figure shows that the rate of change of carpet is greater than that of cardboard. 
# Triplot scaling 1: this figure shows that most species are correlated with the carpet and less so with cardboard.
# Triplot scaling 2: carpet and cardboard has approximately a 110º angle, 
# indicating that there is a strong negative correlation between the two variables.

# Question 7 --------------------------------------------------------------

# The study of the relationships found among the ecological processes 
# in the environment and ecosystems, is referred to as landscape ecology. Cardboard, wood and 
# carpet were the main components of the landscape in this study.
# The carpet covered the 100% of the landscape, whereas the wood and cardboard were only found 
# at specific sites, making up about 5% of the 100% of 
# the landscape. The wood and cardboard formed natural barriers, such as a mountain
# and incline. A stainless steel bowl also created  barrier and made that region inhabitable for the species
# Sites plotted close to one another, shared similar species frequencies.
# Showing that there may be underlying environmental variables playing a role
# in the abundance at the different sites. For example sites without any species were plotted close to 
# one another, possibly due to the same
# environmental factor causing an absence of species to occur at these sites. 
# The wood and cardboard, species richness was low, possibly due to unfavourable environmental conditions 
# for these species to persist.
# The incline and mountain formed by the wood and cardboard, played a role in the distribution 
# of some of the species as it isolated species from the others.
# A microhabitat formed below the incline for some species. The wind played the major role in dispersal of the species.
# The mass of the species, and the strength and direction of the wind influenced dispersal, those that 
# weighed less travelled further. 






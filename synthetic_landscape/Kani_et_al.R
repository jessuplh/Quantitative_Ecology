#Luvuyo Kani
#Practical Assignment _ Sampling (Quantitative Ecology)
#16 August 2018

#Load Packages
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
library(tidyverse)
library(ggplot2)
library(readr)
library(plyr)
library(dplyr)
library(FD)


#Load Data
env <- read_csv("env.csv")
spa <- read_csv("spa.csv")
spe <- read_csv("spe.csv")
spe_mass <- read_csv("spe_mass.csv")


# Question 1 --------------------------------------------------------------
#Dimensions
dim(env)
dim(spa)
dim(spe)


# Question 2 --------------------------------------------------------------
#a
coordinates <- spa

spe_r <- spe %>% #added a column with the number of observations in each of the sites (Species Richness)
  mutate(Species_richness = rowSums (spe, na.rm = FALSE, dims = 1)) %>% 
  select(Species_richness) %>% 
  cbind(coordinates)


ggplot(spe_r, asp = 1, aes(x = x, y = y))+
    geom_point(cex = spe_r$Species_richness, shape = 1, colour = "red")+
  labs(x = "Transacts", y = "Y_coordinates (cm)")+
  ggtitle("Species Richness")+
  theme(plot.title = element_text(hjust = 0.5))
  
#The species richness plot shows that most of the species are found in the range of 150 to 
# 300cm across the y_coordinates of the transact lines. This is because most of the
#dispersal took place at x = 25 and y = 350.  


#b

spe_H <- round(diversity(spe, index = "shannon"),2) %>% 
  cbind(coordinates)
  

ggplot(data=spe_H, aes(x=x, y=y))+
  geom_point(aes(cex = spe_H$.), fill = "red", colour = "Red")+
  labs(x = "Transacts", y = "Diversity (H')")+
  ggtitle("Species Diversity (H')")+
  theme(plot.title = element_text(hjust = 0.5))

#Similarly to the species richness plot, the Shannon-Weaver diversity plot shows that most 
# of the divesity of the sampled shapes lies with the y range of 150 to 300cm across all of
# of the trasnact lines. This again due to the nature of the dispersal as mentioned above.



#c
 spe_S <- round(diversity(spe, index = "simpson"),2) %>% 
   cbind(coordinates)
 
 ggplot(data=spe_S, aes(x=x, y=y)) +
   geom_point(aes(cex = spe_S$.), fill = "red", colour = "Red")+
   labs(x = "Transacts", y = "Diversity (D)")+
   ggtitle("Species Diversity (D)")+
   theme(plot.title = element_text(hjust = 0.5))

#The Simpson's Index plot show a fairly uniform diversity across the site. This is because 
 # the Simpson's measure of diversity is a measure of dominance or concentration, in other
 #words it measures that probability that 2 species drawn from a sample will be the same species #
 #This fairly uniform species diversity across the site reflects low number of species (17) 
 #as well as the fact that the species were intetionally spread out. 
 

# Question 3 --------------------------------------------------------------

 spe_t <- t(spe)
 spe_t_db <- vegdist(spe_t)
spe_t_db_mat <- as.matrix(spe_t_db) %>% 
  na.omit()


 
 coldiss(spe_t_db_mat, byrank = FALSE, diag = TRUE)
 
 
#iii) There is a positive species association between the crumbled species and this could be
 #because they were dispersed from a single point at one time and recieved the same amount 
 #generated wind. There is a negative species association relationship between the crumbled shape
 #species as well as the sugar sweet species. This maybe because they were dispersed at 
 #different dispersion points at different times and there is a significant weight difference between
 #them therefore the generated wind as well as the force of gravity would've had different
 #impacts on their dispersion.
 
 

# Question 4   ------------------------------------------------------------
#1
 spe.pca <- rda(spe, scale = TRUE)
 summary(spe.pca)
 
#Inertia: refers to the sum of the variables (covariance matrix) or correlation matrix diagonals
 #In correlation matrix, inertia equates to the number of variables. In this case inertia is 17.
 
#Eigenvalues: are the amount of variance (importance) explained by each of the new principle
 #components. When correlation matrix is applied, then the eigenvalues will sum up to the
 #number of variables. Eignenvalues can also be expressed as Proportions Explained or 
 #proportions of variation accounted for, by dividing them by the total inertia.
 
# Scaling 2: The eigenvectors(site scores) are scaled to the square root of the eigenvalues
 #This is approriate for interpreting the relationships amongst variables (species)
 

#Species Scores: coordinates of the arrow heads of the variables.
 
#Site scores: coordinates  of the sites in the ordination diagram.
 
 
 #2
 par(mfrow = c(1, 2))
 biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
 biplot(spe.pca, main = "PCA - scaling 2 ('species')")
 
 
 ev <- spe.pca$CA$eig
 evplot(ev)
 
 #The first 6 axes (PCA 1 to 6) must be retained because they collectively account for 75%
 #of the variation within the dataset. 
 
 #PCA - scaling 1 (sites): does not show any positive or negative correlations between the
 #sites, but rather suggests that variation between the sites is random. 
 
 #PCA - scaling 2 (species): the plot suggests that overall most of the species are 
 #associated with each other which probably due to the fact that they are sessile 
 #and do not compete with each other and that its a very small total area.
 #There is a stronger association between the paper species which is probably because 
 # they were dispersed from the same point.
 
 
  
 
 
 
 

# Question 5 --------------------------------------------------------------

apply(spe > 0, 1, sum) #to find out which rows are empty (in the results the "0" are the rows that dont have any species in them
 
spe_1 <- spe[-16,]  #removing row 16
spe_removed <- spe_1[-26:-27,] #removing row 26 to 27

spe_nmds <- metaMDS(spe_removed, distance = "bray")
summary(spe_nmds)

plot(spe_nmds, type = "t", main ="NMDS/Bray-Stress")



# Shepard plot and goodness of fit
dev.new(title = "NMDS - Shepard plot", width = 12, height = 6, noRStudioGD = TRUE)
par(mfrow = c(1, 2))
stressplot(spe_nmds, main = "Shepard plot")
gof <- goodness(spe_nmds)
plot(spe_nmds, type = "t", main = "Goodness of fit")
points(spe_nmds, display = "sites", cex = gof * 300)




# Question 6 --------------------------------------------------------------

#Tranforming the species dataset (Hellinger Tranformation)
#apply the rda functionat the defaultsetting: scaling 2

spe.H <- decostand(spe, "hellinger")
(spe.rda <- rda(spe.H ~ ., env))

#scaling 2 (default)
summary(spe.rda)

#1)#Partitioning of variance: the overall variance is partitioned into constrained and
#unconstrained fractions. The constrained fraction is the amount of variance of
#the Y matrix explained by the explanatory variables. Expressed as a proportion,
#it is equivalent to an R2 in multiple regression; in RDA this quantity is also called
#the bimultivariate redundancy statistic.
#Eigenvalues: are the amount of variance (importance) explained by each of the new principle
#components. When correlation matrix is applied, then the eigenvalues will sum up to the
#number of variables. Eignenvalues can also be expressed as Proportions Explained or 
#proportions of variation accounted for, by dividing them by the total inertia.

#Accumulated constrained eigenvalues: these are cumulative amounts of variance
#expressed as proportions of the total explained variance, as opposed to their
#contribution to the total variance described above.

#Species scores are the coordinates of the tips of the vectors representing the
#response variables in the bi- or triplots. As in PCA, they depend on the scaling
#chosen.
#Site scores (weighted sums of species scores): coordinates of the sites as
#expressed in the space of the response variables Y.

#Site constraints (linear combinations of constraining variables): coordinates of
#the sites in the space of the explanatory variables X. These are the fitted site
#scores.
#Biplot scores for constraining variables: coordinates of the tips of the vectors
#representing the explanatory variables. These coordinates are obtained as follows:
#correlations are computed between the explanatory variables and the fitted
#site scores, and then these correlations are transformed to produce the biplot
#scores.

#Centroids for factor constraints: coordinates of centroids of levels of factor
#variables, i.e. means of the scores of the sites possessing state “1” for a given
#level.



#sum of inertia of CCA and eigenvalues
sum(spe.rda$CCA$eig)

#Unconstaind axis:residual inertia - can not be explained by the original varables
#constrained is the inertia



#ii.	Run the necessary permutation tests. Explain the outcomes. (6)


#All colonical axes
anova(spe.rda, by = "term", permutations = how(nperm = 999)) #p value for each colonial axia

anova(spe.rda, by = "axis", permutations = how(nperm = 999)) #p value for each axia


#Variance inflation factors (VIF)
vif.cca(spe.rda)

#Global test for the RDA result
anova(spe.rda, permutations = how(nperm = 999)) #p= - 0.6


#iii.	Produce the necessary biplots. Explain the findings. (14)

## Triplots of the rda results (wa scores)
## Site scores as weighted averages (vegan's default)
# Scaling 1: distance triplot

dev.new(title="RDA scaling 1 + wa")
par(mfrow = c(2, 2))
plot(spe.rda, scaling=1, 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2 (default): correlation triplot
dev.new(title="RDA scaling 2 + wa")
plot(spe.rda, main="Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
dev.new(title="RDA scaling 1 + lc")
plot(spe.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2
dev.new(title="RDA scaling 2 + lc")
plot(spe.rda, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length=0, lty=1, col="red")




# Question 7 --------------------------------------------------------------

#The sampled area was largely homogenous with the exception of of the cardboard and wood 
#that made the "mountain" as well as the metal bowl. The majority of the area was carpet.
#The dispersion of the species was totally random. There was no ecological interactions 
#between the species mainly because they were all sessile and abiotic. The only major 
#influence of the distribution of species once they were dispersed was the walking around 
#and disturbance from fellow students during sampling. During the dispersal, wind was generated, 
#and thus had different influence on light weight species (paper) and bit more heavy species
#(sweets). 



 
 
 
 
    
  
   

  
    

  

  

  
  

  
  
          
          
  
  

        
  
  
  
  









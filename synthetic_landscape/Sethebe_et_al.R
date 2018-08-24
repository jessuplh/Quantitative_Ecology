# Kholofelo Sethebe
# 3137355
# 24 August 2018
# Quantitative Ecology practical
# BCB743

library(tidyverse)
library(ggpubr)
library(vegan)

# Question 1 (2 marks) ----------------------------------------------------
# Specify the dimensions of the three data sets.

env <- as.tibble(read.csv("~/Quantitative ecology Practical/Practical Assinment/env.csv"))
spa <- as.tibble(read.csv("~/Quantitative ecology Practical/Practical Assinment/spa.csv"))
spe <- as.tibble(read.csv("~/Quantitative ecology Practical/Practical Assinment/spe.csv"))

# dimensions given as rows x columns
dim(env)
dim(spa)
dim(spe)

#Question 2 (6 marks)
#i.	Provide a plot on Cartesian (spatial) coordinates of 
#a) the species richness

rowSums(spe > 0)
spe_rich <- rowSums(spe > 0)

spe_rich <- species richness(spe,index = "species richness")

ggplot(spa, asp = 1, aes(x = x, y = y))+
  geom_point(aes(cex = spe_rich), colour = "red", fill = "red")+
  labs(x = "x-coordinates(cm)", y = "y-coordinates(cm)")+
  ggtitle("Species Richness")+
  theme(plot.title = element_text(hjust = 0.5))


#b) Shannon-Weaver

SW <- diversity(spe,index = "shannon")

ggplot(spa, asp = 1, aes(x = x, y = y))+
  geom_point(aes(cex = SW), colour = "blue", fill = "blue")+
  labs(x = "x-coordinates(cm)", y = "y-coordinates(cm)")+
  ggtitle("Shannon-Weaver")+
  theme(plot.title = element_text(hjust = 0.5))


#c) The Simpson's Index
S <- diversity(spe,index = "simpson")

ggplot(spa, asp = 1, aes(x = x, y = y))+
  geom_point(aes(cex = S), colour = "purple", fill = "purple")+
  labs(x = "x-coordinates(cm)", y = "y-coordinates(cm)")+
  ggtitle("Simpson's Index")+
  theme(plot.title = element_text(hjust = 0.5))

#ii.	Explain each figure using your knowledge of how the 'landscape was created'.



# Question 3 (8 marks) ----------------------------------------------------
#i.	Calculate an association matrix for the species data. (1)

spe.t <- t(spe)
spe.jacc <- vegdist(spe.t)
matrix <- as.tibble(round(as.matrix(spe.jacc),2)) %>% 
  na.omit()

#ii.	Produce a plot that is suitable for showing the kinds of information contained in the association matrix. (3)

dev.new(title = "percentage difference (jaccard, row data", width = 10, height = 5)
coldiss(spe.jacc, byrank = FALSE, diag = TRUE)

#iii.	Briefly describe a few of the more ‘obvious’ positive and negative species associations. If possible, add an ‘ecological’ explanation for why this is so. (4)



# Question 4 (20 marks) ----------------------------------------------------
# i.	Undertake a PCA on the species data. Comprehensively explain the various (all) components of the summary() of the PCA object. (20)

spe.pca <- rda(spe)
summary(spe.pca)

#ii.	Provide plots of the PCA. How many axes must be retained, and why? What patterns become visible from the biplot? Explain your findings in detail. (10)
biplot(spe.pca)
par(mfrow=c(1,2))
biplot(spe.pca,scaling = 1, main="PCA - scaling 1 (sites)")
biplot(spe.pca, main="PCA - scaling 2 (species)")

# Eigenvalues
(ev <- spe.pca$CA$eig)
ev[ev > mean(ev)]

# Same plots using a single function: evplot()
# Plot eigenvalues and % of variance for each axis
evplot(ev)

# Question 5 (10 marks) ----------------------------------------------------
#i.	Undertake an NMDS on the species data, and provide plots of the NMDS. (2)

#removing the rows that have no data for this analysis(removing the N/A)
spe.1 <- spe[-16,]
remove_spe <- spe.1[-26:-27,]
spe.bray <- vegdist(remove_spe,method = "bray", binary = TRUE) #bray_curtis

#nmds on bray
spe.nmds <- metaMDS(spe.bray)
spe.nmds

dev.new(title = "NMDS on species data - percentage difference")
        plot(spe.nmds, type = "t", main = paste("NMDS/percentage difference", round(spe.nmds$stress,3)))

par(mfrow = c(1, 2))
stressplot(spe.nmds, spe.bray)
ordiplot(spe.nmds, type = "t", display = c("sites"),
         main = "NMDS with site scores")

abline(v = 0, h = 0, lty = 3)

dev.new(title = "NMDS - shepard plot", width = 12, height = 6)
par(mfrow = c(1,2))
stressplot(spe.nmds, main = "shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof*300)

#ii.	What patterns become visible from the ordination plot? Explain your findings in detail. (10)



# Question 6 (40 marks) ----------------------------------------------------
#i. Do a constrained analysis on the data, and comprehensively explain all the findings as per the of the summary() of the NMDS object. (20)

#Tranforming the species dataset (Hellinger Tranformation)
#apply the rda functionat the defaultsetting: scaling 2

spe.H <- decostand(spe, "hellinger")
(spe.rda <- rda(spe.H ~ ., env))

#scaling 2 (default)
summary(spe.rda) 

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



# Question 7 (10 marks) ----------------------------------------------------
#Write down your understanding of the ‘ecology’ of the landscape 
#that was sampled and analysed in the various steps, above.


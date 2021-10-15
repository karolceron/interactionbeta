# Prey availability and interaction rewiring drive the spatial and seasonal 
#structure of a predator-prey metaweb

# Karoline Ceron, Diogo B. Provete, Mathias M. Pires, 
# Andrea C. Araujo, Nico Blüthgen, Diego J. Santana

# Load packages
library('betalink')
library('bipartite')
library('BAT')
library('ggplot2')
library('mvabund')

#Interaction beta diversity
#Dry and Wet seasons of same newtork (within-seasons) or two different networks
# from different ecoregions (among ecoregions)

# Weighted matrices
ntw1<-read.table("ntw1.txt", header=TRUE) 
head(ntw1)
ntw2<-read.table("ntw2.txt", header=TRUE) 
head(ntw2)

ntw1<-as.matrix(ntw1)
ntw2<-as.matrix(ntw2)

betalinkr(webs2array(ntw1,ntw2),index = "bray",function.dist="vegdist",binary = FALSE,partitioning="poisot")
#where
# ST changes in species composition 
# OS spatial or temporal rewiring of interactions  
# WN beta diversity of interactions 
# S the species composition dissimilarity


# Anura beta diveristy 

#composition matrix
anura<-read.table("anura.txt", header=TRUE) 
head(anura)

#Calculating anura beta diversity
beta_tax=beta(anura, func = "sorensen")
btotal_tax<-beta_tax$`Btotal`

# Prey beta diversity
prey<-read.table("prey.txt", header=TRUE) 
head(prey)

beta_tax=beta(prey, func = "sorensen")
btotal_tax<-beta_tax$`Btotal`

# Regression among geogrpahic distance vs. beta diversity (prey, predator, 
# interaction and shared species)

#Importing data
beta<-read.table("beta.txt", header = T)
beta$var

#Ploting
labels_Gas <- c(
  "Banura" = expression(Banura), 
  "Bprey" = expression(Bprey), 
  "Brw" = expression(Brw),
  "Bint" = expression(Bint)
)

p1<-ggplot(beta, aes(distance, Beta, shape=var, colour=var, fill=var)) +
  geom_smooth(method="lm") +
  theme_classic() + 
  xlab("Geographic distance") +
  ylab("Beta diversity") +
  ggtitle("") + 
  expand_limits(y=0.25) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  scale_shape(labels = labels_Gas) + 
  scale_colour_discrete(labels = labels_Gas) + 
  scale_fill_discrete(labels = labels_Gas)

p1+theme(axis.text = element_text(colour="black", size=10))


# Variation on prey availability between seasons and across ecoregions
# lines ecoregions, col1 season, cols prey
# Importinga data on prey availability by season and ecoregion
ava<-read.table("season.txt", header=T)

ava_spp <- mvabund(ava[,3:30])
as.factor(ava$season)

#Plotting abundance data and mean variance
par(mar=c(2,10,2,2)) # adjusts the margins
boxplot(ava[,3:30],horizontal = TRUE,las=2, main="Abundance")
plot3<-meanvar.plot(ava_spp~as.factor(ava$season))

#Variation in prey abundances by season and sites
seas<-as.factor(ava$season)
site<-as.factor(ava$site)
plot(ava_spp~site, cex.axis=0.8, cex=0.8)
plot(ava_spp~seas, cex.axis=0.8, cex=0.8)


# Poisson regression
mod1 <- manyglm(ava_spp ~ ava$season*ava$site, family="poisson")
plot(mod1)
anova(mod1)

#There is difference of prey availability between seasons, across ecoregions 
# and between seasons across ecoregions
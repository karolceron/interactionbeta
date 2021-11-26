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

#Poisot
betalinkr(webs2array(ntw1,ntw2),index = "bray",function.dist="vegdist",binary = FALSE,partitioning="poisot")

#or

#Frund
betalinkr(webs2array(ntw1,ntw2),index = "bray",binary = FALSE,partitioning="commondenom")

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

#================================================================
#generates null networks according to prey availability per region and season and computes pairwise interaction beta diversity
#================================================================

require(dplyr)
require(plyr)
require(abind)
require(bipartite)

#....................................
#function to generate null networks
gen.null <- function(mat = full.mat, props = prob){
  
  #number of sampled species in that site
  N <-  ncol(mat)
  #number of consumed itens
  n.itens <- colSums(mat)
  rdiet.mat <- matrix(0,nrow = nrow(mat), ncol = N)
  
  for (i in 1:N){
    #sampling diets according to site-specific proportions
    rand.diet <- sample(1:nrow(mat),n.itens[i], replace = T,prob = props)
    quant <- table(rand.diet)
    rdiet.mat[as.numeric(names(quant)),i] <- quant
  }
  rownames(rdiet.mat) <- rownames(mat)
  colnames(rdiet.mat) <- colnames(mat)
  
  return(rdiet.mat)
}
#...........................

#importing data
env.data <- read.table("season.txt", header = T, as.is = T) #availability
ecoregion <- read.table("ecoregion.txt", header = T, as.is = T) #network info  


list.null <-  list()
list.local <- list()
for (k in 1: nrow(ecoregion)){
  local <- read.table(paste0(ecoregion$network[k],'.txt'), sep="\t", header=T, row.names = 1)
  region <- ecoregion$ecoregion[k]
  season.local <- ecoregion$season[k]
  
  #building full matrix
  all.prey <- names(env.data)[-c(1,2)]
  full.mat <- matrix(0,nrow = length(all.prey),ncol = ncol(local))
  colnames(full.mat) <- colnames(local)
  row.names(full.mat) <- all.prey
  aux <- match(row.names(local),all.prey)
  local <- local[-which(is.na(aux)),]
  aux <- aux[-which(is.na(aux))]
  for(i in 1:ncol(local)){
    full.mat[aux,i] <- local[,i]
  }
  
  #proportion of prey in the environment
  probs <- env.data
  probs[,3:ncol(probs)] <- probs[,3:ncol(probs)]/rowSums(probs[,3:ncol(probs)])
  prob <- filter(probs, site == region & season == season.local)
  
  #generating null networks
  rep = 10 #number of replicates
  null <- replicate(rep,gen.null(mat = full.mat, props = prob[-c(1,2)])) 
  null <- alply(null,3)
  list.null[[k]] <- null
  list.local[[k]] <- full.mat
}

#==============================
#computing pairwise beta
#example using two networks 
#==================================
i <- 8
j <- 20

#parwise beta between  networks i and j
array.local <-  webs2array(list(list.local[[i]],list.local[[j]]))
pw.beta.obs <- betalinkr_multi(array.local, index = 'bray', binary = FALSE) #computing beta

#parwise beta between null networksbased on i and j
array.null <- webs2array(c(web_i = list.null[[i]],web_j = list.null[[j]])) #creating array with all null matrices
pw.beta.null <- betalinkr_multi(array.null, index = 'bray', binary = FALSE) #computing beta

#finding the comparisons bewteen i and j only
aux <- which(unlist(lapply(strsplit(pw.beta.null$j,"j"),length))==2)
pw.beta.null.cross <- pw.beta.null[aux,]
pw.beta.null.cross <- pw.beta.null.cross[1:(rep*rep),]

hist(pw.beta.null.cross$WN, col = "black", border = "black", xlab = "Beta", main = "Beta null distribution")
abline(v = pw.beta.obs$WN, col = "tomato", lwd = 2)

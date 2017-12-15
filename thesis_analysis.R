library(vegan)
library(ggplot2)
library(SoDA)
library(mvabund)
library(multcompView)
library(lme4)
library(MASS)
library(ggmap)
library(bipartite)
library(metacom)

setwd("C:/Users/User/Dropbox/UBC/AVILES LAB/PROJECT/DATA/Analysis")

# Community data set
data.comm<-read.csv("data_comm.csv", row.names = 1) 

# Environmental data set
data.env<-read.csv("data_env.csv", stringsAsFactor = TRUE, row.names = 1) 
data.env$Year <- factor(data.env$Year) #Year as a factor
data.env$Host <- factor(data.env$Host, levels = c("Aglaoctenus", "Kapogea", "A.domingo", "A.eximius")) #Order of hosts

# Species data set
data.species <- read.csv("data_species.csv", stringsAsFactor = TRUE) 

if(all(data.species$species==colnames(data.comm))){ #Same order of species
  data.comm <- data.comm[-which(data.species$type=="generalist")] #Remove GENERALISTS
}

#Transform community data set
comm.hell<-decostand(data.comm, "hellinger") # Hellinger-transformed community data set


# Transect data sets
census.comm <- read.csv("census_comm.csv", row.names = 1)
census.env  <- read.csv("census_env.csv", row.names = 1)


#### FORWARD SELECTION OF MODEL ####

# Forward-selection of environmental variables
corr.full <- rda(comm.hell ~ log10(tangle.vol.cm3)+Host+cld.num+sociality+Year, 
                 scale = FALSE, data = data.env)

step.forward <- ordistep(rda(comm.hell ~ 1, scale = FALSE, data = data.env),
                         scope = formula(corr.full), direction = "forward", pstep = 1000)

step.forward$terms[[3]] #Which terms are selected


#### VARIATION PARTITIONING ####


#Create data sets of no missing GPS cordinates (2015 only)
data.env1 <- data.env[data.env$Year==2015 & !is.na(data.env$longitude),]
row.names(data.env1) <- data.env1$Nest

# Spatial data set
data.space <- geoXY(latitude = data.env1$latitude, longitude = data.env1$longitude, unit = 1)
row.names(data.space) <- row.names(data.env1)

# Modified Hellinger-transformed data set
comm.hell1 <- comm.hell[data.env$Year==2015 & !is.na(data.env$longitude),]

# PCNM
dis <- dist(data.space)
pcnm1 <- pcnm(dis)
data.pcnm <- as.data.frame(pcnm1$vectors)


# Forward selection based on AdjR2
rda.space <- rda(comm.hell1 ~ ., data = data.pcnm) #RDA of all PCNM vectors

fsel.space <- ordiR2step(rda(comm.hell1 ~ 1, data = data.pcnm), scope = formula(rda.space), direction = 'forward')

fsel.space$terms[[3]] #Which vector are selected. Output: PCNM21 + PCNM4 + PCNM7

# Partitioning using spatial terms from ordistep()
mm1 <- model.matrix(~ Host+log10(tangle.vol.cm3), data.env1)[,-1]
mm2 <- model.matrix(~ PCNM21 + PCNM4 + PCNM7, data.pcnm)[,-1]
variation.part <- varpart(comm.hell1, mm1, mm2)
variation.part

indiv.frac <- as.character(round(variation.part$part$indfract[,3], 3))

showvarparts(2, labels = indiv.frac, cex = 1.5,
             bg = c("green","hotpink"), 
             Xnames = "")

# test fraction [a] (environmental) with partial RDA
rda.result <- rda(comm.hell1 ~ mm1 + Condition(mm2))
anova(rda.result, step=200, perm.max=200)

# test fraction [b] (spatial) with partial RDA
rda.result <- rda(comm.hell1 ~ mm2 + Condition(mm1))
anova(rda.result, step=200, perm.max=200)


#### COMMUNITY ANALYSES ####

source("ordination_plot.R")

# Constrained ordination method chosen: RDA

# Forward-selected model for environmental variables (see above):
corr <- rda(comm.hell ~ Host+log10(tangle.vol.cm3)+Condition(Year), data = data.env)
corr
anova(corr)
plot(corr)
plot(corr, type="t", display = c("sp","bp")) #To see the species axes

tiff(filename = "C:/Users/User/Dropbox/UBC/AVILES LAB/PROJECT/DATA/Figures 2017/Fig3_ordination.tiff", width = 800, height = 400)
rda.plot<-ordination_plot(ordi = corr, 
                          group = data.env$Host,
                          species = data.species,
                          polygons=FALSE, 
                          plot.species = TRUE, 
                          sp.bold = TRUE)
dev.off()

# Summary of ordination
corr.sum <- summary(corr)

RsquareAdj(corr)$adj.r.squared

# Kaiser-Guttman criterion to residual axes (Which RDA axes are significant?)
corr$CCA$eig[corr$CCA$eig > mean(corr$CCA$eig)]

# Test significance of terms
anova.cca(corr, by = "terms")

#PERMANOVA
adonis(comm.hell~Host+log10(tangle.vol.cm3),
       method = "euclidean", data = data.env) 

# Which species drive the community differences?
mv.comm <- mvabund(data.comm)
mv.fit <- manyglm(mv.comm ~ log10(tangle.vol.cm3)+Host, 
                  data = data.env, family = "negative.binomial")

plot(mv.fit) #Distribution not skewed, therefore appropriate

#Perform glm on each species
#Takes a long time (15-25 min)
mv.glm.simple <- anova.manyglm(mv.fit, p.uni="none")
mv.glm.allp <- anova.manyglm(mv.fit, p.uni="unadjusted")

#Which associates drive the difference in communitry composition?
colnames(mv.glm.allp$uni.p)[which(mv.glm.allp$uni.p[3,]<0.05)]


#### RAREFIED RICHNESS ####

#Is rarefied richness still significantly influenced by Host and tangle size?
#The size of sample should be smaller than total community size.

rare.min <- 5
data.env2 <- data.env[-which(rowSums(data.comm)<=rare.min),]
data.comm2 <- data.comm[-which(rowSums(data.comm)<=rare.min),]
rare <- rarefy(data.comm2, sample = rare.min)       #Rarefied richness of each nest

#Host
rare.lm.host <- lm(rare ~ Host, data = data.env2)
summary(rare.lm.host)
anova(rare.lm.host)
tuke.rare <- TukeyHSD(aov(rare.lm.host))
multcompLetters(extract_p(tuke.rare)$Host)$Letters

g2<-ggplot(data.env2, aes(Host, rare))+
  geom_boxplot(aes(fill = Host))+ 
  ylab("Rarefied richness")+
  scale_y_continuous(labels=c("   1", "   2", "   3", "   4")) +
  scale_fill_manual(values=c("#FF0000FF","#8000FFFF","#4CBF00FF","#00BAFFFF")) +
  theme_classic() + 
  theme(legend.position="none")+
  theme(axis.title=element_text(size=18))+
  theme(axis.text.x=element_text(size=18, face="italic"))

g2$layers[[length(g2$layers)+1]] <- geom_text(x = 0.7, y = 4, label = "c", size = 14)
g2

#Size
fit.null <- lm(rare ~ 1, data = data.env2)
fit.full <- lm(rare ~ log10(tangle.vol.cm3)*Host, data = data.env2)
step <- stepAIC(fit.null, scope = formula(fit.full), direction="forward")
step$anova # best model is with interaction term

rare.lm.size <- lm(rare ~ log10(tangle.vol.cm3)*Host, data = data.env2)
summary(rare.lm.size)
anova(rare.lm.size)

# Each host LM
for(i in 1:4){
  rare.lm.size <- lm(rare[data.env2$Host==levels(data.env2$Host)[i]] ~ 
                       log10(data.env2$tangle.vol.cm3[data.env2$Host==levels(data.env2$Host)[i]]))
  print(levels(data.env2$Host)[i])
  print(summary(rare.lm.size))
}



# LM line fit for each host
source("plot_lmfit2.R")
anov.fit.rare <- aov(rare.lm.size) 
tuke.rare <- TukeyHSD(anov.fit.rare, which = "Host")
tuke.rare$Host
(letters.rare <- multcompLetters(extract_p(tuke.rare)$Host)$Letters)
coef(rare.lm.size)

g1 <- plot_lmfit(var = rare, int.term = TRUE, ylabel = "Rarefied richness", dat = data.env2)
g1

#### PARTIAL RDA ####

# Effect of Host alone
corr <- rda(comm.hell ~ Host+Condition(log10(tangle.vol.cm3))+Condition(Year), data = data.env)
corr
anova(corr)

# Effect of Size alone
corr <- rda(comm.hell ~ log10(tangle.vol.cm3)+Condition(Host)+Condition(Year), data = data.env)
corr
anova(corr)


#### COHERENCE ####

# Summary of ordination
corr <- rda(comm.hell ~ Host+log10(tangle.vol.cm3)+Condition(Year), data = data.env)
corr.sum <- summary(corr)

#Matrix of first axis (solitary-social)
sp.matrix1 <- as.matrix(t(comm.hell[order(corr.sum$sites[,1]),
                                    order(corr.sum$species[,1])]))
#Matrix of second axis (small-large)
sp.matrix2 <- as.matrix(t(comm.hell[order(corr.sum$sites[,2]),
                                    order(corr.sum$species[,2])]))

## Coherence (Needs to be updated in MS with comm.hell)

#First axis (solitary-social) ALL
Coherence(sp.matrix1, sims = 1000) #Nope

#Second axis (small-large) ALL
Coherence(sp.matrix2, sims = 1000) #Nope


#### SOLITARY VS. SOCIAL ####

# Test for variance in web size
library(car)
data.env$log10tangle <- log10(data.env$tangle.vol.cm3)

tapply(data.env$log10tangle,data.env$Host,FUN = var)

leveneTest(log10tangle ~ Host*sociality, data=data.env)


# Then we run an ANOVA, and post-hoc if necessary (does not work)
levene.dat.aov<-aov(dat.med.res~cat*cat2,df)
summary(levene.dat.aov)
TukeyHSD(levene.dat.aov)

# Test difference in variability btw social and solitary

(bdis <- betadisper(vegdist(comm.hell), group = data.env$Host))
plot(bdis)
anova(bdis)
TukeyHSD(bdis)


#### ACCUMULATION CURVES ####

source("accum_curve.R") #Home-made function

accum_curve(comm = data.comm, env = data.env, group = data.env$Host)

#Accum. richness of each host
for(i in levels(data.env$Host)){
  sp.accum <- specaccum(data.comm[data.env$Host==i,], method = "rarefaction")
  print(paste(i,tail(sp.accum$richness, n=1)))
}


#### AGGRESSIVE/PASSIVE ASSOCIATES ####

fgroups<-read.csv("functional_groups2.csv", stringsAsFactors=FALSE) #open functional_groups.csv
fgroups$Species <- colnames(data.comm)

source("fg_richness_abundance.R")

data.fg.ab<-fg_abundance(data.comm, fgroups)

aggressive <- vector()
passive <- vector()
for(i in 1:length(data.fg.ab[,1])){
  aggressive[i] <- (data.fg.ab$pred[i]+data.fg.ab$kleptoA[i])
  passive[i] <- (data.fg.ab$kleptoP[i]+data.fg.ab$detrit[i]+data.fg.ab$inq[i])
}

agg.pass <- data.frame(aggressive, passive)

aggressive.sum <- vector()
passive.sum <- vector()
for(i in 1:length(unique(data.env$Host))){
  aggressive.sum[i] <- sum(agg.pass$aggressive[which(data.env$Host==unique(data.env$Host)[i])])
  passive.sum[i] <- sum(agg.pass$passive[which(data.env$Host==unique(data.env$Host)[i])])
} 

fg.sums <- data.frame(aggressive.sum, passive.sum)
row.names(fg.sums)<-unique(data.env$Host)
colnames(fg.sums)<-c("aggressive", "passive")

fg.gg <- data.frame(Host = rep(unique(data.env$Host),2),
                    Type = c(rep("aggressive",4), rep("passive", 4)),
                    Total = c(fg.sums[,1], fg.sums[,2]),
                    prop.agg = rep(fg.sums[,1]/(fg.sums[,1]+fg.sums[,2]),2))

#Bar graph                    
g<-ggplot(subset(fg.gg, Type == "aggressive"), aes(Host, prop.agg))
g+geom_bar(stat = "identity", fill = "gray20")+ 
  ylab("Proportion of aggressive associates")+
  theme_classic()


# Test for significance in proportions of aggressive associates and sociality with linear mixed effects model
mix.null <- lmer(prop.aggressive ~  (1 | Host), data = data.env, REML=FALSE)
summary(mix.null)

mix.lm <- lmer(prop.aggressive ~ sociality + (1 | Host), data = data.env, REML=FALSE)
summary(mix.lm)

anova(mix.null, mix.lm)


#### NESTEDNESS WITHIN HOST ####

# NODF (takes a while)
nested.table <- data.frame(Host=rep(NA,4), NODF=rep(NA,4), null.NODF=rep(NA,4), pval=rep(NA,4))

for(i in 1:length(levels(data.env$Host))){
  comm.h <- log(data.comm[data.env$Host==levels(data.env$Host)[i],]+1)
  comm.nodf <- as.numeric(nestednodf(comm.h, weighted = TRUE)$statistic[3])
  
  nulls.nodf <- bipartite::nullmodel(as.matrix(round(comm.h, 2)), N=100, method="swap.web")
  null.nodf <- vector()
  for(j in 1:length(nulls.nodf)){
    null.nodf[j] <- as.numeric(nestednodf(nulls.nodf[[j]], weighted = TRUE)$statistic[3])
  }
  
  tail95.nodf <- sort(null.nodf)[round(length(null.nodf)*0.95,0)]
  pval.nodf <- length(which(null.nodf>comm.nodf))/length(null.nodf)
  pval.nodf
  
  nested.table[i,] <- c(levels(data.env$Host)[i], comm.nodf, mean(null.nodf), pval.nodf)
  
}

nested.table



#### EFFECT OF SIZE ####

for(i in 1:4){
  print(levels(data.env$Host)[i])
  print(paste("min:", min(round(log10(data.env$tangle.vol.cm3[data.env$Host==levels(data.env$Host)[i]]),2)),
        "max:", max(round(log10(data.env$tangle.vol.cm3[data.env$Host==levels(data.env$Host)[i]]),2))))
}

#test and plot each variable

source("plot_lmfit2.R")
source("multiplot.R")

## Abundance
fit.full <- lm(log10(rowSums(data.comm)+1) ~ log10(tangle.vol.cm3)*Host+cld.num, data = data.env)
step <- stepAIC(fit.full, direction="both")
step$anova # best model is without interaction term and without cld

fit.ab <- lm(log10(rowSums(data.comm)+1) ~ log10(tangle.vol.cm3)+Host, data = data.env)
coef(fit.ab)

anova(fit.ab)
summary(fit.ab)
anov.fit.ab <- aov(fit.ab) 
tuke.ab <- TukeyHSD(anov.fit.ab, which = "Host")
tuke.ab$Host
(letters.ab <- multcompLetters(extract_p(tuke.ab)$Host)$Letters)


## Richness
fit.rich <- lm(log10(specnumber(data.comm)+1) ~ log10(tangle.vol.cm3)+Host, data = data.env)
coef(fit.rich)

anova(fit.rich)
summary(fit.rich)
anov.fit.rich <- aov(fit.rich) 
tuke.rich <- TukeyHSD(anov.fit.rich, which = "Host")
tuke.rich$Host
(letters.rich <- multcompLetters(extract_p(tuke.rich)$Host)$Letters)


## Density of individuals per m3
density.m3 <- rowSums(data.comm)/data.env$tangle.vol.cm3*10000
fit.dens <- lm(log10(density.m3+1) ~ log10(tangle.vol.cm3)+Host, data = data.env)
coef(fit.dens)

anova(fit.dens)
summary(fit.dens)
anov.fit.dens <- aov(fit.dens) 
tuke.dens <- TukeyHSD(anov.fit.dens, which = "Host")
tuke.dens$Host
(letters.dens <- multcompLetters(extract_p(tuke.dens)$Host)$Letters)

# GRAPHS
g1 <- plot_lmfit(var = log10(rowSums(data.comm)+1), int.term = TRUE,
           ylabel = "Number of associates (log10-transformed)")
g1$layers[[length(g1$layers)+1]] <- geom_text(x = 2, y = 2.5, label = "a", size = 14)

g3 <- plot_lmfit(var = log10(data.env$density.m3+1), int.term = TRUE,
           ylabel = "Associates/m3 (log10-transformed)")
g3$layers[[length(g3$layers)+1]] <- geom_text(x = 6.3, y = 2.25, label = "b", size = 14)

tiff(filename = "C:/Users/User/Dropbox/UBC/AVILES LAB/PROJECT/DATA/Figures 2017/Fig4.tiff", width = 1000, height = 900)
multiplot(g1,g2, g3, cols = 2)
dev.off()


#### COLONIZATION ####
#Most analyses made using JMP.



#### APPENDIX ####



#### MAP OF WEBS ####

#ggmap

sites.js <- get_stamenmap(bbox = c(-77.635,-1.09,-77.607, -1.06), 
                          source = "stamen", maptype = "terrain-lines")
g <- ggmap(sites.js)
g+geom_point(aes(x = longitude, y = latitude, shape = factor(Host)),#, colour = factor(Year)), 
             size = 3, data = data.env)+
  scale_shape_manual(name  ="Host", values=c(0,6,17,19))+
  #scale_colour_manual(name  ="Year", values = c("gray30", "black"))+ 
  ylab("Latitude")+
  xlab("Longitude")+
  theme_bw()


#### TRANSECT DATA ####

species.prop <- data.species

species.prop$species <- factor(species.prop$species, 
                               levels = c(as.character(species.prop$species[order(species.prop$prop.transects, decreasing = TRUE)]), "All others"))

transect.gg <- species.prop[species.prop$prop.transects>0,]
transect.gg <- rbind(transect.gg, c("All others", NA, 0, NA))
transect.gg$species <- factor(transect.gg$species)
transect.gg$prop.transects <- as.numeric(transect.gg$prop.transects)*100


g<-ggplot(transect.gg, aes(species, prop.transects))
g+geom_bar(fill = "grey20", stat = "identity")+
  ylab("Total percentage of encountered arthropods")+
  ylim(0,10)+
  theme_classic()+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12))


#### web size between host ####

#test differences

fit.siz <- aov(log10(tangle.vol.cm3) ~ Host, data = data.env)
fit.siz
summary(fit.siz)
tuke <- TukeyHSD(fit.siz) 

letters <- multcompLetters(extract_p(tuke)$Host)$Letters

# Boxplot

g<-ggplot(data.env, aes(Host, log10(tangle.vol.cm3)))
g+geom_boxplot(aes(fill = Host))+ #, outlier.shape = NA)+
  geom_text(x = 1, y = max(log10(data.env$tangle.vol.cm3))+1, label = "ab", size = 10)+
  geom_text(x = 2, y = max(log10(data.env$tangle.vol.cm3))+1, label = "a", size = 10)+
  geom_text(x = 3, y = max(log10(data.env$tangle.vol.cm3))+1, label = "b", size = 10)+
  geom_text(x = 4, y = max(log10(data.env$tangle.vol.cm3))+1, label = "c", size = 10)+
  ylab("Tangle volume (log10-transformed)")+
  ylim(c(0,max(log10(data.env$tangle.vol.cm3))+1))+
  scale_fill_discrete(guide=FALSE)+
  theme_classic()


#### Height from ground ####

g<-ggplot(subset(data.env, sociality == "solitary"), aes(Host, log(h.from.ground)))
g+geom_boxplot(aes(fill = Host))+
  ylab("Height from ground (cm, log-transformed)")+
  scale_fill_discrete(guide=FALSE)+
  theme_classic()

#### Number of hosts per associates ####
occ <- vector()
# Occurrence in each host
for(i in 1:length(colnames(data.comm))){
  occ[i]<-length(unique(data.env$Host[data.comm[i]>0]))
  names(occ)[i]<-colnames(data.comm)[i]
}

length(which(occ==4))/length(colnames(data.comm))
length(which(occ==3))

# Plot means new vs old webs
census.means <- read.csv("census_means.csv")
census.means$Source <- factor(census.means$Source, 
                              levels = c("transplant", "jungle"))

census.means$delta <- rep(census.means$Mean.prop[census.means$Source=="jungle"]-
                            census.means$Mean.prop[census.means$Source=="transplant"], 2)

census.means$Species<- gsub("Faiditus sp.1", "Faiditus sp.1 *", census.means$Species)
census.means$Species<- gsub("Ranzovius", "Ranzovius *", census.means$Species)
census.means$Species<- gsub("Faiditus sp.6", "Faiditus sp.6 *", census.means$Species)

census.means$Species <- factor(census.means$Species, 
                               levels = as.character(census.means$Species[
                                 order(census.means$delta[1:13], decreasing = TRUE)]))


#Color palette, red = highest decrease, green = highest increase
delta.col <- rainbow(13, end = 2/6)[rank(census.means$delta[order(census.means$delta[1:13], decreasing = TRUE)])]

#Plot asin-transformed mean proportions and CIs
g<-ggplot(census.means, aes(Source, Mean.prop.asin, colour=Species, group = Species))
g+geom_line() +
  geom_errorbar(width=.1, aes(ymin=Mean.prop.asin-ci.prop.asin, ymax=Mean.prop.asin+ci.prop.asin)) +
  geom_point(shape=21, size=3, fill="white")+
  scale_color_manual(values = delta.col)+
  xlab(" ")+
  ylab("Mean proportion in webs (arcsin-transformed)")+
  scale_x_discrete(labels= c("New webs", "Old webs"))+
  theme_classic()


#### NOTES ####

# 1. Araneidae sp.2 in analysis is named Araneidae sp.1 in MS







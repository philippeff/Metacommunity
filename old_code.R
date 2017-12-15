# OLD CODE AND DATA FRAMING THESIS PROJECT 


#### VARIATION PARTITIONING ####

#visualize each pcnm vector
pcnm.vector <- 3 #vector to plot
ordisurf(data.space, scores(pcnm1, choi=pcnm.vector), bubble = 4, main = "PCNM 3")

#RDA of all PCNM vectors
rda.space <- rda(comm.hell1 ~ ., data = data.pcnm)

#Funky way of graphing
grid.newpage()
grid.rect()
indiv.frac <- round(variation.part$part$indfract[,3], 3)
venn <- draw.pairwise.venn(area1 = indiv.frac[1]+indiv.frac[2], area2 = indiv.frac[2]+indiv.frac[3], cross.area = indiv.frac[2], 
                           category = c("Environmental", "Spatial"), lwd = 1, 
                           fontfamily = "sans", cat.fontfamily = "sans",
                           fill = c("light green", "light pink"), cex = c(4,1.5,1.5),
                           cat.pos = c(355, 25), cat.dist = c(-0.06, 0.03), cat.cex = c(2,1.5))
grid.text(label = paste("Residuals = ", indiv.frac[4]), 0.8, 0.05, gp = gpar(cex = 1.5))

#### RDA ####

#Nested anova for sociality nested within Host
nested.npmanova(comm.hell ~ sociality+Host, data = data.env, permutations = 1000)
#Signif.

#Nested anova for cld nested within Host
nested.npmanova(comm.hell ~ Host+CLD, data = data.env, permutations = 1000)
#NS

# Specialists graph

# Graph abundances of indicative species
comm.gg <- data.frame(Nest = rep(data.env$Nest, 10), 
                      Host = rep(data.env$Host, 10), 
                      Associate = c(rep(colnames(data.comm)[4],150),
                                    rep(colnames(data.comm)[5],150),
                                    rep(colnames(data.comm)[11],150),
                                    rep(colnames(data.comm)[13],150),
                                    rep(colnames(data.comm)[15],150),
                                    rep(colnames(data.comm)[16],150),
                                    rep(colnames(data.comm)[17],150),
                                    rep(colnames(data.comm)[18],150),
                                    rep(colnames(data.comm)[19],150),
                                    rep(colnames(data.comm)[20],150)),
                      Abundance.log = log10(c(data.comm[,4],
                                              data.comm[,5],
                                              data.comm[,11],
                                              data.comm[,13],
                                              data.comm[,15],
                                              data.comm[,16],
                                              data.comm[,17],
                                              data.comm[,18],
                                              data.comm[,19],
                                              data.comm[,20])+1))

#Wrap by associate                            
g<-ggplot(comm.gg, aes(Host, Abundance.log))
g+geom_bar(fill = "grey20", stat = "identity", aes(group = Host))+ #, outlier.shape = NA)+
  ylab("Mean abundance (log10-transformed)")+
  facet_wrap(~Associate, ncol = 5)+
  theme_bw()+
  theme(axis.title.y  = element_text(size = 17),
        axis.title.x  = element_text(size = 17))+
  theme(axis.text.x  = element_text(angle=80, vjust=0.5, hjust = 0.5, size=12))


#Wrap by host                            
g<-ggplot(comm.gg, aes(Associate, Abundance.log))
g+geom_bar(fill = "grey70", stat = "identity", aes(group = Host))+
  ylab("Mean abundance (log10-transformed)")+
  facet_wrap(~Host, ncol = 1)+
  theme_bw()+
  theme(axis.title.y  = element_text(size = 17),
        axis.title.x  = element_text(size = 17))+
  theme(axis.text.x  = element_text(angle=60, vjust = 1.1, hjust = 1.25, size=12))

#Test difference in variance OF WEB SIZE
var.test(log(data.env$tangle.vol.cm3[data.env$sociality=="social"]),
         log(data.env$tangle.vol.cm3[data.env$sociality=="solitary"]))

### richness ####

g2 <- plot_lmfit(var = log10(specnumber(data.comm)+1), int.term = TRUE,
                 ylabel = "Associate richness (log10-transformed)")
g2$layers[[length(g2$layers)+1]] <- geom_text(x = 2, y = 1.07, label = "c", size = 14)

#### PATTERNS (COHERENCE) ####
corr <- rda(comm.hell ~ Host+log10(tangle.vol.cm3)+Condition(Year), data = data.env)
corr
anova(corr)
corr.sum <- summary(corr)

pattern <- data.comm[order(corr.sum$sites[,1]),order(corr.sum$species[,1])]

#Visualize

# Look for modularity or nestedness
visweb(comm.hell, type = "diagonal") #Nope
visweb(comm.hell, type = "nested") #Nope


#Mean (normalised, between 0-1) number of checkerboard combinations across all species.
# quantitative index of community organization (coherence) 
comm.cscore <- networklevel(as.matrix(data.comm), index = "C score")

nulls.cs <- nullmodel(as.matrix(data.comm), N=100, method="swap.web")
null.cs <- unlist(sapply(nulls.cs, networklevel, index="C score"))

tail95.cs <- sort(null.cs[1,])[round(length(null.cs[1,])*0.95,0)]

#P-value
(pval.cs <- length(which(null.cs[1,]>comm.cscore[1]))/length(null.cs[1,]))

hist(null.cs[1,], xlim = c(0.5,0.8), main = "C score")
abline(v=comm.cscore[1], col = "red", lwd = 2)


### STRUCTURE ####

# Summary of ordination
corr
corr.sum <- summary(corr)

#Matrix of first axis (solitary-social)
sp.matrix1 <- as.matrix(t(comm.hell[order(corr.sum$sites[,1]),
                                    order(corr.sum$species[,1])]))
#Matrix of second axis (small-large)
sp.matrix2 <- as.matrix(t(comm.hell[order(corr.sum$sites[,2]),
                                    order(corr.sum$species[,2])]))

## Coherence (Needs to be updated in MS with comm.hell)

#First axis (solitary-social) ALL
Coherence(sp.matrix1, sims = 100) #Nope

#Second axis (small-large) ALL
Coherence(sp.matrix2, sims = 100) #Nope


# Test significance in proportion differences

likelihood.test(x = census.means$Mean.prop[which(census.means$Source=="jungle")], 
                y = census.means$Mean.prop[which(census.means$Source=="transplant")], 
                conservative = FALSE)

chisq.test(census.means$Mean[which(census.means$Source=="jungle")], 
           census.means$Mean[which(census.means$Source=="transplant")])

data(InsectSprays)

likelihood.test(InsectSprays$count>7,InsectSprays$spray)

###### DATAFRAMING #####

data15<-read.csv(file.choose(), stringsAsFactor = TRUE) # open data2015_full.csv
data14<-read.csv(file.choose(), stringsAsFactor = TRUE) # open data2014_full.csv 

data<-read.csv(file.choose(), stringsAsFactor = TRUE) # open data associates.csv in 2015 data folder

# Select only species that were observed at least once
# Include only species observed >5 times

data <- data %>%
  select(Nest, Host, (which(colSums(data[4:(length(data)-1)])>5)+3)) %>%
  arrange(Nest)


write.csv(data, file="data2015.csv", row.names = FALSE)



# Put nest measurements in data.csv (July 28, 2015)

nest<-read.csv(file.choose(), stringsAsFactor = TRUE) # open data associates.csv in 2015 data folder

nest <- nest %>%
  arrange(Nest)

#data1 <- data.frame(data[3:length(data)], nest, stringsAsFactors = TRUE)

write.csv(data1, file="data2015_full.csv", row.names = FALSE)


# Exclude Sumaco in data14

data14 <- data14[-grep("Sumaco", as.character(data14$Site)),]
data14$Site <- factor(data14$Site)

# Exclude ROAD nests of eximius in data14

data14 <- data14[-grep("ROAD", as.character(data14$Nest)),]
data14$Nest <- factor(data14$Nest)

data14 <- data14 %>%
  select((which(colSums(data14[4:(which(colnames(data14)=="Total.Associates")-1)])>5)+3), 
         Nest, Host, which(colnames(data14)=="Collect.Date"):which(colnames(data14)=="Notes")) %>%
  arrange(Nest)

data14$Nest <- factor(gsub("funnel", "age", data14$Nest))
data14$Nest <- factor(gsub("grid", "mecy", data14$Nest))
data14$Host <- factor(gsub("funnel", "agelenopsis", data14$Host))
data14$Host <- factor(gsub("grid", "mecynogea", data14$Host))

data14 <- data14  %>%
  arrange(Nest)

write.csv(data14, file="data2014_full.csv", row.names = FALSE)

# merge both years

data.both <- merge(data, data14, all = TRUE)

data.both <- data.both  %>%
  arrange(Nest)

write.csv(data.both, file="data_both.csv", row.names = FALSE)
# Then worked on it manually in Excel

# Remove empty nests (2014age16 2014dom25)

data.comm <- data.both[1:(which(colnames(data.both)=="Unidentified")-1)]
which(rowSums(data.comm)==0)
data.both <- data.both[-which(rowSums(data.comm)==0),]
data.both$Nest <- factor(data.both$Nest)

write.csv(data.both, file="data_both.csv", row.names = FALSE)


#Add CLD for all mecynogea and missing agenlenopsis 2014
data.both$CLD[which(data.both$Host=="agelenopsis"&data.both$Year==2014)] <- "some"
data.both$CLD[which(data.both$Host=="mecynogea")] <- "sparse"

data.both$CLD <- gsub("Dense", "dense", data.both$CLD)
data.both$CLD <- gsub("Many", "many", data.both$CLD)
data.both$CLD <- gsub("Some", "some", data.both$CLD)
data.both$CLD <- gsub("Sparse", "sparse", data.both$CLD)
data.both$CLD <- factor(data.both$CLD)

write.csv(data.both, file="data_both.csv", row.names = FALSE)


# Modify data_both.csv with Sam's data (Sept 10, 2015)

data.both1<-read.csv("data_both_add.csv", stringsAsFactor = TRUE) 

data.both <- rbind(data.both, data.both1)

data.both <- data.both  %>%
  arrange(Nest)

write.csv(data.both, file="data_both.csv", row.names = FALSE)


data.both<-read.csv("data_both.csv", stringsAsFactor = TRUE) # open data.both.csv 

#Exclude Araneus
data <- data.both[-grep("araneus", as.character(data.both$Host)),]
data$Host <- factor(data$Host)
#Year as a factor
data$Year <- factor(data$Year)

data.comm <- data[1:(which(colnames(data)=="Unidentified")-1)]
data.env <- data[which(colnames(data)=="Nest"):length(data)]

data.env$Simpson <- diversity(data.comm, index = "simpson") # Simpson evenness index
data.env$Shannon <- diversity(data.comm, index = "shannon") # Shannon entropy index
data.env$abundance <- rowSums(data.comm)
data.env$richness <- specnumber(data.comm) #richness of each nest

#Mysmenidae instead of Theridiosomatidae
colnames(data)<- gsub("Theridiosomatidae", "Mysmenidae", colnames(data))

write.csv(data.env, file="data_env.csv", row.names = FALSE)
write.csv(data.comm, file="data_comm.csv", row.names = FALSE)


### Modify quality in data.env

data.env$Quality <- gsub("Excellent", "excellent", data.env$Quality)
data.env$Quality <- gsub("Good", "good", data.env$Quality)
data.env$Quality <- gsub("Fair", "fair", data.env$Quality)
data.env$Quality <- gsub("Poor", "poor", data.env$Quality)

data.env$density.m3 <- data.env$abundance/data.env$tangle.vol.cm3*10000


write.csv(data.env, file="data_env.csv", row.names = FALSE)


# Add sociality (Jan 25, 2016)
sociality <- vector()
for(i in 1:length(data.env$Host)){
  if(data.env$Host[i]=="eximius"|data.env$Host[i]=="domingo"){
    sociality[i] <- "social"
  } else if(data.env$Host[i]=="agelenopsis"|data.env$Host[i]=="mecynogea"){
    sociality[i] <- "solitary"
  }
}

data.env$sociality <- sociality

write.csv(data.env, file="data_env.csv", row.names = FALSE)


# Change CLD to ordinal (Feb 5, 2016)
levels(data.env$CLD)

cld.num <- vector()
for(i in 1:length(data.env$CLD)){
  if(data.env$CLD[i]=="sparse"){
    cld.num[i] <- 1
  } else if(data.env$CLD[i]=="some"){
    cld.num[i] <- 2
  } else if(data.env$CLD[i]=="many"){
    cld.num[i] <- 3
  } else if(data.env$CLD[i]=="dense"){
    cld.num[i] <- 4
  }
}

data.env$cld.num <- cld.num

write.csv(data.env, file="data_env.csv", row.names = FALSE)


##### COLONIZATION DATAFRAMING ####

census<-read.csv("census1.csv", stringsAsFactors=TRUE) # open census1.csv (for analysis)
census.size<-read.csv("census_day24.csv", stringsAsFactors=TRUE) # open census_day24.csv

census24.comm <- census.size[1:16]
row.names(census24.comm) <- census.size$NEST
colnames(census24.comm)[6] <- "Mysmenidae"
census24.comm <- census24.comm[-c(10,14,15)]
census24.env  <- census.size[16:length(colnames(census.size))]

colnames(census24.comm) 
colnames(data.comm) 

data.comm.census <- data.comm %>%
  transmute(Araneidae = Araneidae06+Araneidae12,
            Theridiidae = Theridiidae07+Theridiidae01+Theridiidae03+Theridiidae06,
            Argyrodes = Argyrodes.02+Argyrodes.06,
            Neospintharus = Neospintharus01+Neospintharus06,
            Rhomphaea01 = Rhomphaea01,
            Mysmenidae = Mysmenidae01+Mysmenidae02,
            Uloboridae = Uloboridae01,
            Salticidae = Salticidae01,
            Mimetidae = Mimetidae01,
            Hemiptera01 = Hemiptera01,
            Hemiptera04 = Hemiptera04,
            Ants = Crematogaster01+Wasmannia01,
            Lathridiidae = Lathridiidae)

row.names(data.comm.census) <- data.env$Nest
data.comm.census <- data.comm.census[which(data.env$Host=="eximius"),]

colnames(data.comm.census) == colnames(census24.comm)
write.csv(data.comm.census, "communities_forcensus.csv")

census.comm <- rbind(data.comm.census, census24.comm)
write.csv(census.comm, "communities_nestandcensus.csv")

colnames(census24.env) 
colnames(data.env) 


census <- read.csv("communities_nestandcensus.csv", row.names = 1)
census.prop <- read.csv("census_tab.csv")

census.comm <- census[1:13]
census.env <- census[14:15]

census.comm[which(census.env$source=="transplant"),
            which(colnames(census.comm)=="Theridiidae")] <-
  c(0,2,6,0,0,0,2,1,0,0,2,6,0,0,2,1,3,0,1,8,1,0,2,0,4)

census.comm[which(census.env$source=="transplant"),
            which(colnames(census.comm)=="Neospintharus")] <-
  c(0,4,6,2,0,0,2,2,2,0,0,0,0,0,1,2,1,7,3,4,1,1,6,3,1)

census.comm[which(census.env$source=="forest"),
            which(colnames(census.comm)=="Theridiidae")] <-
  data.comm[which(data.env$Host=="eximius"), which(colnames(data.comm)=="Theridiidae01")]

census.comm[which(census.env$source=="forest"),
            which(colnames(census.comm)=="Neospintharus")] <-
  data.comm[which(data.env$Host=="eximius"), which(colnames(data.comm)=="Neospintharus01")]

colnames(census.comm) <- gsub("Theridiidae", "Theridiidae01", colnames(census.comm))
colnames(census.comm) <- gsub("Neospintharus", "Neospintharus01", colnames(census.comm))

write.csv(census.comm, "census_comm.csv")
write.csv(census.env, "census_env.csv")

census.tab <- as.data.frame(t(prop.table(as.matrix(census.sums), 1)))

census.props <- census.comm
for(i in 1:length(census.props[,1])){
  for(j in 1:length(census.props[1,])){
    census.props[i,j] <- census.props[i,j]/sum(census.props[i,])
  }
}
census.props$nest <- row.names(census.props)

census.props <- census.props %>%
  select(Theridiidae01, Neospintharus01, Hemiptera01, nest) %>%
  gather("associate", "proportion", 1:3)
census.props$host.source <- rep(c(rep("forest", 42),rep("transplant", 25)), 3)

write.csv(census.props, "census_props.csv")

census.prop <- census.prop[which(census.prop$associate=="Theridiidae"|
                                   census.prop$associate=="Neospintharus"|
                                   census.prop$associate=="Hemiptera01"),]
census.prop$associate <- factor(census.prop$associate, 
                                levels = c("Theridiidae", "Neospintharus", "Hemiptera01"))


grays<- gray(0:2/2)

barchart(proportion~host.source|associate, data=census.prop,
         layout=c(3,1), col=grays, between=list(x = 0.2, y = 1), box.ratio = 10,
         scales=list(cex=c(1,1),rot=0,alternating=3, x=list(cex=1.2,relation="free",rot=45)),
         ylab="Proportion", xlab="Associate", ylim=c(0,0.4))


#Change names and order of hosts in data.env
data.env$Host <- gsub("agelenopsis", "Aglaoctenus", data.env$Host)
data.env$Host <- gsub("mecynogea", "Kapogea", data.env$Host)
data.env$Host <- gsub("eximius", "A.eximius", data.env$Host)
data.env$Host <- gsub("domingo", "A.domingo", data.env$Host)
data.env$Host <- factor(data.env$Host, levels = c("Aglaoctenus", "Kapogea", "A.domingo", "A.eximius"))

write.csv(data.env, "data_env.csv")

#Change associate names
colnames(data.comm) <- c("Argyrodes.sp1",
                         "Neospintharus.sp1",
                         "Mimetus",
                         "Faiditus.sp6",
                         "Faiditus.sp5",
                         "Faiditus.sp4",
                         "Emesinae",
                         "Araneidae.sp2",
                         "Latridiidae",
                         "Tetramorium",
                         "Rhomphaea",
                         "Hypaeus",
                         "Faiditus.sp1",
                         "Faiditus.sp3",
                         "Faiditus.sp2",
                         "Mysmenopsis.sp1",
                         "Mysmenopsis.sp2",
                         "Philoponella",
                         "Ranzovius",
                         "Crematogaster",
                         "Araneidae.sp1",
                         "Wasmannia")

write.csv(data.comm, "data_comm.csv")

## Latridiidae is actually Paratenetus (Tenebrionidae)

colnames(data.comm)[colnames(data.comm)=="Latridiidae"] <- "Paratenetus"

colnames(data.comm)

write.csv(data.comm, "data_comm.csv")











#### OLD CODE ####

source("null_model_analysis_WNODF_QuaBiMo.R")
source("null_model_analysis_NODF_QuaBiMo.R")


null.nodf<-null_model_analysis_NODF_QuaBiMo(host.comm, 
                                            type = "NODF", 
                                            observed_value = nodf, N_null_webs = 100, 
                                            null_model = "swap.web")
mean.null.NODF[i] <- as.numeric(null.nodf[2])
P.value.NODF[i] <- 2*pnorm(-abs(as.numeric(null.nodf[1])))

ind # Output

ind$NestedConfig$DegreeMatrix # Ordered matrix

ind$Bin_t1$NODF$Measure # NODF score
ind$Bin_t1$NODF$pvalue #represents the probability of attaining a more nested matrix than the one under consideration.
ind$Bin_t1$NODF$NormalisedTemperature # If above 1.0 matrix is nested (when calculating NODF)


#Visualize nestedness

# Look for modularity
visweb(data.comm, type = "diagonal")

# Within host
host.comm <- data.comm[which(data$Host==host),]

visweb(host.comm, type = "nested") #type = "diagonal")

#### NESTEDNESS ANALYSIS USING FALCON PACKAGE

source("InteractiveMode.R") # Works with binary data only (presence/absence)


# Manual entry:
ind<-PERFORM_NESTED_TEST(MAT = host.matrix, 
                         bintest = 1, # binary? 0 = quant 1 = binary
                         sortVar = 1, # sort matrix?
                         functhand = "NODF",
                         nullmodels = 1, # 1=SS, 2=FF, 3=CC, 4=DD, 5=EE, c()=all 
                         EnsembleNumber = c(),
                         plotON = 1)

###EXAMPLE PERMANOVA


## create a design matrix of the contrasts for "treat"
contrasts(treat) <- cbind(c(0,1,0),c(0,0,1))
Treat <- model.matrix(~ data.env$Host)[, -1]

imp.in.t1 <- Imp * ifelse(treat == "t1", 1, 0)
imp.in.t2 <- Imp * Treat[, 1]
imp.in.t3 <- Imp * Treat[, 2]

## specify the orthogonal contrasts for "treat"
contrasts(data.env$Host) <- cbind(c(1, -1, 0), c(1, 0, -1))

## specify the design matrix of the orthogonal
## contrasts for "treat"
Treat.ortho <- model.matrix(~ data.env$Host)[, -1]

## create a factor for each of the orthogonal "treat" contrasts
treat1vs2 <- Treat.ortho[, 1]
treat1vs3 <- Treat.ortho[, 2]
treat1vs4 <- Treat.ortho[, 3]


## do the pm-manova with the full model
fm1 <- adonis(spdf ~ treat * imp, method = "euclidean", perm = 999)

## do the pm-manova with the orthogonal contrasts for imp and treat'
## and the interaction contrasts of interest
fm2 <- adonis(comm.hell ~ treat1vs2 + treat1vs3 +
                treat1vs4,
              method = "euclidean", perm = 999)
fm2


### COMPARE COMMUNITIES
cs<-capscale(comm.hell
             ~log(tangle.vol.cm3)*Host+Condition(Year), distance = "bray", data = data.env) 

cs

plot(cs)

# ordiellipse draws 95 % confidence ellipses around class centroids. 
# If these confidence ellipses do no overlap, the classes probably are significantly 
# different at level P ??? 0.05.

ordiplot(corr)

corr.ell <- ordiellipse(corr, groups = data.env$Host, kind="se",
                        scaling = 2, conf=0.95, lwd=1, col="black")
summary(corr.ell)


corr.hull <- ordihull(corr, groups = data.env$Host, label = TRUE)
summary(corr.hull)

### RAREFACTION CURVES

data.rare <- data %>%
  group_by(Host) %>% 
  select(which(colnames(data) %in% colnames(data.comm))) %>%
  summarise_each(funs(sum)) 

data.rare <- as.data.frame(as.matrix(data.rare[which(colnames(data) %in% colnames(data.comm))]))
row.names(data.rare) <- unique(data.env$Host)

rarecurve(data.rare, label = TRUE)

for(i in 1:length(data.rare[which(colnames(data) %in% colnames(data.comm))])){
  data.rare[which(colnames(data) %in% colnames(data.comm))]
}

rownames(data.comm) <- data.env$Nest
rarecurve(data.comm, label = TRUE)

#### ANALYSE ORDINATION PLOT

#Difference btw centroids 

fm1<-adonis(comm.hell
            ~log(tangle.vol.cm3)+Host, data = data.env, method = "euclidean") 
fm1

fm2<-adonis(comm.hell~ log(tangle.vol.cm3) + Year, data = data.env, method = "euclidean") 
fm2



#Mean distance btw centroids

comm.dist<-vegdist(comm.hell)
comm.dist.mean <- meandist(comm.dist, grouping = data.env$Host) #of a vegdist() object + betadisper()
plot(comm.dist.mean)
TukeyHSD(comm.dist.mean)
summary(comm.dist.mean)


corr.mrp<-mrpp(comm.hell, grouping = data.env$Host)
summary(corr.mrp)











#Rarefied richness is still significantly influenced by Host and tangle size
data.env2 <- data.env[-which(rowSums(data.comm)==1),]
data.comm2 <- data.comm[-which(rowSums(data.comm)==1),]
rare <- rarefy(data.comm2, sample = 2)
rare.lm <- lm(rare ~ Host, data = data.env2)

plot(x = data.env$Host[-which(rowSums(data.comm)==1)], y = rare)
text(x=4, y=1.7, main = "*")
summary(rare.lm)
anova(rare.lm)
TukeyHSD(aov(rare.lm))

g<-ggplot(data.env2, aes(Host, rare))
g+geom_boxplot(aes(fill = Host))+ 
  geom_text(x = 4, y = 1.7, label = "*", size = 15)+
  ylab("Rarefied richness")+
  theme_classic() + 
  theme(legend.position="none")


### CONNECTANCE ###

rare <- rarefy(data.env$connectance, sample = 2)
rare.lm <- lm(rare ~ tangle.vol.cm3+Host, data = data.env)
summary(rare.lm)
anova(rare.lm)

g<-ggplot(data.env, aes(x = richness, y = connectance))
g+geom_point(aes(colour=Host), size=2)+
  geom_smooth(aes(group=Host, colour=Host), fullrange=FALSE, 
              method="lm", size=1, fill = NA)+
  theme_classic()

# Function to draw polygons in NMDS plot

polygon<-function(fac){
  # vector holding the colors
  cols <- ncolors     #cols <- c("blue", "green", "red", "purple", "yellow", "orange")
  # empty plot
  ordiplot(m, type = "n")
  # Add points colored by Environmental Variable (Host)
  points(m, col = cols[fac], pch = 16)
  #Add trajectory arrows for each nest
  #ordiarrows(m, groups = factor(census$NEST))
  # add legend
  legend("topright", title = "Species", legend=levels(fac), col=cols, pch = 16)
  # add stress
  legend("bottomright", bty = "n", cex = 0.8, legend=paste("Stress:", round(m$stress,4)))
  # Draw polygons for each day
  colors<-vector()
  for(i in unique(fac)){
    colors[which(fac==i)]<-cols[which(unique(fac)==i)]
  }
  #Plot convex hulls with colors based on treatment
  for(i in unique(fac)) {
    ordihull(m$point[which(fac==i),],draw="polygon",
             groups=fac[fac==i],
             col=colors[which(fac==i)],label=F) 
  } 
}





data.comm1 <- data.comm[order(data.env$Host, -data.env$tangle.vol.cm3),]
data.env1 <- data.env[order(data.env$Host, -data.env$tangle.vol.cm3),]


source("nested_analysis.R")

nest.summ <-  nested_analysis(community.matrix = data.comm, 
                              save.nested = FALSE,
                              sort.by.size = TRUE,
                              env = data.env)

View(nest.summ)

nest.summ.w <-  nested_analysis(community.matrix = data.comm, 
                                weighted = TRUE, 
                                save.nested = TRUE)



library(metacom)

comm.hell<-decostand(data.comm, "hellinger")

corr <- rda(comm.hell ~ log(tangle.vol.cm3)+Host+cld.num+Condition(Year), 
            scale = FALSE, data = data.env)

# Summary of ordination
corr.sum <- summary(corr)

#Community data of social only
soc.comm <- data.comm[which(data.env$sociality=="social"),]
#Community data of social only
sol.comm <- data.comm[which(data.env$sociality=="solitary"),]

#Matrix of second axis (social)
soc.sites <- corr.sum$sites[which(data.env$sociality=="social"),]
soc.matrix2 <- as.matrix(t(soc.comm[order(soc.sites[,2]),
                                    order(corr.sum$species[,2])]))

#Matrix of second axis (solitary)
sol.sites <- corr.sum$sites[which(data.env$sociality=="solitary"),]
sol.matrix2 <- as.matrix(t(sol.comm[order(sol.sites[,2]),
                                    order(corr.sum$species[,2])]))

sol.matrix2 <- sol.matrix2[rowSums(sol.matrix2)>0,]

#Second axis SOCIAL
Coherence(soc.matrix2, sims = 100)

#Second axis SOLITARY
Coherence(sol.matrix2, sims = 100)



xy.dist <- dist(data.xy, method = "euclidean")
comm.dist <- vegdist(comm.hell[!is.na(data.env$longitude),], method = "bray")

mantel(comm.dist, xy.dist, method ="spearman")


part.all <- varpart(comm.hell1, data.env1, data.space)
part.all
plot(part.all, digits = 2)


# CONNECTIVITY?

d <- dist(cbind(data.env$longitude[!is.na(data.env$longitude)],
                data.env$latitude[!is.na(data.env$longitude)]))

dists <- as.matrix(d)

nnb <- vector()
for(i in 1:length(dists[,1])){
  nnb[i] <- (min(dists[i,-i])*1000)
}
range(nnb) 

alpha <- 0.5
edis <- as.matrix(exp(-alpha * d))
diag(edis) <- 0
A <- rep(1, 138)
edis <- sweep(edis, 2, A, "*")
S <- rowSums(edis)
S



# Modularity

# The cluster coefficient for a network is the average cluster coefficients of its members,
# i.e. simply the number of realised links divided by the number of possible links. 
# Introduced by Watts & Strogatz (1998)
comm.cluscoef <- networklevel(as.matrix(data.comm), index = "cluster coefficient")

nulls.cc <- nullmodel(as.matrix(data.comm), N=1000, method="swap.web")
null.cc <- unlist(sapply(nulls.cc, networklevel, index="cluster coefficient"))

tail95.cc <- sort(null.cc[1,])[round(length(null.cc[1,])*0.95,0)]
pval.cc <- length(which(null.cc[1,]>comm.cluscoef[1]))/length(null.cc[1,])

hist(null.cc[1,], xlim = c(0,0.2), main = "Cluster coefficient") #Overall network
abline(v=comm.cluscoef[1], col = "red", lwd = 2)

hist(null.cc[2,], xlim = c(0.5,0.8), main = "Cluster coefficient") #Highel Level (species)
abline(v=comm.cluscoef[2], col = "red", lwd = 2)

# Asymmetry (higher vs. lower trophic level) of specialisation
# positive values indicate a higher specialisation of the higher trophic level.
networklevel(as.matrix(data.comm), index = "SA", SAmethod = "log")

# Explaining dependence asymmetry is also a measure of specialisation, 
# across both trophic levels. Proposed by Bascompte et al. (2006)
# Positive values indicate higher dependence in the higher trophic level
networklevel(as.matrix(data.comm), index = "ISA")


#New variable: nearest neighbour
d <- as.matrix(dist(cbind(data.env$longitude[!is.na(data.env$longitude)],
                          data.env$latitude[!is.na(data.env$longitude)])))
nnb <- vector()
for(i in 1:length(d[,1])){
  nnb[i] <- (min(d[i,-i])*1000)
}
data.env1$near.nb <- nnb

# Older version of ggplot2
packageurl <- "https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_2.0.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

### Variation partitioning
comm.hell<-decostand(data.comm, "hellinger")

#Create datasets of no missing GPS cordinates (9 missing)
data.env1 <- data.env[!is.na(data.env$longitude),] #Select Host, tangle size and cld.num
row.names(data.env1) <- data.env$Nest[!is.na(data.env$longitude)]
#data.env1$Host <- as.numeric(data.env1$Host) #varpart() takes only numeric

data.space <- data.frame(longitude = data.env$longitude, latitude = data.env$latitude)
data.space <- geoXY(latitude = data.env$latitude, longitude = data.env$longitude, unit = 1)

data.space <- data.space[!is.na(data.env$longitude),]
row.names(data.space) <- row.names(data.env1)

comm.hell1 <- comm.hell[!is.na(data.env$longitude),]

#Figure of variation partitioning (empty)
showvarparts(2)

#Spatial component using polynomials
#Centered xy coordinates
xy.c <- scale(data.space, scale = FALSE)

space.poly <- poly(as.matrix(xy.c), degree = 3, raw=TRUE) #raw=FALSE is orthogonal
colnames(space.poly) <- c("X", "X2", "X3", "Y", "XY", "X2Y", "Y2", "XY2", "Y3")
# RDA of all 9 polynomials (unrealistic)
rda.space <- rda(comm.hell1 ~ ., data = as.data.frame(space.poly))
rda.space

#Forward selection: Which terms of polynomial to keep?
# evaluates variables according to two criteria: 
# 1. if their inclusion into the model leads to significant increase of explained variance
# 2. if the AIC of the new model is lower than AIC of the more simple model. 
ordistep(rda(comm.hell1 ~ 1, data = as.data.frame(space.poly)), scope = formula(rda.space), direction = 'forward')
#Output: keep X2Y + X3 + XY

# Or use selection based on AdjR2?
ordiR2step(rda(comm.hell1 ~ 1, data = as.data.frame(space.poly)), scope = formula(rda.space), direction = 'forward')
#Output: keep X2Y + X3

data.space <- as.data.frame(space.poly)

# Partitioning using spatial terms from ordistep()
mm1 <- model.matrix(~ log(tangle.vol.cm3)+Host+cld.num, data.env1)[,-1]
mm2 <- model.matrix(~ X2Y + X3 + XY, data.space)[,-1]
variation.part <- varpart(comm.hell1, mm1, mm2)
variation.part

indiv.frac <- as.character(round(variation.part$part$indfract[,3], 3))

showvarparts(2, labels = indiv.frac, cex = 1.5,
             bg = c("green","hotpink"), 
             Xnames = "")

text(-0.3, 0.8, "Environmental", cex = 1.5)
text(1.5, 0.8, "Spatial", cex = 1.5)

## Colonization long time scale

census.props <- read.csv("census_props.csv")

census.props$associate <- gsub("Theridiidae01", "Faiditus.sp1", census.props$associate)
census.props$associate <- gsub("Neospintharus01", "Faiditus.sp6", census.props$associate)
census.props$associate <- gsub("Hemiptera01", "Ranzovius", census.props$associate)
census.props$associate <- factor(census.props$associate, 
                                 levels = c("Faiditus.sp1", "Faiditus.sp6", "Ranzovius"))
census.props$host.source <- gsub("forest", "Old webs", census.props$host.source)
census.props$host.source <- gsub("transplant", "New webs", census.props$host.source)
census.props$host.source <- factor(census.props$host.source,
                                   levels = c("New webs", "Old webs"))

# Mosaic plot of communities
census.sums <- rbind(jungle = colSums(census.comm[which(census.env$source=="forest"),]),
                     transplants = colSums(census.comm[which(census.env$source=="transplant"),]))
census.sums <- as.data.frame(census.sums[,order(colSums(census.sums), decreasing = T)])

mosaicplot(census.sums, main = "Compositional difference \nbetween forest and census",
           las = 1, color = rainbow(13), off = 10)

# Differences in proportions of 3 specialized associates
bp <- with(census.props, 
           bargraph.CI(x.factor=associate, group=host.source, response=proportion,
                       lc=FALSE, xlab="Species", ylab = "Average proportion in community", ylim = c(0,0.6),
                       legend=TRUE, x.leg=6, y.leg = 0.6, cex.leg=1, cex.names=1, cex.lab = 1,
                       ci.fun=function(x) {c(mean(x) - 1.96*se(x), mean(x) + 1.96*se(x))}))

text(x=bp$xvals[2], y=bp$CI[4]+0.05, "*", cex = 2.5)
text(x=bp$xvals[4], y=bp$CI[8]+0.05, "*", cex = 2.5)
text(x=bp$xvals[6], y=bp$CI[12]+0.05, "*", cex = 2.5)


census.prop <- census.comm

for(i in 1:length(census.prop[,1])){
  for(j in 1:length(census.prop[1,])){
    census.prop[i,j] <- census.prop[i,j]/rowSums(census.matrix)[i]
  }
}

for(i in 1:length(census.prop[1,])){
  print(colnames(census.prop)[i])
  print(t.test(census.prop[census.env$source=="forest", i], census.prop[census.env$source=="transplant", i]))
}

census.comm <- read.csv("census_comm.csv", row.names = 1)
census.env  <- read.csv("census_env.csv", row.names = 1)

# Plot means new vs old webs

census.means <- read.csv("census_means.csv")
census.means$Source <- factor(census.means$Source, 
                              levels = c("transplant", "jungle"))

#Or order them by proportion in jungle?
census.means$Species <- factor(census.means$Species, 
                               levels = as.character(census.means$Species[
                                 order(census.means$Mean.prop[
                                   census.means$Source=="jungle"], 
                                   decreasing = TRUE)]))

delta.col <- rainbow(13, end = 2/6)[rank(census.means$delta[order(census.means$Mean.prop[1:13], decreasing = TRUE)])]


#Plot mean proportions and CIs
g<-ggplot(census.means, aes(Source, Mean.prop, colour=Species, group = Species))
g+geom_line() +
  geom_errorbar(width=.1, aes(ymin=Mean.prop-ci.prop, ymax=Mean.prop+ci.prop)) +
  geom_point(shape=21, size=3, fill="white")+
  xlab(" ")+
  ylab("Mean proportion in webs")+
  scale_x_discrete(labels= c("New webs", "Old webs"))+
  theme_classic()



#Plot mean proportions and CIs
g<-ggplot(census.means, aes(Source, Mean.prop, colour=Species, group = Species))
g+geom_line() +
  geom_errorbar(width=.1, aes(ymin=Mean.prop-ci.prop, ymax=Mean.prop+ci.prop)) +
  geom_point(shape=21, size=3, fill="white")+
  xlab(" ")+
  ylab("Mean proportion in webs")+
  scale_x_discrete(labels= c("New webs", "Old webs"))+
  theme_classic()

#Plot asin-transformed mean proportions and CIs
g<-ggplot(census.means, aes(Source, Mean.prop.asin, colour=Species, group = Species))
g+geom_line() +
  geom_errorbar(width=.1, aes(ymin=Mean.prop.asin-ci.prop.asin, ymax=Mean.prop.asin+ci.prop.asin)) +
  geom_point(shape=21, size=3, fill="white")+
  xlab(" ")+
  ylab("Mean proportion in webs (arcsin-transformed)")+
  scale_x_discrete(labels= c("New webs", "Old webs"))+
  theme_classic()

#Different color for specialized
g<-ggplot(census.means, aes(Source, Mean.prop, colour=Specialized, group = Species))
g+geom_line() +
  geom_errorbar(width=.1, aes(ymin=Mean.prop-ci.prop, ymax=Mean.prop+ci.prop)) +
  geom_point(aes(shape=Species), size=5)+
  scale_shape_manual(name  ="Species",
                     values=c(0,11,2,5,15,16,6,7,9,10,1,17,18))+
  scale_colour_manual(values=c("grey70", "grey30"))+
  xlab(" ")+
  ylab("Mean proportion")+
  scale_x_discrete(labels= c("New webs", "Old webs"))+
  theme_classic()

census.means$Species[(census.means$Mean.prop[census.means$Source=="jungle"]-census.means$ci.prop[census.means$Source=="jungle"])>
                       (census.means$Mean.prop[census.means$Source=="transplant"]+census.means$ci.prop[census.means$Source=="transplant"])]


#Plot on ordination
m<-metaMDS(census.comm, distance = "bray") #Gower is good on normalized comm data, BUT NOT FOR COMMUNITIES (it's simmetrical)
plot(m, type="t", display = c("species")) #To see the species NMDS axes

ordiplot(m, type = "n")
points(m, col = rainbow(2)[factor(census.env$source)], pch = 16)
legend("topright", title = "Source", legend=row.names(census.sums), 
       cex = 0.7, col=rainbow(2), pch = 16)
legend("bottomright", bty = "n", cex = 0.8, legend=paste("Stress:", round(m$stress,4)))


# SURVIVAL ANALYSIS
library(survival)
census.24 <- read.csv("census_time.csv")

census.24 <- census.24 %>%
  gather("species", "abundance", 1:13)

presence <- vector()
for(i in 1:length(census.24$abundance)){
  if(census.24$abundance[i]>0){
    presence[i] <- 1
  } else{
    presence[i] <- 0
  }
}
census.24$presence <- presence

surv <- Surv(census.24$DAY[census.24$species=="Faiditus.all"], 
             census.24$presence[census.24$species=="Faiditus.all"])
surv.fit <- survfit(surv ~ 1)

plot(surv.fit, fun = "event")

surv <- Surv(census.24$DAY, census.24$presence)
surv
surv.fit <- survfit(surv ~ census.24$species)
surv.fit

plot(surv.fit, fun = "event")


# Mean first day of colonization
census.24 <- read.csv("census_time.csv")

spec <- vector()
first.col <- vector()
for(i in 1:13){
  first.col.app <- vector()
  for(j in 1:length(unique(census.24$NEST))){
    first.col.app[j] <- census.24$DAY[which(census.24[
      census.24$NEST==as.character(unique(census.24$NEST)[j]), i]>0)[1]]
  }
  first.col <- append(first.col, first.col.app) 
  spec <- append(spec, rep(colnames(census.24)[i], 25))
  
}

colonize <- data.frame(spec, first.col)
colonize$first.col[is.na(colonize$first.col)] <- 24

write.csv(colonize, "first_colo.csv")

# Colonization odds used in Altermatt et al (2008)
# Colonization odds are the number of available empty habitats that were colonized 
# divided by those that were not colonized. 
# For example, if there were 100 available habitats and 20 of these got colonized, 
# colonization odds were as follows: 20/80=0.25

census.24 <- read.csv("census_time.csv")

col.odds <- vector()

for(i in 1:13){
  num.col <- 0
  for(j in as.character(unique(census.24$NEST))){
    if(sum(census.24[which(census.24$NEST==j),i])>0){
      num.col <- num.col+1
    } 
  }
  col.odds[i] <- num.col/(length(unique(census.24$NEST))-num.col)
}

spec <- colnames(census.24)[c(1,13)]



# TREND IN COLONIZATION? 

census.24 <- read.csv("census_time.csv")

rda(census.24[1:13]~DAY, data = census.24)

#### SPECIES DIFFERENTIATION ####

source("ordination_plot.R")

# RDA axes explain more variance and species plot looks better
#Hellinger transformation to reduce the importance of large abundances
comm.hell<-decostand(data.comm, "hellinger")

#Only environmental variables
corr1 <- rda(comm.hell ~ log(tangle.vol.cm3)+Host+cld.num+sociality+Year, 
             scale = FALSE, data = data.env)

step.forward <- ordistep(rda(comm.hell ~ 1, scale = FALSE, data = data.env),
                         scope = formula(corr1), direction = "forward", pstep = 1000)

#With spatial ones (?)
corr2 <- rda(comm.hell1 ~ log(tangle.vol.cm3)+Host+sociality+basket.height+cld.num+
               data.space$X2+data.space$X2Y+data.space$XY+near.nb+Year, 
             scale = FALSE, data = data.env1)


step.forward <- ordistep(rda(comm.hell1 ~ 1, scale = FALSE, data = data.env1),
                         scope = formula(corr2), direction = "forward", pstep = 1000)

# Forward-selected model (plus covariate Year):
corr <- rda(comm.hell ~ Host+log(tangle.vol.cm3)+Condition(Year), data = data.env)
corr
anova(corr)


#Nested anova for sociality nested within Host
nested.npmanova(comm.hell ~ sociality+Host, data = data.env, permutations = 1000)
#Signif.

#Nested anova for cld nested within Host
nested.npmanova(comm.hell ~ Host+CLD, data = data.env, permutations = 100)
#NS

#With spatial (but not needed to answer Question 2)
corr <- rda(comm.hell1 ~ Host+log(tangle.vol.cm3)+data.space$X2Y+Condition(Year), data = data.env1)
corr

# Test overall model fit
anova(corr)

plot(corr)

source("ordination_plot.R")

rda.plot<-ordination_plot(ordi = corr, 
                          group = data.env$Host,
                          polygons=FALSE, 
                          plot.species = TRUE, 
                          sp.bold = TRUE)


plot(corr, type="t", display = c("sp","bp")) #To see the species axes



# Summary of ordination
corr.sum <- summary(corr)

RsquareAdj(corr)$adj.r.squared

# Test overall model fit
anova(corr, strata = data.env$Host)

# Kaiser-Guttman criterion to residual axes (Axes above the mean are significant)
corr$CCA$eig[corr$CCA$eig > mean(corr$CCA$eig)]

# Test significance of axes
anova.cca(corr, by = "axis")

# Test significance of terms
anova.cca(corr, by = "terms")
anova(rda.full, by="margin")

#PERMANOVA
adonis(comm.hell~Host+log(tangle.vol.cm3),
       method = "euclidean", data = data.env) 

# Which species drive the community differences?

# "negative binomial" doesn't look skewed, so use it in glm

mv.comm <- mvabund(data.comm)
mv.fit <- manyglm(mv.comm ~ log(tangle.vol.cm3)+Host, 
                  data = data.env, family = "negative.binomial")

plot(mv.fit)

#Perform glm on each species
#Takes a long time (15-25 min)
mv.glm.simple <- anova.manyglm(mv.fit, p.uni="none")
mv.glm.allp <- anova.manyglm(mv.fit, p.uni="unadjusted")


# Graph abundances of indicative species
comm.gg <- data.frame(Nest = rep(data.env$Nest, 10), 
                      Host = rep(data.env$Host, 10), 
                      Associate = c(rep(colnames(data.comm)[4],150),
                                    rep(colnames(data.comm)[5],150),
                                    rep(colnames(data.comm)[11],150),
                                    rep(colnames(data.comm)[13],150),
                                    rep(colnames(data.comm)[15],150),
                                    rep(colnames(data.comm)[16],150),
                                    rep(colnames(data.comm)[17],150),
                                    rep(colnames(data.comm)[18],150),
                                    rep(colnames(data.comm)[19],150),
                                    rep(colnames(data.comm)[20],150)),
                      Abundance.log = log10(c(data.comm[,4],
                                              data.comm[,5],
                                              data.comm[,11],
                                              data.comm[,13],
                                              data.comm[,15],
                                              data.comm[,16],
                                              data.comm[,17],
                                              data.comm[,18],
                                              data.comm[,19],
                                              data.comm[,20])+1))

#Wrap by associate                            
g<-ggplot(comm.gg, aes(Host, Abundance.log))
g+geom_bar(fill = "grey20", stat = "identity", aes(group = Host))+ #, outlier.shape = NA)+
  ylab("Mean abundance (log10-transformed)")+
  facet_wrap(~Associate, ncol = 5)+
  theme_bw()+
  theme(axis.title.y  = element_text(size = 17),
        axis.title.x  = element_text(size = 17))+
  theme(axis.text.x  = element_text(angle=80, vjust=0.5, hjust = 0.5, size=12))


#Wrap by host                            
g<-ggplot(comm.gg, aes(Associate, Abundance.log))
g+geom_bar(fill = "grey70", stat = "identity", aes(group = Host))+
  ylab("Mean abundance (log10-transformed)")+
  facet_wrap(~Host, ncol = 1)+
  theme_bw()+
  theme(axis.title.y  = element_text(size = 17),
        axis.title.x  = element_text(size = 17))+
  theme(axis.text.x  = element_text(angle=60, vjust = 1.1, hjust = 1.25, size=12))



#TABLE of effects of env. variables on each species?
#For only some species - mvabund significant ones?
faid.rda <- rda(comm.hell$Faiditus.sp2 ~ Host+log(tangle.vol.cm3)+cld.num, data = data.env)
faid.rda

anova(faid.rda, by = "terms")




### Matt Barbour code:

rda.0 <- rda(comm.hell ~ 1, data = data.env)
rda.full <- rda(comm.hell ~ log(tangle.vol.cm3)+Host+cld.num+Year, data = data.env)
RsquareAdj(rda.full)
anova(rda.full, by="margin")

rda.R2step <- ordiR2step(rda.0, scope=formula(rda.full))
RsquareAdj(rda.R2step)
plot(rda.R2step, display=c("sp","bp"))
rda.R2step.anova <- anova(rda.R2step, by="margin", step = 1000)
anova(rda.R2step, step = 1000)

rda.R2step.tab <- data.frame(Response = c("Community composition",rep("",1)), 
                             Variable = c("Host","log_size"), 
                             coef.mean = rep("",2), coef.se = rep("",2), P = rda.R2step.anova$"Pr(>F)"[-3], 
                             partial.r2 = round((rda.R2step.anova$Var/sum(rda.R2step.anova$Var))[-3],3)) 

#maps
map(regions = "Ecuador", xlim=c(-77.635,-77.61), ylim=c(-1.09,-1.06), 
    fill = TRUE, col = "grey90") 
map.scale()
points(x = data.env$longitude, y = data.env$latitude, 
       pch = as.numeric(data.env$Host)+19, cex = 1.5, bg = "white")
legend("topright", title = "Host", 
       legend=levels(data.env$Host), 
       pch = unique(as.numeric(data.env$Host))+19, cex = 1, pt.cex = 2)

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


hist(null.nodf, xlim = c(0,50), main = "WNODF")
abline(v=comm.nodf, col = "red", lwd = 2)

# NTC
comm.temp <- as.numeric(nestedtemp(data.comm)$statistic)
comm.temp

nulls.temp <- nullmodel(as.matrix(data.comm), N=100, method="swap.web")
null.temp <- vector()
for(i in 1:length(nulls.temp)){
  null.temp[i] <- as.numeric(nestedtemp(nulls.temp[[i]])$statistic)
}

tail95.temp <- sort(null.temp[1,])[round(length(null.temp[1,])*0.95,0)]
pval.temp <- length(which(null.temp[1,]>comm.temp[1]))/length(null.temp[1,])

hist(null.temp, xlim = c(0,40), main = "Temperature")
abline(v=comm.temp, col = "red", lwd = 2)

# BINMATNEST

comm.temp <- nestedness(as.matrix(data.comm))

# Determine two orders to compare:
nest.order<-comm.temp$pack.order.row    #Nestedness function order
size.order<-order(data.env$tangle.vol.cm3, decreasing = TRUE)   #Nest size order (decreasing top-down)
corr <- cor.test(nest.order,size.order)

corr$estimate
corr$p.value


#This function calculates the nestedness (when rows ordered by size) of each host
#Takes some time (20min)
source("nested_analysis2.R")

nest.summ <-  nested_nodf_size(community.matrix = data.comm, 
                               env = data.env,
                               save.nested = FALSE)

View(nest.summ)


#PARTIAL RDA of Rarefied communities
median(rowSums(data.comm))
rare.min <- 28
data.env2 <- data.env#[-which(rowSums(data.comm)<=rare.min),]
data.comm2 <- data.comm#[-which(rowSums(data.comm)<=rare.min),]
rare.comm <- rrarefy(data.comm2, sample = rare.min)       #Rarefied communities of each nest

#Host
corr <- rda(rare.comm ~ Host+Condition(log10(tangle.vol.cm3))+Condition(Year), data = data.env2)
corr
anova(corr)

#Size
corr <- rda(rare.comm ~ log10(tangle.vol.cm3)+Condition(Host)+Condition(Year), data = data.env2)
corr
anova(corr)

# Test for variance in web size
mix.null <- lmer(log10(data.env$tangle.vol.cm3) ~  (1 | Host), data = data.env, REML=FALSE)
summary(mix.null)

mix.lm <- lmer(log10(data.env$tangle.vol.cm3) ~ sociality + (1 | Host), data = data.env, REML=FALSE)
summary(mix.lm)

anova(mix.null, mix.lm)

#Test difference in variance OF WEB SIZE solitary vs. social
var.test(log(data.env$tangle.vol.cm3[data.env$sociality=="social"]),
         log(data.env$tangle.vol.cm3[data.env$sociality=="solitary"]))

var.test(log(data.env$tangle.vol.cm3)~sociality, data = data.env)

var(log10(data.env$tangle.vol.cm3[data.env$Host=="A.eximius"]))


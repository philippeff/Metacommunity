#### VARIATION PARTITIONING ####

comm.hell<-decostand(data.comm, "hellinger")


#Create datasets of no missing GPS cordinates (9 missing)
data.env1 <- data.env[data.env$Year==2015,]
data.env1 <- data.env1[!is.na(data.env1$longitude),] #Select Host, tangle size and cld.num
row.names(data.env1) <- data.env1$Nest

data.space <- geoXY(latitude = data.env1$latitude, longitude = data.env1$longitude, unit = 1)
row.names(data.space) <- row.names(data.env1)

comm.hell1 <- comm.hell[data.env$Year==2015 & !is.na(data.env$longitude),]


# PCNM

dis <- dist(data.space)
pcnm1 <- pcnm(dis)
data.pcnm <- as.data.frame(pcnm1$vectors)

#RDA of all PCNM vectors
rda.space <- rda(comm.hell1 ~ ., data = data.pcnm)

# Or use selection based on AdjR2
ordiR2step(rda(comm.hell1 ~ 1, data = data.pcnm), 
           scope = formula(rda.space), direction = 'forward')
#Output: keep PCNM4 + PCNM21

#visualize each pcnm vector
ordisurf(data.space, scores(pcnm1, choi=4), bubble = 4, main = "PCNM 4")
ordisurf(data.space, scores(pcnm1, choi=21), bubble = 4, main = "PCNM 21")

# Partitioning using spatial terms from ordistep()
mm1 <- model.matrix(~ log(tangle.vol.cm3)+Host, data.env1)[,-1]
mm2 <- model.matrix(~ PCNM4 + PCNM21, data.pcnm)[,-1]
variation.part <- varpart(comm.hell1, mm1, mm2)
variation.part

indiv.frac <- as.character(round(variation.part$part$indfract[,3], 3))

showvarparts(2, labels = indiv.frac, cex = 1.5,
             bg = c("green","hotpink"), 
             Xnames = "")

#Funky way of graphing
indiv.frac1<-round(as.numeric(indiv.frac),3)
indiv.frac1[2]<-0.02

grid.newpage()
grid.rect()
venn <- draw.pairwise.venn(area1 = indiv.frac1[1]+indiv.frac1[2], 
                           area2 = indiv.frac1[2]+indiv.frac1[3],
                           cross.area = indiv.frac1[2], 
                           category = c("Environmental", "Spatial"), lwd = 1, 
                           fontfamily = "sans", cat.fontfamily = "sans",
                           fill = c("light green", "light pink"), cex = c(4,0,1.5),
                           cat.pos = c(355, 25), cat.dist = c(-0.06, 0.03), cat.cex = c(2,1.5))
grid.text(label = paste("Residuals =", indiv.frac1[4]), 0.8, 0.05, gp = gpar(cex = 1.5))
grid.text(label = "0.007", 0.77, 0.5, gp = gpar(cex = 1.5))

# test fraction [a] (environmental)
rda.result <- rda(comm.hell1 ~ mm1 + Condition(mm2))
anova(rda.result, step=200, perm.max=200)

# test fraction [b] (spatial)
rda.result <- rda(comm.hell1 ~ mm2 + Condition(mm1))
anova(rda.result, step=200, perm.max=200)


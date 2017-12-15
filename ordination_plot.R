
### ORDINATION PLOT 

ordination_plot <- function(ordi, group, species, polygons=FALSE, plot.species = FALSE, 
                            sp.bold = FALSE, label.site = NULL, env = data.env){
  
  require(vegan)
  
  if(length(grep("rda", class(ordi)))==0 & 
       length(grep("cca", class(ordi)))==0 &
       length(grep("capscale", class(ordi)))==0){
    stop("argument is not of class \"rda\" or \"cca\" or \"capscale\"")
  }
  
  if(missing(ordi)){
    stop("Need to include ordination object as argument")
  }
  
  if(missing(group)){
    stop("Must specify argument \"group\" as a column in env data of same length as comm data")
  } 
  
  sum.ordi <- base::summary(ordi)
  species <- species[-which(species$type=="generalist"),]
  
  perc.ca1 <- round(100*sum.ordi$cont$importance[2,1], 2) #Percent of inertia explained by CA1
  perc.ca2 <- round(100*sum.ordi$cont$importance[2,2], 2) #Percent of inertia explained by CA2
  
  #ncolors<-c("#FF0000FF","#8000FFFF","#80FF00FF","#00FFFFFF") #Determine the colors (original)
  ncolors<-c("#FF0000FF","#8000FFFF","#4CBF00FF","#00BAFFFF") #Determine the colors 
  
  if(plot.species == TRUE){
    par(mfrow=c(1,2))
  }
  
  
  xl <- c(min(sum.ordi$sites[,1])-0.25,max(sum.ordi$sites[,1]))
  yl <- c(min(sum.ordi$sites[,2]),max(sum.ordi$sites[,2]))
  
  plot(ordi, type ="n", 
       xlab = paste("RDA 1: ", perc.ca1, "%", sep = ""),
       xlim = xl,
       ylab = paste("RDA 2: ", perc.ca2, "%", sep = ""),
       ylim = yl)
  
  points(ordi, col = ncolors[group], pch = 16)
  
  # Draw polygons
  if(polygons==TRUE){
    for(i in unique(group)) {
      ordihull(ordi,draw="polygon", show.groups = i, groups=group,
               border=ncolors[which(levels(group)==i)],label=FALSE) 
    }
  }
  
  legend("topleft", title = "Host", legend=levels(group), bty = "n", col=ncolors, pch = 16, cex = 0.7)

  # label each nest
  if(!is.null(label.site)){
    if(length(label.site)!=length(group)){
      stop("Argument label.site is not the same length as number of sites")
    } else{
      text(ordi, labels = label.site, cex = 0.5) 
    }
  }
  
  intersect(rownames(sum.ordi$species),as.character(species$species[which(species$mvabund.pval<0.05)]))

  
  
  # Ordination plot of species
  if(plot.species==TRUE){
    if(sp.bold==TRUE){
      mvab.sp <- which(species$mvabund.pval<0.05) #Which species drive the difference in mvabund()?
      plot(ordi, type = "n", display = c("sp"), 
           xlim = xl,xlab = paste("RDA 1: ", perc.ca1, "%", sep = ""),
           ylim = yl, ylab = "")
      points(sum.ordi$centroids[,c(1,2)], pch = "+", col = ncolors, cex = 3)
      text(x = sum.ordi$species[mvab.sp,1], 
           y = sum.ordi$species[mvab.sp,2], 
           labels = row.names(sum.ordi$species)[mvab.sp], cex = 0.7, font = 3, col = "red") 
      text(x = sum.ordi$species[-mvab.sp,1], 
           y = sum.ordi$species[-mvab.sp,2], 
           labels = row.names(sum.ordi$species)[-mvab.sp], cex = 0.7, font = 3)
      arrows(x0 = 0, y0 = 0, x1 = sum.ordi$biplot[4,1], y1 = sum.ordi$biplot[4,2],
             length = 0.2, angle = 20, col = "black", lwd = 2)
      text(x = sum.ordi$biplot[4,1]-0.2, y = sum.ordi$biplot[4,2]+0.02,labels = "Patch size")

      
    }else{
      plot(ordi, type="t", display = c("species"), 
           xlim = xl,
           ylim = yl)
    }
  }
  
  par(mfrow=c(1,1))
  
}
  
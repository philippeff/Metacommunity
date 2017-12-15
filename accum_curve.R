

accum_curve <- function(comm = data.comm, env = data.env, group = data.env$Host){
  
  require(vegan)
  
  ncolors<-c("#FF0000FF","#8000FFFF","#4CBF00FF","#00BAFFFF") #Determine the number of colors 
  
  plot(specaccum(comm[which(group=="A.domingo"),], method = "rarefaction")$richness, 
       type ="n", 
       xlab = "Accumulated number of webs", 
       xlim = c(0, max(specaccum(comm[which(group=="A.domingo"),])$sites)),
       ylab = "Accumulated number of species",
       ylim = c(0, max(specaccum(comm[which(group=="A.domingo"),])$richness))
  )
  
  for(i in levels(group)){
    sp.accum<-specaccum(comm[which(group==i),], method = "rarefaction")
    lines(sp.accum$sites, sp.accum$richness, 
          col = ncolors[which(levels(group) %in% i)],
          lwd = 2)
    legend("bottomright", title = "Host", legend=levels(group), 
           col=ncolors, pch = 15)
  }
  
}
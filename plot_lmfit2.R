

### Function to plot variable

plot_lmfit <- function(var, ylabel = "Variable", int.term = TRUE, dat = data.env){
  
  require(ggplot2)
  require(multcompView)
  
  fit1 <- lm(var ~ log10(tangle.vol.cm3)+Host, data = dat)
  anov.fit1 <- aov(fit1) 
  tuke <- TukeyHSD(anov.fit1, which = "Host")
  coeffs <- coef(fit1)

  letters <- multcompLetters(extract_p(tuke)$Host)$Letters
  letters <- c(letters[4], letters[1], letters[2], letters[3])
  labs <- paste(names(letters), " ", letters)

  #print(tuke$Host)
  
  labs[1]<-gsub("  ", "     ", labs[1])
  labs[2]<-gsub("  ", "          ", labs[2])
  labs[3]<-gsub("  ", "       ", labs[3])
  labs[4]<-gsub("  ", "        ", labs[4])

  #Graph
  
  cols <- c("#FF0000FF","#8000FFFF","#4CBF00FF","#00BAFFFF")
  
  if(int.term == TRUE){
    g<-ggplot(dat, aes(log10(tangle.vol.cm3), var))+#, environment = environment())+
      geom_point(aes(colour=Host), size=2)+
      geom_smooth(fullrange=FALSE, aes(colour=Host), 
                  method="lm", size=1, fill = "grey75")+
      scale_colour_manual(values=cols, labels = labs) +
      #geom_abline(intercept = 0) +
      xlab("Tangle volume in cm3 (log10-transformed)")+
      ylab(ylabel)+
      theme_classic()+
      theme(legend.position="none")
      #theme(axis.title=element_text(size=18))
  }else{
    g<-ggplot(dat, aes(log10(tangle.vol.cm3), var), environment = environment())+
      geom_point(aes(colour=Host), size=2)+
      scale_colour_manual(values=cols, labels = labs) +
      geom_abline(intercept = coeffs[1], slope = coeffs[2], 
                  size = 1, colour = cols[1]) +
      geom_abline(intercept = coeffs[1]+coeffs[3], slope = coeffs[2], 
                  size = 1, colour = cols[2]) +
      geom_abline(intercept = coeffs[1]+coeffs[4], slope = coeffs[2], 
                  size = 1, colour = cols[3]) +
      geom_abline(intercept = coeffs[1]+coeffs[5], slope = coeffs[2], 
                  size = 1, colour = cols[4]) +
      xlab("Tangle volume in cm3 (log10-transformed)")+
      ylab(ylabel)+
      theme_classic()
  }
  
  
  return(g)
}

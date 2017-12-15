nested_nodf_size <- function(community.matrix, 
                            env = data.env,
                            save.nested = FALSE){
  
  require(vegan)

  #needs to source "PERFORM_NESTED_TEST.R" from FALCON package
  source("PERFORM_NESTED_TEST.R") #FALCON package
  
  Host= vector() 
  host.NODF = vector()
  mean.null.NODF= vector()  
  P.value.NODF= vector() 

  for(i in 1:length(unique(env$Host))){
    host <- as.character(unique(env$Host)[i]) # "eximius" "domingo" "mecynogea" "agelenopsis"
    Host[i]<-host
    host.comm <- community.matrix[which(env$Host==host),]
    host.env <- env[which(env$Host==host),]
    host.comm <- host.comm[order(host.env$tangle.vol.cm3, decreasing = TRUE),
                           order(colSums(host.comm), decreasing =  TRUE)]
    host.env <- host.env[order(host.env$tangle.vol.cm3, decreasing = TRUE),]
    #remove row/col sums that are 0
    host.comm <- host.comm[which(rowSums(host.comm)>0), which(colSums(host.comm)>0)]
    
    null.nodf <-PERFORM_NESTED_TEST(MAT = as.matrix(host.comm),
                                    bintest = 1, # binary? 0 = quant 1 = binary
                                    sortVar = 0, # sort matrix?
                                    functhand = "NODF",
                                    nullmodels = c(), # 1=SS, 2=FF, 3=CC, 4=DD, 5=EE, c()=all 
                                    EnsembleNumber = c(),
                                    plotON = 0)
    
    host.NODF[i] <- as.numeric(nestednodf(host.comm, order = FALSE, weighted = FALSE)$statistic[3]) #FALCON uses ths function, so its OK
    
    mean.null.NODF[i] <- null.nodf$Bin_t1$NODF$Mean
    P.value.NODF[i] <- null.nodf$Bin_t1$NODF$pvalue
    

    #Save resulting nested matrix for each host?
    if(save.nested == TRUE){
      write.csv(as.data.frame(host.comm), 
                paste("nested_", host, ".csv", sep = "")) # ordered matrix, Rows and columns have to be >0
    }
    
    print(paste(host, ":"))
    
    print(paste("Mean NODF (SS): ", null.nodf$Bin_t1$NODF$Measure, 
    "Mean null of SS: ", null.nodf$Bin_t1$NODF$Mean, 
    "pval = ", null.nodf$Bin_t1$NODF$pvalue))
    
    print(paste("Mean NODF (FF): ", null.nodf$Bin_t2$NODF$Measure, 
    "Mean null of FF: ", null.nodf$Bin_t2$NODF$Mean, 
    "pval = ", null.nodf$Bin_t2$NODF$pvalue))
    
    print(paste("Mean NODF (CC): ", null.nodf$Bin_t3$NODF$Measure, 
    "Mean null of CC: ", null.nodf$Bin_t3$NODF$Mean, 
    "pval = ", null.nodf$Bin_t3$NODF$pvalue))
    
    print(paste("Mean NODF (DD): ", null.nodf$Bin_t4$NODF$Measure, 
    "Mean null of DD: ", null.nodf$Bin_t4$NODF$Mean, 
    "pval = ", null.nodf$Bin_t4$NODF$pvalue))
    
    print(paste("Mean NODF (EE): ", null.nodf$Bin_t5$NODF$Measure, 
    "Mean null of EE: ", null.nodf$Bin_t5$NODF$Mean, 
    "pval = ", null.nodf$Bin_t5$NODF$pvalue))

  }
  return(data.frame(Host, host.NODF, mean.null.NODF, P.value.NODF))
  
}

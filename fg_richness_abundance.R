
fg_richness <- function(community.matrix, fg.data){
  # RICHNESS of each functional group
  pred <- rep(0,length = length(community.matrix[,1]))
  kleptoA <- rep(0,length = length(community.matrix[,1]))
  kleptoP <- rep(0,length = length(community.matrix[,1]))
  inq <- rep(0,length = length(community.matrix[,1]))
  detrit <- rep(0,length = length(community.matrix[,1]))
  
  for(m in 1:length(community.matrix[,1])){
    list<-colnames(community.matrix)[which(community.matrix[m,]>0)] # returns which species are present at site m
    for(i in 1:length(list)){
      if(fg.data$functional.group[which(fg.data$Species==list[i])]=="kleptoparasite P"){
        kleptoP[m] <- kleptoP[m] + 1
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="kleptoparasite A"){
        kleptoA[m] <- kleptoA[m] + 1
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="predator"){
        pred[m] <- pred[m] + 1
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="inquiline"){
        inq[m] <- inq[m] + 1
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="detritivore"){
        detrit[m] <- detrit[m] + 1
      }
    }
  }
  
  return(data.frame(pred, kleptoA, kleptoP, detrit, inq))
  
}


# ABUNDANCE of each functional group
fg_abundance <- function(community.matrix, fg.data){
  pred <- rep(0,length = length(community.matrix[,1]))
  kleptoA <- rep(0,length = length(community.matrix[,1]))
  kleptoP <- rep(0,length = length(community.matrix[,1]))
  inq <- rep(0,length = length(community.matrix[,1]))
  detrit <- rep(0,length = length(community.matrix[,1]))
  
  for(m in 1:length(community.matrix[,1])){
    list<-colnames(community.matrix)[which(community.matrix[m,]>0)] # returns which species are present at site m
    for(i in 1:length(list)){
      if(fg.data$functional.group[which(fg.data$Species==list[i])]=="kleptoparasite P"){
        kleptoP[m] <- kleptoP[m] + community.matrix[m,which(colnames(community.matrix)==list[i])]
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="kleptoparasite A"){
        kleptoA[m] <- kleptoA[m] + community.matrix[m,which(colnames(community.matrix)==list[i])]
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="predator"){
        pred[m] <- pred[m] + community.matrix[m,which(colnames(community.matrix)==list[i])]
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="inquiline"){
        inq[m] <- inq[m] + community.matrix[m,which(colnames(community.matrix)==list[i])]
      } else if(fg.data$functional.group[which(fg.data$Species==list[i])]=="detritivore"){
        detrit[m] <- detrit[m] + community.matrix[m,which(colnames(community.matrix)==list[i])]
      }
    }
  }
  
  return(data.frame(pred, kleptoA, kleptoP, detrit, inq))
}



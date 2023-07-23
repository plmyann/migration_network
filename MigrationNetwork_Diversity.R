library(dplyr)
library(data.table)
library(taxize)
library(FD)
library(picante)


setwd("./MigrationNetwork2020/Rspace")
######################################################################################################
######################################diversity calculation###################################

####get all the species list

splist <- sapply(sporder, function(od){
  
  cell.dt <- readRDS(paste(od,'rds',sep = '.'))
  
  splist <- colnames(cell.dt[,colnames(cell.dt)!= "CellID" & colnames(cell.dt)!= "NewID"])
  
  splist <- gsub(".", " ", splist, fixed = T)
  
  return(splist)
  
})

specieslist <- do.call(c, splist)

specieslist <- data.table(Spanish_Common_Name = specieslist)

specieslist <- merge(specieslist, spelist, by.x = "Spanish_Common_Name", by.y = "Scientific.name")

colnames(specieslist)[colnames(specieslist)=="English.name"] <- "English_Common_Name"
###Function to calculate the diversity

traits <- Otraits[,c(8:19,24:31,36)]

traits$`Diet-Ver` <- rowSums(traits[,4:7])

traits[,4:7] <- NULL

traits <- traits[!is.na(traits$Scientific),]

traits <- traits[with(traits,order(Scientific)),]

traits <- traits[traits$Scientific %in% specieslist$Spanish_Common_Name,]

specieslist <- specieslist[specieslist$Spanish_Common_Name %in% traits$Scientific,]

traits <- traits[with(traits,order(Scientific)),]

rownames(traits) <- traits$Scientific

######calculate the pairwise distances between using the Gower distance.
##fuzzy variable as the proportion
tabF <- traits[,c(3:8,18,9:15
                  #,16,17
)]
tabFp <- prep.fuzzy(tabF, c(7,7), labels = c("Diet", "ForStrat"))
##quantity and dichotomous variable
tabQ <- data.frame(bodymass = traits[,17], row.names = rownames(traits)) ##body mass
#tabD <- data.frame(pelagic = traits[,16], row.names = rownames(traits)) ##pelagic or not
tabC <- data.frame(clutchsize = traits[,19], row.names = rownames(traits))
###calculate the distance dissmissy
ktabl <- ktab.list.df(list(tabFp, tabQ, tabC))
distrait <- dist.ktab(ktabl, c("F", "Q", "Q"))
#distrait <- gowdis(tabF)
######build the dendrogram using UPGMA clustering approach.
dendro <- hclust(distrait, method = "average")

########build a function to calculate taxonomic, functional and phylogenetic diversity
TFP.cal <- function(df){
  
  options(warn = -1)
  
  df <- df[,df[1,]!=0, drop = F]
  
  newid <- unique(df$NewID)
  
  cellid <- unique(df$CellID)
  
  df <- df[,which(colnames(df)!= "CellID" & colnames(df)!= "NewID"), drop =F]
  
  colnames(df) <- gsub(".", " ", colnames(df), fixed = T)
  
  ###############Diversity calculation###########
  
  phylo.uc <- rep(0,100)
  
  set.seed(2020)
  
  tmpp <- sample(phylotree, 100)
  
  for (i in 1:100) {
    
    tmp <- tmpp[[i]]
    
    tmp$`tip.label` <- sub("_"," ", tmp$`tip.label`)
    
    tmp <- keep.tip(tmp, tmp$`tip.label`[tmp$`tip.label` %in% colnames(df)])
    
    phylo.uc[i] <- pd(df,tmp)$PD[1]
    
  }
  
  DI <- data.frame(CellID = cellid, NewID = newid, uTD = length(df), 
                   uFD = treedive(comm = df, tree = dendro, match.force = T, verbose = F)[1],
                   uPD = mean(phylo.uc))
  
  return(DI)
  
}

splist <- sapply(sporder, function(od)readRDS(paste(od,'rds',sep = '.')), USE.NAMES = T)

tmp <- lapply(1:26, function(x){
  
  dt <- data.frame(splist[[x]])
  
  colnames(dt) <- gsub(".", " ", colnames(dt), fixed = T)
   
  specieslist <- data.table(Spanish_Common_Name = colnames(dt))
  
  specieslist <- merge(specieslist, spelist, by.x = "Spanish_Common_Name", by.y = "Scientific.name")
  
  colnames(specieslist)[colnames(specieslist)=="English.name"] <- "English_Common_Name"
  
  colnames(dt)[3:ncol(dt)] <- specieslist$Spanish_Common_Nam 
  
  colss <- c("CellID","NewID",colnames(dt)[3:ncol(dt)][colnames(dt)[3:ncol(dt)] %in% traits$Scientific])
  
  dt <- dt[,colss]
   
  return(dt)
  
})

### for each order

splist <- tmp

div <- data.frame(splist[[i]])

div1 <- div[,which(colnames(div)!= "CellID" & colnames(div)!= "NewID"), drop = F]

div <- div[rowSums(div1)!=0,]

rm(div1)

divlist <- lapply(1:nrow(div), function(i)TFP.cal(div[i,]))

divlist <- bind_rows(divlist)

#for all merged orders

div <- Reduce(function(...) merge(...,by = c("CellID","NewID"), all=T), splist)

div1 <- div[,which(colnames(div)!= "CellID" & colnames(div)!= "NewID"), drop = F]

div <- div[rowSums(div1)!=0,]

rm(div1)

prefixes <- unique(sub("\\..*", "", grep("[0-9]", colnames(div),value = T)))

div1 <- sapply(prefixes, function(x)rowSums(div[,startsWith(colnames(div), x)]))

div2 <- div %>% select(!starts_with(prefixes))

div <- cbind(div1,div2)

save.image("MNdiv.RData")


######resolve the name issues (a, a1. etc) for further calculation

for(t in 1:length(splist)){
  
  prefixes <- unique(sub("\\..*", "", grep("[0-9]", colnames(splist[[t]]),value = T)))
  
  div1 <- sapply(prefixes, function(x)rowSums(splist[[t]][,startsWith(colnames(splist[[t]]), x)]))
  
  div2 <- splist[[t]] %>% select(!starts_with(prefixes))
  
  if(length(div1) == 0) splist[[t]] <- div2 else splist[[t]] <- cbind(div1,div2)
  
  rm(div1);rm(div2)
  
}


#######add each unprotected cell to see the conservation outputs

#######read the unprotected key nodes

keyUnPA <- read.csv("./Datainput/keynode/unprotnodeFULL.csv")

colnames(keyUnPA)[colnames(keyUnPA)=="ID"] <- "NewID"

keyUnPA$ord <- sub("^Bet","",keyUnPA$ord) %>% sub(pattern = ".csv$",replacement = "")  

keyUnPA <- keyUnPA[keyUnPA$bet>0.001,]
######select the specific order
keyUnPA_suborder <- keyUnPA[keyUnPA$ord==sporder[i],]

######################read the key node cells that were protected  
##################################################################
keyPA <- read.csv("./Datainput/keynode/protkeynode.csv")

colnames(keyPA)[colnames(keyPA)=="ID"] <- "NewID"
###delete the irrelevant charaters of the order names
keyPA$ord <- sub("Bet","",keyPA$ord) %>% sub(pattern = ".csv$",replacement = "")  

######select the specific order
keyPA_suborder <- keyPA[keyPA$ord==sporder[i],]

###########merge the species list for each protected key node cells
keyPA_suborder <- merge(keyPA_suborder, splist[[i]], all.x = T, all.y = F, by = "NewID")

##############remove irrelevant columns
keyPA_suborder$X <- keyPA_suborder$bet <- keyPA_suborder$ord <- NULL

########remove the species not listed for this cells

keyPA_suborder[nrow(keyPA_suborder)+1,] <- colSums(keyPA_suborder)

#######remove the empty columns
keyPA_suborder <- keyPA_suborder[,keyPA_suborder[nrow(keyPA_suborder),]!=0]

######calcualte the diversity for protected key node as background diversity
div0 <- TFP.cal(keyPA_suborder[nrow(keyPA_suborder),])

keyPA_suborder <- keyPA_suborder[1:(nrow(keyPA_suborder)-1),] 

div0$CellID <- div0$NewID <- 0

#######read the diversity for each cell

divcell <- readRDS(paste("./Rspace/Div_EachCell/o",i,".rds",sep = ""))

#divcell <- setDT(divcell)  

divcell <- divcell[divcell$NewID %in% keyUnPA_suborder$NewID,]

divclass <- c("uTD","uFD","uPD")

for(ttt in divclass){
  
  div1 <- div0
  
  divcell <- divcell[order(-divcell[,ttt]),] 
  
  ###########add each cell to see the diversity change
  
  keyPA_id <- data.frame(NewID = keyPA_suborder$NewID, CellID = keyPA_suborder$CellID)
  
  for(t in 1:nrow(divcell)){
    
    keyPA_id[nrow(keyPA_id)+1,] <- cbind(divcell$NewID[t],divcell$CellID[t])
    
    keyPA_subb <- merge(keyPA_id, splist[[i]],all.x = T, all.y = F, by = c("NewID","CellID"))
    
    keyPA_subb[nrow(keyPA_subb)+1,] <- colSums(keyPA_subb)
    
    keyPA_subb <- keyPA_subb[,keyPA_subb[nrow(keyPA_subb),]!=0]
    
    div1[nrow(div1)+1,] <- TFP.cal(keyPA_subb[nrow(keyPA_subb),])
    
    keyPA_subb <- keyPA_subb[1:(nrow(keyPA_suborder)-1),] 
    
  }
  
  write.csv(div1,paste("Basedon_",ttt,".csv",sep = ""))
  
}


########for the betweenness and random ones

colnames(rd999)[colnames(rd999)=="ID"] <- "NewID"

colnames(rd999)[5:length(rd999)] <- paste("R",1:999,sep = "")

rd999_sub <- rd999[rd999$ord == sporder[i],]

rd999_sub <- rd999_sub[rd999_sub$NewID %in% keyUnPA_suborder$NewID,]

divclass2 <- colnames(rd999_sub)[4:ncol(rd999_sub)]

for(ttt in divclass2){
  
  div1 <- div0
  
  rd999_sub <- rd999_sub[order(-rd999_sub[,ttt]),] 
  
  ###########add each cell to see the diversity change
  
  keyPA_id <- data.frame(NewID = keyPA_suborder$NewID)
  
  for(t in 1:nrow(rd999_sub)){
    
    keyPA_id[nrow(keyPA_id)+1,] <- cbind(rd999_sub$NewID[t])
    
    keyPA_subb <- merge(keyPA_id, splist[[i]],all.x = T, all.y = F, by = "NewID")
    
    keyPA_subb[nrow(keyPA_subb)+1,] <- colSums(keyPA_subb)
    
    keyPA_subb <- keyPA_subb[,keyPA_subb[nrow(keyPA_subb),]!=0]
    
    div1[nrow(div1)+1,] <- TFP.cal(keyPA_subb[keyPA_subb$NewID==rd999_sub$NewID[t],])
    
    keyPA_subb <- keyPA_subb[1:(nrow(keyPA_suborder)-1),] 
    
  }
  
  write.csv(div1,paste("Basedon_",ttt,".csv",sep = ""))
  
}











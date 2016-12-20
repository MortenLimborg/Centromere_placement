##############################################
###### Place centromeres  with haploids ######
##############################################

rm(list=ls())

setwd("C:/Users/.../") #Write here the path to the folder where you have the input data file

#Read in the code for the function 'plot.pseudo.y()' below and then run the function on individual LGs to test if it's working. 
#Specifying plot=TRUE will create the accompanying plots seen in the paper
plot.pseudo.y <- function(LG, plot){
  
  ##Read in map with phased genotypes:
  Map <- as.matrix(read.csv(file = 'datafile_phased_haplotypes.csv', sep = ',')) #Read input file (e.g. 'datafile_phased_haplotypes.csv') with map info and haploid genotypes 
  is.matrix(Map)
  rownames(Map) <- Map[,1]                #Change column 1 into rownames
  Map <- Map[,2:ncol(Map)]        
  
  sample.size <- nrow(Map)-5
  
  #Making subset matrix for LG_X:
  LG_X <- NULL
  for (i in 1:length(Map[1,])){
    if (Map[1,i] == LG){  
      LG_X <- cbind(LG_X, Map[,i])
    }
  }
  
  ##Make new matrix where a & b phases are recoded to 1 and 0:
  XO.vec.kw <- matrix(nrow=nrow(Map), ncol=ncol(LG_X))  #Make Cross-over (XO) vector (vec) using kernel-window (kw) for focal LG
  XO.vec.kw[1:5,] <- LG_X[1:5,]
  rownames(XO.vec.kw) <- rownames(LG_X)
  
  for (j in 6:nrow(XO.vec.kw)){                 #This loop starts looking for data in row 6, which is why you need the top 5 columns in the input file
    for (i in 1:ncol(XO.vec.kw)){
      if (as.character(LG_X[j,i]) == "a"){
        XO.vec.kw[j,i] <- 0  
      }
      if (as.character(LG_X[j,i]) == "b"){
        XO.vec.kw[j,i] <- 1  
      }
    }
  }
  
  #Impute first and last genotype if missing:
  for (j in 6:nrow(XO.vec.kw)){
    first.ten.phases <- as.numeric(XO.vec.kw[j,1:10])
    start.phase <- round(mean(first.ten.phases, na.rm=T), digits=0)
    if (as.character(LG_X[j,1]) != start.phase){
      XO.vec.kw[j,1] <- start.phase  
    }
    
    last.ten.phases <-  as.numeric(XO.vec.kw[j,(ncol(XO.vec.kw)-9):ncol(XO.vec.kw)])
    end.phase <- round(mean(last.ten.phases, na.rm=T), digits=0)
    if (as.character(LG_X[j,ncol(XO.vec.kw)]) != end.phase){
      XO.vec.kw[j,ncol(XO.vec.kw)] <- end.phase  
    }
  } 
      
  #Make sub-matrix with only genotypes and change NA to "-":
  XO.phased.genot <- XO.vec.kw[6:nrow(Map),]
  
  for (j in 1:nrow(XO.phased.genot)){
    for (i in 1:ncol(XO.phased.genot)){
      if (is.na(XO.phased.genot[j,i])){
        XO.phased.genot[j,i] <- "-"  
      }
    }
  }
  
  #Sort matrix by first marker and count number of "0"
  XO.phased.genot <- as.data.frame(XO.phased.genot)
  colnames(XO.phased.genot) <- paste("marker",c(1:ncol(XO.phased.genot)), sep="") 
  XO.phased.genot <- XO.phased.genot[order(XO.phased.genot$marker1),] #Sorted by first marker
  XO.phased.genot.right <- XO.phased.genot[order(XO.phased.genot[,ncol(XO.phased.genot)]),] #Sorted by last marker
  
  #Find number of individuals having 0 at first marker
  no.zero <- sum(XO.phased.genot$marker1 == 0)
  no.zero.right <- sum(XO.phased.genot[,ncol(XO.phased.genot)] == 0)
  
  #Convert XO.phased.genot back to matrix in order to edit numeric contents
  XO.phased.genot <- as.matrix(XO.phased.genot)
  XO.phased.genot.right <- as.matrix(XO.phased.genot.right)
  #is.matrix(XO.phased.genot.right)
  
  
  ##Make new temporary matrix for the ~50% offspring that needs to get phases switched (i.e. having 0 at marker1)
  
  #From left:
  Matrix.convert <- NULL
  converted.genot <- NULL #Temporary vector for each individual to be used within loop
  
  for (j in (no.zero+1):nrow(XO.phased.genot)){  #Loop over individuals having the genotype 1 at the first marker
    for (i in 1:ncol(XO.phased.genot)){ #For individual j: Loop over all markers in LG_X
      converted.genot <- cbind(converted.genot, abs(as.numeric(XO.phased.genot[j,i])-1))
      
    }
    Matrix.convert <- rbind(Matrix.convert, converted.genot)
    converted.genot <- NULL
  }
  
  #From right:
  Matrix.convert.right <- NULL
  converted.genot.right <- NULL #Temporary vector for each individual to be used within loop
  
  for (j in (no.zero.right+1):nrow(XO.phased.genot.right)){  #Loop over individuals having the genotype 1 at the last marker
    for (i in 1:ncol(XO.phased.genot.right)){ #For individual j: Loop over all markers in LG_X
      converted.genot.right <- cbind(converted.genot.right, abs(as.numeric(XO.phased.genot.right[j,i])-1))
      
    }
    Matrix.convert.right <- rbind(Matrix.convert.right, converted.genot.right)
    converted.genot.right <- NULL
  }
  
  
  ##Merge synced phased offspring into final matrix for counting cross-over frequencies
  #From left:
  phased.matrix <- rbind(XO.phased.genot[1:no.zero,], Matrix.convert)
  rownames(phased.matrix) <- rownames(XO.phased.genot)
  
  #Make vector for pseudo y values for LG_X from left
  pseudo.y.LG_X <- NULL
  
  for (i in 1:ncol(phased.matrix)){
    pseudo.y <- sum(as.numeric(phased.matrix[,i]), na.rm = TRUE)/sample.size
    pseudo.y.LG_X <- cbind(pseudo.y.LG_X, pseudo.y)
  }
  
  #From right:
  phased.matrix.right <- rbind(XO.phased.genot.right[1:no.zero.right,], Matrix.convert.right)
  rownames(phased.matrix.right) <- rownames(XO.phased.genot.right)
  
  #Make vector for pseudo y values for LG_X from right
  pseudo.y.LG_X.right <- NULL
  
  for (i in 1:ncol(phased.matrix.right)){
    pseudo.y.right <- sum(as.numeric(phased.matrix.right[,i]), na.rm = TRUE)/sample.size
    pseudo.y.LG_X.right <- cbind(pseudo.y.LG_X.right, pseudo.y.right)
  }
    
  #Plot RFm values from left and right:
  if (plot == "Y"){
  plot(NULL, pch=19, col="red", xlab="cM", ylab="RFm", xlim=c(0,360), ylim=c(0,1.2), las=1, main=paste("LG",LG))
  lines(x=c(-10,260), y=c(0.5, 0.5), lty=2)  
  lines(x=c(-10,260), y=c(1, 1), lty=2, col="grey")  
  
  points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X, pch=19, col="red")
  points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X.right, pch=19, col="blue")
  }
  
  #Making output matrix with RFm values
  output <- rbind(as.numeric(LG_X[2,]), as.character(LG_X[4,]), as.numeric(LG_X[5,]), 
                  pseudo.y.LG_X, pseudo.y.LG_X.right)
  
  #Rounding values above 0.5 to 0.500 to facilitate comparison with y:
  output[4,][which(output[4,] > 0.5)] <- 0.5
  output[5,][which(output[5,] > 0.5)] <- 0.5
  
  colnames(output) <- colnames(XO.phased.genot)
  rownames(output) <- c(paste("LG", LG, "_cM", sep=""), paste("LG_arm", sep=""), paste("y (Ho)", sep=""),
                        paste("LG", LG, "_Left", sep=""), paste("LG", LG, "_Right", sep=""))
  return(output)
}



###########################################
### Prepping plot for all LG's with RFm ###
###########################################

#Define LG names:
LGs <- as.character(unique(Map[1,2:ncol(Map)]))

#Define which LGs needs to have axes drawn on plot
LGs_Y_axis <- LGs[c(1,7,13,19,25)]  
LGs_X_axis <- LGs[c(25:30)]

#Define uniform length of x-axes based on the longest LG on the map (Alternatively, separate x-axes can be made for each LG)
xmax <- 150

#######################
### PLOTTING SCRIPT ###
#######################
pdf(file="RFm_plots_speciesname_map.pdf")
par(mfrow=c(5,6),oma=c(4,3,3,2), mar=c(1.5,1,1,0))

for (j in LGs){
  plot(x=NULL, y=NULL, xlim=c(0,xmax), axes=F,ann=F,
       ylim=c(0,1))
  axis(1, las=1, at=c(0,50,100,150,200,250), labels=rep("",6), col="gray80")
  mtext(text=paste("LG", j), side=3, line=-1.5, cex=0.6)
  lines(x=c(-10,260), y=c(0.5, 0.5), lty=2) #Maximum RFm value expected with compelte interference
  box(col="gray80")

  if (j %in% LGs_Y_axis) {     #true when i is 1, 5, 9, 13 (%% is remainder when i divided by 4)
    axis(2, las=1, at=seq(0,1.4,0.2), col="gray80")
    mtext(text=expression(paste("RF"[m])) , side=2, line=3, cex=0.6, las=0)
  }
  if (j %in% LGs_X_axis){
    axis(1, las=1, col="gray80")  #reduce the number of axis labels
    mtext(text="cM", side=1, line=2, cex=0.6)
  }

  #Plotting RFm
  pseudo.y.matrix <- plot.pseudo.y(LG=j, plot="N")
  points(x=pseudo.y.matrix[1,], y=as.numeric(pseudo.y.matrix[4,]), pch=19, col="#0000FF33", cex=0.6) #Plotting RFm from left
  points(x=pseudo.y.matrix[1,], y=as.numeric(pseudo.y.matrix[5,]), pch=19, col="#FF000033", cex=0.6) #Plotting RFm from right

}  

dev.off()

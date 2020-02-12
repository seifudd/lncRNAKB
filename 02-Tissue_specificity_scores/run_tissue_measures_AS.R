#Used librarys
library("gplots") #function rich.colors
library("preprocessCore") #function normalize.quantiles

#Color used for graphes
my.col <- colorRampPalette(c("#FFFFFF", "black", "blue", "#FA8072","#00A2FF", "#00CC00", "#E0E0E0"))(7) #1:Backgroundcolor for all graphs, 2: Foregroundcolor for all graphs (E6E6E6), 3: Fill for histograms, 4: Red, for boxplots, 5: Blue, for boxplots, 6: Green, for boxplots, 7: Light gray



orgExpression <- read.table("./TPM_by_tissue_filtered_final_v6.txt", sep=",", header=TRUE, row.names=1)
#orgExpression_old <- read.table("./TPM_by_tissue_filtered.txt", sep=",", header=TRUE, row.names=1)
tissuesNames <- colnames(orgExpression)
nTissues <- 31

fPlotExpression <- function(x, dataName, fileName, names) #names=tissuesPrintNames #dataName=string to be printed on the x-axis
{
  dev.new(height=9, width=12)
  par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
  palette(rev(rich.colors(ncol(x))))
  
  plot(density(x[,1],n=1000), main = "Expression values among different tissues", xlab=dataName,col=(1), lwd=3)
  for(i in c(2:length(names)))
  {	
    lines(density(x[,i],n = 1000), col=(i), lwd=3)
  }
  legend("topright", names, col=(1:length(names)), lty="solid", lwd=3)
  
  dev.copy2pdf(device=quartz, file=paste(folder, organism, "Expression", expDataSource, "", fileName, add, ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
  #dev.off()
  
  return()
}	
###***###***###

###+++###
#Function requires data frame to be normalized
#1. All 0 are set to NA, to exclude them from quatile normalization
#2. Data are quantile normalized
#3. 0 values (the one set to NA) are set back to 0
fQN <- function(x) #
{
  x[x==0] <- NA
  x_m <- as.matrix(x)
  x <- normalize.quantiles(x_m)
  x[is.na(x)] <- 0
  return(data.frame(x))
}	
###***###***###

###+++###	
#Function require a vector with expression of one gene in different tissues.
#Mean is calculated taking in account tissues with 0 expression. 2+0+4=2
fmean <- function(x)
{
  if(!all(is.na(x)))
  {
    res <- mean(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}
###***###***###	

###+++###	
#Function require a vector with expression of one gene in different tissues.
#Max is calculated taking in account tissues with 0 expression. 2+0+4=2
fmax <- function(x)
{
  if(!all(is.na(x)))
  {
    res <- max(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}
###***###***###

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fGini <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x!=0))
      {
        res <- gini(x)*(length(x)/(length(x)-1))
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }	 		
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}
###***###***###

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fTsi <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x!=0))
      {
        res <- max(x) / sum(x)
      } else {
        res <- 0
      } 	
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }	
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Function requires setting of a treshold (rpkm)	
fCounts <- function(x, rpkm)
{
  if(all(!is.na(x)))
  {
    res <- length(which(x > rpkm))	
    if (res > 0)
    {
      res <- (1 - res/length(x))*(length(x)/(length(x)-1))  #Modification: To bring to normalized scale
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}
###***###***###		

###+++###
#Function require a data frame with expression data, and give back a vector with EEi values for each gene
#If expression for one tissue is not known, gene specificity for this gene is NA
fEe <- function(x)
{
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE))
    x <- rbind(x, c=colSums(x, na.rm=TRUE))	
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)])
    
    res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max)
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale
  } else {
    res <- NA
    print("No data avalable.")
  }
  return(res)
}
###***###***###

###+++###
#Function require a data frame with expression data, and give back a vector with PEM scores
#If expression for one tissue is not known, gene specificity for this gene is NA
fPem <- function(x)
{
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
    x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score
    
    x[x<1] <- 1
    x <- log10(x)
    
    x<- abs(x)				
    res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max) #choose only the maximal score for each gene
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
    #mylist <- list(res, x)
  } else {
    res <- NA
    print("No data avalable.")
  }
  #return(mylist)
  return(res)
}
###***###***###

###+++###
#Hg entropy
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fHg <- function(x)
{	
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x) !=0)
      {
        p <- x / sum(x)
        res <- -sum(p*log2(p), na.rm=TRUE)
        res <- 1 - (res/log2(length(p))) #Modification: To bring to normalized scale
      } else {
        res <- 0
      } 		
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}
###***###***###	

###+++###
#Z-score
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fZ <- function(x)
{	
  if(all(!is.na(x)))
  {
    res <-  apply(scale(t(x), center=TRUE, scale=TRUE),2,max)/((length(x[1,])-1)/sqrt(length(x[1,])))
    res_two <-  (scale(t(x), center=TRUE, scale=TRUE))/((length(x[1,])-1)/sqrt(length(x[1,])))
    res[is.na(res)] <- 0
    res_two[is.na(res_two)] <- 0
    #mylist <- list(res, res_two)
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
  #return(mylist)
}
###***###***###	

###+++###
#SPM score from TISGED
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fSpm <- function(x)
{	
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {	
      if(sum(x) !=0)
      {
        spm <- x^2/(x%*%x)
        res <- max(spm) #Modification:To bring to normalized scale. Choose max
      } else {
        res <- 0
      }	 		
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}






orgExpression$Tau <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, fTau)
#orgExpression$Gini <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, fGini)
#orgExpression$Tsi <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, fTsi)
#orgExpression$Counts <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, function(x){x <- fCounts(x, 1)})
#orgExpression$Hg <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, fHg)
#orgExpression$Zscore <- fZ(orgExpression[,tissuesNames[1:nTissues]])
#orgExpression$Spm <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, fSpm)
#orgExpression$Ee <- fEe(orgExpression[,tissuesNames[1:nTissues]])
orgExpression$Pem <- fPem(orgExpression[,tissuesNames[1:nTissues]])
orgExpression_old$Pem <- fPem(orgExpression_old[,tissuesNames[1:nTissues]])
test <- fPem(orgExpression[,tissuesNames[1:nTissues]])
#test_two <- fZ(orgExpression[,tissuesNames[1:nTissues]])

#orgExpression$Mean <- apply(orgExpression[,tissuesNames[1:nTissues]], 1, fmean)
#orgExpression$Max <- apply(orgExpression[, tissuesNames[1:nTissues]], 1, fmax)



dev.new(height=9, width=12)
par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
palette(rev(rich.colors(10)))
#palette(rev(blues9))

plot(density(orgExpression[,"Tau"],n=1000), main = "Tissue Specificity", xlab="Tissue specificity",col=(1), lwd=4, lty=1
     ,ylim=c(0,8), xlim=c(-0.1,1.1)
)
#lines(density(orgExpression[,"Gini"],n = 1000), col=(2), lwd=4, lty=2)
#lines(density(orgExpression[,"Tsi"],n = 1000), col=(3), lwd=4, lty=1)
#lines(density(orgExpression[,"Counts"],n = 1000), col=(1), lwd=4, lty=1)
#lines(density(orgExpression[,"Ee"],n = 1000), col=(5), lwd=4, lty=1)
#lines(density(orgExpression[,"Hg"],n = 1000), col=(6), lwd=4, lty=2)
#lines(density(orgExpression[,"Zscore"],n = 1000), col=(5), lwd=4, lty=1)
#lines(density(orgExpression[,"Spm"],n = 1000), col=(8), lwd=4, lty=2)
lines(density(orgExpression[,"Pem"],n = 1000), col=(5), lwd=4, lty=1)

legend("topright",c("Tau", "PEM"),col=(c(1,5)), lwd=4, lty=1, bty="n", seg.len=4)
dev.copy(png, filename="specificity_comparison_tau_PEM_v6.png", width =12, height=9, unit="in", res=500)
# dev.copy2pdf(device=quartz, file="./specificity_comparison_figure.pdf
dev.off()
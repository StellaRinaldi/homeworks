############################################################
# Homework 2 - RNAseq data analysis
# Course: Computational Genomics 2025/2026
# Date: 19/11/2025
#
# Description:
# This script implements three functions for RNAseq data analysis:
# 1. ...
# 2. ...
# 3. ...
#
# Group members:
# Lucas Cappelletti, 
# Soukaina Elboulrhaiti, 
# Maria Camilla Rigo,
# Stella Rinaldi,
# Filippo Tiberio
############################################################


# FIRST FUNCTION
MvAplot <- function(exprData, pdffilename, pcolor="black", lcolor="red") {
  reference <- exprData[,2]
  reference_name <- exprData$X[1]
  samples <- (length(exprData))
  
  pdf(file = pdffilename)
   for(i in 3:samples){
    current <- exprData[,i]
    
    M<-(log2(reference)-log2(current))[[1]]
    A<-((log2(reference)+log2(current))/2)[[1]]
    
    max_max <- max(M[!is.infinite(M)], na.rm=T)
    max_min <- max(abs(M[!is.infinite(M)]), na.rm=T)
    
    ylim <- max(c(max_max, max_min), na.rm=T)
    
    plot(x = A, y = M,
         pch = ".", col = pcolor,
         ylim=c(-ylim, ylim),
         main=paste("MvA plot of sample 1 (", reference_name, ") vs sample ",i, " (", colnames(exprData)[i], ")"),
         xlab="M",
         ylab="A",
    )
    abline(h = 0, col = lcolor)
  }
  
  dev.off()
}

# SECOND FUNCTION
TMMnorm <- function(exprData, annot, Mtrim=0.02, Atrim = c(0,8)) {
  #scaling by the sequencing depth and 
  SD<-colSums(exprData[,2:ncol(exprData)], na.rm = T) #ricontrollare
  data<-sweep(exprData[,2:ncol(exprData)],2,SD, FUN='/')*10^6 
  rownames(data)<-exprData[,1]
  normData<-data
  #SCALING FACTORS
  ni<-ncol(data)
  # initialization of SF vector
  SF<-rep(0,times=ni)
  
  for (i in 2:ni){
    # computing A and M with an offset
    offset=0.0001
    M<-log2(data[,1]+offset)-log2(data[,i]+offset)
    A<-(log2(data[,1]+offset)+log2(data[,i]+offset))/2
    
    indA<-A>Atrim[1]&A<Atrim[2]
    SF[i]<-mean(M[indA], trim=Mtrim)
    SF[i]<-2^(SF[i])
    normData[,i]<-normData[,i]*SF[i]
  }
  SF[1]<-2^(SF[1])
  # Scaling for the length 
  normData <- sweep((normData),1,annot$Length, FUN="/")*(10^3)
  # let's add again the column with the names
  normData<-cbind(exprData[,1],normData)
  # returning the list
  l<-list(normData,SF)
  return(l)
}


DATA <- read.table("raw_count.txt", sep="\t", row.names=1, header=TRUE)
annotations <- read.table("gene_annot.txt", sep="\t", row.names=2, header=TRUE, quote = "\"")

MvAplot(DATA, "MvA_Plot.pdf")
normalized <- TMMnorm(DATA, annotations)
MvAplot(normalized[[2]], "TEST_normalized.pdf")


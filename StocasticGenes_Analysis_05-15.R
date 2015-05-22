##Code for finding genes that are highly noisy but have low expression change between two sets. 
#Uses results from the MultiRankSeq analysis done in the MRSR_Analysis_05-15.R program
##See http://dx.doi.org/10.1155/2014/248090 for paper publication
##Christian Meyer 
##Contact christian.t.meyer@vanderbilt.edu with questions
##Quaranta Lab May 2015
##Analysis of the GSE58135 data set from Varley et al PMID: 24929677
#Dataset consists of 
# 28 Cell Line samples (abbreviated CL in code)
# 42 TNBC samples
# 42 ER/HER2+ samples (abbr. ER in code)
# 30 Unmatched Adjacent ER/HER2+ samples (abbr. adjER in code)
# 5  No cancer mammary cell samples (abbr. NC in code)
# 21 Unmatched Adjacent TNBC samples (abbr. adjTNBC in code)

#NOTE: changing directories is done in the Lib_RNAseq.R file!
#Run Lib_RNAseq.R to compile necessary libraries
#NOTE: Code written for running in Linux.  

##################################################################################
##################################################################################
#Functions for code:

##Functions for getting Filenames and reading in Data
#Read in data fileName is the name of each file type
fileName <- function(file_names) {
  toKeep <- gregexpr("_", file_names)  #Just keep the GSE... part of file name as identifier
  name <- substr(file_names,1,toKeep[[1]][2]-1)
  return(name)
}

#Read in files in directory between init and fin
getData <- function(init,fin,file_names,sample_names) {
  list_dataset = list(FPKM=NULL, CI=NULL)
  for(i in init:fin) {
    if (i == init) {
      temp_dataset <- read.delim(file_names[i])
      #Change depending on what kind of name you want for the rows
      #rnames <- paste0(temp_dataset[[4]],'_',temp_dataset[[5]])
      rnames <- temp_dataset[[4]]
      list_dataset[[1]] <- cbind(temp_dataset[[10]])
      list_dataset[[2]] <- cbind(temp_dataset[[12]] - temp_dataset[[11]])
      rownames(list_dataset[[1]]) <- rnames  
      rownames(list_dataset[[2]]) <- rnames
    }
    if (i != init){
      temp_dataset <- read.delim(file_names[i])
      list_dataset[[1]] <- cbind(list_dataset[[1]],temp_dataset[[10]])
      list_dataset[[2]] <- cbind(list_dataset[[2]],temp_dataset[[12]] - temp_dataset[[11]])      
    }
  }
  colnames(list_dataset[[1]]) <- sample_names[init:fin]
  colnames(list_dataset[[2]]) <- sample_names[init:fin]
  return(list_dataset)
}

#Get the gene short names from the dataset.
geneShort <- function(file_names){
  temp_dataset <- read.delim(file_names[3])
  gene_shortName = temp_dataset[[5]]
  return(gene_shortName)
}


#Function to generate vector describing stocastic noise of each gene in a sample type
#Randomly select samples out of the columns of dataset to calculate noise (median absolute deviation)
#Give the name of the type of cell (dataType) to the columns
#Return confidence interval of noise
#Pass in the number of times this function has been called for this sample
stocFun = function(list_dataset,dataType,k){
  indx = NULL
  colnm = NULL
  #Compute the rowSds for the dataset
  ret[[1]] = rSD = rowMads(list_dataset[[1]],na.rm=TRUE)
  ret[[2]] = rCI = rowMads(list_dataset[[2]], na.rm=T)  
  #Keep the column name
  colnm = paste0(dataType, "_", k)
  rownames(ret[[1]]) <- row.names(list_dataset[[1]])
  rownames(ret[[2]]) <- row.names(list_dataset[[1]])
  colnames(ret[[1]]) <- colnm
  colnames(ret[[2]]) <- colnm
  return(ret)
}


#Function to cluster the data by cell type and display dendrogram
clusterPlot <- function(dataset,groupCodes,colorCodes,title,subtitle,mag) {
  #Cluster all the datasets to estimate how well they fit into each cluster by gene expression
  len = dist(t(dataset),method = "euclidean") #add one to avoid zero problem (infinite distance from zeros)
  hc = hclust(len, method = "complete")
  dend = as.dendrogram(hc)
  #Change labels colors using dendextend library
  labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
  # reduced label size to fit with input magnification (mag)
  # save to jpg file or output to screen
  X11()
  par(cex=mag, mar=c(7, 1, 5, 10))
  plot(dend, xlab="", ylab="", main="", sub="", axes=FALSE, horiz = TRUE)
  par(cex=1)
  title(main=title,sub=subtitle)
  x = par("usr")
  #Location of the subtitle
  #text(x[1],x[4],subtitle)
  #text(locator(1),subtitle)
  axis(1)
}



##################################################################################
##################################################################################
#Only run if havn't read in data in another program
if(!exists("file_names")&&!exists("sFileNM")&&!exists("Tot_Data")){
  #First get files and data from current directory and put it all into a lists of lists Tot_Data
  #Get current directory file names
  file_names <- dir();
  #Get shortened names of each file to use as headers
  sFileNM <- lapply(file_names,fileName)
  
  #Declare List to hold the RNAseq data information for each sample type
  #Each component of the list is also a list of three values
  #List[[1]] is the FPKM values  $FPKM
  #List[[2]] is the difference between the high and low FPKM confidence limits $CI
  Tot_Data = list(CL=NULL,TNBC=NULL,ER=NULL,adjER=NULL,NC=NULL,adjTNBC=NULL)
  NumSamp = c(28,42,42,30,5,21)
  colorCodes <- c(CL = "red", TNBC = "green", ER = "blue", adjER = "black", NC = "pink", adjTNBC = "purple")
  
  #Load in the data for each cell line
  for (i in 1:length(NumSamp)){
    if(i==1){
      Tot_Data[[i]] = getData(1,NumSamp[i],file_names,sFileNM)
    }else{
      Tot_Data[[i]]=getData((sum(NumSamp[1:i-1])+1),sum(NumSamp[1:i]),file_names,sFileNM)
    }
  }
  #Grab the short names of the genes
  gene_shortNames = geneShort(file_names)
}

#Return an error message if MultiRankSeqReport is not in folder
if(file.exists(....)){
  stop('No MRSR is in the currently specified folder')
}

##################################################################################
##################################################################################
#Find the median absolute distance (NOISE) for each gene within the sample set and cluster 
#Sample types by their respective "noise" levels

#Select which cell lines to analyze and subsequently cluster
#to_Clust = c("CL","TNBC","ER","adjER","NC","adjTNBC")
to_Clust = c("CL","TNBC","adjTNBC","adjER")

#Color codes for the hierarchical clustering
colorCodes <- c(CL = "red", TNBC = "green", ER = "blue", adjER = "black", NC = "pink", adjTNBC = "purple")
for (i in 1:length(to_Clust)){
  if(i == 1){
    temp_data = stocFun(Tot_Data[[to_Clust[i]]],to_Clust[i],1)
    GeneNoise = cbind(temp_data[[1]])
    GNerror = cbind(temp_data[[2]])
    groupCodes = c(rep(to_Clust[i],1))
  }else{
    temp_data = stocFun(Tot_Data[[to_Clust[i]]],to_Clust[i],1)
    GeneNoise = cbind(GeneNoise,temp_data[[1]])
    GNerror = cbind(GeneNoise,temp_data[[2]])
    groupCodes = c(groupCodes,rep(to_Clust[i],1))
  }
}

#Select maginfication of the sample name in cluster plot and title
mag = 1
title = "Clustering the expression noise across samples, \n no weights"
sub_title = "All samples used"


clusterPlot(GeneNoise,groupCodes,colorCodes,title,sub_title,mag)

#Box plot of the distribution of noise between samples...
X11()
boxplot(GeneNoise, main = title)

rm(GeneNoise,temp_data,GNerror)
##################################################################################
##################################################################################
#Now to find the genes with significant noise, but low fold change as found in the MRSR program

#Read in results from MRSR folder...
#find log fold change in gene noise between two sample types
#Pick the same two cell lines that the MRSR report has been run for...
to_Clust = c("TNBC","adjTNBC")

#Find the log fold change in each gene between types
for (i in 1:length(to_Clust)){
  if(i == 1){
    temp_data = stocFun(Tot_Data[[to_Clust[i]]],to_Clust[i],1)
    GeneNoise = cbind(temp_data[[1]])
    GNerror = cbind(temp_data[[2]])
    groupCodes = c(rep(to_Clust[i],1))
  }else{
    temp_data = stocFun(Tot_Data[[to_Clust[i]]],to_Clust[i],1)
    GeneNoise = cbind(GeneNoise,temp_data[[1]])
    GNerror = cbind(GeneNoise,temp_data[[2]])
    groupCodes = c(groupCodes,rep(to_Clust[i],1))
  }
}

logFoldChange_Noise = foldChanges(GeneNoise[,1],GeneNoise[,2],to_Clust[1],to_Clust[2],log=T)
#calculate the uncertainty in the fold change
Noise_Uncertainty = log(GNerror[,1]/GeneNoise[,1]+GNerror[,2]/GeneNoise[,2], base = 2)


#Get the fold change in expression from the MultiRankSeq report
MRSR <- read.csv("~/Documents/Lab/WGCNA_analysis/Report_5-16/MultiRankSeqReport_TNBCvadjTNBC/MRSR_TNBCvadjTNBC_Large.html.diff.csv")

logFoldChange_Expres = 


























##################################################################################
##################################################################################
#Find how robust the noise is to downsampling
persample = .5 #Percent of samples to be used
numboot = 20 #Number of times to bootstrap

#Select which cell lines to cluster
#to_Clust = c("CL","TNBC","ER","adjER","NC","adjTNBC")
to_Clust = c("CL","TNBC","adjTNBC", "adjER")
for (i in 1:length(to_Clust)){
  for (j in 1:numboot){
    num = round(dim(Tot_Data[[to_Clust]][[1]])[2]*persample) #Number of samples to downsize to
    if(num<=2){
      num = 3
    }
    #Random samples to take
    rnd = sample(1:dim(Tot_Data[[to_Clust]][[1]])[2], num)
    if(j == 1){
      temp_data = stocFun(Tot_Data[[to_Clust[i]]],to_Clust[i],1)
      GeneNoise = cbind(temp_data[[1]])
      GNerror = cbind(temp_data[[2]])
      groupCodes = c(rep(to_Clust[i],1))
    }else{
      temp_data = stocFun(Tot_Data[[to_Clust[i]]],to_Clust[i],1)
      GeneNoise = cbind(GeneNoise,temp_data[[1]])
      GNerror = cbind(GeneNoise,temp_data[[2]])
      groupCodes = c(groupCodes,rep(to_Clust[i],1))
    }
  }
}

#Select maginfication of the sample name in cluster plot and title
mag = .4
str = paste0(numsample*100, "% downsampled for each set")
title = "Bootstrapped measure of noise"
sub_title = str

#Plot the noise subgroups to show robust clustering distance to increased distance
clusterPlot(GeneNoise,groupCodes,colorCodes,title,sub_title,mag)



X11()
#to look at the noise of each boot strapped sample
temp_data = matrix(0,dim(GeneNoise)[1],length(to_Clust))
for (j in 1:length(to_Clust)){
  temp_data[,j] = rowMeans(GeneNoise[,((j-1)*numDataSets+1):(j*numDataSets)],na.rm=T)
}
boxplot(temp_data,main = title,sub=sub_title)
#For space rm GeneNoise and temp_data
rm(GeneNoise)
}





#################################################################################################
#extra code
# for (i in 1:length(to_Clust)){
#   if(i == 1){
#     GeneNoise2 = cbind(Tot_Data[[to_Clust[i]]][[1]])
#   }else{
#     GeneNoise2 = cbind(GeneNoise2,Tot_Data[[to_Clust[i]]][[1]])
#   }
# }
# 
# #Threshold the row values
# GeneNoise = thresh(GeneNoise,GeneNoise2,1)
# rm(GeneNoise2)


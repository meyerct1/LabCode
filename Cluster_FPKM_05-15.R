##Code for clustering samples by FPKM value
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

#Analysis of the clustering of samples in dataset by FPKM values
#Run Lib_RNAseq.R to compile necessary libraries
#NOTE: changing directories is done in the Lib_RNAseq.R file!
#NOTE: Code written for running in Linux.  If using a windows or mac change the X11() function
#in the clustPlot function to avoid errors

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


#Function to cluster the data by cell type and display dendrogram
#This uses the dendexend library to color code the leafs on the plot
#Plot is in a new window by calling X11() which opens a new graphical device in Linux
clustPlot <- function(dataset,groupCodes,colorCodes,title,subtitle,mag,graphname,tosave) {
  #Cluster all the datasets to estimate how well they fit into each cluster by gene expression
  len = dist(t(dataset),method = "euclidean") #add one to avoid zero problem (infinite distance from zeros)
  hc = hclust(len, method = "complete")
  dend = as.dendrogram(hc)
  #Change labels colors using dendextend library
  labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
  # reduced label size to fit with input magnification (mag)
  # save to jpg file or output to screen
  
  if (tosave == 1) {
    jpeg(graphname)
    par(cex=mag, mar=c(7, 1, 5, 10))
    plot(dend, xlab="", ylab="", main="", sub="", axes=FALSE, horiz = TRUE)
    par(cex=1)
    title(xlab="", ylab="", main=title)
    text("topleft",subtitle)
    axis(1)
    dev.off()
  } 
  else{
    X11()
    par(cex=mag, mar=c(7, 1, 5, 10))
    plot(dend, xlab="", ylab="", main="", sub="", axes=FALSE, horiz = TRUE)
    par(cex=1)
    title(main=title,sub=subtitle)
    x = par("usr")
    #text(x[1],x[4],subtitle)
    #text(locator(1),subtitle)
    axis(1)
  }
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



##################################################################################
##################################################################################
#Cluster with the FPKM values

tosave = 0  #Option to save the plots 1 = save 0 = output plots in new window
#Name to save the plot by if tosave=1
graphname = "All.jpg"

#Pick which cell lines to cluster.  
#to_Clust = c("CL","TNBC","ER","adjER","NC","adjTNBC")
to_Clust = c("TNBC","adjTNBC","NC","CL")
#Color codes for the hierarchical clustering
colorCodes <- c(CL = "red", TNBC = "green", ER = "blue", adjER = "black", NC = "pink", adjTNBC = "purple")

for (i in 1:length(to_Clust)){
  if(i == 1){
    foo = cbind(Tot_Data[[to_Clust[i]]][[1]])
    groupCodes = c(rep(to_Clust[i],dim(Tot_Data[[to_Clust[i]]][[1]])[2]))
  }else{
    foo = cbind(foo,Tot_Data[[to_Clust[i]]][[1]])
    groupCodes = c(groupCodes,rep(to_Clust[i],dim(Tot_Data[[to_Clust[i]]][[1]])[2]))
  }
}

#Select maginfication of the sample name in cluster plot and title
mag = 1-.1*length(to_Clust)
title = "Clustering sample type by FPKM values"
#sub_title = "CL=red \n TNBC=green \n ER=blue \n adjER=black \n NC=pink \n adjTNBC=purple"
sub_title = "TNBC=green \n adjTNBC=purple \n CL=red \n NC=pink"
#NOTE the output plot takes a right click to place sub_title
clustPlot(foo,groupCodes,colorCodes,title,sub_title,mag,graphname,tosave)

#Remove foo to save space
rm(foo)
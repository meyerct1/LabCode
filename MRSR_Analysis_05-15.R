##Code for running Yan Guo's Multi Rank Seq Package
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

#Run Lib_RNAseq.R to compile necessary libraries
#NOTE: Code written for running in Linux.  
#NOTE: changing directories is done in the Lib_RNAseq.R file!

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
#Now run a MultiSeqRank report generator code developed by Yan Guo et al. http://dx.doi.org/10.1155/2014/248090
#to find significantly up or down regulated genes
#This produces a csv file, html report and a folder with images
#For this analysis differential expression was considered siginificant if the rank for all three
#methods (baySeq,edgeR,DESeq2)
#NOTE: THE OUTPUT DIRECTORY MUST BE THE SAME AS THE CURRENT WORKING DIRECTORY!


#Do a random subset of the TNBC v adjTNBC as the MultiRankSeqReport can't do all of them
#Time for sample of 20 is ~___hours

setwd("~/Desktop/MRSR/")
tic
#Number of samples
numsample = 20
rndT = sample(1:dim(Tot_Data$TNBC[[1]])[2],numsample)
rndA = sample(1:dim(Tot_Data$adjTNBC[[1]])[2],numsample)
TNBC_sub = cbind(Tot_Data$TNBC[[1]][,rndT],Tot_Data$adjTNBC[[1]][,rndA])
group_sub = c(rep(0,numsample),rep(1,numsample))
#Requires integer values so multiply by 10000 and round to retain all significant digits of FPKM values
resTNBC_adjTNBC <- MultiRankSeqReport(output = "~/Desktop/MRSR/MRSR_TNBCvadjTNBC_Large.html",rawCounts = round(10000*TNBC_sub),group = group_sub)


numsample = 20
rndT = sample(1:dim(Tot_Data$TNBC[[1]])[2],numsample)
rndA = sample(1:dim(Tot_Data$adjER[[1]])[2],numsample)
TNBC_sub = cbind(Tot_Data$TNBC[[1]][,rndT],Tot_Data$adjER[[1]][,rndA])
group_sub = c(rep(0,numsample),rep(1,numsample))
#Requires integer values so multiply by 10000 and round to retain all significant digits of FPKM values
resTNBC_adjER <- MultiRankSeqReport(output = "~/Desktop/MRSR/MRSR_TNBCvadjER_Large.html",rawCounts = round(10000*TNBC_sub),group = group_sub)

numsample = 20
rndT = sample(1:dim(Tot_Data$CL[[1]])[2],numsample)
rndA = sample(1:dim(Tot_Data$TNBC[[1]])[2],numsample)
TNBC_sub = cbind(Tot_Data$CL[[1]][,rndT],Tot_Data$TNBC[[1]][,rndA])
group_sub = c(rep(0,numsample),rep(1,numsample))
#Requires integer values so multiply by 10000 and round to retain all significant digits of FPKM values
resTNBC_adjER <- MultiRankSeqReport(output = "~/Desktop/MRSR/MRSR_CLvTNBC_Large.html",rawCounts = round(10000*TNBC_sub),group = group_sub)
toc


#Write out list of significantly up or down regulated genes...
#Read in MRSR from folder
res <- read.csv("~/Documents/Lab/WGCNA_analysis/Report_5-16/MultiRankSeqReport_TNBCvadjTNBC/MRSR_TNBCvadjTNBC_Large.html.diff.csv")

#get index of all genes which had a rank less than 200 for all methods
idx = resTNBC[[12]] <= 200 & resTNBC[[11]] <= 200 & resTNBC[[13]] <= 200
#idx = resTNBC[[14]] <= 200
sum(idx==TRUE) #Display the number which should coorelate with the 3rd venn diagram in results
rnames = rownames(Tot_Data[[1]][[1]])
true_loc = which(idx==TRUE)
loc = rep(0,length(true_loc))
for (i in 1:length(true_loc)){
  loc[i]= which(resTNBC[[1]][true_loc[i]]==rnames)
}

#Write out names of significant genes to file MRSR_sigGenes for analysis with WebGestalt...
write(substr(rnames[loc],1,15),"~/Documents/Lab/WGCNA_analysis/Report_5-16/MRSR_sigGenes.txt",sep = "\n")
MRSRgenes = rnames[loc]

#Remove TNBC_sub to save space
rm(TNBC_sub)

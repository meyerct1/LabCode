#Library for RNAseq analysis
#05-15
#Christian Meyer
#Contact christian.t.meyer@vanderbilt.edu
#Functions and libraries to load initially

#location of Data:
dataLoc = "~/Documents/Lab/WGCNA_analysis/GSE58135_RAW/FPKM"
#location of MRSR report
MRSRreport = ""
#location of 
##Initial library loading and workspace configuration
# Please adapt the following path. Note that we use / instead of \ in the path.
setwd("~/Documents/Lab/WGCNA_analysis/GSE58135_RAW/FPKM")


# read in the R libraries
library(MASS) # standard, no need to install
library(class) # standard, no need to install
library(cluster)
library(impute)# install it for imputing missing value
library(Hmisc) # install it for the C-index calculations
library(survival)
library(dendextend)
library(WGCNA)
library(MultiRankSeq)
library(R.utils)
library(matrixStats)
library(pracma)
library(extremevalues)
allowWGCNAThreads()
options(stringsAsFactors=F)



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

#Functions to measure time in R similary to matlab execution
tic <- 1
class(tic) <- "tic"

toc <- 1
class(toc) <- "toc"

print.tic <- function(x,...) {
  if (!exists("proc.time"))
    stop("cannot measure time")
  gc(FALSE)
  assign(".temp.tictime", proc.time(), envir = .GlobalEnv)
}

print.toc <- function(x,...) {
  if (!exists(".temp.tictime", envir = .GlobalEnv))
    stop("Did you tic?")
  time <- get(".temp.tictime", envir = .GlobalEnv)
  rm(".temp.tictime", envir = .GlobalEnv)
  print(res <- structure(proc.time() - time,
                         class = "proc_time"), ...)
  invisible(res)
}
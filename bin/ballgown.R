#!/usr/bin/Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.rstudio.com/")

BiocManager::install("ballgown")

# args <-commandArgs(TRUE)
# id_tuple <-dir(args[1])
# folderName<- paste(id_tuple,".tablemaker ")

data_directory = system.file('extdata', package='ballgown') # automatically finds ballgown's installation directory
# examine data_directory:
#data_directory

bg = ballgown(dataDir='.', samplePattern='brain', meas='all')
bg
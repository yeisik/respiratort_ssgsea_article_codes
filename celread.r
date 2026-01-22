# install the core bioconductor packages, if not already installed
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.16")

# install additional bioconductor libraries, if not already installed
#BiocManager::install("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
#BiocManager::install("affy")  # Methods for Affymetrix Oligonucleotide Arrays
#BiocManager::install("hgu133a2.db", type = "source")  # GSE1297: Platform_title = [HG-U133A]
#BiocManager::install("gpl14604hs133av2entrezcdf")
#install.packages("gpl14604hs133av2entrezcdf", type = "source", repos = NULL)
#install.packages("hgu133a2hsentrezgcdf", type = "source", repos = NULL)


#setwd("C:/Users/heroc/Desktop/CEL_FILE_READER/data")

#
# change path.
#
setwd("CEL_File_Reader//GSE73072//GSE73072_RAW")
cels = list.files("CEL_File_Reader//GSE73072//GSE73072_RAW")
print(cels[1:10])

library(affy)
library(hgu133a2.db)


#raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133a2cdf")

#library(gpl14604hs133av2entrezcdf)
#raw.data = ReadAffy(verbose = FALSE, filenames = cels,cdfname="gpl14604hs133av2entrezcdf")

#library(hgu133a2hsentrezgcdf)
raw.data = ReadAffy(verbose = FALSE, filenames = cels,cdfname="hgu133a2hsentrezgcdf")

# change path.
setwd("CEL_File_Reader//")

data.rma.norm = rma(raw.data)

rma = exprs(data.rma.norm)

# USE UPDATED CDF:
#write.table(rma, file = "GSE_ALL_hgu133a2hsentrezgcdf.txt", quote = FALSE, sep = "\t")
# Library 
library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)

## prerequisite     
# set your working directory (ex: setwd("C:/MyFolder/PRM/") )
# where your PRM_ExportSkyline.csv has been stored
setwd("C:/MyFolder/PRM/")

## Import files
dat <- read.csv("PRM_ExportSkyline.csv", 
                header = T, stringsAsFactors = F, dec = ",", sep = ";") %>%
  mutate_all(funs(str_replace(., ",", "."))) %>%  
  mutate_all(funs(str_replace(., "\\#N/A", "NA"))) %>%
  mutate_if(str_detect(colnames(.), ".Total.Area|.Average.Mass.Error.PPM|.Library.Dot.Product|.Detection.Q.Value"),  as.numeric) %>%
  mutate(Feature = paste(Peptide.Sequence, Precursor.Charge, Precursor.Mz, sep = "_"))

### Peptides info
dat_info <- dat %>%
  dplyr::filter(Protein.Name != "Biognosys standards") %>%
  dplyr::select(Feature)

### Dot Product file
dat_DotP <- dat %>%
  dplyr::filter(Protein.Name != "Biognosys standards") %>%
  dplyr::select(contains(".Library.Dot.Product")) %>%
  rename_all(. %>%  gsub(".Library.Dot.Product", "", .) ) 

### Mass error ppm file
dat_ppm <- dat %>%
  dplyr::filter(Protein.Name != "Biognosys standards") %>%
  dplyr::select(contains(".Average.Mass.Error.PPM")) %>%
  rename_all(. %>%  gsub(".Average.Mass.Error.PPM", "", .) ) 
dat_ppm[is.na(dat_ppm)] <- 100

## Discretize data
dat_ML_disc <- dat_ppm
for(x in 1:nrow(dat_ppm)) {
  for (y in 1:ncol(dat_ppm)) {
    dat_ML_disc[x,y] <- ifelse( ((dat_DotP[x,y] > 0.85 & dat_ppm[x,y] < 10) | (dat_DotP[x,y] > 0.75 & dat_ppm[x,y] < 3 )), "TRUE", "FALSE")
  }
}

## Export data
Samples <- colnames(dat_ML_disc)
x <- dat_ML_disc %>%
  transpose %>%
  as.data.frame 
colnames(x) <- dat_info$Feature
rownames(x) <- Samples
x=as.data.frame(cbind(rownames(x),x))
colnames(x)[1]="Sample"
write.table(x, "PRM_InputML_discrete.csv", sep = "\t", quote=FALSE, row.names = FALSE)


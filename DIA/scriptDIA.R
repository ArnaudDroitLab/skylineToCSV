# Libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)

## prerequisite     
# Store all your Skyline exports in a directory called "ExportSkyline"
# set your working directory where the ExportSkyline have been created (ex: setwd("C:/MyFolder/DIA") )
setwd("C:/MyFolder/DIA/")

## create the output folder of processed skyline exports
dir.create(file.path(".", "InputML"))


## import files from ExportSkyline folder----
all_files <- list.files(path = "./ExportSkyline/", pattern = ".csv", full.names = T)

dat = list()
dat_info = list()
dat_DotP = list()
dat_Qval = list()
dat_ppm = list()


## parse content
for (i in 1:length(all_files)) {
  ### Input file  
  dat[[i]] <- read.csv(all_files[[i]], stringsAsFactors = F, header = T, dec = ".", sep = ";")  %>%
    mutate_all(funs(str_replace(., ",", "."))) %>%  
    mutate_all(funs(str_replace(., "\\#N/A", "NA"))) %>%
    mutate_if(str_detect(colnames(.), ".Total.Area|.Average.Mass.Error.PPM|.Library.Dot.Product|.Detection.Q.Value"),  as.numeric) %>%
    mutate_if(is.numeric, abs) %>%
    dplyr::filter(Protein.Name != "Decoys", 
                  Missed.Cleavages == 0, 
                  !str_detect(Peptide.Modified.Sequence.Unimod.Ids, "unimod"), 
                  !str_detect(Peptide.Sequence, "M"), 
                  !str_detect(Peptide.Sequence, "C"), 
                  nchar(Peptide.Sequence) >= 8) %>%
    mutate(Feature = paste(Peptide.Sequence, Precursor.Charge, Precursor.Mz, sep = "_"))
  
  ### Peptides info
  dat_info[[i]] <- dat[[i]] %>%
    dplyr::filter(Protein.Name != "Biognosys standards") %>%
    dplyr::select(Feature)
  
  ### Dot Product file
  dat_DotP[[i]] <- dat[[i]] %>%
    dplyr::filter(Protein.Name != "Biognosys standards") %>%
    dplyr::select(contains(".Library.Dot.Product")) %>%
    rename_all(. %>%  gsub(".Library.Dot.Product", "", .) ) 
  
  ### Detection QValue file
  dat_Qval[[i]] <- dat[[i]] %>%
    dplyr::filter(Protein.Name != "Biognosys standards") %>%
    dplyr::select(contains(".Detection.Q.Value")) %>%
    rename_all(. %>%  gsub(".Detection.Q.Value", "", .) )
  
  
  ### Mass error ppm file
  dat_ppm[[i]] <- dat[[i]] %>%
    dplyr::filter(Protein.Name != "Biognosys standards") %>%
    dplyr::select(contains(".Average.Mass.Error.PPM")) %>%
    rename_all(. %>%  gsub(".Average.Mass.Error.PPM", "", .) ) 
  dat_ppm[[i]][is.na(dat_ppm[[i]])] <- 100
}

## For ML ---- 
dat_ML_disc = list()
for (i in 1:length(all_files)) {
  dat_ML_disc[[i]] <- dat_ppm[[i]]
  for(x in 1:nrow(dat_ML_disc[[i]])) {
    for (y in 1:ncol(dat_ML_disc[[i]])) {
      dat_ML_disc[[i]][x,y] <- ifelse((dat_DotP[[i]][x,y] > 0.75 & dat_Qval[[i]][x,y] < 0.01), "TRUE", "FALSE")
    }
  }
  dat_ML_disc[[i]] <- cbind(dat_info[[i]], dat_ML_disc[[i]]) 
}

df_disc <- dat_ML_disc[[1]]
for (i in 2:length(all_files)) {
  df_disc <- merge(df_disc, dat_ML_disc[[i]], by = "Feature", all = TRUE)    
}

## Export
Feature <- df_disc$Feature
Samples <- colnames(df_disc)[2:ncol(df_disc)]
x <- df_disc %>% 
  dplyr::select(-Feature) %>%
  transpose %>%
  mutate_each(funs(replace(., is.na(.), "FALSE"))) %>%
  mutate(Sample = Samples, 
         class = gsub("_.*", "", Sample)) %>%
  dplyr::select(Sample, class, 3:ncol(.)-2)
colnames(x)[3:ncol(x)] <- Feature
write.table(x, "InputML/DIA_InputML_discrete.csv", sep = "\t", quote=FALSE, row.names = FALSE)

# skylineToCSV
Convert skyline outputs to CSV files compatible with BioDiscML (and also most of machine learning libraries)

There are 2 folders: 
- DIA 
- PRM

DIA contains a R script, scriptDIA.R, that retreive the exported files from Skyline in the folder ExportSkyline to create a machine learning CSV compatible file DIA_InputML_discrete.csv.

PRM contains a R script, scriptPRM.R, that retreive the exported file PRM_ExportSkyline.csv from Skyline to create a machine learning CSV compatible file PRM_InputML_discrete.csv.

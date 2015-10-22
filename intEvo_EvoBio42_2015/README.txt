This data package includes raw data, as well as R code and other information required to reproduce the analyses in this paper:

Haber, A. 2015. The Evolution of Morphological Integration in the Ruminant Skull. Evolutionary Biology 42:99-114.

Raw data and specimen information are provided as both tab-limited text file and R object. Phylogenetic trees are provided as both nexus file and R object (of class phyllo; see package ape). Other information is provided as R objects.

A_NkmSymm.Rdata / A_NkmSymm.tab
The data matrix of landmark configurations for all specimens included in the study, in the format of N specimens by k*m variables, where k is number of landmarks (43) and m is number of dimensions (3). Columns are ordered as {x1,y1,z1,x2,y2,z2,etc.}. See table S1 for landmark names and definitions. Configurations are bilateral and symmetric.

specimenInfo.Rdata / specimenInfo.tab
A matrix of information for each specimen in the data matrix. Row names are as in the data matrix, and compose of the museum abbreviation and the specimen number as appears in the museum record. Columns are: mf, male or female; tx, taxon code name as appears in the phylogenetic trees; species.FV, the species it belongs to in the FV phylogeny; species.CT, the species it belongs to in the CT phylogeny; species.MR, the species it belongs to in the MR phylogeny.

LtI.Rdata / LtI.nexus
The three phylogenetic hypotheses used in this study including branch lengths. The R object is a phylo object suitable for ape package.


ildefL.Rdata
Definition tables for calculating interlandmark distances. This is a list of two elements: "IL32", first two columns provide the two landmarks defining each interlandmark distance in the IL32 dataset, and the third column provide the postulated module to which that distance belong (see table S2); "ILtes" a mart of two columns only, providing the two landmarks defining each of the 107 interlandmark distances for the tessellation-based dataset.

VerL.Rdata
A list with two elements: "IL32", is a list of 51 measurement (error) covariance matrices, one for each taxon, based on the IL32 dataset;  "ILtes", is a list of 51 measurement (error) covariance matrices, one for each taxon, based on the ILtes dataset. See supplemental information I for details on how they were generated.

prtcl-EvoOfInt_EvoBio.R
An annotated R protocol for running all the analyses in this study. Each part in the code is run as a stand-along, in the sense that it assumes an empty space, calling for whatever is needed and ends with its own output and workspace image file. However, it requires objects that were generated in previous parts and assumes that they are present in the same folder.

functions-EvoOfInt_EvoBio.R
An annotated R source file with all the functions required to run the R protocol provided here. These funciton were written specifically for this study by the author. This file is sourced when running the protocol. Please cite this data package if using them for your own needs.

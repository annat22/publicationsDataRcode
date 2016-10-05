Raw data, R code, and other information required to reproduce the analyses in:

Haber, A., and I. Dworkin. 2016. Dis-integrating the fly: A mutational perspective on phenotypic integration and covariation. Evolution in press
BioRxiv doi: http://dx.doi.org/10.1101/023333.

inf_gt.Rdata - R object of the same name with information about the genotypes. 25-by-6 data frame. Row names are genotype ID. Column names are: WT, wild type; allele, indicating the mutant allele; gene; gene abbr; pathway; N, sample size

inf_spec.Rdata - R object of the same name with information about the specimens. 1880-by-4 data frame. rownames are specimen ID. Column names are: WT, wild type; genotype, genotype ID as in inf_gt; repID, replicate ID; centsize, centroid size

XallData.Rdata/XallData.tab - R object named Xall and a tab-delimited text file with landmark coordinates after Procrustes superimposition and projection to tangent space. 1880-by-96 numeric matrix. Columns ordered as “x1, y1, x2, y2…”. Rows ordered by specimen ID.


JxcTPS.tab - a tab delimited table with LORY data, using TPS as interpolation function. The original Procrustes data was not regressed on size

JxcEBS.tab - a tab delimited table with LORY data, using EBSS as interpolation function. The original Procrustes data was not regressed on size

JxcTPSsr.tab - a tab delimited table with LORY data, using TPS as interpolation function. The original Procrustes data was regressed on size before running in LORY

JxcEBSsr.tab - a tab delimited table with LORY data, using EBSS as interpolation function. The original Procrustes data was regressed on size before running in LORY

JxcCo.tab - a tab delimited table with coordinates for the LORY data. These are the coordinates of the evaluation points at which the Jacobians are calculated.

functions-HaberDworkin2015.R - an R source file with custom-made functions that are needed to run the protocol

protocol-HaberDworkin2015.R - an R source file with the protocol to run all the analyses and produce all the figures and tables presented in Haber & Dworkin 2015


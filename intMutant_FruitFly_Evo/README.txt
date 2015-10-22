{\rtf1\ansi\ansicpg1252\cocoartf1265\cocoasubrtf210
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\margl1440\margr1440\vieww22940\viewh12540\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural

\f0\fs24 \cf0 \ul \ulc0 Input:\
\
\ulnone inf_gt.Rdata - R object of the same name with information about the genotypes. 25-by-6 data frame. Row names are genotype ID. Column names are: WT, wild type; allele, indicating the mutant allele; gene; gene abbr; pathway; N, sample size\
\
inf_spec.Rdata - R object of the same name with information about the specimens. 1880-by-4 data frame. rownames are specimen ID. Column names are: WT, wild type; genotype, genotype ID as in inf_gt; repID, replicate ID; centsize, centroid size
\fs26 \
\
XallData.Rdata - R object named Xall with landmark coordinates after Procrustes superimposition. 1880-by-96 numeric matrix. Columns ordered as \'93x1, y1, x2, y2\'85\'94\
\
JxcTPS.dat - a tab limited table with LORY data, using TPS as interpolation function. Data was not regressed on size\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural
\cf0 JxcEBS.dat - a tab limited table with LORY data, using EBSS as interpolation function. Data was not regressed on size\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural
\cf0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural
\cf0 JxcTPSsr.dat - a tab limited table with LORY data, using TPS as interpolation function. Data was regressed on size before running in LORY\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural
\cf0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural
\cf0 JxcEBSsr.dat - a tab limited table with LORY data, using EBSS as interpolation function. Data was regressed on size before running in LORY\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural
\cf0 \
functions-HaberDworkin2015.R - an R source file with custom-made functions that are needed to run the protocol\
\
prtcl-HaberDworkin2015.R - an R source file with the protocol to run all the analyses and produce all the figures and tables presented in Haber & Dworkin 2015\
\
}
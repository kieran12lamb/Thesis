args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library("dplyr")
library("tidyr")
library("Biostrings")
library("ggplot2")
library("cowplot")
library("pals")
source("../sig2aa/functions.R")

trinucleotideContextDF <- ParseComsicSignatures(filePath = args[1], args[3])
marginalFrequencies <- c("T" = 0.25, "A" = 0.25, "G" = 0.25, "C" = 0.25)
frameProb <- c(1/3, 1/3, 1/3)
aminoAcidHits <-ConvertTrinucleotideContextToAimnoAcidContext(trinucleotideContextDF,
                                                              marginalFrequencies, 
                                                              frameProb)
aminoAcidSig <- MapAminoAcidHits(aminoAcidHits)
PlotInAndOutAminoAcidSignatures(aminoAcidSig, filePath = args[2])
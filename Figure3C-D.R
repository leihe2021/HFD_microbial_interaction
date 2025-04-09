#----------------------------------------------------------------------------------#
# Figure 3C and D 
#----------------------------------------------------------------------------------#
source("scripts/figs_scripts/utils.R")
library(tidyverse)
library(data.table)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)

select_gene <- read_rds("LocusZoom_selectGenes.rds")
brca_gwas <- fread("GCST90011804_buildGRCh37.tsv.gz")
input <- read_rds("input.rds")
snps <- read_rds("sig_snp2gene.rds")

rg = 1e5
pcutoff = 1e-5
chr = 3
ens_db = "EnsDb.Hsapiens.v75"
# select SLC4A7 for Roseburia 
pdf(file = "figure_3C_locusZoom_SLC4A7.pdf", 
    width = 7.12, height = 3.62)
plot_LocusZoom(res_coloc = input, gwas = brca_gwas, gene = "SLC4A7", sig.snps = snps)
dev.off()

# gene = RASIP1 for Rumi
pdf(file = "locusZoom_RASIP1.pdf", 
    width = 7.12, height = 3.62)
plot_LocusZoom(res_coloc = input, gwas = brca_gwas, gene = "RASIP1", sig.snps = snps)
dev.off()

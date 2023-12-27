library(tidyverse)
library(msa) # for sequence alignment
library(ape) # for neighbor joining
library(phangorn) # for MP, ML, Boostrapping
library(Biostrings)
library(seqinr)
library(phytools)

setwd("/home/james/Documents/BFX_SL40115")

### Functions ------------------------------------------------------------------

clean_aa_names <- function(vector_obj) {
  names <- vector_obj
  new_names <- c()
  for (x in 1:length(names)) {
    string <- str_split(names, "\\[")[[x]][2]
    string <- sub("]", "", string)
    print(string)
    new_names <- c(new_names, string)
  }
  return(new_names)
}

clean_dna_names <- function(vector_obj) {
  names <- vector_obj
  new_names <- c()
  for (x in 1:length(names)) {
    genus <- str_split(names, " ")[[x]][2]
    species <- str_split(names, " ")[[x]][3]
    string <- paste(genus, species)
    print(string)
    new_names <- c(new_names, string)
  }
  return(new_names)
}

accessions <- function(vector_obj) {
  names <- vector_obj
  new_names <- c()
  for (x in 1:length(names)) {
    string <- str_split(names, " ")[[x]][1]
    # string <- sub("]", "", string)
    print(string)
    new_names <- c(new_names, string)
  }
  return(new_names)
}

do_nj_bootstrap <- function(alignment, protein, branch_lengths, replicates) {
  set.seed(123)
  phydat_obj <- as.phyDat(alignment)
  tree_nj <- nj(dist.hamming(phydat_obj))
  NJtrees <- bootstrap.phyDat(phydat_obj,
                              FUN=function(x)NJ(dist.hamming(x)), bs=replicates)
  png(filename = str_glue("output/{protein}_nj_boostrap.png"), width = 2400, height = 2400, res = 150)
  plotBS(tree_nj, NJtrees, "phylogram")
  dev.off()
  return(tree_nj)
}

do_mp_bootstrap <- function(alignment, protein, branch_lengths, replicates) {
  set.seed(123)
  phydat_obj <- as.phyDat(alignment)
  treeMP <- pratchet(phydat_obj)
  treeMP <- acctran(treeMP, phydat_obj)
  BStrees <- bootstrap.phyDat(phydat_obj, pratchet, bs = replicates, multicore = TRUE)
  
  png(filename = str_glue("output/{protein}_mp_boostrap.png"), width = 2400, height = 2400, res = 150)
  plotBS(treeMP, BStrees, "phylogram") 
  add.scale.bar()
  dev.off()
  return(treeMP)
}

make_cophylo_plot <- function(tree1, tree2, file_name) {
  cophylo_plot <- cophylo(tree1, tree2)
  png(str_glue("output/coplot_{file_name}.png"), width = 2400, height = 2400)
  plot(cophylo_plot)
  dev.off()
}

## Initial script ----------------------------------------------------------------------

# load amino acid sequences
nadh_aa_seqs <-readAAStringSet("sequences/protein/nadh_dehydrogenase.txt")
cyt_aa_seqs <-readAAStringSet("sequences/protein/cyt_b.txt")
orn_aa_seqs <-readAAStringSet("sequences/protein/ornithine_dehydrogenase.txt")

# load DNA sequences
nadh_dna_seqs <-readDNAStringSet("sequences/gene/nadh_dehydrogenase.txt")
cyt_dna_seqs <-readDNAStringSet("sequences/gene/cyt_b.txt")
orn_dna_seqs <-readDNAStringSet("sequences/gene/ornithine_dehydrogenase.txt")

# list of outgroups
outgroups <- c("Malurus cyaneus", "Climacteris rufus", "Epthianura aurifrons") # incomplete

# cleap up names for phylogram
nadh_aa_seqs@ranges@NAMES <- clean_aa_names(nadh_aa_seqs@ranges@NAMES)
cyt_aa_seqs@ranges@NAMES <- clean_aa_names(cyt_aa_seqs@ranges@NAMES)
orn_aa_seqs@ranges@NAMES <- clean_aa_names(orn_aa_seqs@ranges@NAMES)

nadh_dna_seqs@ranges@NAMES <- clean_dna_names(nadh_dna_seqs@ranges@NAMES)
cyt_dna_seqs@ranges@NAMES <- clean_dna_names(cyt_dna_seqs@ranges@NAMES)
orn_dna_seqs@ranges@NAMES <- clean_dna_names(orn_dna_seqs@ranges@NAMES)

### AA NJ trees ----

nadh_aa_nj_tree <- do_nj_bootstrap(alignment = msa(nadh_aa_seqs),
                                   protein = "NADH",
                                   branch_lengths = TRUE,
                                   replicates = 500)
cyt_aa_nj_tree <- do_nj_bootstrap(alignment = msa(cyt_aa_seqs),
                                  protein = "Cytochrome_B",
                                  branch_lengths = TRUE,
                                  replicates = 500)
orn_aa_nj_tree <- do_nj_bootstrap(alignment = msa(orn_aa_seqs),
                                  protein = "Ornithine_decarboxylase",
                                  branch_lengths = TRUE,
                                  replicates = 500)

### AA MP trees ----
nadh_aa_mp_tree <- do_mp_bootstrap(alignment = msa(nadh_aa_seqs),
                                   protein = "NADH",
                                   branch_lengths = TRUE,
                                   replicates = 500)
cyt_aa_mp_tree <- do_mp_bootstrap(alignment = msa(cyt_aa_seqs),
                                  protein = "Cytochrome_B",
                                  branch_lengths = TRUE,
                                  replicates = 500)
orn_aa_mp_tree <- do_mp_bootstrap(alignment = msa(orn_aa_seqs),
                                  protein = "Ornithine_decarboxylase",
                                  branch_lengths = TRUE,
                                  replicates = 500)

# AA tree metrics
treedist(nadh_aa_mp_tree, nadh_aa_nj_tree)
treedist(cyt_aa_mp_tree, cyt_aa_nj_tree)
treedist(orn_aa_mp_tree, orn_aa_nj_tree)

## Modified script -------------------------------------------------------------
# The following code was written following the previous exploratory analysis

# make co-plots but only for NADH as good tree metrics

make_cophylo_plot(nadh_aa_mp_tree, nadh_dna_mp_tree, file_name = "nadh_mp_aa_vs_dna")
make_cophylo_plot(nadh_aa_nj_tree, nadh_aa_mp_tree, file_name = "nadh_aa_nj_vs_mp")

# DNA plot, only for NADH as only the NADH protein plot good enough for comparison
nadh_dna_tree <- do_mp_bootstrap(alignment = msa(nadh_dna_seqs),
                                 protein = "NADH_DNA",
                                 branch_lengths = TRUE,
                                 replicates = 500)

# remove any sequences not present in both DNA and AA datasets
# (these are automatically ignore in the co-plots)
common_tips <- intersect(names(nadh_aa_seqs), names(nadh_dna_seqs))
# tree metrics (NADH DNA vs AA)
treedist(keep.tip(nadh_aa_mp_tree, common_tips), keep.tip(nadh_dna_tree, common_tips))
# make a comparisn plot
make_cophylo_plot(nadh_aa_mp_tree, nadh_dna_tree, file_name = "nadh_dna_vs_aa")

# protein super alignment
# not enough data to do a gene super alignment
combined_alignment <- msa(c(nadh_aa_seqs, cyt_aa_seqs, orn_aa_seqs))
aa_super_alignment_mp_tree <- do_mp_bootstrap(alignment = combined_alignment,
                                              protein = "super_alignment_mp_for_dna",
                                              branch_lengths = TRUE,
                                              replicates = 75)

# finally co plot NADH AA with super alignment of AA
make_cophylo_plot(nadh_aa_mp_tree, aa_super_alignment_mp_tree, file_name = "nadh_aa_vs_super_aa")
common_tips <- intersect(names(nadh_aa_seqs), names(cyt_aa_seqs))
common_tips <- intersect(common_tips, names(orn_aa_seqs))
treedist(keep.tip(nadh_aa_mp_tree, common_tips), keep.tip(aa_super_alignment_mp_tree, common_tips))


## Extra functions -------------------------------------------------------------
# Writing alignments to file

# saveWidth <- getOption("width")
# options(width=100)
# sink("output/ORN_DNA.txt")
# print(msa(orn_dna_seqs), show="complete", halfNrow=-1)
# sink()
# options(width=saveWidth)



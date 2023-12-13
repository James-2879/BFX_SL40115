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

make_nj_tree <- function(alignment, protein, branch_lengths) {
  
  alignment <- msaConvert(alignment, type="seqinr::alignment")
  distance_matrix <- seqinr::dist.alignment(alignment, 
                                            matrix = "identity")
  tree <- nj(distance_matrix)
  
  if (!branch_lengths) {
    png(filename = str_glue("output/{protein}_nj_no_lengths.png"), width = 2400, height = 2400, res = 150)
    plot.phylo(tree, main="Phylogenetic Tree",
               use.edge.length = F)
    mtext(text = "Rooted, no branch lengths")
  } else if (branch_lengths) {
    png(filename = str_glue("output/{protein}_nj_lengths.png"), width = 2400, height = 2400, res = 150)
    plot.phylo(tree, main="Phylogenetic Tree", 
               use.edge.length = T)
    mtext(text = "Rooted, with branch lengths")
  }
  dev.off()
}

make_mp_tree <- function(alignment, protein, branch_lengths) {
  
  phydat_obj <- as.phyDat(alignment)
  tree_nj <- nj(dist.hamming(phydat_obj))
  
  parsimony(tree_nj, phydat_obj) # check parsimony
  best_tree <- optim.parsimony(tree_nj, phydat_obj) # optimize tree
  parsimony(best_tree, phydat_obj) # check parsimony again
  
  if (!branch_lengths) {
    png(filename = str_glue("output/{protein}_mp_no_lengths.png"), width = 2400, height = 2400, res = 150)
    plot.phylo(best_tree, main="Phylogenetic Tree",
               use.edge.length = F)
    mtext(text = "Rooted, no branch lengths")
  } else if (branch_lengths) {
    png(filename = str_glue("output/{protein}_mp_lengths.png"), width = 2400, height = 2400, res = 150)
    plot.phylo(best_tree, main="Phylogenetic Tree", 
               use.edge.length = T)
    mtext(text = "Rooted, with branch lengths")
  }
  dev.off()
}

do_nj_bootstrap <- function(alignment, protein, branch_lengths, replicates) {
  set.seed(123)
  phydat_obj <- as.phyDat(alignment)
  tree_nj <- nj(dist.hamming(phydat_obj))
  NJtrees <- bootstrap.phyDat(phydat_obj,
                              FUN=function(x)NJ(dist.hamming(x)), bs=replicates)
  # treeNJ <- plotBS(tree_nj, NJtrees, "phylogram")
  
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

translate_dna_sequences <- function(sequences, protein) {
  translated_sequences <- Biostrings::translate(sequences, if.fuzzy.codon = "X")
  con <- file(str_glue("sequences/translated/{protein}.txt"), "w")
  for (i in seq_along(translated_sequences)) {
    writeLines(c(paste0(">", names(translated_sequences)[i]), as.character(translated_sequences[i]), "\n"), con)
  }
  close(con)
}

make_super_alignment <- function(alignments) {
  phydat_alignments <- c()
  for (alignment in alignments) {
    converted_alignment <- as.phyDat(alignment)
    phydat_alignments <- c(phydat_alignments, converted_alignment)
  }
  super_alignment <- do.call(cbind, phydat_alignments)
  return(super_alignment)
}

make_cophylo_plot <- function(tree1, tree2, file_name) {
  cophylo_plot <- cophylo(tree1, tree2)
  png(str_glue("output/coplot_{file_name}.png"), width = 2400, height = 2400)
  plot(cophylo_plot)
  dev.off()
}

## Script ----------------------------------------------------------------------

# load amino acid sequences
nadh_aa_seqs <-readAAStringSet("sequences/protein/nadh_dehydrogenase.txt")
cyt_aa_seqs <-readAAStringSet("sequences/protein/cyt_b.txt")
orn_aa_seqs <-readAAStringSet("sequences/protein/ornithine_dehydrogenase.txt")

# load DNA sequences
nadh_dna_seqs <-readDNAStringSet("sequences/gene/nadh_dehydrogenase.txt")
cyt_dna_seqs <-readDNAStringSet("sequences/gene/cyt_b.txt")
orn_dna_seqs <-readDNAStringSet("sequences/gene/ornithine_dehydrogenase.txt")

# translate DNA sequences
translate_dna_sequences(nadh_dna_seqs, protein = "nadh_dehydrogenase")
translate_dna_sequences(cyt_dna_seqs, protein = "cyt_b")
translate_dna_sequences(orn_dna_seqs, protein = "ornithine_dehydrogenase")

# load translated DNA sequences
nadh_translated_dna_seqs <-readAAStringSet("sequences/translated/nadh_dehydrogenase.txt")
cyt_translated_dna_seqs <-readAAStringSet("sequences/translated/cyt_b.txt")
orn_translated_dna_seqs <-readAAStringSet("sequences/translated/ornithine_dehydrogenase.txt")

# list of outgroups
outgroups <- c("Malurus cyaneus", "Climacteris rufus", "Epthianura aurifrons") # incomplete

# cleap up names for phylogram
nadh_aa_seqs@ranges@NAMES <- clean_aa_names(nadh_aa_seqs@ranges@NAMES)
cyt_aa_seqs@ranges@NAMES <- clean_aa_names(cyt_aa_seqs@ranges@NAMES)
orn_aa_seqs@ranges@NAMES <- clean_aa_names(orn_aa_seqs@ranges@NAMES)

nadh_dna_seqs@ranges@NAMES <- clean_dna_names(nadh_dna_seqs@ranges@NAMES)
cyt_dna_seqs@ranges@NAMES <- clean_dna_names(cyt_dna_seqs@ranges@NAMES)
orn_dna_seqs@ranges@NAMES <- clean_dna_names(orn_dna_seqs@ranges@NAMES)

nadh_translated_dna_seqs@ranges@NAMES <- clean_dna_names(nadh_translated_dna_seqs@ranges@NAMES)
cyt_translated_dna_seqs@ranges@NAMES <- clean_dna_names(cyt_translated_dna_seqs@ranges@NAMES)
orn_translated_dna_seqs@ranges@NAMES <- clean_dna_names(orn_translated_dna_seqs@ranges@NAMES)

### User -------------------

# neighbor joining trees
alignment <- msa(nadh_aa_seqs)
nadh_aa_nj_tree <- do_nj_bootstrap(alignment = alignment,
             protein = "NADH",
             branch_lengths = TRUE,
             replicates = 500)

alignment <- msa(cyt_aa_seqs)
cyt_aa_nj_tree <- do_nj_bootstrap(alignment = alignment,
             protein = "Cytochrome_B",
             branch_lengths = TRUE,
             replicates = 500)

alignment <- msa(orn_aa_seqs)
orn_aa_nj_tree <- do_nj_bootstrap(alignment = alignment,
             protein = "Ornithine_decarboxylase",
             branch_lengths = TRUE,
             replicates = 500)

# maximum parsimony trees
alignment <- msa(nadh_aa_seqs)
nadh_aa_mp_tree <- do_mp_bootstrap(alignment = alignment,
             protein = "NADH",
             branch_lengths = TRUE,
             replicates = 500)

alignment <- msa(cyt_aa_seqs)
cyt_aa_mp_tree <- do_mp_bootstrap(alignment = alignment,
             protein = "Cytochrome_B",
             branch_lengths = TRUE,
             replicates = 500)

alignment <- msa(orn_aa_seqs)
orn_aa_mp_tree <- do_mp_bootstrap(alignment = alignment,
             protein = "Ornithine_decarboxylase",
             branch_lengths = TRUE,
             replicates = 500)

# neighbor joining trees
alignment <- msa(nadh_translated_dna_seqs)
nadh_dna_nj_tree <- do_nj_bootstrap(alignment = alignment,
                protein = "NADH_in_silico",
                branch_lengths = TRUE,
                replicates = 500)

alignment <- msa(cyt_translated_dna_seqs)
cyt_dna_nj_tree <- do_nj_bootstrap(alignment = alignment,
                protein = "Cytochrome_B_in_silico",
                branch_lengths = TRUE,
                replicates = 500)

alignment <- msa(orn_translated_dna_seqs)
orn_dna_nj_tree <- do_nj_bootstrap(alignment = alignment,
                protein = "Ornithine_decarboxylase_in_silico",
                branch_lengths = TRUE,
                replicates = 500)

# maximum parsimony trees
alignment <- msa(nadh_translated_dna_seqs)
nadh_dna_mp_tree <- do_mp_bootstrap(alignment = alignment,
                protein = "NADH_in_silico",
                branch_lengths = TRUE,
                replicates = 500)

alignment <- msa(cyt_translated_dna_seqs)
cyt_dna_mp_tree <- do_mp_bootstrap(alignment = alignment,
                protein = "Cytochrome_B_in_silico",
                branch_lengths = TRUE,
                replicates = 500)

alignment <- msa(orn_translated_dna_seqs)
orn_dna_mp_tree <- do_mp_bootstrap(alignment = alignment,
                protein = "Ornithine_decarboxylase_in_silico",
                branch_lengths = TRUE,
                replicates = 500)

# coplots

make_cophylo_plot(nadh_aa_mp_tree, nadh_dna_mp_tree, file_name = "nadh_mp_aa_vs_dna")
make_cophylo_plot(nadh_aa_nj_tree, nadh_aa_mp_tree, file_name = "nadh_aa_nj_vs_mp")

make_cophylo_plot(cyt_aa_mp_tree, cyt_dna_mp_tree, file_name = "cyt_mp_aa_vs_dna")
make_cophylo_plot(cyt_aa_nj_tree, cyt_aa_mp_tree, file_name = "cyt_aa_nj_vs_mp")

make_cophylo_plot(orn_aa_mp_tree, orn_dna_mp_tree, file_name = "orn_mp_aa_vs_dna")
make_cophylo_plot(orn_aa_nj_tree, orn_aa_mp_tree, file_name = "orn_aa_nj_vs_mp")


treedist(cyt_aa_mp_tree, cyt_aa_nj_tree)

matrix1 <- cophenetic.phylo(nadh_aa_mp_tree)
matrix2 <- cophenetic.phylo(nadh_dna_mp_tree)

super_alignment <- make_super_alignment(c(msa(nadh_translated_dna_seqs), msa(cyt_translated_dna_seqs), msa(orn_translated_dna_seqs),
                                          msa(nadh_aa_seqs), msa(cyt_aa_seqs), msa(orn_aa_seqs)))


# do one for the dna sequences as well
 combined <- c(nadh_aa_seqs, cyt_aa_seqs, orn_aa_seqs)
 combined_alignment <- msa(combined)

 orn_dna_mp_tree <- do_mp_bootstrap(alignment = combined_alignment,
                                    protein = "super_alignment",
                                    branch_lengths = TRUE,
                                    replicates = 50)

#  -------



nadh_protein_df <- map_dfr(.x = nadh_aa_seq_names,
              .f = function(x) {
                org_string <- str_split(x, "\\[")[[1]][2]
                org_string <- sub("]", "", org_string)
                acc_string <- str_split(x, " ")[[1]][1]
                vec <- data.frame(organism = org_string, nadh_protein_accession = acc_string)
              })


nadh_aa_seq_names <- nadh_aa_seqs@ranges@NAMES
nadh_protein_df <- map_dfr(.x = nadh_aa_seq_names,
                           .f = function(x) {
                             org_string <- str_split(x, "\\[")[[1]][2]
                             org_string <- sub("]", "", org_string)
                             acc_string <- str_split(x, " ")[[1]][1]
                             vec <- data.frame(organism = org_string, nadh_protein_accession = acc_string)
                           })

cyt_aa_seq_names <-  cyt_aa_seqs@ranges@NAMES
cyt_protein_df <- map_dfr(.x = cyt_aa_seq_names,
                           .f = function(x) {
                             org_string <- str_split(x, "\\[")[[1]][2]
                             org_string <- sub("]", "", org_string)
                             acc_string <- str_split(x, " ")[[1]][1]
                             vec <- data.frame(organism = org_string, cyt_protein_accession = acc_string)
                           })

orn_aa_seq_names <-  orn_aa_seqs@ranges@NAMES
orn_protein_df <- map_dfr(.x = orn_aa_seq_names,
                           .f = function(x) {
                             org_string <- str_split(x, "\\[")[[1]][2]
                             org_string <- sub("]", "", org_string)
                             acc_string <- str_split(x, " ")[[1]][1]
                             vec <- data.frame(organism = org_string, orn_protein_accession = acc_string)
                           })
 
nadh_dna_seq_names <-  nadh_dna_seqs@ranges@NAMES
nadh_dna_df <- map_dfr(.x = nadh_dna_seq_names,
                           .f = function(x) {
                             org_string <- str_split(x, "\\[")[[1]][2]
                             org_string <- sub("]", "", org_string)
                             acc_string <- str_split(x, " ")[[1]][1]
                             vec <- data.frame(organism = org_string, nadh_dna_accession = acc_string)
                           })

cyt_dna_seq_names <-  cyt_dna_seqs@ranges@NAMES
cyt_dna_df <- map_dfr(.x = cyt_dna_seq_names,
                       .f = function(x) {
                         org_string <- str_split(x, "\\[")[[1]][2]
                         org_string <- sub("]", "", org_string)
                         acc_string <- str_split(x, " ")[[1]][1]
                         vec <- data.frame(organism = org_string, cyt_dna_accession = acc_string)
                       })

orn_dna_seq_names <-  orn_dna_seqs@ranges@NAMES
orn_dna_df <- map_dfr(.x = orn_dna_seq_names,
                       .f = function(x) {
                         org_string <- str_split(x, "\\[")[[1]][2]
                         org_string <- sub("]", "", org_string)
                         acc_string <- str_split(x, " ")[[1]][1]
                         vec <- data.frame(organism = org_string, orn_dna_accession = acc_string)
                       })



### Next steps -----------------------------------------------------------------

# maybe look into different subsitution matrices also
# Examine the bootstrap support values on the branches of your phylogenetic tree. Support values indicate the proportion of bootstrap replicates that support a particular branch. Higher values (e.g., 70% or above) are generally considered more reliable.

# Can concatenate protein sequences BUT this is less good if they have different rates of evolution

####


 
 
 
### Printing the alignment -----------------------------------------------------
# print(nadh_alignment, show="complete")

# msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
#                showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
# 
# msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis",
#                showNames="none", showLogo="none", askForOverwrite=FALSE)



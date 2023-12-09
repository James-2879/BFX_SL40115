library(tidyverse)
library(msa) # for sequence alignment
library(ape) # for neighbor joining

setwd("/home/james/Documents/phylogenetics")

### Functions ------------------------------------------------------------------

clean_names <- function(vector_obj) {
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

## Script ----------------------------------------------------------------------

nadh_seqs <-readAAStringSet("/home/james/Documents/phylogenetics/sequences/protein/nadh_dehydrogenase.txt")
cyt_seqs <-readAAStringSet("/home/james/Documents/phylogenetics/sequences/protein/cyt_b.txt")
orn_seqs <-readAAStringSet("/home/james/Documents/phylogenetics/sequences/protein/ornithine_dehydrogenase.txt")

# make names ncier for phylogram
nadh_seqs@ranges@NAMES <- clean_names(nadh_seqs@ranges@NAMES)
cyt_seqs@ranges@NAMES <- clean_names(cyt_seqs@ranges@NAMES)
orn_seqs@ranges@NAMES <- clean_names(orn_seqs@ranges@NAMES)

align_and_tree <- function(sequences, rooted, branch_lengths) {
  
  alignment <- msa(sequences)
  alignment <- msaConvert(alignment, type="seqinr::alignment")
  
  ### Printing the alignment -----------------------------------------------------
  # print(nadh_alignment, show="complete")
  
  # msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
  #                showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
  # 
  # msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis",
  #                showNames="none", showLogo="none", askForOverwrite=FALSE)
  
  ### Making trees ---------------------------------------------------------------
  
  distance_matrix <- seqinr::dist.alignment(alignment, 
                                            matrix = "identity")
  
  tree <- nj(distance_matrix)
  
  ## agh need outgroups
  
  if (!rooted & !branch_lengths) {
  plot.phylo(tree, main="Phylogenetic Tree",
             type = "unrooted",
             use.edge.length = F)
  mtext(text = "UNrooted, no branch lengths")
  } else if (!rooted & branch_lengths) {
    plot.phylo(tree, main="Phylogenetic Tree",
               type = "unrooted",
               use.edge.length = T)
    mtext(text = "UNrooted, with branch lengths")
  } else if (rooted & !branch_lengths) {
  plot.phylo(tree, main="Phylogenetic Tree",
             use.edge.length = F)
  mtext(text = "Rooted, no branch lenths")
  } else if (rooted & branch_lengths) {
  plot.phylo(tree, main="Phylogenetic Tree", 
             use.edge.length = T)
  mtext(text = "Rooted, with branch lenths")
  }
}

align_and_tree(sequences = cyt_seqs,
               rooted = TRUE,
               branch_lengths = TRUE)















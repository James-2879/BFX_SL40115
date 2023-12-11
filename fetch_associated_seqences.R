library(tidyverse)
library(rentrez) # NCBI queries
library(Biostrings)

setwd("/home/james/Documents/BFX_SL40115")

### Functions ------------------------------------------------------------------

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

entrez_gene_from_protein <- function(accessions) {
  transcripts <- c()
  print("Script will pause for 1 second between requests.")
  for (accession in accessions) {
    result <- tryCatch(
      expr = {
        print(paste0("Running accession: ", accession))
        search_result <- entrez_search(db = "protein", term = accession, retmax = 1)
        linked_seq_ids <- entrez_link(dbfrom="protein", id=search_result$ids, db="nuccore")
        
        linked_transcripts <- linked_seq_ids$links$protein_nuccore
        transcript <- entrez_fetch(db = "nuccore", id = linked_transcripts, rettype = "fasta")
        transcripts <- c(transcripts, transcript)
        remove(search_result, linked_seq_ids, linked_transcripts, transcript)
        Sys.sleep(1) # prevent 429 errors (too many requests)
      },
      error = function(err) {
        print(paste0("Failed for accession: ", accession))
      }
    )
  }
  return(transcripts)
}

write_sequences <- function(vector, file) {
  collapsed_vector <- paste(vector, collapse = "\n\n")
  writeLines(collapsed_vector, file)
}

## Script ----------------------------------------------------------------------

# load protein sequences
nadh_seqs <-readAAStringSet("/home/james/Documents/BFX_SL40115/sequences/protein/nadh_dehydrogenase.txt")
cyt_seqs <-readAAStringSet("/home/james/Documents/BFX_SL40115/sequences/protein/cyt_b.txt")
orn_seqs <-readAAStringSet("/home/james/Documents/BFX_SL40115/sequences/protein/ornithine_dehydrogenase.txt")

# pull accessions
nadh_accessions <- accessions(nadh_seqs@ranges@NAMES)
cyt_accessions <- accessions(cyt_seqs@ranges@NAMES)
orn_accessions <- accessions(orn_seqs@ranges@NAMES)

# fetch associated sequences
nadh_sequences <- entrez_gene_from_protein(nadh_accessions)
cyt_sequences <- entrez_gene_from_protein(cyt_accessions)
orn_sequences <- entrez_gene_from_protein(orn_accessions)

# write sequences to file
write_sequences(nadh_sequences, "sequences/gene/nadh_dehydrogenase.txt")
write_sequences(cyt_sequences, "sequences/gene/cyt_b.txt")
write_sequences(orn_sequences, "sequences/gene/ornithine_dehydrogenase.txt")









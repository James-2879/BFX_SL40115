library(tidyverse)
library(Biostrings)

setwd("/home/james/Documents/BFX_SL40115")

# get amino acid sequences headers
nadh_aa_seq_names <-readAAStringSet("sequences/protein/nadh_dehydrogenase.txt") %>% names()
cyt_aa_seq_names <-readAAStringSet("sequences/protein/cyt_b.txt") %>% names()
orn_aa_seq_names <-readAAStringSet("sequences/protein/ornithine_dehydrogenase.txt") %>% names()

# get DNA sequence headers
nadh_dna_seq_names <-readDNAStringSet("sequences/gene/nadh_dehydrogenase.txt") %>% names()
cyt_dna_seq_names <-readDNAStringSet("sequences/gene/cyt_b.txt") %>% names()
orn_dna_seq_names <-readDNAStringSet("sequences/gene/ornithine_dehydrogenase.txt") %>% names()

nadh_protein_df <- map_dfr(.x = nadh_aa_seq_names,
                           .f = function(x) {
                             org_string <- str_split(x, "\\[")[[1]][2]
                             org_string <- sub("]", "", org_string)
                             acc_string <- str_split(x, " ")[[1]][1]
                             vec <- data.frame(organism = org_string, nadh_protein_accession = acc_string)
                           })


# nadh_aa_seq_names <- nadh_aa_seqs@ranges@NAMES
nadh_protein_df <- map_dfr(.x = nadh_aa_seq_names,
                           .f = function(x) {
                             org_string <- str_split(x, "\\[")[[1]][2]
                             org_string <- sub("]", "", org_string)
                             acc_string <- str_split(x, " ")[[1]][1]
                             vec <- data.frame(organism = org_string, nadh_protein_accession = acc_string)
                           })

# cyt_aa_seq_names <-  cyt_aa_seqs@ranges@NAMES
cyt_protein_df <- map_dfr(.x = cyt_aa_seq_names,
                          .f = function(x) {
                            org_string <- str_split(x, "\\[")[[1]][2]
                            org_string <- sub("]", "", org_string)
                            acc_string <- str_split(x, " ")[[1]][1]
                            vec <- data.frame(organism = org_string, cyt_protein_accession = acc_string)
                          })

# orn_aa_seq_names <-  orn_aa_seqs@ranges@NAMES
orn_protein_df <- map_dfr(.x = orn_aa_seq_names,
                          .f = function(x) {
                            org_string <- str_split(x, "\\[")[[1]][2]
                            org_string <- sub("]", "", org_string)
                            acc_string <- str_split(x, " ")[[1]][1]
                            vec <- data.frame(organism = org_string, orn_protein_accession = acc_string)
                          })

# nadh_dna_seq_names <-  nadh_dna_seqs@ranges@NAMES
nadh_dna_df <- map_dfr(.x = nadh_dna_seq_names,
                       .f = function(x) {
                         org_string_1 <- str_split(x, " ")[[1]][2]
                         org_string_2 <- str_split(x, " ")[[1]][3]
                         org_string <- paste(org_string_1, org_string_2)
                         acc_string <- str_split(x, " ")[[1]][1]
                         vec <- data.frame(organism = org_string, nadh_dna_accession = acc_string)
                       })

# cyt_dna_seq_names <-  cyt_dna_seqs@ranges@NAMES
cyt_dna_df <- map_dfr(.x = cyt_dna_seq_names,
                      .f = function(x) {
                        org_string_1 <- str_split(x, " ")[[1]][2]
                        org_string_2 <- str_split(x, " ")[[1]][3]
                        org_string <- paste(org_string_1, org_string_2)
                        acc_string <- str_split(x, " ")[[1]][1]
                        vec <- data.frame(organism = org_string, cyt_dna_accession = acc_string)
                      })

# orn_dna_seq_names <-  orn_dna_seqs@ranges@NAMES
orn_dna_df <- map_dfr(.x = orn_dna_seq_names,
                      .f = function(x) {
                        org_string_1 <- str_split(x, " ")[[1]][2]
                        org_string_2 <- str_split(x, " ")[[1]][3]
                        org_string <- paste(org_string_1, org_string_2)
                        acc_string <- str_split(x, " ")[[1]][1]
                        vec <- data.frame(organism = org_string, orn_dna_accession = acc_string)
                      })


remove(nadh_aa_seq_names, cyt_aa_seq_names, orn_aa_seq_names,
       nadh_dna_seq_names, cyt_dna_seq_names, orn_dna_seq_names)

accessions <- nadh_protein_df %>% 
  left_join(cyt_protein_df) %>% 
  left_join(orn_protein_df) %>% 
  left_join(nadh_dna_df) %>% 
  left_join(cyt_dna_df) %>% 
  left_join(orn_dna_df) %>% 
  rename(organism = "Organism",
         nadh_protein_accession = "MT-ND2 (Protein)",
         cyt_protein_accession = "MT-CYB (Protein)",
         orn_protein_accession = "ODC1 (Protein)",
         nadh_dna_accession = "MT-ND2 (Gene)",
         cyt_dna_accession = "MT-CYB (Gene)",
         orn_dna_accession = "ODC1 (Gene)")

remove(nadh_protein_df, cyt_protein_df, orn_protein_df,
       nadh_dna_df, cyt_dna_df, orn_dna_df)





#!/usr/bin/env Rscript
#Validation dataset generation script
#Loading the required packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

#Importing external arguments
ext_args <- commandArgs(trailingOnly = T)
SCAFF_LEN_MIN <- as.numeric(ext_args[1])
SCAFF_LEN_BREAKS <- as.numeric(unlist(strsplit(ext_args[2], ",")))
TOTAL_SEQ_N <- as.numeric(ext_args[3])
SEED <- as.numeric(ext_args[4])
SPONGE_FRAG_N <- as.numeric(ext_args[5])
CONTAM_FRAG_N <- as.numeric(ext_args[6])

#Defining functions
#Importing FASTA files from directory (function)
import_data <- function(directory) {
 files_in_directory <- grep(".fa$", list.files(path = directory), value = T)
 imported_sequence_list <- list()
 for(i in 1 : length(files_in_directory)) {
  imported_sequences <- readDNAStringSet(paste(directory, files_in_directory[i], sep = "/"))
  names(imported_sequences) <- sub(" .*", "", names(imported_sequences))
  names(imported_sequences) <- paste(names(imported_sequences), files_in_directory[i], sep = "_")
  imported_sequence_list[[i]] <- imported_sequences
 }
 imported_sequence_list <- do.call(c, imported_sequence_list)
 return(imported_sequence_list)
}

#Importing analysed sponge scaffolds from directory. Defining scaffold length bins
analysed_sponge_scaffolds <- import_data("analysed_sponge_genomes")
analysed_sponge_scaffolds <- analysed_sponge_scaffolds[width(analysed_sponge_scaffolds) >= SCAFF_LEN_MIN]
analysed_sponge_scaffold_length <- width(analysed_sponge_scaffolds)
analysed_sponge_scaffold_length_bins <- cut(analysed_sponge_scaffold_length, breaks = SCAFF_LEN_BREAKS, include.lowest = T)
length_bin_proportion <- rev(table(analysed_sponge_scaffold_length_bins) / length(analysed_sponge_scaffold_length_bins))
rm(analysed_sponge_scaffolds, analysed_sponge_scaffold_length_bins)

#Importing validated sponge chromosomes from directory
sponge_chromosomes <- import_data("validation_sponge_chromosomes")
sponge_chromosomes <- sponge_chromosomes[width(sponge_chromosomes) >= SCAFF_LEN_MIN]
names(sponge_chromosomes) <- paste(names(sponge_chromosomes), "sponge_chromosome", sep = "_")

#Importing spike-in contaminant scaffolds from directory
contaminant_scaffolds <- import_data("validation_contaminant_scaffolds")
contaminant_scaffolds <- contaminant_scaffolds[width(contaminant_scaffolds) >= SCAFF_LEN_MIN]
names(contaminant_scaffolds) <- paste(names(contaminant_scaffolds), "contaminant_scaffold", sep = "_")

#Joining sponge and contaminant sequences
validation_sequences <- c(sponge_chromosomes, contaminant_scaffolds)
rm(sponge_chromosomes, contaminant_scaffolds)

#Fragmenting sequences
length_bin_lower_bounds <- rev(SCAFF_LEN_BREAKS[-length(SCAFF_LEN_BREAKS)])
length_bin_upper_bounds <- rev(SCAFF_LEN_BREAKS[-1])
length_bin_number <- floor(length_bin_proportion * TOTAL_SEQ_N)

validation_dataset <- list()
for(i in 1 : length(length_bin_number)) {
 sample_N <- length_bin_number[i]
 length_bin_lower_bound <- length_bin_lower_bounds[i]
 length_bin_upper_bound <- length_bin_upper_bounds[i]
 for(j in 1 : sample_N) {
  length_sample_max <-  suppressWarnings(max(width(validation_sequences)[width(validation_sequences) >= length_bin_lower_bound & width(validation_sequences) <= length_bin_upper_bound]))
  if(is.infinite(length_sample_max) == F) {
   set.seed(SEED + width(validation_sequences[length(validation_sequences)]))
   length_sample <- sample(length_bin_lower_bound : length_sample_max, size = 1)
   candidate_scaffs_to_frag <- validation_sequences[width(validation_sequences) >= length_sample]
   set.seed(SEED + i + j)
   scaff_to_frag <- sample(candidate_scaffs_to_frag, size = 1)
   scaff_to_frag_len <- width(scaff_to_frag)
   set.seed(SEED + i + j)
   start_sample <- sample(1 : (scaff_to_frag_len - length_sample + 1), size = 1)
   end_sample <- start_sample + length_sample - 1
   scaff_frag <- DNAStringSet(scaff_to_frag[[1]][start_sample : end_sample])
   names(scaff_frag) <- paste(names(scaff_to_frag), "frag", start_sample, end_sample, sep = "_")
   if((grepl("sponge_chromosome", names(scaff_to_frag), fixed = T) & str_count(names(scaff_to_frag), "frag") < SPONGE_FRAG_N) | (grepl("contaminant_scaffold", names(scaff_to_frag), fixed = T) & str_count(names(scaff_to_frag), "frag") < CONTAM_FRAG_N)) {
    frag_scaff_five_prime <- c()
    frag_scaff_three_prime <- c()
    if(start_sample > 1) {
     frag_scaff_five_prime <- DNAStringSet(scaff_to_frag[[1]][1 : start_sample - 1])
     names(frag_scaff_five_prime) <- paste(names(scaff_to_frag), "frag", "1", start_sample - 1, sep = "_")
    }
    if(end_sample < scaff_to_frag_len) {
     frag_scaff_three_prime <- DNAStringSet(scaff_to_frag[[1]][(end_sample + 1) : scaff_to_frag_len])
     names(frag_scaff_three_prime) <- paste(names(scaff_to_frag), "frag", end_sample + 1, scaff_to_frag_len, sep = "_")
    }
    validation_sequences <- c(validation_sequences, frag_scaff_five_prime, frag_scaff_three_prime)
   }
   validation_dataset <- c(validation_dataset, scaff_frag)
   validation_sequences <- validation_sequences[names(validation_sequences) != names(scaff_to_frag)]
  } else {
   break
  }
 }
}
validation_dataset <- do.call(c, validation_dataset)

#Removing sequences with more than 10% ambiguous bases
n_counts <- vcountPattern("N", validation_dataset, fixed = T)
n_fraction <- n_counts / width(validation_dataset)
validation_dataset <- validation_dataset[n_fraction <= 0.1]

#Writing the final validation dataset to file
writeXStringSet(validation_dataset, file = "validation_dataset.fa", format = "fasta")

#Plotting the distribution of sequence length in the validation dataset and sponge genome assemblies analysed for contamination (Supplementary Figure 1)
sequence_length_dt <- data.table(category = c(rep("Analyzed sponges", length(analysed_sponge_scaffold_length)), rep("Validation dataset", length(validation_dataset))), scaff_length = c(analysed_sponge_scaffold_length, width(validation_dataset)))
supp_fig1 <- ggplot(data = sequence_length_dt, aes(x = scaff_length, fill = category)) +
 geom_density(aes(y = after_stat(density / sum(density))), alpha = 0.5, adjust = 2) +
 theme_minimal() +
 scale_fill_manual(values = c("#DC0000FF", "#4FA3BCFF")) +
 theme(legend.position = "bottom") +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 7)) +
 theme(legend.key.size = unit(0.3, "cm")) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70"))  +
 scale_x_log10(breaks = 10^c(3, 4, 5, 6, 7, 8), labels = scales::label_log()) +
 xlab("Sequence length") +
 ylab("Relative density") +
 theme(axis.title = element_text(size = 7)) +
 theme(axis.text = element_text(size = 6)) 

ggsave(supp_fig1, file = "Supplementary_figure1.tiff", height = 3, width = 4, dpi = 600, bg = "white")
# script with an example run for wisard

library(devtools)
setwd("/home/abc/wisard/") # change depending on your system
install("wisard")
library(wisard)
# reading an xml file
xml_file <- system.file("extdata", "example.xml", package = 'wisard')
blast_result <- read_blast_xml(xml_file, grange_output = T, keep_sequences = F)
blast_result$alignments



# calculating the BLAST sum score of a set of alignments:
sum_score(blast_result$alignments$raw_score, blast_result$metadata)

# filter records:
blast_result$alignments <- filter_blast(blast_result$alignments,
                                        min_len = 10,
                                        max_e = 1,
                                        min_identity = 0)

# run wis
wis_results <- run_get_wis(blast_result$alignments, is_circular = T,
                           use_strand = T,
                           score_col = "raw_score")



# run greedy
greedy_results <- greedy_algorithm(blast_result$alignments,
                                   max_score = "raw_score",
                                   is_circular = T,
                                   use_strand = T)
length(blast_result$alignments)
length(greedy_results)

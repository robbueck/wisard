# run wisard on HIV annotation, can be deleted I think

options(warn=2)
library(GenomicRanges)
library(checkmate)
library(XML)
library(GenomicFeatures)
library(ggbio)

# library(ggplot2)
options(error=traceback)
source("../R/WIS_functions.R")
source("../R/greedy.R")
source("../R/read_blast.R")

dirname <- "/data/meyer/rob/phages/genome_alignments/method_comparison/other_data/HIV"
dirname <- "/home/abc/wisard/HIV/"
setwd(dirname)
gff_files <- list.files(path=".", pattern="*gff", full.names=TRUE, recursive=FALSE)

run_test <- function (in_file){
  print(in_file)
  #data <- makeTxDbFromGFF(in_file, format = "gff3")
  data <- read.table("HIV1_KC503852.gff", quote = "", sep = "\t", comment.char = "#", header = F, skip = 7)
  colnames(data) <- c("genome", "source", "type", "start", "end", "something", "strand", "phase", "comments")
  granges_data <- GRanges(seqnames = data$genome, ranges = IRanges(start = data$start, end = data$end), strand = data$strand, mcols = data$comments)
  #granges_data <- transcripts(data)
  if (length(granges_data) < 1){ return(0)}
  granges_data$score <- width(granges_data)
  isCircular(granges_data) <- F
  length(granges_data)
  seqlengths((granges_data))
  autoplot(granges_data)



  test_wis <- run_get_wis(granges_data, is_circular = F, use_strand = F, score_col = "score")
  test_greedy <- greedy_algorithm(granges_data, max_score = "score", is_circular = F, use_strand = F)
  length(test_wis)
  length(test_greedy)
  (wis_sum <- sum(width(test_wis)))
  (greedy_sum <- sum(width(test_greedy)))
  return((wis_sum - greedy_sum) / greedy_sum) * 100


}

outstats <- lapply(gff_files, run_test)


in_file <- "SIV_U79412.gff"
data <- makeTxDbFromGFF(in_file, format = "gff3")
data <- read.table(in_file, quote = "", sep = "\t", comment.char = "#", header = F, skip = 7)
colnames(data) <- c("genome", "source", "type", "start", "end", "something", "strand", "phase", "comments")
granges_data <- GRanges(seqnames = data$genome, ranges = IRanges(start = data$start, end = data$end), strand = data$strand, mcols = data$comments)
granges_data <- transcripts(data)
granges_exons <- exons(data)
granges_data$score <- width(granges_data)
isCircular(granges_data) <- F
length(granges_data)
seqlengths((granges_data))
granges_data$type <- "transcripts"
granges_exons$type <- "exons"
test_data <- c(granges_data, granges_exons)
autoplot(test_data, facets = type~seqnames)



test_wis <- run_get_wis(granges_data, is_circular = F, use_strand = F, score_col = "score")
test_greedy <- greedy_algorithm(granges_data, max_score = "score", is_circular = F, use_strand = F)
length(test_wis)
length(test_greedy)
(wis_sum <- sum(width(test_wis)))
(greedy_sum <- sum(width(test_greedy)))
return((wis_sum - greedy_sum) / greedy_sum) * 100

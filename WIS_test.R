options(warn=0)
library(GenomicRanges)
library(checkmate)
library(XML)
library(ggplot2)
options(error=traceback)
source("/data/meyer/rob/phages/robert_phage/wisard/R/WIS_functions.R")
source("/data/meyer/rob/phages/robert_phage/wisard/R/greedy.R")
source("/data/meyer/rob/phages/robert_phage/wisard/R/read_blast.R")
library(ggbio)


autoplot(out_data_1$alignments, aes(color = strand, fill = strand), facets = strand~seqnames)
autoplot(out_data_2$alignments, aes(color = strand, fill = strand), facets = strand~seqnames)



#create a test range object
test_case_plus <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 500), end = seq(700, 32108, 500)), strand = "+", score = sample(10:1000, 63, replace = T), Hsp_hit.frame = sample(1:3, 63, replace = T))
test_case_plus <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 1000), end = seq(700, 32108, 1000)), strand = "+", score = sample(10:1000, 32, replace = T), Hsp_hit.frame = sample(1:3, 32, replace = T))
# test_case_plus_one <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 200), end = seq(700, 32108, 200)), strand = "+", score = sample(10:1000, 158, replace = T), Hsp_hit.frame = 1)
# test_case_plus_two <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 200), end = seq(700, 32108, 200)), strand = "+", score = sample(10:1000, 158, replace = T), Hsp_hit.frame = 2)
# test_case_plus_three <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 200), end = seq(700, 32108, 200)), strand = "+", score = sample(10:1000, 158, replace = T), Hsp_hit.frame = 3)
test_case_minus <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 1140800, 50), end = seq(700, 1210800, 50)), strand = "-", score = sample(10:100000, 24203, replace = T), Hsp_hit.frame = sample(-1:-3, 24203, replace = T))
# test_case_minus_one <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 200), end = seq(700, 32108, 200)), strand = "-", score = sample(10:1000, 158, replace = T), Hsp_hit.frame = -1)
# test_case_minus_two <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 200), end = seq(700, 32108, 200)), strand = "-", score = sample(10:1000, 158, replace = T), Hsp_hit.frame = -2)
# test_case_minus_three <- GRanges(seqnames = seqnames(blast_data$alignments)[1], ranges = IRanges(seq(1, 31408, 200), end = seq(700, 32108, 200)), strand = "-", score = sample(10:1000, 158, replace = T), Hsp_hit.frame = -3)
#test_case_save <- c(test_case_plus_one, test_case_plus_two, test_case_plus_three, test_case_minus_one, test_case_minus_two, test_case_minus_three)
test_case_save <- c(test_case_plus, test_case_minus)
names(test_case_save) <- 1:length(test_case_save)
isCircular(test_case_save) <- T
seqlengths(test_case_save) <- seqlengths(blast_data$alignments)
test_case_save$score <- sample(10:1000, length(test_case_save), replace = T)
test_case <- test_case_save

test_case <- test_case[which(test_case$score > 100)]

autoplot(test_case, aes(color = strand, fill = strand))
source("/data/meyer/rob/phages/robert_phage/wisard/R/WIS_functions.R")
# blast_data$alignments
# test_case
#test_case <- test_case[1:150]
test_wis <- run_get_wis(test_case, is_circular = T, use_strand = F,
                                           score_col = "score", overlap = 0,
                                           score_f = function(x, y, n) { x + y })
# test_wis$alignments <- sort(test_wis$alignments, by = ~start)
# test_wis$alignments <- test_wis$alignments[1:30]

autoplot(test_wis, aes(color = strand, fill = strand))
# , facets = Hsp_hit.frame ~ seqnames
length(test_case)
length(test_wis)
overlaps1 <- countOverlaps(test_case, test_case, ignore.strand = T)
overlaps2 <- countOverlaps(test_wis, test_wis, ignore.strand = T)
cov_test <- coverage(test_case)
cov_test > 0

# check all overlaps
test_s <- test_case_plus
test_q <- test_case_minus
test <- test_s
test$query <- test_q

source("/data/meyer/rob/phages/robert_phage/wisard/R/read_blast.R")

setwd("/data/meyer/rob/phages/robert_phage/genome_alignments/method_comparison/bacteria_comparison/")
blast_data <-read_blast_xml("GCF_000829055_GCF_000831645_ab_tblastx.xml", keep_sequences = F)

###################################################
# random test case
source("/data/meyer/rob/phages/robert_phage/wisard/R/WIS_functions.R")

starts <- sample(1:100000, 300, replace = T)
widths <- sample(1:9000, 300, replace = T)
ends <- starts + widths
score <- sample(1:1000, 30, replace = T)
score <- widths
test_case <- GRanges(seqnames = "bla", ranges = IRanges(start = starts, end = ends), strand = "+", score = score)
names(test_case) <- 1:length(test_case)
isCircular(test_case) <- T
seqlengths(test_case) <- 99000
autoplot(test_case, aes(color = score, fill = score))
test_wis <- run_get_wis(test_case, is_circular = T, use_strand = F, score_col = "score", overhang = 1000, score_f = function(x, y, n) { x + y })
test_greedy <- greedy_algorithm(test_case, max_score = "score", overhang = 1000, is_circular = T, use_strand = F)
(sum_wis <- sum(test_wis$score))
(sum_greedy <- sum(test_greedy$score))
sum_wis/sum_greedy
mean( coverage( test_case )[which( coverage( test_case ) > 0 )])
mean( coverage( test_case ))
autoplot(test_wis, aes(color = score, fill = score))
autoplot(test_greedy, aes(color = score, fill = score))


################################################################
# functino test of second_range argument
source("/data/meyer/rob/phages/robert_phage/wisard/R/WIS_functions.R")
source("/data/meyer/rob/phages/robert_phage/wisard/R/greedy.R")

dirname <- '/data/meyer/rob/phages/robert_phage/genome_alignments/method_comparison/bacteria_comparison/'
setwd(dirname)
blast_files <- list.files(path=dirname, pattern="*_ab_blastn.xml", full.names=TRUE, recursive=FALSE)
sizes <- file.info(blast_files)$size
in_file <- blast_files[which.min(sizes)]
in_file <- "GCF_003606365_GCF_900560965_ab_blastn.xml"
blast_data <- read_blast_xml(in_file, grange_output = T, keep_sequences = F)
granges_data <- blast_data$alignments
granges_data$score <- granges_data$Hsp_align.len
granges_wis <- two_ranges_run_get_wis(granges_data,
            is_circular = F,
            overhang = 0,
            use_strand = F,
            score_col = "score",
            overlap = overlap,
            score_f = function(x, y, n) {x + y},
            second_range = "query_ranges"
)
granges_greedy <-  two_ranges_greedy(granges_data,
                                   is_circular = F,
                                   overhang = 0,
                                   use_strand = F,
                                   max_score = "score",
                                   overlap = 0, second_range = "query_ranges")
better_granges_greedy <-  greedy_algorithm(granges_data,
                                     is_circular = F,
                                     overhang = 0,
                                     use_strand = F,
                                     max_score = "score",
                                     overlap = 0, second_range = "query_ranges")

one_granges_greedy <- greedy_algorithm(granges_data,
                                       is_circular = F,
                                       overhang = 0,
                                       use_strand = F,
                                       max_score = "score",
                                       overlap = 0)



one_granges_wis <- run_get_wis(granges_data,
                                      is_circular = F,
                                      overhang = 0,
                                      use_strand = F,
                                      score_col = "score",
                                      overlap = overlap,
                                      score_f = function(x, y, n) {x + y})
granges_wis_greedy <- greedy_algorithm(one_granges_wis,
                                       is_circular = F,
                                       overhang = 0,
                                       use_strand = F,
                                       max_score = "score",
                                       overlap = 0, second_range = "query_ranges")

granges_wis_invert <- run_get_wis_two_ranges(granges_data,
                                      is_circular = F,
                                      overhang = 0,
                                      use_strand = F,
                                      score_col = "score",
                                      overlap = overlap,
                                      score_f = function(x, y, n) {x + y},
                                      second_range = "query_ranges")


(sum_greedy <- sum(granges_greedy$score))
(sum_wis_invert <- sum(granges_wis_invert$score))
(sum_wis <- sum(granges_wis$score))
(sum_better_greedy <- sum(better_granges_greedy$score))
(sum_one_greedy <- sum(one_granges_greedy$score))
(sum_one_wis <- sum(one_granges_wis$score))
(sum_granges_wis_greedy <- sum(granges_wis_greedy$score))
autoplot(granges_data$query_ranges)
autoplot(one_granges_greedy$query_ranges)
autoplot(one_granges_wis$query_ranges)

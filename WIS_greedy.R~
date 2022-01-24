options(warn=2)
library(GenomicRanges)
library(checkmate)
library(XML)

# library(ggplot2)
options(error=traceback)
source("/data/meyer/rob/phages/robert_phage/wisard/R/WIS_functions.R")
source("/data/meyer/rob/phages/robert_phage/wisard/R/greedy.R")
source("/data/meyer/rob/phages/robert_phage/wisard/R/read_blast.R")


library(parallel)
library(data.table)



# function to run everything, input: iput filename and filter criterias.
run_everything <- function(in_file, len_filter=0, e_filter=10,
                           ident_filter=0.0, is_circular_g = F,
                           overhang = 0, overlap = 0, is_stranded = F, framed) {
  # read input, calculate score for each seq and apply filter
  blast_data <- read_blast_xml(in_file, grange_output = T, keep_sequences = F)
  print(length(blast_data))
  if (length(blast_data) < 1) {
    return(list())
  }
  blast_data$alignments <- filter_blast(blast_data$alignments,
                                        min_len = len_filter,
                                        max_e = e_filter,
                                        min_identity = ident_filter)

  # create second  range for query compatibility


  # check if anything passed the filter
  if (length(blast_data$alignments) < 1) {
    warning("Dataset empty after filtering, continue anyway")
    out_list <- list(final = blast_data,
                     raw = blast_data,
                     empty = T,
                     wis_time = 0,
                     greedy_time = 0)
    return(out_list)
  } else {
    # create score to maximize over
    mdata <- blast_data$metadata
    # blast_data$alignments$score <- (blast_data$alignments$raw_score *
    #                                   mdata$Statistics_lambda - 1 -
    #                                   log(mdata$Statistics_kappa *
    #                                         mdata$query_len *
    #                                         mdata$Statistics_db.len))
    start_g <- Sys.time()
    out_greedy <- list( alignments = greedy_algorithm(blast_data$alignments,
                                                      is_circular = is_circular_g,
                                                      overhang = overhang,
                                                      use_strand = is_stranded,
                                                      max_score = "raw_score",
                                                      overlap = overlap),
                        mdata = mdata)
    end_g <- Sys.time()
    cat("Greedy algorithm needed: ", (end_g - start_g), "\n")
    out_list <- list(greedy = out_greedy,
                     raw_in = blast_data,
                     empty = F,
                     greedy_time = (end_g - start_g))
    return(out_list)
  }
}



#################################### RUN #####################################
start_total <- Sys.time()
dirname <- '/data/meyer/rob/phages/robert_phage/genome_alignments/method_comparison/memory_comp_bact_greedy'
setwd(dirname)
circular <- T
blast_overhang <- 0
threads <- 3
outstats_name <- "outstats_greedy_comp.txt"

# filter criteria
len <- 20
e <- 1
identity <- 0
allowed_overlap <- 0
stranded <- F
#frame_column <-"Hsp_hit.frame"



# create data.table to collect stats for all analyzed files
outstats <- data.frame()
# loop through all non-empty files, calculate statistics and write wis to file
blast_files <- list.files(path=dirname, pattern="*365_ab_blastn.xml", full.names=TRUE, recursive=FALSE)
sizes <- file.info(blast_files)$size
blast_files <-  blast_files[which(sizes != 0)]




outstats <- mclapply(blast_files, function(x) {
  print(x)
  default_time <- Sys.time()
  out <- run_everything(x, len_filter = len, e_filter = e, ident_filter = identity, is_circular_g = circular, overlap = allowed_overlap, is_stranded = stranded)
  outfile <- gsub(".xml", "max_raw_wis.tab", x)
  # outfile <- gsub(".xml", "_wis.tab", gsub("vs_aeruginosa_pao1_", "", x))

  #add statistics to main stat output
  gname <- gsub(".xml", "", basename(x))
  default_time <- default_time - default_time
  if (length(out) == 0 || out$empty) {
    stats <- data.frame(genome=gname,
                        wis_file=outfile,
                        cov_raw=0,
                        cov_greedy=0,
                        N_raw=0,
                        N_greedy=0,
                        mean_pident_greedy=0,
                        mean_len_greedy=0,
                        total_length_greedy=0,
                        mean_e_greedy=0,
                        max_e_greedy=0,
                        sum_score_greedy=0,
                        average_raw_greedy=0,
                        time_greedy=default_time
    )
    return(stats)
  }

  mdata <- out$final$mdata
  greedy_alis <- out$greedy$alignments


  stats <- data.table(genome=gname,
                      wis_file=outfile,
                      cov_raw=mean(coverage(out$raw$alignments)),
                      cov_greedy=mean(coverage(greedy_alis)),
                      N_raw=length(out$raw$alignments),
                      N_greedy=length(greedy_alis),
                      mean_pident_greedy=mean((greedy_alis$Hsp_identity / greedy_alis$Hsp_align.len)),
                      mean_len_greedy=mean(greedy_alis$Hsp_align.len),
                      total_length_greedy=sum(width(greedy_alis)),
                      mean_e_greedy=mean(greedy_alis$Hsp_evalue),
                      max_e_greedy=max(greedy_alis$Hsp_evalue),
                      sum_score_greedy=sum_score(greedy_alis$raw_score, out$final$mdata),
                      average_raw_greedy=mean(greedy_alis$raw_score),
                      time_greedy=out$greedy_time
  )

  # print(outfile)
  # out_pdf <- gsub(".xml", "_wis_stats.pdf", gsub("vs_aeruginosa_pao1_", "", x))
  # out_ggplot <- gsub(".xml", "_wis_gen_plot.png", gsub("vs_aeruginosa_pao1_", "", x))
  out_table <- data.frame(seqnames=seqnames(greedy_alis),
                          start=start(greedy_alis),
                          end=end(greedy_alis),
                          strand=strand(greedy_alis),
                          mcols(greedy_alis),
                          query_len=mdata$query_len)

  write.table(out_table, file = outfile, append = F, sep = "\t")

  # library(ggplot2)
  # genome_plot <- ggplot(data=out_table, aes(x=start, y=Hsp_query.from, group=strand)) +
  #     geom_point(aes(color=strand), size=out_table$Hsp_align.len/62644)
  # ggsave(out_ggplot)
  # pdf(out_pdf)
  # # percent identity and length distribution and mean
  # hist(final_alis$Hsp_identity, breaks = 60, freq = F, xlim = c(0, 100), ylim = c(0, 0.1))
  # abline(v = mean(final_alis$Hsp_identity),
  #        lty = 1,
  #        lwd = 2)
  # abline(v = median(final_alis$Hsp_identity),
  #        lty = 2,
  #        lwd = 2)
  # legend(x = "topright", # location of legend within plot area
  #        c("Mean", "Median"),
  #        lty = c(1, 2),
  #        lwd = c(2, 2))
  #
  # hist(final_alis$Hsp_align.len, breaks = 100, freq = F, xlim = c(0, 15000), ylim = c(0,0.005))
  # abline(v = mean(final_alis$Hsp_align.len),
  #        lty = 1,
  #        lwd = 2)
  # abline(v = median(final_alis$Hsp_align.len),
  #        lty = 2,
  #        lwd = 2)
  # legend(x = "topright", # location of legend within plot area
  #        c("Mean", "Median"),
  #        lty = c(1, 2),
  #        lwd = c(2, 2))
  #
  # hist(final_alis$Hsp_evalue, breaks = 100, freq = F)
  # abline(v = mean(final_alis$Hsp_evalue),
  #        lty = 1,
  #        lwd = 2)
  # abline(v = median(final_alis$Hsp_evalue),
  #        lty = 2,
  #        lwd = 2)
  # legend(x = "topright", # location of legend within plot area
  #        c("Mean", "Median"),
  #        lty = c(1, 2),
  #        lwd = c(2, 2))
  # dev.off()
  return(stats)
  #})
}, mc.cores = threads)
print("finished, writing stats to file")

outstats <- rbindlist(outstats, fill = T)
# as.data.table(do.call(rbind, outstats))
write.table(outstats, file = outstats_name, append = F, sep = "\t")
end_total <- Sys.time()
print("total amount of time:")
print(end_total - start_total)
#

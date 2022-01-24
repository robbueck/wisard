
# function to combine two data frames with unequal set of columns. missing
# columns in one frame are filled with NA. empty data frames allowed
rbind.all.columns <- function(x, y) {
  if ((ncol(x) > 0) & (ncol(y) > 0)) {
    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))
    x[, c(as.character(y.diff))] <- NA
    y[, c(as.character(x.diff))] <- NA
    return(rbind(x, y))
  } else {
    return(rbind(x, y))
  }
}

#' Parse BLAST XML output
#'
#' @description Parses a BLAST xml-file to a
#'   \linkS4class{GRangesList} or \link[base]{data.frame} object.
#'
#' @param x A xml file created with NCBI-BLAST option -outfmt 5 or AB-BLAST
#'   option -mformat=7.
#' @param grange_output If set TRUE, output is a \linkS4class{GRangesList}
#'   object, otherwise a \link[base]{data.frame} object.
#' @param keep_sequences If FALSE, sequence and alignment columns are discarded.
#' @return A list with two elements accessible via: \cr \code{$metadata}:
#'   BLAST-metadata such as run parameters, sequence information and
#'   statistics.\cr \code{$alignments}: Information about the alignments.\cr If
#'   \code{grange_output} is set FALSE, output is a \link[base]{data.frame}
#'   object, else a \linkS4class{GRangesList} object with start, end, strand,
#'   sequence name and length taken from the subject. Information about the
#'   query is presered in the metadata columns.
#'
#' @details This function parses BLAST-xml output files using
#'   \link[XML]{xmlParse}. Setting \code{keep_sequences} FALSE will discard all
#'   sequence and alignment columns to save memory.
#'
#' @export
#'


read_blast_xml <- function(x, grange_output = T, keep_sequences = F) {

  Check <- checkmate::makeAssertCollection()
  checkmate::assertFileExists(x, add = Check)
  checkmate::assertLogical(grange_output, add = Check)
  checkmate::assertLogical(keep_sequences, add = Check)

  checkmate::reportAssertions(Check)

  xml_file <- xmlParse(x)

  if (xmlSize(xmlRoot(xml_file)[[9]]) > 1) {
    stop("more than one sequence in query, either change you blast search
         or this function")
  }

  # get metadata
  blast_stats <- data.frame(t(xmlSApply(getNodeSet(xml_file,
                                                   "//Statistics"
  )[[1]],
  function(p) xmlSApply(p, xmlValue))),
  row.names=NULL)
  blast_stats <- as.data.frame(lapply(blast_stats,
                                      function(x) as.numeric(as.character(x))))
  blast_stats$program <- xmlValue(getNodeSet(xml_file,
                                             "//BlastOutput_program")[[1]])
  blast_stats$query_name <- xmlValue(getNodeSet(xml_file,
                                                "//BlastOutput_query-def")[[1]])
  blast_stats$query_ID <- xmlValue(getNodeSet(xml_file,
                                              "//BlastOutput_query-ID")[[1]])
  blast_stats$query_len <- as.numeric(xmlValue(getNodeSet(xml_file,
                                                          "//BlastOutput_query-len")[[1]]))
  blast_params <- data.frame(t(xmlSApply(getNodeSet(xml_file,
                                                    "//Parameters")[[1]],
                                         function(p) xmlSApply(p, xmlValue))),
                             row.names=NULL)

  # convert numbers to numeric, different format for each blast-program
  if (blast_stats$program  == "blastn") {
    blast_params[,1:5] <- as.data.frame(lapply(blast_params[1:5],
                                               function(x) as.numeric(as.character(x))))
  } else if (blast_stats$program %in% c("tblastx", "blastp", "tblastn")) {
    blast_params[,2:4] <- as.data.frame(lapply(blast_params[2:4],
                                               function(x) as.numeric(as.character(x))))
  } else if (blast_stats$program == "blastp") {
    stop('invalid blasttype\n', metadata_frame$program)
  }
  metadata_frame <- merge(blast_stats, blast_params)
  # get rid of - for compatibility reasons
  colnames(metadata_frame) <- gsub(".text", "", colnames(metadata_frame))
  colnames(metadata_frame) <- gsub("-", ".", colnames(metadata_frame))


  # get all alignments for each sequence in the database
  xmlallhits <- getNodeSet(xml_file, "//Iteration_hits")[[1]]
  n_seqs <- xmlSize(xmlallhits)

  seq_frame <- data.frame()

  for (s in 1: n_seqs) {
    xmltop <- xmlallhits[[s]]
    seq_len <- as.numeric(xmlValue(getNodeSet(xmltop, "Hit_len")))
    seq_name <- xmlValue(getNodeSet(xmltop, "Hit_id"))
    try_out <- tryCatch(
      expr = {
        new_frame <- xmlToDataFrame(getNodeSet(xmltop, "Hit_hsps")[[1]])
        new_frame$seq_len <- seq_len
        new_frame$seq_ID <- seq_name
        if (!keep_sequences) {
          new_frame$Hsp_qseq <- NULL
          new_frame$Hsp_hseq <- NULL
          new_frame$Hsp_midline <- NULL
        }
        seq_frame <- rbind.all.columns(seq_frame, new_frame)
        T
      },
      error=function(cond) {
        return(F)
      }
    )
  }

  # if (!try_out) {
  #   print("No hits found in BLAST-result, returning empty data.frame")
  #   return(data.frame())
  # }
  # get rid of - for compatibility
  colnames(seq_frame) <- gsub("-", ".", colnames(seq_frame))
  seq_frame$query_id <- rep(blast_stats$query_ID, nrow(seq_frame))
  seq_frame$query_len <- rep(blast_stats$query_len, nrow(seq_frame))

  # check if all columns are there, else create missing ones
  required_colnames <- c('Hsp_num', 'Hsp_bit.score', 'Hsp_score',
                         'Hsp_evalue', 'Hsp_query.from', 'Hsp_query.to',
                         'query_id', 'query_len', 'Hsp_hit.from', 'Hsp_hit.to',
                         'Hsp_query.frame', 'Hsp_hit.frame', 'Hsp_identity',
                         'Hsp_positive', 'Hsp_gaps', 'Hsp_align.len',
                         'seq_len', 'seq_ID')

  missing_colnames <- setdiff(required_colnames, colnames(seq_frame))
  for (m in missing_colnames) {
    seq_frame[[m]] <- 0
  }
  seq_frame <- seq_frame[,c('Hsp_num', 'Hsp_bit.score', 'Hsp_score',
                            'Hsp_evalue', 'Hsp_query.from', 'Hsp_query.to',
                            'Hsp_hit.from', 'Hsp_hit.from', 'Hsp_hit.to',
                            'Hsp_query.frame', 'Hsp_hit.frame', 'Hsp_identity',
                            'Hsp_positive', 'Hsp_gaps', 'Hsp_align.len',
                            'seq_len', 'query_len', 'seq_ID',  'query_id')]

  seq_frame[,1:17] <- as.data.frame(lapply(seq_frame[1:17],
                                           function(x) as.numeric(as.character(x))))

  # check that start is always smaller than end
  change_sstrand <- which(seq_frame$Hsp_hit.from > seq_frame$Hsp_hit.to)
  if (length(change_sstrand) > 0) {
    # print("changing sstart")
    seq_frame$X <- seq_frame$Hsp_hit.from
    seq_frame[change_sstrand,]$Hsp_hit.from <- seq_frame[change_sstrand,]$Hsp_hit.to
    seq_frame[change_sstrand,]$Hsp_hit.to <- seq_frame[change_sstrand,]$X
    seq_frame$X <- NULL
  }
  change_qstrand <- which(seq_frame$Hsp_query.from > seq_frame$Hsp_query.to)
  if (length(change_qstrand) > 0) {
    # print("changing qstart")
    seq_frame$X <- seq_frame$Hsp_query.from
    seq_frame[change_qstrand,]$Hsp_query.from <- seq_frame[change_qstrand,]$Hsp_query.to
    seq_frame[change_qstrand,]$Hsp_query.to <- seq_frame[change_qstrand,]$X
    seq_frame$X <- NULL
  }


  seq_frame$sstrand <- "+"
  # check which hits are reverse complement on one genome. for blastn,
  # mark always subject strand as "-" for WIS
  if (metadata_frame$program == "blastn") {
    print("file is from blastn")
    minus_strand <- which(seq_frame$Hsp_query.frame < 0)
    # then from ab-blast, strand info in on query frame field
    if (length(minus_strand) > 0) {
      seq_frame[minus_strand,]$sstrand <- "-"
      # to avoid later confution, move information to hit-strand
      seq_frame$Hsp_hit.frame <- seq_frame$Hsp_query.frame
      seq_frame$Hsp_query.frame <- 1
    } else if (length(minus_strand) < 1) {
      # then possibly from ncbi-blast, strand info in hit/subject frame field
      minus_strand <- which(seq_frame$Hsp_hit.frame < 0)
      if (length(minus_strand) > 0) {
        # if not, everything seems to be on one strand
        seq_frame[minus_strand,]$sstrand <- "-"
      }
    }
    # for tblastx, it matters on which genome the minus strand is aligned
  } else if (metadata_frame$program == "tblastx") {
    print("file is from tblastx")
    minus_sstrand <- which(seq_frame$Hsp_hit.frame < 0)
    seq_frame[minus_sstrand,]$sstrand <- "-"
  } else if ( metadata_frame$program == "blastp") {
    print("file is from blastp")
    wrong_sstrand <- which(seq_frame$Hsp_hit.frame != 0)
    wrong_qstrand <- which(seq_frame$Hsp_query.frame != 0)
    if (length(wrong_qstrand) > 0 | length(wrong_sstrand) > 0) {
      stop("File from ", metadata_frame$program, "contains invalid values
           for strand, should be 0")
    }

    # no strandedness expected here, nothing to do
    } else if (metadata_frame$program == "blastx") {
      print("file is from blastx")
      wrong_sstrand <- which(seq_frame$Hsp_hit.frame != 0)
      if (length(wrong_sstrand) > 0) {
        stop("File from ", metadata_frame$program, "contains invalid values
             for hit strand, should be 0")
      }
      # hits are proteins, so always of frame 0, nothing to do either for now
      } else if (metadata_frame$program == "tblastn") {
        print("file is from tblastn")
        wrong_qstrand <- which(seq_frame$Hsp_query.frame != 0)
        if (length(wrong_qstrand) > 0) {
          stop("File from ", metadata_frame$program, "contains invalid values
               for query strand, should be 0")
        }
        minus_sstrand <- which(seq_frame$Hsp_hit.frame < 0)
        seq_frame[minus_sstrand,]$sstrand <- "-"
        } else {
          stop('invalid blasttype\n', metadata_frame$program)
        }

  seq_frame$raw_score <- seq_frame$Hsp_score
  seq_frame$ID <- 1: nrow(seq_frame)
  seq_frame$filter_pass <- TRUE                       # create filter column

  #create granges output
  if (grange_output) {
    seq_granges <- makeGRangesFromDataFrame(seq_frame, keep.extra.columns = T,
                                            seqnames.field = 'seq_ID',
                                            start.field = 'Hsp_hit.from',
                                            end.field = 'Hsp_hit.to',
                                            strand.field = 'sstrand')
    # create metadata column with granges object for query intervals
    query_ranges <- makeGRangesFromDataFrame(seq_frame, keep.extra.columns = F,
                                             seqnames.field = 'query_id',
                                             start.field = 'Hsp_query.from',
                                             end.field = 'Hsp_query.to',
                                             ignore.strand = T)
    names(seq_granges) <- seq_frame$ID
    names(query_ranges) <- seq_frame$ID
    # print(unique(seq_frame$seq_len))

    seqlengths(seq_granges) <- unique(seq_frame$seq_len)
    seqlengths(query_ranges) <- unique(seq_frame$query_len)
    seq_granges$query_ranges <- query_ranges
    out_data <- list(alignments = seq_granges, metadata = metadata_frame)
  } else {
    out_data <- list(alignments = seq_frame, metadata = metadata_frame)
  }
  return(out_data)
  }


#' Simple BLAST-results filter
#'
#' @description Filter BLAST-results by different criteria. Removes every
#'   alignment not fulfilling the criteria
#'
#' @param raw_granges A \linkS4class{GRangesList} as
#'   returned from [read_blast_xml].
#' @param grange_output If set TRUE, output it a
#'   \linkS4class{GRangesList} object, otherwise a
#'   \link[base]{data.frame} object.
#' @param keep_sequences If FALSE, sequence and alignment columns are discarded.
#' @param min_len,max_e,min_identity Integer or float to filter the
#'   BLAST-results by minimum length, maximum e-value or minimum percent idenity
#'
#' @export
#'
filter_blast <- function(raw_granges, min_len = 0, max_e = 100,
                         min_identity = 0) {
  # filter alignments according to filter criteria
  filter_len <- which ((raw_granges$Hsp_align.len < min_len))
  raw_granges[filter_len]$filter_pass <- rep(FALSE, length(filter_len))
  filter_e <- which ((raw_granges$Hsp_evalue > max_e))
  raw_granges[filter_e]$filter_pass <- rep(FALSE, length(filter_e))
  filter_identity <- which ((raw_granges$Hsp_identity < min_identity))
  raw_granges[filter_identity]$filter_pass <- rep(FALSE, length(filter_identity))

  print("Number of alignments filtered out")
  n_filter_bad <- length(which(raw_granges$filter_pass == F))
  print(n_filter_bad)
  print("Number of alignments that passed filter")
  n_filter_good <- length(raw_granges) - n_filter_bad
  print(n_filter_good)
  raw_granges <- raw_granges[which(raw_granges$filter_pass == TRUE)]
  if (n_filter_good < 1) {
    warning("No alignemnts that match filter criteria.")
    return(raw_granges[0])
  }
  return(raw_granges)
}

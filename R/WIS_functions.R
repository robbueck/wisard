#' Log factorial
#'
#' @description Calculates the natural logarithm of a factorial for high numbers
#'
#' @param n Integer or vector of integers

#' @return log(n!) as float. If input is a vector of integers, output is a
#' vecor with log(n!) for each value in n
#'
#' @export
#'
log_factorial <- function(n) {
  if (length(n) == 1) {
    if (n%%1 != 0) {
      stop("n has to be an integer")
    } else if (n < 40) {
      return(log(factorial(n)))
    } else {
      return(n*log(n)-n)
    }
  } else if (length(n) > 1) {
    n_logfact <- sapply(n, log_factorial)
    return(n_logfact)
  } else {
    stop("function log_factorial not supported for", n)
  }
}



#' BLAST sum-score
#'
#' @description Caclulated BLAST sum-score from raw-scores and BLAST search
#' statistics
#'
#' @param b Vector of integers or floats containing the raw scores.
#' @param metadata BLAST metadata as returned by [read_blast_xml].
#' A \link[base]{data.frame} object with one row containing the columns
#' Statistics_lambda, Statistics_kappa, query_len and Statistics_db.le
#' @return Sum score as float.
#'
#' @details This function calculates the BLAST sum-score as described  in
#' Equation 4-17 here:
#' \url{http://etutorials.org/Misc/blast/Part+II+Theory/Chapter+4.+Sequence+Similarity/4.7+Sum+Statistics+and+Sum+Scores/}
#'
#' @export
#'
sum_score <- function(b, metadata) {
  l <- metadata$Statistics_lambda
  lnkmn <- log(metadata$Statistics_kappa * metadata$query_len *
                 metadata$Statistics_db.len)
  r <- length(b)
  # for testing, keep the approximation,
  # for final version, take the original one
  #return(l * sum(b) - (r * lnkmn) + log_factorial(r))
  return(sum(b * l - lnkmn - 1) + (r * log(r)))
}


#' WIS-algorithm for a set of intervals
#'
#' @description Performs the WIS algorithm on a set of weighted intervals given
#'   as a \linkS4class{GRangesList} object. Better to use wrapper function [get_wis]
#'   to run it.
#'
#' @param granges_list A set of intervals as a \linkS4class{GRangesList} object.
#'   Requires an interval-score to maximize over as a column specified in
#'   \code{raw_score}.
#' @param overlap A non-negative integer. The allowed overlap between two
#'   intervals.
#' @param max_score Name of a column in \code{granges_list} containing integers or
#'   floats to be used as an interval-score to maximize over.
#' @param score_funct Function to update the score of a set each time an element
#'   is added. Can take up to three arguments: x (set-score of set so far) y
#'   (interval-score of added item) and n (number of items in set including the
#'   one added). Returnes set-score, default: x + y.
#' @param max_funct Function to calculate a set-score to maximize over. This
#'   set-score is not used to calculate itself, but is calculated using the
#'   set-score from \code{score_funct}. Takes up to three arguments: x (set-score of
#'   set so far), y (interval-score of added item), n (number of items in set
#'   including the one added). Returns set-score, default: x + y. For initial
#'   calculation y = 0 and n = 2. If in doubt about its use, set equal to
#'   \code{score_funct}.
#' @return A \linkS4class{GRangesList} object with intervals compatible
#'   according to the specified options.
#'
#' @details This function calculates a set of non-overlapping intervals while
#'   maximizing over a given set-score. The behaviour of this set-score when an
#'   element is added to a list is specified by the two input parameters
#'   \code{score_funct} and \code{max_funct}. For a simple sum or product set-score, set both
#'   functions to x + y or x * y.
#'
#'   The \code{max_score} option can be used, if a set-score contains non linear
#'   elements, such as log(n), e.g. in the BLAST scum-score:
#'
#'   \code{score_funct = function (x, y, n) (x + y)}
#'
#'   \code{max_funct = function (x, y, n) (x + y) + log(n!)}
#'
#'   If n is included in one of the formulars in a non-linear way, it is not
#'   guaranteed anymore that the final set is the best possible set.
#'
#'   Strand information is ignored here, has to be taken into account before
#'   calling this fucntion.
#'
#' @export


 # @param second_range Optional, name of a column in \code{granges_list} containing
 #  another \linkS4class{GRangesList} object. If set, for intervals in
 #  granges_list to be compatible, they also have to be compatible in
 #  \code{second_range}. (better remove this and solve problem differently, no optimal
 #  solution here).
get_wis <- function(granges_list, overlap=0, max_score = "score",
                    score_funct = function(x, y, n) x + y) {
  # check arguments
  Check <- checkmate::makeAssertCollection()
  checkmate::assertClass(granges_list, classes = "GRanges", add = Check)
  checkmate::assertInt(overlap, na.ok = F, lower = 0, add = Check)
  checkmate::assertSubset(max_score,
               choices = colnames(mcols(granges_list)),
               empty.ok = F,
               add = Check)
  checkmate::assertFunction(score_funct, nargs = 3, add = Check)
  # checkmate::assertFunction(max_funct, nargs = 3, add = Check)

  checkmate::reportAssertions(Check)
  # cat("get_wis was called with following parameters:\noverlap\t", overlap, "\nmax_score\t", max_score, "\n")


  #sort everything by end position and score
  granges_list <- sort(granges_list, by = ~end + -eval(as.name(max_score)))

  # check if granges_list contains more than one element,
  # otherwise return original list
  if (length(granges_list) > 1) {
    granges_list$ID <- 1:length(granges_list)        # prevent things from getting mixed up
    n_alis <- length(granges_list)
    wis <- data.frame(matrix(ncol = 3, nrow =n_alis))
    colnames(wis) <- c("score", "r", "maxscore")

    # initialisation
    # wis$r <- 2                                            # number of alignments in WIS currently ending at the specific alignment if another one is added
    wis$score <- mcols(granges_list)[[max_score]]
    # wis$maxscore <- max_funct( x = wis$score, y = 0, n = wis$r )

    # neg <- which(wis$maxsumscore < 0)
    P <- rep(0, n_alis)                                   # list with pointers for each alignment to the previous one in its max-score set
    for ( i in 2: n_alis ) {
      # get all compatible previous intervals
      no_overlaps <- which(countOverlaps(granges_list[1:(i-1)],
                                         granges_list[i],
                                         ignore.strand = T,
                                         minoverlap = (overlap + 1 )) == 0)

      if ( length(no_overlaps) > 0) {
        # if no non-overlapping range before, this is the first in set, just
        # continue
        # find maximum compatible set and store pointer
        j_alis <- granges_list[no_overlaps,]
        j_scores <- wis[no_overlaps,]$score
        max_ID <- j_alis[which.max(j_scores)]$ID

        P[i] <- max_ID

        # update r and score
        # wis[i,]$r <- (wis[max_ID,]$r + 1 )

        # wis[i,]$maxscore <- max_funct( x = wis[max_ID,]$score,
        #                                y = wis[i,]$score,
        #                                n = wis[i,]$r )
        wis[i,]$score <- score_funct( x = wis[max_ID,]$score,
                                      y = wis[i,]$score,
                                      n = wis[i,]$r )
      }
    }
    wis$r <- wis$r - 1                            # nothing more is added

    # wis$maxscore <- max_funct( x = wis$score, y = 0, n = wis$r )

    # backtracing
    K <- which.max(wis$score)
    cat("max score is: ", max(wis$score), "\n")
    final_set_ids <- c()
    while ( K > 0 ) {
      final_set_ids <- c(final_set_ids, K)
      K = P[K]
    }

    return(list( alignments = granges_list[final_set_ids], max_score = max(wis$score) ))
  } else if (length(granges_list) == 1) {
    # input only contains one interval
    return(list( alignments = granges_list,
                 max_score = score_funct(x = mcols(granges_list)[[max_score]][1],
                                       y = 0,
                                       n = 1) ))
  } else {
    # input is empty
    warning("Input does not contain any intervals, returning empty set")
    return(list( alignments = granges_list, max_score = -Inf ))
  }
}



#' Get best start for WIS
#'
#' @description For a circular scale, this function finds the best interval(s)
#'   to start the WIS alsogrithm with.
#'
#' @param granges_list A \linkS4class{GRangesList} object
#'   containing the intervals.
#' @param overhang  A non-negative integer specifying for a linearized circular
#'   scale how much of the end is overlapping with the start.
#' @return Returns the index of the optimal interval to start a WIS-algorithm with. If no
#'   such interval exists, returns an integer vector with indices for the best
#'   possible interval and all intervalls overlapping its starting position.
#' @details Takes a circular  \linkS4class{GRangesList}
#'   object and finds the best intervals to start with. Ideally, this is an
#'   inverval, whose start position is not overlapped by any other interval. If
#'   no such interval exists, the interval with the least overlaps at its start
#'   position and each interval overlapping this position are returned.
#'
#'   Strand information is ignored here, has to be taken into account before
#'   calling this fucntion.
#'
#'   Also checks if overhang is long enough or should be extended, e.g. in
#'   BLAST-searches by copying a longer sequence from the beginning to the end.
#'
#' @export

get_best_start <- function(granges_list, overhang = 0) {
  Check <- checkmate::makeAssertCollection()
  checkmate::assertClass(granges_list, classes = "GRanges", add = Check)
  checkmate::assertInt(overhang, na.ok = F, lower = 0, add = Check)
  checkmate::assertLogical(isCircular(granges_list), add = Check)

  checkmate::reportAssertions(Check)
  if (length(granges_list) < 1) {
    print("warning, GRanges object empty, return 0")
    return(0)
  }

  new_start <- c(0)

  # check if pos 1 of the genome is covered by anything
  circ_test <- GRanges(seqnames = seqnames(granges_list)[1],
                       ranges = IRanges(1, end = 1))
  olap_start <- which(countOverlaps(granges_list, circ_test,
                                    ignore.strand = T) > 0)
  # if pos 1 is covered by something
  if (length(olap_start) > 0){
    # as soon as it's 1, the first alignemt could be overlapping with the last one, keep it 0!
    cat("rearranging genome\n")
    # count overlaps at start postitions
    granges_list_olaps <- granges_list
    end(granges_list_olaps) <- start(granges_list)          # granges object only covering the start positions
    olaps <- countOverlaps(granges_list_olaps, granges_list, ignore.strand = T)
    # get position of best starting point
    possible_starts <- c(which.min(olaps))
    new_start_pos <- start(granges_list)[possible_starts]
    olap_test <- GRanges(seqnames = seqnames(granges_list)[1],
                         ranges = IRanges(new_start_pos,
                                          end = new_start_pos))

    # get indices for everything overlapping the start point of optimal start
    new_start <- which(countOverlaps(granges_list,
                                     olap_test,
                                     ignore.strand = T) > 0)

    # check if any alignments overlap position 1 and are longer than overhang
    # they might not be good alignments, overhang should be increased then
    problematic <- which(
      end(granges_list[olap_start]) == seqlengths(granges_list) + overhang)
    if (length(problematic) > 0) {
      better_len <- max(end(granges_list[which(start(granges_list) == 1)]))
      cat("Following alignments span more than the ", overhang,
          " added bases at the end.\nBetter values could be at least ",
          better_len, "\n")
    }
  } else {
    # best start is first alignment
    new_start <- c(1)
  }
  # check output
  if(length(new_start) < 1 || new_start == 0) {
    stop("Granges object not empty, but no start position found,
         something went wrong")
  }
  return(new_start)
}




#' Run get_wis function
#'
#' @description Runs the [get_wis] function on a set of weighted intervals given
#'   as a \linkS4class{GRangesList} object by taking into
#'   account several parameters.
#'
#' @param granges_list A set of intervals as a
#'   \linkS4class{GRangesList} object. Requires a
#'   interval-score to maximize over as a column specified in \code{score_col}.
#' @param is_circular Boolean, set TRUE if scale is circular.
#' @param overhang A non-negative integer. A non-negative integer specifying for
#'   a linearized circular scale how much of the end is overlapping with the
#'   start. Ignored if is_circular is FALSE.
#' @param use_strand Boolean, if TRUE, overlapping intervals on different strand
#'   are considered compatible.
#' @param score_col Name of a column in \code{granges_list} containing integers or
#'   floats to be used as a score to maximize over.
#' @param score_f Function to update the score of a set each time an element is
#'   added. Passed to \code{score_funct} in [get_wis], see there for detailed
#'   description.
#' @param max_f Function to calculate a set-score to maximize over. Passed to
#'   \code{max_funct} in [get_wis], see there for detailed description.
#' @param overlap A non-negative integer. The allowed overlap between two
#'   intervals.
#' @param frame_col Optional, column in \code{granges_list} containing frame information. If
#'   set, overlapping intervals on different strand are considered compatible.
#' @param second_range Optional, name of a column in granges_list containing
#'   another \linkS4class{GRangesList} object. If set, for intervals in
#'   granges_list to be compatible, they also have to be compatible in
#'   second_range.

#'
#' @return A \linkS4class{GRangesList} object with intervals
#'   compatible according to the specified options.
#'
#' @details This function runs the [get_wis] function on a set of intervals and allows to distinguish between intervals on different stradns or frames. For each of these settings,[get_wis] is run seperately. For detailed description of the output, see [get_wis].
#'
#' @export
#'

# @param second_range Optional, name of a column in \code{granges_list} containing
  # another \linkS4class{GRanges} object. If set, for
  # intervals in granges_list to be compatible, they also have to be compatible
  # in second_range. (take this out of here, does not really work like this)

run_get_wis <- function(granges_list, is_circular, overhang = 0, use_strand = F,
                        frame_col = NULL, score_col = "score",
                        score_f = function(x, y, n) x + y,
                        # max_f = function(x, y, n) x + y,
                        overlap = 0,
                        second_range = NULL)
{
  Check <- checkmate::makeAssertCollection()
  checkmate::assertClass(granges_list, classes = "GRanges", add = Check)
  checkmate::assertInt(length(granges_list), na.ok = F, lower = 1, add = Check)
  checkmate::assertLogical(is_circular, add = Check)
  checkmate::assertInt(overhang, na.ok = F, lower = 0, add = Check)
  checkmate::assertInt(overlap, na.ok = F, lower = 0, add = Check)
  checkmate::assertLogical(use_strand, add = Check)
  checkmate::assertSubset(frame_col, choices = colnames(mcols(granges_list)), empty.ok = T,
               add = Check)
  checkmate::assertFunction(score_f, nargs = 3, add = Check)
  # checkmate::assertFunction(max_f, nargs = 3, add = Check)

  checkmate::reportAssertions(Check)
  # cat("run_get_wis was called with following parameters:\nis_circular\t", is_circular, "\n overhang\t", overhang,
  #     "\nuse_strand\t", use_strand, "\nframe_col\t", frame_col, "\nscore_col\t", score_col, "\noverlap\t", overlap, "\n")
  start_l <- Sys.time()

  if(length(isCircular(granges_list)) > 1) {
    warning("BLAST-Hits consist of more than one sequence,behaviour
            of algorithm is uncertain in that case")
  }
  isCircular(granges_list) <- rep(is_circular, length(isCircular(granges_list)))

  # remove all alignments in overhang
  if (is_circular) {
    # print(seqlengths(granges_list))
    # not here, should be done before
    # seqlengths(granges_list) <- seqlengths(granges_list) - overhang
    granges_list <- granges_list[start(granges_list) < seqlengths(granges_list)]

  }



  if (is.null(frame_col)) {
    # not use frame
    ###############
    if (use_strand) {
      # use strand but not frame
      ##########################
      print("strand yes, frame no")
      # split based on strand to analyze each strand specifically
      # otherwise output is not right anymore
      granges_plus <- granges_list[strand(granges_list) == "+"]
      granges_minus <- granges_list[strand(granges_list) == "-"]

      # get start positions to calculate WIS for
      if (is_circular) {
        start_pos_plus <- get_best_start(granges_plus, overhang = overhang)
        start_pos_minus <- get_best_start(granges_minus, overhang = overhang)
      } else {
        start_pos_plus <- 1
        start_pos_minus <- 1
      }

      # run get_wis for all starting positions and find the best result
      # what about *?

      # for plus strand
      if (start_pos_plus[1] == 0) {
        print("no intervals on plus strand, continue anyway")
        best_data_plus <- granges_plus
      } else {
        print("analyzing + strand")
        score <- -Inf
        for (s in start_pos_plus) {
          cat("calculate wis for start alignment ", s, "\n")
          # rearrange for current start
          granges_plus<- c(granges_plus[s:length(granges_plus)], granges_plus[0:(s-1)])
          from_get_wis <- get_wis(granges_plus,
                                  max_score = score_col,
                                  score_funct = score_f,
                                  # max_funct = max_f,
                                  overlap = overlap)
          new_data_plus <- from_get_wis$alignments
          sum_score <- from_get_wis$max_score
          if (sum_score > score) {
            score <- sum_score
            best_data_plus <- new_data_plus
          }
        }
      }

      # for minus strand
      if (start_pos_minus[1] == 0) {
        print("no intervals on minus strand, continue anyway")
        best_data_minus <- granges_minus
      } else {
        print("analyzing - strand")
        score <- -Inf
        for (s in start_pos_minus) {
          cat("calculate wis for start alignment ", s, "\n")
          # rearrange for current start
          granges_minus<- c(granges_minus[s:length(granges_minus)], granges_minus[0:(s-1)])
          from_get_wis <- get_wis(granges_minus,
                                  max_score = score_col,
                                  score_funct = score_f,
                                  # max_funct = max_f,
                                  overlap = overlap)
          new_data_minus <- from_get_wis$alignments
          sum_score <- from_get_wis$max_score
          if (sum_score > score) {
            score <- sum_score
            best_data_minus <- new_data_minus
          }
        }
      }
      final_data <- c(best_data_minus, best_data_plus)


    } else {
      # use neither strand nor frame
      ##############################
      print("strand no, frame no")
      if(is_circular) {
        start_pos <- get_best_start(granges_list, overhang = overhang)
      } else {
        start_pos <- 1
      }

      if (start_pos[1] == 0) {
        stop("no intervals found to compute WIS")
      } else {
        final_data <- list()
        score <- -Inf
        for (s in start_pos) {
          # rearrange for current start
          cat("calculate wis for start pos: ", s, "\n")
          granges_list <- c(granges_list[s:length(granges_list)], granges_list[0:(s-1)])
          from_get_wis <- get_wis(granges_list,
                                  max_score = score_col,
                                  score_funct = score_f,
                                  # max_funct = max_f,
                                  overlap = overlap)
          new_data <- from_get_wis$alignments
          sum_score <- from_get_wis$max_score
          if (sum_score > score) {
            score <- sum_score
            final_data <- new_data
          }
        }
      }
    }
  } else {
    # use frame
    ###########
    print("frame yes")
    # if (!is.null(second_range)) {
    #   warning("Arguments second_range and frame_col are not compatible,
    #           at least one of them has to be NULL.\nArgument second_range
    #           will be ignored")
    # }
    if (use_strand) {
      frames <- levels(factor(mcols(granges_list)[[frame_col]]))
    } else {
      # here strand information is lost (but should be kept in the strand field
      # and be able to get restored?)
      frames <- levels(factor(abs(mcols(granges_list)[[frame_col]])))
    }
    for ( f in as.numeric(frames) ) {
      cat("calculating WIS for frame ", f, "\n")
      if (use_strand) {
        current_frame_granges <- granges_list[mcols(granges_list)[[frame_col]] == f]
      } else {
        current_frame_granges <- granges_list[abs(mcols(granges_list)[[frame_col]]) == f]
      }
      if ( is_circular ) {
        start_pos <- get_best_start(current_frame_granges, overhang = overhang)
      } else { start_pos <- 1 }
      if ( start_pos[1] == 0 ) {
        cat("no intervals in frame ", f, ", continue anyway\n")
        next
      } else {
        score <- -Inf
        for ( s in start_pos ) {
          cat("calculate wis for start alignment ", s, "\n")
          # rearrange for current start
          current_frame_granges <- c(current_frame_granges[s:length(current_frame_granges)],
                                  current_frame_granges[0:(s-1)])
          from_get_wis <- get_wis(current_frame_granges,
                                  max_score = score_col,
                                  score_funct = score_f,
                                  # max_funct = max_f,
                                  overlap = overlap)
          new_data <- from_get_wis$alignments
          sum_score <- from_get_wis$max_score
          if( sum_score > score ) {
            score <- sum_score
            best_data <- new_data
          }
        }
      }
      # collect everything from all frames
      if (!exists("final_data")) {
        final_data <- best_data
      } else {
        final_data <- c(final_data, best_data)
      }
    }
    }
  if(!exists("final_data")) {
    warning("No intervals found to return")
  }
  return(final_data)

}



# run second range function, switches the main range with a second indicated
# range, runs get_wis() and switches the ranges back
two_ranges_run_get_wis <- function(granges_list, is_circular = F, overhang = 0, use_strand = F,
                             frame_col = NULL, score_col = "score",
                             score_f = function(x, y, n) x + y,
                             overlap = 0,
                             second_range = NULL) {
  print("running get wis on first range")
  granges_list <- run_get_wis(granges_list, is_circular = is_circular, overhang = overhang, use_strand = use_strand,
                              frame_col = frame_col, score_col = score_col,
                              score_f = score_f,
                              overlap = overlap,
                              second_range = second_range)
  if (!is.null(second_range)) {
    print("running get wis for given second range: ")
    print(second_range)
    granges_list$secondID <- 1:length(granges_list)
    second_list <- mcols(granges_list)[[second_range]]
    mcols(second_list) <- mcols(granges_list)
    second_list <- run_get_wis(second_list, is_circular = is_circular, overhang = overhang, use_strand = use_strand,
                               frame_col = frame_col, score_col = score_col,
                               score_f = score_f,
                               overlap = overlap,
                               second_range = second_range)
    granges_list <- granges_list[second_list$secondID]
    granges_list$secondID <- NULL
  }
  return(granges_list)
}

# run second range function, switches the main range with a second indicated
# range, runs get_wis() and switches the ranges back
run_get_wis_two_ranges <- function(granges_list, is_circular = F, overhang = 0, use_strand = F,
                                   frame_col = NULL, score_col = "score",
                                   score_f = function(x, y, n) x + y,
                                   overlap = 0,
                                   second_range = NULL) {
  print("running get wis on first range")
  if (!is.null(second_range)) {
    print("running get wis for given second range: ")
    print(second_range)
    granges_list$secondID <- 1:length(granges_list)
    second_list <- mcols(granges_list)[[second_range]]
    mcols(second_list) <- mcols(granges_list)
    second_list <- run_get_wis(second_list, is_circular = is_circular, overhang = overhang, use_strand = use_strand,
                               frame_col = frame_col, score_col = score_col,
                               score_f = score_f,
                               overlap = overlap,
                               second_range = second_range)
    granges_list <- granges_list[second_list$secondID]
    granges_list$secondID <- NULL
  }
  granges_list <- run_get_wis(granges_list, is_circular = is_circular, overhang = overhang, use_strand = use_strand,
                              frame_col = frame_col, score_col = score_col,
                              score_f = score_f,
                              overlap = overlap,
                              second_range = second_range)
  return(granges_list)
}

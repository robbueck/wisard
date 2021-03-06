% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WIS_functions.R
\name{get_wis}
\alias{get_wis}
\title{WIS-algorithm for a set of intervals}
\usage{
get_wis(
  granges_list,
  overlap = 0,
  max_score = "score",
  score_funct = function(x, y, n) x + y,
  max_funct = function(x, y, n) x + y
)
}
\arguments{
\item{granges_list}{A set of intervals as a \linkS4class{GRangesList} object.
Requires an interval-score to maximize over as a column specified in
\code{raw_score}.}

\item{overlap}{A non-negative integer. The allowed overlap between two
intervals.}

\item{max_score}{Name of a column in \code{granges_list} containing integers or
floats to be used as an interval-score to maximize over.}

\item{score_funct}{Function to update the score of a set each time an element
is added. Can take up to three arguments: x (set-score of set so far) y
(interval-score of added item) and n (number of items in set including the
one added). Returnes set-score, default: x + y.}

\item{max_funct}{Function to calculate a set-score to maximize over. This
set-score is not used to calculate itself, but is calculated using the
set-score from \code{score_funct}. Takes up to three arguments: x (set-score of
set so far), y (interval-score of added item), n (number of items in set
including the one added). Returns set-score, default: x + y. For initial
calculation y = 0 and n = 2. If in doubt about its use, set equal to
\code{score_funct}.}
}
\value{
A \linkS4class{GRangesList} object with intervals compatible
according to the specified options.
}
\description{
Performs the WIS algorithm on a set of weighted intervals given
as a \linkS4class{GRangesList} object.
}
\details{
This function calculates a set of non-overlapping intervals while
maximizing over a given set-score. The behaviour of this set-score when an
element is added to a list is specified by the two input parameters
\code{score_funct} and \code{max_funct}. For a simple sum or product set-score, set both
functions to x + y or x * y.

The \code{max_score} option can be used, if a set-score contains non linear
elements, such as log(n), e.g. in the BLAST scum-score:

\code{score_funct = function (x, y, n) (x + y)}

\code{max_funct = function (x, y, n) (x + y) + log(n!)}

If n is included in one of the formulars in a non-linear way, it is not
guaranteed anymore that the final set is the best possible set.

Strand information is ignored here, has to be taken into account before
calling this fucntion.
}

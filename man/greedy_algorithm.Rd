% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedy.R
\name{greedy_algorithm}
\alias{greedy_algorithm}
\title{Greedy algorithm for a set of overlapping intervals.}
\usage{
greedy_algorithm(
  granges_list,
  max_score = "raw_score",
  overlap = 0,
  is_circular = F,
  overhang = 0,
  use_strand = T,
  second_range = NULL
)
}
\arguments{
\item{granges_list}{A set of intervals as a \linkS4class{GRangesList} object.
Requires a score column to be specified in \code{max_score}.}

\item{max_score}{Name of a column in \code{granges_list} containing integers
or floats to be used as a score to maximize over.}

\item{overlap}{A non-negative integer specifying the allowed overlap between
two intervals.}

\item{is_circular}{Boolean. If set true, intervals are assumed to be on a
circular scale.}

\item{overhang}{A non-negative integer how much of the end is overlapping
with the start. specifying for a linearized circular scale  Ignored if
is_circular is FALSE.}

\item{use_strand}{Boolean. If TRUE, strand information is used to detect
overlaps}

\item{second_range}{Optional, name of a column in granges_list containing
another \linkS4class{GRangesList} object. If set, for intervals in
granges_list to be compatible, they also have to be compatible in
second_range.}
}
\value{
A \linkS4class{GRangesList} object with intervals compatible
according to the specified options.
}
\description{
A function to create a set of non-overlapping intervals using a
greedy algorithm by prefering intervals with higeher scores.
}
\details{
The algorithm seeks to maximize the sum of a score given in
max_score, but does not necessarily return the best possible solution.
Created mainly for comparison with the WIS-algorithm function \link{get_wis}.

When is_circular is TRUE, intervals exceeding the end of the scale given by
(seqlength(granges_list) - overhang) are wrapped to the beginning.
}
\seealso{
\url{https://doi.org/10.1093/bioinformatics/bth324} for description
of the algorithm
}

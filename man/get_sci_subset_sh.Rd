% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sci_unscaled.R
\name{get_sci_subset_sh}
\alias{get_sci_subset_sh}
\title{Obtain SCI subsets for the animal's year of life}
\usage{
get_sci_subset_sh(
  df,
  members_l,
  focals_l,
  females_l,
  interactions_l,
  min_res_days,
  directional,
  legacy_sci
)
}
\arguments{
\item{df}{One row individual-year-of-life data}

\item{members_l}{A subset of members table produced by the function 'subset_members'}

\item{focals_l}{A subset of focals produced by the function 'subset_focals'}

\item{females_l}{A subset of female counts produced by the function 'subset_females'}

\item{interactions_l}{A subset of interactions data produced by the function 'subset_interactions'}

\item{min_res_days}{The minimum number of residence days needed to be included. Defaults to 60 days.}

\item{directional}{Logical value indicating whether to preserve directionality}
}
\value{
The input row with an additional list column containing the subset
}
\description{
Obtain SCI subsets for the animal's year of life
}

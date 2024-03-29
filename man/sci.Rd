% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sci_unscaled.R
\name{sci}
\alias{sci}
\title{Calculate Social Connectedness Index variables from individual-year-of-life data}
\usage{
sci(
  my_iyol,
  members_l,
  focals_l,
  females_l,
  interactions_l,
  min_res_days = 60,
  parallel = FALSE,
  ncores = NULL,
  directional = FALSE,
  legacy_sci = FALSE
)

sci(
  my_iyol,
  members_l,
  focals_l,
  females_l,
  interactions_l,
  min_res_days = 60,
  parallel = FALSE,
  ncores = NULL,
  directional = FALSE,
  legacy_sci = FALSE
)
}
\arguments{
\item{my_iyol}{Individual-year-of-life data.}

\item{members_l}{A subset of members table produced by the function 'subset_members'}

\item{focals_l}{A subset of focals produced by the function 'subset_focals'}

\item{females_l}{A subset of female counts produced by the function 'subset_females'}

\item{interactions_l}{A subset of interactions data produced by the function 'subset_interactions'}

\item{min_res_days}{The minimum number of coresidence days needed for dyad to be included. Defaults to 60 days.}

\item{parallel}{Logical value indicating whether to process in parallel}

\item{ncores}{Integer value indicating how many cores to use in parallel processing}

\item{directional}{Logical value indicating whether to preserve directionality}
}
\value{
The input data with an additional list columsn containing the full subset and variables.

The input data with an additional list columsn containing the full subset and variables.
}
\description{
Calculate Social Connectedness Index variables from individual-year-of-life data

Calculate Social Connectedness Index variables from individual-year-of-life data
}

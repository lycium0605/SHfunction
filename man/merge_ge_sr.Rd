% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SR-related.R
\name{merge_ge_sr}
\alias{merge_ge_sr}
\title{merge_ge_sr}
\usage{
merge_ge_sr(ge_meta, sci, agi, dsi, dsi_sum, type = "wt_avg")
}
\arguments{
\item{ge_meta}{should contain sname,dart_date,treatment,group,sex,rank,age}

\item{sci}{sci input from Ramboseli::sci}

\item{agi}{agi input from Ramboseli::agi}

\item{dsi}{dsi_pop input from Ramboseli::dsi_pop}

\item{dsi_sum}{dsi_pop_sumamry input from Ramboseli::dsi_pop_summary}

\item{type}{"wt_avg"(default),"dir_avg", or "grp_only". calculating weighted average by days_present across group, direct average, or keep only the group observed on the darting date}
}
\value{
a data frame with merged meta info
}
\description{
Merge the GE meta info with SR info
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotCNstatus.chr}
\alias{plotCNstatus.chr}
\title{plot copy number status estimated from getCNstatus chr by chr}
\usage{
plotCNstatus.chr(CNstatus, BAF, LRR, sampleName)
}
\arguments{
\item{CNstatus}{output from getCNstatus()}

\item{BAF}{three columns (chr, position, BAF values)}

\item{LRR}{three columns (chr, position, log2(read depth ratio) values)}

\item{sampleName}{string; output prefix}
}
\description{
plot copy number status estimated from getCNstatus chr by chr
}


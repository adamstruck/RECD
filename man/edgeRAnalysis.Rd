\name{edgeRAnalysis}
\alias{edgeRAnalysis}
\title{edgeRAnalysis}
\usage{
  edgeRAnalysis(countsTable, designTable, output_file)
}
\arguments{
  \item{countsTable}{}

  \item{designTable}{}

  \item{output_file}{path to and name of file to output
  results to.}
}
\value{
  TSV file
}
\description{
  Try to detect repeats that have gotten longer than the
  typical inner mate distance
}
\details{
  Uses the edgeR package to try and detect differential
  read class frequencies (facing versus flanking) between
  tumor and normal samples.
}
\author{
  Adam Struck - Intern
}


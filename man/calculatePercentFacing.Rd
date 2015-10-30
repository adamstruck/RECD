\name{calculatePercentFacing}
\alias{calculatePercentFacing}
\title{calculatePercentFacing}
\usage{
  calculatePercentFacing(facing_file, flanking_file,
    output_dir)
}
\arguments{
  \item{facing_file}{.facing file}

  \item{flanking_file}{.flanking file}

  \item{output_dir}{path to output results}
}
\value{
  TSV file
}
\description{
  Calculate the ratio of facing to flanking read pairs for
  each repeat.
}
\details{
  ratio = [ facing / ( flanking + facing ) ] * 100
}
\author{
  Adam Struck - Intern
}


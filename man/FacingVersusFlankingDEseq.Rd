\name{FacingVersusFlankingDEseq}
\alias{FacingVersusFlankingDEseq}
\title{FacingVersusFlankingDEseq}
\usage{
  FacingVersusFlankingDEseq(facing_file_list,
    flanking_file_list, pheno_table, output_file)
}
\arguments{
  \item{facing_file_list}{list of all .facing sample files}

  \item{flanking_file_list}{list of all .flanking sample
  files}

  \item{pheno_table}{ExpressionPlot type pheno table}

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
  Uses the DESeq package to try and detect differential
  read class frequencies (facing versus flanking) between
  tumor and normal samples.
}
\author{
  Adam Struck - Intern
}


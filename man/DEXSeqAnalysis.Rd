\name{DEXSeqAnalysis}
\alias{DEXSeqAnalysis}
\title{DEXSeqAnalysis}
\usage{
  DEXSeqAnalysis(counts_table, fdata_table, pheno_table,
    condition, grouping, output_file)
}
\arguments{
  \item{countsTable}{}

  \item{fdata}{}

  \item{pheno_table}{}

  \item{condition}{condition to base the differential
  analysis on.}

  \item{grouping}{list specifying a specific variable name
  and level to subset the data by. Only samples that have
  this trait (level) will be examined. [ Example: grouping
  = list(name="Tissue.Group.Oncology", level="Colon") ]}

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
  Uses the DEXSeq package to try and detect differential
  read class frequencies (facing versus flanking) between
  tumor and normal samples.
}
\author{
  Adam Struck - Intern
}


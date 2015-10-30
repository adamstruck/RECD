\name{makeDEXSeqTables}
\alias{makeDEXSeqTables}
\title{makeDEXSeqTables}
\usage{
  makeDEXSeqTables(facing_file_list, flanking_file_list,
    output_dir = "~", prefix)
}
\arguments{
  \item{facing_file_list}{list of all .facing sample files}

  \item{flanking_file_list}{list of all .flanking sample
  files}

  \item{output_dir}{path to and name of file to output
  results to.}
}
\value{
  TSV file
}
\description{
  Makes the necessary input counts table and fdata table
  for running DEXSeq
}
\author{
  Adam Struck - Intern
}


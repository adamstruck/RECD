\name{mergeSampleData}
\alias{mergeSampleData}
\title{mergeSampleData}
\usage{
  mergeSampleData(input_files, attribute, output_file)
}
\arguments{
  \item{file_list}{character vector of sample files to
  merge}

  \item{attribute}{feature to merge from files}

  \item{output_file}{path to and name of desired output
  file}
}
\value{
  Summary of implied repeat lengths for all samples
}
\description{
  mergeFiles take in a list of files as its input and
  merges them into one file.
}
\details{
  Description of output file format: rows = repeat.id
  columns = attribute (such as mean.irl) from each input
  file; column names are derived from the input filename
  (all text preceding the first period in the filename).
  For CGP data, this should be the SAM identification.
}
\author{
  Adam Struck - Intern
}


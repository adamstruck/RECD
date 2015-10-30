\name{makeEdgeRTables}
\alias{makeEdgeRTables}
\title{makeEdgeRTables}
\usage{
  makeEdgeRTables(facing_file_list, flanking_file_list,
    pheno_table, min_total_counts = 50, output_dir = "~",
    prefix)
}
\arguments{
  \item{facing_file_list}{list of all .facing sample files}

  \item{flanking_file_list}{list of all .flanking sample
  files}

  \item{pheno_table}{ExpressionPlot type pheno table}

  \item{min_total_counts}{miniumum total number of read
  pairs associated with a given repeat across all samples}

  \item{output_dir}{path to and name of file to output
  results to.}
}
\value{
  TSV file
}
\description{
  Creates the necessary counts and design tables for
  running edgeR.
}
\author{
  Adam Struck - Intern
}


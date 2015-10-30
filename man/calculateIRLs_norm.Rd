\name{calculateIRLs_norm}
\alias{calculateIRLs_norm}
\title{calculateIRLs_norm}
\usage{
  calculateIRLs_norm(flanking_file, imd_dist,
    median_imd = NULL, upper_limit = 0.95, margin = 37.5,
    min_coverage = 10, output_dir)
}
\arguments{
  \item{flanking_file}{_flanking read pair file}

  \item{imd_dist}{.insert_len_dist file (i.e. the inner
  mate distance distribution for the sample)}

  \item{median_imd}{the median inner mate distance of your
  sample; can be determined empirically if imd_dist is
  supplied}

  \item{upper_limit}{the value used in mapATR %
  processLastAlignedBase [default= 0.95]}

  \item{margin}{the value used in mapATR %
  processLastAlignedBase [default= 37.5]}

  \item{min_coverage}{minimum number of flanking read
  paired that map to a given repeat required for IRL
  calculation [default= 10]}

  \item{output_dir}{path to output directory}
}
\value{
  _implied_repeat_lengths file
}
\description{
  For repeats that have read pairs flanking them, the
  implied expansion and length of that repeat is
  determined.
}
\details{
  Summary stats for the implied repeat exansions/lengths
  are also reported.
}
\author{
  Adam Struck - Intern
}


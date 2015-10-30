\name{processLastAlignedBase}
\alias{processLastAlignedBase}
\title{processLastAlignedBase}
\usage{
  processLastAlignedBase(con, table_name, input_file,
    imd_dist = NULL, upper_limit = 0.95,
    margin = (read_length/2), read_length = 75,
    min_repeat_length = (read_length/2), output_dir = ".",
    overwrite = FALSE, append_chr = TRUE)
}
\arguments{
  \item{con}{database connection}

  \item{table_name}{name of table in database}

  \item{input_file}{last_aligned_base.tsv file}

  \item{imd_dist}{.insert_len_dist file (i.e. the inner
  mate distance distribution for the sample)}

  \item{upper_limit}{values between 0 and 1 will be used to
  set prob in quantile() and then used to calculate the
  maximum allowed distance from the repeat. Values greater
  than 1 will be applied as the maximum. Facing reads whose
  ends are farther from the repeat than the value set by
  this parameter will be not be used. [default= 0.95]}

  \item{read_length}{read length [default= 75]}

  \item{output_dir}{output directory to write results
  [default= "."]}

  \item{overwrite}{overwrite output files if they exist
  [default= FALSE]}

  \item{append_chr}{appends "chr" to chromosome names of
  inner_mate_ranges file [default= TRUE]}

  \item{max_gap}{maximum distance between the last aligned
  base and the start or end of the repeat [default = 500]}

  \item{margin}{maximum bp overlap between each read and a
  given repeat repeat [default= read_length/2]}

  \item{min_repeat_length}{minimum repeat length in
  reference genome [default= read_length/2]}
}
\value{
  _facing, and _fully_repetitive tsv files
}
\description{
  Sorts read pairs into two categories: facing, and fully
  repetitive.
}
\details{
  Description of output files: facing: one read flanks or
  is anchored outside the repeat (minimum = margin) while
  its pair aligned inside the repeat (or didn't meet
  minimum anchor requirements) fully_repetitive: both reads
  in pair aligned inside the repeat (or didn't meet minimum
  anchor requirements)
}
\author{
  Adam Struck - Intern
}


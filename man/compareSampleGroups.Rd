\name{compareSampleGroups}
\alias{compareSampleGroups}
\title{compareSampleGroups}
\usage{
  compareSampleGroups(merged_sample_data, pheno_table,
    grouping = list(list(name = "plot.metaclass")),
    test = "ttest", matched_patient_tissues = FALSE,
    output_file)
}
\arguments{
  \item{merged_sample_data}{A data frame of numerical data
  that you wish to group and compare using a student's
  t-test.}

  \item{pheno_table}{A data frame giving the group labels
  (phenotypic data) for each point in ‘merged_sample_data’.
  It should have the same number of rows as
  ‘merged_sample_data’ is long. If ‘merged_sample_data’ is
  named and ‘pdata’ has rownames then they are checked for
  equality.}

  \item{grouping}{A list specifying the levels of grouping
  your data based on the ‘pdata’. Each entry of the list
  has several componenets: $name The name of a field from
  ‘pdata’ $level a specific level of the specified $name to
  filter the data by. Only samples with this level will be
  considered.  ONLY ONE $NAME CAN BE SPECIFIED WITHOUT A
  CORRESPONDING LEVEL.}

  \item{test}{"ttest" or "cor"; perform either a student's
  t test on the supplied data or calculate the Spearmans'
  correlation between patient tumor/normal samples.}

  \item{matched_patient_tissues}{logical, if TRUE only
  samples with matched tumor and normal samples will be
  considered.}

  \item{output_file}{Desired name of output file}
}
\value{
  TSV file containing t-test statistics for each repeat
}
\description{
  Perform a Student's t-test or a Spearmans's correlation
  test on different groups of your data.
}
\details{
  See below:
}
\author{
  Adam Struck - Intern
}


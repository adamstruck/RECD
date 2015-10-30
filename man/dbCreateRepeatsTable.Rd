\name{dbCreateRepeatsTable}
\alias{dbCreateRepeatsTable}
\title{dbCreateRepeatsTable}
\usage{
  dbCreateRepeatsTable(con, file, table_name)
}
\arguments{
  \item{con}{database connection}

  \item{file}{repeats data to be stored (e.g. output from
  parseTRFoutput.py)}

  \item{table_name}{name of table to create}
}
\value{
  SQL/SQLite database table
}
\description{
  Creates a table in your SQLite database.
}
\details{
  .. content for \details{} ..
}
\author{
  Adam Struck
}


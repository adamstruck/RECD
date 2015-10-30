\name{downloadSatellog}
\alias{downloadSatellog}
\title{downloadSatellog}
\usage{
  downloadSatellog(URL = "http://satellog.bcgsc.ca/db_files/satellog_db_tables.tar.gz",
    DIR, DBNAME = "satellog_db", USR = "root")
}
\arguments{
  \item{URL}{url for database download}

  \item{DIR}{desired output path}

  \item{DBNAME}{desired name of database that will be
  created to store the data}

  \item{USR}{user name for database access}
}
\value{
  untarred Satellog table files and MySQL database
}
\description{
  This package is used to download all tables from Satellog
  and store them in a local MySQL database
}
\details{
  .. content for \details{} ..
}
\author{
  Adam Struck
}


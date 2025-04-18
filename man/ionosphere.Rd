% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kknn-package.R
\docType{data}
\name{ionosphere}
\alias{ionosphere}
\title{Johns Hopkins University Ionosphere Database}
\format{
A data frame with 351 observations on the following 35 variables.
The first 34 continuous covariables are used for the prediction. The 35th
attribute is either \code{g} ("good") or \code{b} ("bad") according to the
definition summarized above. This is a binary classification task.
}
\source{
% \itemize{\item Vince Sigillito (vgs@aplcen.apl.jhu.edu), Space
Physics Group, Applied Physics Laboratory, Johns Hopkins University, Johns
Hopkins Road, Laurel, MD 20723

The data have been taken from the UCI Machine Learning Database Repository
\cr \url{https://archive.ics.uci.edu/}\cr and were
converted to R format by \email{klaus.schliep@gmail.com }.
}
\description{
This radar data was collected by a system in Goose Bay, Labrador.  This
system consists of a phased array of 16 high-frequency antennas with a total
transmitted power on the order of 6.4 kilowatts.  See the paper for more
details.  The targets were free electrons in the ionosphere. "Good" radar
returns are those showing evidence of some type of structure in the
ionosphere.  "Bad" returns are those that do not; their signals pass through
the ionosphere.
}
\details{
Received signals were processed using an autocorrelation function whose
arguments are the time of a pulse and the pulse number.  There were 17 pulse
numbers for the Goose Bay system.  Instances in this database are described
by 2 attributes per pulse number, corresponding to the complex values
returned by the function resulting from the complex electromagnetic signal.
}
\examples{

data(ionosphere)

}
\keyword{datasets}

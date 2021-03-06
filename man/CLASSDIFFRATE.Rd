% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMFASTICAALG.R
\name{CLASSDIFFRATE}
\alias{CLASSDIFFRATE}
\title{CLASSDIFFRATE}
\usage{
CLASSDIFFRATE(factor1, factor2)
}
\arguments{
\item{factor1}{First factor.}

\item{factor2}{Second factor of the same length as the First factor.}
}
\description{
A function calculates classification difference rate between two factors. It
is used in interpreting info stored in EMFASTICA object.
}
\examples{
## An example evaluates the classification difference rate
## between two classification results in the form of factors
fac1<-factor(c(1,4,2,3,1,1,3,3,1,2,2,1))
fac2<-factor(c(3,1,2,2,1,2,4,3,2,3,1,1))
CLASSDIFFRATE(fac1, fac2)

}
\keyword{~kwd1}
\keyword{~kwd2}

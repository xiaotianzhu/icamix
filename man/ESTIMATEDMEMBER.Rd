% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMFASTICAALG.R
\name{ESTIMATEDMEMBER}
\alias{ESTIMATEDMEMBER}
\title{ESTIMATEDMEMBER}
\usage{
ESTIMATEDMEMBER(rstICAMIX)
}
\arguments{
\item{rstICAMIX}{An EMFASTICAALG object.}
}
\value{
%% ~Describe the value returned %% If it is a LIST, use
\item{estimatedmember }{A factor with levels representing the estimated
classes.} %% \item{comp2 }{Description of 'comp2'} %% ...
}
\description{
A function calculates estimated class membership from an EMFASTICAALG
object.
}
\examples{

## An Example that runs the NSMM-ICA algorithm on Cohen's tone data
data(tonedata, package="mixtools")

b <- EMFASTICAALG(tonedata, 2, h=0, tol=1e-8)
estimatedgroup <- ESTIMATEDMEMBER(b) # estimated species info


}
\keyword{~kwd1}
\keyword{~kwd2}

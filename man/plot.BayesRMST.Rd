\name{plot.BayesRMST}
\alias{plot.BayesRMST}
\title{Plots for BayesRMST output.}
\description{
This function is used to obtain survival, rmst and omega plots.
}
\usage{
\method{plot}{BayesRMST}(x,
                           survival = FALSE,
                           omega = FALSE,
                           rmst = FALSE,
                           CI = TRUE,
                           ...)
##
}
\arguments{
\item{x}{Object of class inheriting from "BayesRMST".}
\item{survival}{If TRUE then plot survival curves.}
\item{omega}{If TRUE then plot omega curves.}
\item{rmst}{If TRUE then plot RMST curves.}
\item{CI}{If TRUE then plot CI for the survival curves.}
\item{...}{Other arguments.}
}
\seealso{
\code{\link{survRMST}}.
}
\examples{
\dontrun{
##

plot(out) # where out is the output from BayesRMST class
plot(out, survival = TRUE)
plot(out, omega = TRUE)
plot(out, rmst = TRUE)

##
}
}
\keyword{spT}

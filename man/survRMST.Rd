\name{survRMST}
\alias{survRMST}
\title{Bayesian survival model and RMST estimation.}
\description{
This function is used to obtain optimal time for RMST using Bayesian approach.
}
\usage{
survRMST(data, modelType = c("KM"), family = "weibull",
    iter = 2000, chains = 1, cores = 1, control = list(adapt_delta = 0.95), seed = 1234,
    prior_expert_list = NULL, eta = 0.95)
}
\arguments{
\item{data}{Data frame containing the variables of survival data with columns in order: "time", "event/status" and "treatment".}
\item{modelType}{Choices are: "BayesPara", "BayesMspline" and "BayesNonPara".}
\item{family}{Only applicable if "BayesPara" is used, default family is "weibull" distribution.}
\item{iter}{...}
\item{chains}{...}
\item{cores}{...}
\item{control}{...}
\item{seed}{...}
\item{prio_expert_list}{... list(c(5, 6), c(10, 12), c(8.5, 10))}
\item{eta}{...}
}
\value{
\item{rmst}{...}
\item{model_fit}{...}
\item{decision}{...}
}
}
\references{

 Bakar, K. S. et al. (2026).

}
\seealso{
\code{\link{plot.BayesRMST}}
\examples{
\dontshow{

}
\donttest{
##
######################################################

## Datasets for evaluation of cetuximab in head and neck cancer

data(cetux)
y <- cetux[,c(2:4)]
names(y) <- c("X", "delta", "A")
y$A <- as.numeric(y$A)-1

# Bayesian M-spline

out_ms <- survRMST(data=y, modelType="BayesMspline")
print(out_ms)
plot(out_ms)
plot(out_ms, omega=TRUE)

# Bayesian non-parametric (DDP)

out_dp <- survRMST(data=y, modelType="BayesNonPara")
print(out_dp)
plot(out_dp)

# Bayesian parametric

out_para <- survRMST(data=y, modelType="BayesPara")
print(out_para)
plot(out_para)


######################################################
##
}
}
\keyword{spT}

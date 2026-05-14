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
\item{family}{Only applicable if "BayesPara" is used, default family is "weibull" distribution, see \link{brmsfamily}.}
\item{iter}{Number of MCMC iterations, see \link{stan}.}
\item{chains}{Number of chains for MCMC samples, see \link{stan}.}
\item{cores}{Number of cores, see \link{stan}.}
\item{control}{A named list of parameters to control the sampler's behavior, see \link{stan}.}
\item{seed}{Seed for random number generation.}
\item{prior_expert_list}{A list argument for expert knowledge on tau, e.g., list(c(5, 6), c(10, 12), c(8.5, 10))}
\item{eta}{Decision threshold for calculating Bayesian power, default eta = 0.95}
}
\value{
\item{rmst}{Optimal RMST table.}
\item{model_fit}{Model fit object.}
\item{decision}{Decision for Bayesian Power.}
}
\references{

 Bakar, K. S. et al. (2026).

}
\seealso{
\code{\link{plot.BayesRMST}}
}
\examples{
\dontrun{

##
######################################################

## Datasets for evaluation of cetuximab in head and neck cancer

data(cetux, package="survextrap")
head(cetux)
y <- cetux[,c("years","d","treat")] # time, event, trt
y$treat = as.numeric(y$treat)-1
head(y)

## North Central Cancer Treatment Group: Lung Cancer Data

data(cancer, package="survival")
x <- cancer[,c("time","status","sex")] # time, event, trt
x$sex = x$sex-1; x$status = x$status-1
head(x)

## Acute Myelogenous Leukemia survival data

head(aml) # time, event, trt
aml$x = as.numeric(aml$x)-1
head(aml)

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

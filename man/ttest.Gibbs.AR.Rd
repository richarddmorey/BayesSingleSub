\name{ttest.Gibbs.AR}
\alias{ttest.Gibbs.AR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain Bayesian t test and posterior distributions for single case data
}
\description{
This function computes Bayes factors for the mean difference between two phases of a single subject data sequence, 
using Gibbs sampling. Posterior samples of parameters are also provided.
}
\usage{
ttest.Gibbs.AR(before, after, iterations = 1000, areaNull = c(-0.2, 0.2), 
              leftSided = TRUE, treat = NULL, r.scale = 1, alphaTheta = 1, 
		          betaTheta = 5, sdMet = 0.3, progress = TRUE, 
              return.chains = FALSE, return.onesided = FALSE)
}
\arguments{
  \item{before}{
A vector of observations, in time order, taken in Phase 1 (e.g., before the treatment).
}
  \item{after}{
A vector of observations, in time order, taken in Phase 2 (e.g., after the treatment). 
}
  \item{iterations}{
Number of Gibbs sampler iterations to perform.
}
\item{areaNull}{
Only used if \code{return.chains=TRUE}. Bounds for the interval null hypothesis for delta.
}
\item{leftSided}{
Only used if \code{return.onesided=TRUE}. Should the one sided Bayes factor be left sided?
}
\item{treat}{
Vector of dummy coding, indicating Phase 1 and Phase 2; default is -.5 for Phase 1 and .5 for Phase 2.
}
 \item{r.scale}{
Prior scale for delta (see Details below).
}
\item{alphaTheta}{
The alpha parameter of the beta prior on theta (see Details below).
}
 \item{betaTheta}{
The beta parameter of the beta prior on theta (see Details below).
}
 \item{sdMet}{
Scale for the Metropolis-Hastings sampling of theta (see Details below).
}
 \item{progress}{
Report progress with a text progress bar?
}
\item{return.chains}{
Return posterior samples of parameters and area null Bayes factor?.
}
\item{return.onesided}{
Return one sided Bayes factor?
}
}
\details{
This function computes Bayes factors for the mean difference between two data sequences from a single subject,
 using Gibbs sampling. The Bayes factors compare null hypotheses of no true mean difference against 
alternative hypotheses of a true mean difference. A Bayes factor larger than 1 supports the null hypothesis,
 a Bayes factor smaller than 1 supports the alternative hypothesis. Auto-correlation of the errors is modeled by a first
 order auto-regressive process.

Posterior distributions of the model parameters can also be obtained. Model parameters of interest include mu0 (overall mean), 
delta (standardized mean difference), sig2 (variance of the random shocks), and rho (auto-correlation).

A Cauchy prior is placed on the standardized mean difference delta. The \code{r.scale} argument controls the scale of
 this Cauchy prior, with \code{r.scale=1} yielding a standard Cauchy prior. A noninformative Jeffreys prior is 
placed on the variance of the random shocks of the auto-regressive process. A beta prior is placed on the auto-correlation
 rho. The \code{alphaTheta} and \code{betaTheta} arguments control the form of this beta prior.
 
 Missing data are sampled from the likelihood function, conditioned at the observed data, at each iteration of the Gibbs sampler.
}
\value{
A list containing the following:
\item{logbf}{An MCMC estimate of the log two sided point null Bayes factor, computed using the Savage-Dickey method (Morey, Rouder, Pratte, and Speckman, submitted).}
\item{chains}{Only returned if \code{return.chains=TRUE}. An object of type MCMC containing the chains for each parameter.}
\item{acc}{Only returned if \code{return.chains=TRUE}. The Metropolis-Hastings acceptance rate.}
\item{logbfArea}{Only returned if \code{return.chains=TRUE}. An MCMC estimate of the log two sided interval null Bayes factor.}
\item{logbfOnesided}{Only returned if \code{return.onesided=TRUE}. An MCMC estimate of the log one sided point null Bayes factor.}
}
\references{
De Vries, R. M. \& Morey, R. D. (submitted). Bayesian hypothesis testing Single-Subject Data. Psychological Methods.

R code guide: \url{http://drsmorey.org/research/rdmorey/}
}
\author{
Richard D. Morey and Rivka de Vries
}
\note{
For a more accurate method of computing the Bayes factor, see \code{\link{ttest.MCGQ.AR}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ttest.MCGQ.AR}}, \code{\link{trendtest.Gibbs.AR}}, \code{\link{trendtest.MC.AR}}
}
\examples{
## Define data
data = c(87.5, 82.5, 53.4, 72.3, 94.2, 96.6, 57.4, 78.1, 47.2,
 80.7, 82.1, 73.7, 49.3, 79.3, 73.3, 57.3, 31.7, 50.4, 77.8,
 67, 40.5, 1.6, 38.6, 3.2, 24.1)

## Obtain log Bayes factor
logBF = ttest.Gibbs.AR(data[1:10], data[11:25])

## Obtain log Bayes factor, chains, and log interval null Bayes factor
output = ttest.Gibbs.AR(data[1:10], data[11:25], return.chains = TRUE, areaNull = c(-0.2, 0.2))

## Look at the posterior distribution of the mean
plot(output$chains[,1])

## Obtain summary statistics of posterior distributions
summary(output$chains)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{models}
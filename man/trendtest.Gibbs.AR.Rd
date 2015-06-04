\name{trendtest.Gibbs.AR}
\alias{trendtest.Gibbs.AR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain Bayesian trend test and posterior distributions for single case data
}
\description{
This function computes Bayes factors for the trend and intercept differences between two phases of a single subject data sequence, 
using Gibbs sampling. Posterior samples of parameters are also provided.
}
\usage{
trendtest.Gibbs.AR(before, after, iterations = 1000, 
                  intArea = c(-0.2,0.2), slpArea = c(-0.2, 0.2), 
                  leftSidedInt = TRUE, leftSidedSlp = TRUE, 
                  r.scaleInt = 1, r.scaleSlp = 1, 
		  alphaTheta = 1, betaTheta = 5, sdMet = 0.3, 
                  progress = TRUE, return.chains = FALSE,
                  return.onesided = FALSE)
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
\item{intArea}{
Only used if \code{return.chains=TRUE}. Bounds for the interval null hypothesis for the intercept difference.
}
\item{slpArea}{
Only used if \code{return.chains=TRUE}. Bounds for the interval null hypothesis for the trend difference.
}
\item{leftSidedInt}{
Only used if \code{return.onesided = TRUE}. Should the one sided Bayes factor for the intercept difference be one sided?
}
\item{leftSidedSlp}{
Only used if \code{return.onesided = TRUE}. Should the one sided Bayes factor for the trend difference be one sided?
}
 \item{r.scaleInt}{
Prior scale for the intercept difference (see Details below).
}
\item{r.scaleSlp}{
Prior scale for the trend difference (see Details below).
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
Return posterior samples of parameters?.
}
\item{return.onesided}{
Return one sided Bayes factors?
}
}
\details{
This function computes Bayes factors for the differences in trend and intercept between two data sequences from a single subject,
 using Gibbs sampling. The Bayes factor for trend difference compares the null hypothesis of no true trend difference against the 
alternative hypothesis of a true trend difference. The Bayes factor for intercept difference compares the null hypothesis of no 
true intercept difference against the alternative hypothesis of a true intercept difference. Also, a joined Bayes factor for the 
trend and intercept combined is provided. Bayes factors larger than 1 support the null hypothesis, Bayes factors smaller than 1 support
 the alternative hypothesis. Auto-correlation of the errors is modeled by a first order auto-regressive process.

Posterior distributions of the model parameters can also be obtained. Model parameters of interest include mu0 (overall mean), 
sig*delta (difference between intercepts), beta0 (overall trend), sig*beta1 (difference between trends), sig2 (variance of the random 
shocks), and rho (auto-correlation).

Cauchy priors are placed on the standardized trend and intercept differences. The \code{r.scaleInt} and \code{r.scaleSlp} arguments 
control the scales of these Cauchy priors, with \code{r.scaleInt = 1} and \code{r.scaleSlp = 1} yielding standard Cauchy priors. 
A noninformative Jeffreys prior is placed on the variance of the random shocks of the auto-regressive process. A beta prior is 
placed on the auto-correlation theta. The \code{alphaTheta} and \code{betaTheta} arguments control the form of this beta prior.

Missing data are sampled from the likelihood function, conditioned at the observed data, at each iteration of the Gibbs sampler.
}
\value{
A list containing the following:
\item{logbf}{MCMC estimates of the log two sided point null Bayes factors, computed using the Savage-Dickey method (Morey, Rouder, Pratte, and Speckman, submitted).}
\item{chains}{Only returned if \code{return.chains=TRUE}. An object of type MCMC containing the chains for each parameter.}
\item{acc}{Only returned if \code{return.chains=TRUE}. The Metropolis-Hastings acceptance rate.}
\item{logbfArea}{Only returned if \code{return.chains=TRUE}. MCMC estimates of the log two sided interval null Bayes factors.}
\item{logbfOnesided}{Only returned if \code{return.onesided=TRUE}. MCMC estimates of the log one sided point null Bayes factors.}
}
\references{
De Vries, R. M. \& Morey, R. D. (submitted). Bayesian hypothesis testing Single-Subject Data. Psychological Methods.

R code guide: \url{http://drsmorey.org/research/rdmorey/}
}
\author{
Richard D. Morey and Rivka de Vries
}
\note{
For a more accurate method of computing the Bayes factor, see \code{\link{trendtest.MC.AR}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trendtest.MC.AR}}, \code{\link{ttest.Gibbs.AR}}, \code{\link{ttest.MCGQ.AR}}
}
\examples{
## Define data
data = c(87.5, 82.5, 53.4, 72.3, 94.2, 96.6, 57.4, 78.1, 47.2,
 80.7, 82.1, 73.7, 49.3, 79.3, 73.3, 57.3, 31.7, 50.4, 77.8,
 67, 40.5, 1.6, 38.6, 3.2, 24.1)

## Obtain log Bayes factors
logBFs = trendtest.Gibbs.AR(data[1:10], data[11:25])

## Obtain log Bayes factors, chains, and log interval null Bayes factors
output = trendtest.Gibbs.AR(data[1:10], data[11:25], 
        return.chains = TRUE, intArea = c(-0.2,0.2), 
        slpArea = c(-0.2, 0.2))

## Look at the posterior distribution of the mean
plot(output$chains[,1])

## Obtain summary statistics of posterior distributions
summary(output$chains)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{models}
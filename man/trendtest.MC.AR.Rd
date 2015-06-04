\name{trendtest.MC.AR}
\alias{trendtest.MC.AR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain Bayesian trend test or single case data
}
\description{
This function computes Bayes factors for the trend and intercept differences between two phases of a single subject data sequence, 
using Monte Carlo integration.
}
\usage{
trendtest.MC.AR(before, after, iterations = 1000, 
                r.scaleInt = 1, r.scaleSlp = 1, 
                alphaTheta = 1, betaTheta = 5, 
                progress = TRUE)
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
 \item{progress}{
Report progress with a text progress bar?
}
}
\details{
This function computes Bayes factors for the differences in trend and intercept between two data sequences from a single subject,
 using monte carlo integration. The Bayes factor for trend difference compares the null hypothesis of no true trend difference against 
the alternative hypothesis of a true trend difference. The Bayes factor for intercept difference compares the null hypothesis of no 
true intercept difference against the alternative hypothesis of a true intercept difference. Also, a joined Bayes factor for the 
trend and intercept combined is provided. Bayes factors larger than 1 support the null hypothesis, Bayes factors smaller than 1 support
 the alternative hypothesis. Auto-correlation of the errors is modeled by a first order auto-regressive process. 

Cauchy priors are placed on the standardized trend and intercept differences. The \code{r.scaleInt} and \code{r.scaleSlp} arguments 
control the scales of these Cauchy priors, with \code{r.scaleInt = 1} and \code{r.scaleSlp = 1} yielding standard Cauchy priors. 
A noninformative Jeffreys prior is placed on the variance of the random shocks of the auto-regressive process. A beta prior is 
placed on the auto-correlation theta. The \code{alphaTheta} and \code{betaTheta} arguments control the form of this beta prior.

Missing data are handled by removing the locations of the missing data from the design matrix and error covariance matrix. 
}
\value{
A matrix containing the Monte carlo estimates of the log Bayes factors.
}
\references{
De Vries, R. M. \& Morey, R. D. (submitted). Bayesian hypothesis testing Single-Subject Data. Psychological Methods.

R code guide: \url{http://drsmorey.org/research/rdmorey/}
}
\author{
Richard D. Morey and Rivka de Vries
}
\note{
To obtain posterior distributions and interval null Bayes factors, see \code{\link{trendtest.Gibbs.AR}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trendtest.Gibbs.AR}}, \code{\link{ttest.Gibbs.AR}}, \code{\link{ttest.MCGQ.AR}}
}
\examples{
## Define data
data = c(87.5, 82.5, 53.4, 72.3, 94.2, 96.6, 57.4, 78.1, 47.2,
 80.7, 82.1, 73.7, 49.3, 79.3, 73.3, 57.3, 31.7, 50.4, 77.8,
 67, 40.5, 1.6, 38.6, 3.2, 24.1)

## Obtain log Bayes factors
logBFs = trendtest.MC.AR(data[1:10], data[11:25])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{models}
\name{BayesSingleSub-package}
\alias{BayesSingleSub-package}
\alias{BayesSingleSub}
\docType{package}
\title{
Functions to obtain Bayes factor hypothesis tests and posterior samples of parameters for single case data.
}
\description{
This package contains functions to compute Bayes factors for the mean difference, trend difference, and intercept 
difference in single case data sequences and data sequences from groups of subjects. Some functions also provide 
posterior samples of parameters and interval null Bayes factors.
}
\details{
\tabular{ll}{
Package: \tab BayesSingleSub\cr
Type: \tab Package\cr
Version: \tab 0.6.3\cr
Date: \tab 2014-09-16\cr
License: \tab GPL 2.0\cr
LazyLoad: \tab yes\cr
}
The following functions are currently implemented, with more to follow:

ttest.Gibbs.AR(): Compute Bayes factor for mean difference based on single subject, using Gibbs sampling, and obtain posterior samples of parameters
ttest.MCGQ.AR(): Compute Bayes factor for mean difference based on single subject, using monte carlo or gaussian quadrature integration
trendtest.Gibbs.AR(): Compute Bayes factors for trend and intercept differences based on single subject, using Gibbs sampling, and obtain posterior samples of parameters
trendtest.MC.AR(): Compute Bayes factors for trend and intercept differences based on single subject, using monte carlo integration
bayesMLint(): compute Bayes factors for mean difference based on group of subjects, using Gibbs sampling, and obtain posterior samples of parameters
bayesMLtrend(): compute Bayes factors for trend and intercept differences based on group of subjects, using Gibbs sampling, and obtain posterior samples of parameters

}
\author{
Richard D. Morey and Rivka de Vries

Maintainer: Richard D. Morey <richarddmorey@gmail.com>
}
\references{
De Vries, R. M., Hartogs, B. M. A., \& Morey, R.D. (submitted). A tutorial on computing Bayes factors for single-subject designs.
De Vries, R. M. \& Morey, R. D. (2013). Bayesian hypothesis testing for Single-Subject Data. Psychological Methods, 18, 165-185.
De Vries, R.M., Morey, R.D., \& Tendeiro, J.N. In preparation.

R code guide: \url{http://supp.apa.org/psycarticles/supplemental/a0031037/newsupplementalfilefromauthor122012.pdf}
}
\keyword{htest}
\examples{
## See specific functions for examples
}

\name{bayesMLtrend}
\alias{bayesMLtrend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain Bayesian test for intercept and trend change and posterior distributions for repeated measurement data
}
\description{
This function computes Bayes factors regarding the intercept and trend difference between two phases, for subjects measured repeatedly over time, using Gibbs sampling. Posterior samples of parameters are also provided.
}
\usage{
bayesMLtrend(y, subject, occ, int, intpoint = 0, 
            testPatt=rep(0, length(unique(subject))), 
            iter=80000, n.burnin = 100, 
            thin = 20, r1=.5, r2=.5, r3=.1, r4=.1, beta_theta=5, trunc = 1,
            aInt = -.0001, bInt=.0001, aTrend = -.0001, bTrend = .0001,
            mu_mu_sig2_z=0, sig2_mu_sig2_z=1, return.chains = TRUE, 
            return.BFnrchanged = TRUE, raftery.check = FALSE)
}
\arguments{
  \item{y}{
A vector of observations, in time order per subject, with all subjects combined in one vector.
}
  \item{subject}{
A numeric vector indicating subject number. Should be in increasing order.
}
  \item{occ}{
A numeric vector indicating measurement occasion; time intervals between measurement occasions are taken into account.
}
  \item{int}{
A numeric, dichotomous vector indicating whether the observation was made in Phase 1 or in Phase 2. 
}
  \item{intpoint}{
A numeric vector indicating the intervention point for each subject; should be the distance in time units between the first measurement occasion of the subject and the intervention point for the subject (see Details below).
}
  \item{testPatt}{
A specific pattern of change over subjects to test (see Details below). 
}
  \item{iter}{
Number of Gibbs sampler iterations to perform. 
}
  \item{n.burnin}{
Number of burn in iterations. 
}
  \item{thin}{
Thinning rate. 
}
  \item{r1}{
r1 scale  (see Details below). 
}
  \item{r2}{
r2 scale  (see Details below). 
}
  \item{r3}{
r3 scale  (see Details below). 
}
  \item{r4}{
r4 scale  (see Details below). 
}
 \item{beta_theta}{
The beta parameter of the Beta prior on theta (see Details below).
}
 \item{trunc}{
Truncation of the right side of the Beta prior on theta (see Details below).
}
 \item{aInt}{
Lower bound for the interval null hypothesis for the individual standardized intercept difference delta_i.
}
 \item{bInt}{
Upper bound for the interval null hypothesis for the individual standardized intercept difference delta_i.
}
 \item{aTrend}{
Lower bound for the interval null hypothesis for the individual standardized trend difference beta1_i.
}
 \item{bTrend}{
Upper bound for the interval null hypothesis for the individual standardized trend difference beta1_i.
}
 \item{mu_mu_sig2_z}{
Mean of the Normal prior on mu_sig2_e.
}
 \item{sig2_mu_sig2_z}{
Variance of the Normal prior on mu_sig2_e.
}
\item{return.chains}{
Return posterior samples of parameters?.
}
\item{return.BFnrchanged}{
Return Bayes factor that at least X subjects have changed? (See Details below).
}
\item{raftery.check}{
Compute raftery diagnostics?.
}
}
\details{
Based on a data set consisting of a group of subjects measured repeatedly during two phases of a design (e.g., before and after an intervention), this function can compute 
(1) logbfSubjectInt: a Bayes factor for intercept change per subject: compares the null hypothesis of negligible intercept change to the alternative hypothesis of non-negligible intercept change for each subject.

(2) logbfSubjectTrend: a Bayes factor for trend change per subject: compares the null hypothesis of negligible trend change to the alternative hypothesis of non-negligible trend change for each subject.

(3) logbfSubjectIT: a Bayes factor for joint intercept and trend change per subject: compares the null hypothesis of negligible intercept and trend change to the alternative hypothesis of non-negligible intercept and trend change for each subject.

(4) logbfPattern: for intercept, trend, and joint intercept and trend change, respectively, a Bayes factor for a specific pattern of change over subjects. Here a null hypotheses can be specified which states that, for instance, all subjects from a control group have only negligibly changed (coded 0), while all subjects from an experimental group have non-negligibly changed (coded 1), and

(5) logbfAtleastX: for intercept, trend, and joint intercept and trend change, respectively, Bayes factors comparing the null hypothesis that fewer than X subjects have non-negligibly changed to the alternative hypothesis that at least X subjects have non-negligibly changed, where X ranges from 1 to N, with N the number of subjects.

Intercept change represents level change at intervention point. Bayes factors larger than 1 support the null hypothesis, Bayes factors smaller than 1 support the alternative hypothesis. Auto-correlation of the errors is modeled by a first order auto-regressive process. JAGS (Plummer, 2003) is used for Gibbs sampling of parameters. The Bayes factors are estimated using a procedure analogous to conditional marginal density estimation (Gelfand & Smith, 1990; Morey et al., 2011). The null hypotheses of no change for each subject are interval null hypotheses, but by making the interval very small around zero, the interval null hypotheses closely resemble point null hypotheses.

Posterior distributions of the model parameters can also be obtained. Model parameters of interest include mu_i (individual overall mean), delta*sig2z_i (individual intercept difference), beta0_i (individual overall trend), beta1*sig2z_i (individual trend difference), sig2z_i (individual variance of random shocks), mu_delta*sig2z (population mean of the individual intercept differences), mu_beta1*sig2z (population mean of the individual trend differences), musig2z (location parameter for log normal prior on individual variances of random shocks), sig2sig2z (scale parameter for log normal prior on individual variances of random shocks), and rho (auto-correlation)).

Cauchy priors are placed on the standardized individual intercept and trend differences delta_i and beta1_i and on their population means mu_delta and mu_beta1. The \code{r1} argument represents the r scale of the Cauchy prior for mu_delta. The sum of the \code{r1} and code{r2} arguments sets the r scale of the Cauchy priors for the delta_i parameters. The \code{r3} argument represents the r scale of the Cauchy prior for mu_beta1 and the sum of the \code{r3} and code{r4} arguments sets the r scale of the Cauchy priors for the beta1_i parameters. A flat prior is placed on the individual overall means and trends, mu_i and beta0_i, and a log Normal(musig2z, sig2sig2z) prior on the individual variances of the random shocks. A Normal(mu_mu_sig2_z, sig2_mu_sig2_z) prior is put on the location parameter musig2z and an Inverse Gamma(.5,.5) prior is placed on the scale parameter sig2sig2z. The \code{mu_mu_sig2_z} and  \code{sig2_mu_sig2_z} arguments control the mean and variance, respectively, of the Normal(mu_mu_sig2_z, sig2_mu_sig2_z) prior. A Beta prior is placed on the auto-correlation rho, where the \code{beta_theta} argument controls the b parameter of this Beta prior.

Initial values of the chains are randomly drawn from the prior distributions by JAGS.

The \code{intpoint} argument indicates the intervention point for each subject. The intervention point should be somewhere between the last measurement occasion from Phase 1 and the first measurement occasion from Phase 2, and should not be equal to one of these time points. If set to 0 (or another scalar), the intervention point for each subject is set to the average of the last measurement occasion from Phase 1 and the first measurement occasion from Phase 2.

The \code{trunc} argument can be used to truncate the right side of the Beta prior on the auto-correlation rho. It prevents the MCMC chain to obtain rho values close to 1, which might happen when time intervals between observations are large. 

The \code{return.BFnrchanged} argument controls whether the Bayes factor logbfAtleastX should be computed. The computation requires a loop going through 2^N iterations, which can take a lot of time when N is large. Due to memory limits this Bayes factor cannot be computed for N > 39.
}
\value{
A list containing the following:
\item{logbfSubjectInt}{A matrix containing the MCMC estimates of the log Bayes factors logbfSubjectInt for all subjects.}
\item{logbfSubjectTrend}{A matrix containing the MCMC estimates of the log Bayes factors logbfSubjectTrend for all subjects.}
\item{logbfSubjectIT}{A matrix containing the MCMC estimates of the log Bayes factors logbfSubjectIT for all subjects.}


\item{logbfPattern}{MCMC estimates of the log Bayes factors logbfPattern.}
\item{logbfAtleastX}{Only returned if \code{return.BFnrchanged=TRUE} and number of subjects N < 40. A matrix containing the MCMC estimates of the log Bayes factors logbfAtleastX, for X = 1 to N.}
\item{chains}{Only returned if \code{return.chains=TRUE}. A matrix containing the MCMC chains for the model parameters.}
\item{converged}{Only returned if \code{raftery.check = TRUE}. An indicator value for convergence: equals 1 if the lengths of the thinned chains are larger than the lengths required according to raftery diagnostic, 0 otherwise.}

}
\references{
Gelfand, A. and Smith, A. F. M. (1990). Sampling based approaches to calculating marginal densities. Journal of the American Statistical Association, 85, pp. 398-409.
Morey, R. D., Rouder, J. N., Pratte, M. S., and Speckman, P. L. (2011). Using MCMC chain outputs to efficiently estimate Bayes factors. Journal of Mathematical Psychology, 55, pp. 368-378.
Plummer, M. (2003). JAGS: a program for analysis of Bayesian graphical models
using Gibbs sampling. In: Proceedings of the 3rd international workshop on
distributed statistical computing.
}
\author{
Rivka M. de Vries and Richard D. Morey
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bayesMLint}}
}
\examples{
## Define data
y = rnorm(30) # observations
sub = c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5) # subject 
occ = c(0:5,0:5,0:5,0:5,0:5) # measurement occasion
int = c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2) # phase indicator


## Obtain log Bayes factors and chains
output = bayesMLtrend(y, sub, occ, int, iter=1000, n.burnin = 100, thin = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{models}
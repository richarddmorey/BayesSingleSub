\name{bayesMLint}
\alias{bayesMLint}
\title{
Obtain Bayesian test for mean change and posterior distributions for repeated measurement data
}
\description{
This function computes Bayes factors regarding the mean difference between two phases, for subjects measured repeatedly over time, 
using Gibbs sampling. Posterior samples of parameters are also provided.
}
\usage{
bayesMLint(y, subject, occ, int, 
          testPatt=rep(0, length(unique(subject))), 
          iter=10000, thin=1, r1=.5, r2=.5, 
          beta_theta=5, trunc = 1,
          a = -.0001, b=.0001, thres = 1e-12,
          sdmetlogsig2z=1.5, sdmetmusig2z=.9, sdmettheta=.3,
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
  \item{testPatt}{
A specific pattern of change over subjects to test (see Details below). 
}
  \item{iter}{
Number of Gibbs sampler iterations to perform. 
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
 \item{beta_theta}{
The beta parameter of the Beta prior on theta (see Details below).
}
 \item{trunc}{
Truncation of the right side of the Beta prior on theta (see Details below).
}
 \item{a}{
Lower bound for the interval null hypothesis for the individual standardized mean difference delta_i.
}
 \item{b}{
Upper bound for the interval null hypothesis for the individual standardized mean difference delta_i.
}
 \item{thres}{
Threshold for inverse check (see Details below).
}
 \item{sdmetlogsig2z}{
Scale for the Metropolis-Hastings sampling of log sig2_e (see Details below).
}
 \item{sdmetmusig2z}{
Scale for the Metropolis-Hastings sampling of mu_sig2_e (see Details below).
}
 \item{sdmettheta}{
Scale for the Metropolis-Hastings sampling of theta (see Details below).
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
(1) logbfSubject: a Bayes factor for mean change per subject: compares the null hypothesis of negligible change to the alternative hypothesis of non-negligible change for each subject. 

(2) logbfPattern: a Bayes factor for a specific pattern of mean change over subjects: here a null hypotheses can be specified which states that, for instance, all subjects from a control group have only negligibly changed (coded 0), while all subjects from an experimental group have non-negligibly changed (coded 1), and

(3) logbfAtleastX: Bayes factors comparing the null hypothesis that fewer than X subjects have non-negligibly changed to the alternative hypothesis that at least X subjects have non-negligibly changed, where X ranges from 1 to N, with N the number of subjects.

Bayes factors larger than 1 support the null hypothesis, Bayes factors smaller than 1 support the alternative hypothesis. Auto-correlation of the errors is modeled by a first order auto-regressive process. Gibbs sampling is used for sampling of parameters. The Bayes factors are estimated using a procedure analogous to conditional marginal density estimation (Gelfand & Smith, 1990; Morey et al., 2011). The null hypotheses of no change for each subject are interval null hypotheses, but by making the interval very small around zero, the interval null hypotheses closely resemble point null hypotheses.

Posterior distributions of the model parameters can also be obtained. Model parameters of interest include mu_i (individual overall mean), delta*sig2z_i (individual mean difference), sig2z_i (individual variance of random shocks), mu_delta*sig2z (population mean of the individual mean differences), musig2z (location parameter for log normal prior on individual variances of random shocks), sig2sig2z (scale parameter for log normal prior on individual variances of random shocks), and rho (auto-correlation)).

Cauchy priors are placed on the standardized individual mean differences delta_i and on their population mean mu_delta. The \code{r1} argument represents the r scale of the Cauchy prior for mu_delta. The sum of the \code{r1} and code{r2} arguments sets the r scale of the Cauchy priors for the delta_i parameters. A flat prior is placed on the individual overall means mu_i, and a log Normal(musig2z, sig2sig2z) prior on the individual variances of the random shocks. A Normal(mu_mu_sig2_z, sig2_mu_sig2_z) prior is put on the location parameter musig2z and an Inverse Gamma(.5,.5) prior is placed on the scale parameter sig2sig2z. The \code{mu_mu_sig2_z} and  \code{sig2_mu_sig2_z} arguments control the mean and variance, respectively, of the Normal(mu_mu_sig2_z, sig2_mu_sig2_z) prior. A Beta prior is placed on the auto-correlation rho, where the \code{beta_theta} argument controls the b parameter of this Beta prior.

The \code{trunc} argument can be used to truncate the right side of the Beta prior on the auto-correlation rho. It prevents the MCMC chain to obtain rho values close to 1, which might happen when time intervals between observations are large. 

At each iteration of the Gibbs sampler, the inverse of the auto-correlation matrix is estimated using Levinson Durbin recursion. The \code{Thres} argument sets the required precision of the inverse. When precision is lower than this threshold in one or more iterations, a warning is given. In case a negative value is assigned to the \code{Thres} argument, no check is performed. 

The \code{sdmetlogsig2z}, \code{sdmetmusig2z}, and \code{sdmettheta} arguments set the metropolis scales for the logsig2_z, mu_sig2_z, and theta parameters, respectively. The acceptance rates should be between .25 and .5 for most efficient estimation. Acceptance rates are returned when the thinng rate is 1.

The \code{return.BFnrchanged} argument controls whether the Bayes factor logbfAtleastX should be computed. The computation requires a loop going through 2^N iterations, which can take a lot of time when N is large. Due to memory limits this Bayes factor cannot be computed for N > 39.
}
\value{
A list containing the following:
\item{logbfSubject}{A matrix containing the MCMC estimates of the log Bayes factors logbfSubject for all subjects.}
\item{logbfPattern}{An MCMC estimate of the log Bayes factor logbfPattern.}
\item{logbfAtleastX}{Only returned if \code{return.BFnrchanged=TRUE} and number of subjects N < 40. A matrix containing the MCMC estimates of the log Bayes factors logbfAtleastX, for X = 1 to N.}
\item{chains}{Only returned if \code{return.chains=TRUE}. A matrix containing the MCMC chains for the model parameters.}
\item{converged}{Only returned if \code{raftery.check = TRUE}. An indicator value for convergence: equals 1 if the lengths of the thinned chains are larger than the lengths required according to raftery diagnostic, 0 otherwise.}
\item{accRates}{Only returned if \code{thin=1}. The Metropolis-Hastings acceptance rates.}
}
\references{
Gelfand, A. and Smith, A. F. M. (1990). Sampling based approaches to calculating marginal densities. Journal of the American Statistical Association, 85, pp. 398-409.
Morey, R. D., Rouder, J. N., Pratte, M. S., and Speckman, P. L. (2011). Using MCMC chain outputs to efficiently estimate Bayes factors. Journal of Mathematical Psychology, 55, pp. 368-378.
}
\author{
Rivka M. de Vries and Richard D. Morey
}
\seealso{
\code{\link{bayesMLtrend}}
}
\examples{
## Define data
y = rnorm(30) # observations
sub = c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5) # subject 
occ = c(0:5,0:5,0:5,0:5,0:5) # measurement occasion
int = c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2) # phase indicator


## Obtain log Bayes factors and chains
output = bayesMLint(y, sub, occ, int, iter=1000)
}
\keyword{htest}
\keyword{models}

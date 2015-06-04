
void sampleMiss(double *y, double B0, double B1, double *t, double sig2_e, double *Psi,
                int Nobs, int Nmiss, int *miss, double *ySample);

SEXP RMCTwoSampleAR(SEXP yR, SEXP NR, SEXP distMatR, SEXP tR, SEXP rscaleR, SEXP alphaThetaR, SEXP betaThetaR, SEXP iterationsR);

double MCTwoSampleAR(double *y, int N, int *distMat, double *t, double rscale, double alphaTheta, double betaTheta, int iterations);

SEXP RgibbsTwoSampleAR(SEXP yR, SEXP NR, SEXP tR, SEXP rscaleR, SEXP alphaThetaR, SEXP betaThetaR,
                      SEXP loAreaR, SEXP upAreaR, SEXP iterationsR, SEXP leftSidedR, SEXP sdmetR, SEXP missR, SEXP NmissR,
                      SEXP progressR,SEXP pBar, SEXP rho);

void gibbsTwoSampleAR(double *y, int N, double *t, double rscale, double alphaTheta, double betaTheta,
                      double loArea, double upArea, int iterations, int leftSided, double sdmet, double *chains,
                      double *postp, int progress, int *miss, int Nmiss, SEXP pBar, SEXP rho);
                      
double sampThetaAR(double theta, double mu, double delta, double sig2, double g, double *y, int N, double *t, double alphaTheta, double betaTheta , double sdmet);

double thetaLogLikeAR(double theta, double mu, double delta, double sig2, double g, double *y, int N, double *t, double alphaTheta, double betaTheta);




void sampleMisstrend(double *y, double B0, double B1, double B2, double B3, double *x1, double *x2, double *x3, double sig2_e, double *Psi,
                int Nobs, int Nmiss, int *miss, double *ySample);

double thetaLogLikeAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta);

double sampThetaAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta, double sdmet);

void gibbsTwoSampleAR_trend(double *y, int N, double *X, int p, double rscaleInt, double rscaleSlp,
                            double alphaTheta, double betaTheta, double loInt, double upInt, double loSlp, double upSlp,
                            int iterations, int leftSidedInt, int leftSidedSlp, double sdmet, double *chains, double *postp, double *postdens, int progress, int *miss, int Nmiss,
                            SEXP pBar, SEXP rho);

SEXP RgibbsTwoSampleAR_trend(SEXP yR, SEXP NR, SEXP XR, SEXP pR, SEXP rscaleIntR, SEXP rscaleSlpR, SEXP alphaThetaR, SEXP betaThetaR, 
                      SEXP loIntR, SEXP upIntR, SEXP loSlpR, SEXP upSlpR, SEXP iterationsR, SEXP leftSidedIntR, SEXP leftSidedSlpR,
                      SEXP sdmetR, SEXP missR, SEXP NmissR, SEXP progressR, SEXP pBar, SEXP rho);

SEXP MCAR_trend(SEXP yR, SEXP NR, SEXP NmissR, SEXP distMatR, SEXP alphaThetaR, SEXP betaThetaR, SEXP rsqIntR, SEXP rsqSlpR, SEXP X0R, SEXP X1R, SEXP iterationsR, SEXP progressR, SEXP pBar, SEXP rho);

void MCmargLogLikeAR_trend(double theta,  int Nmiss, int *distMat, double g1, double g2, double *y, int N, double alphaTheta, double betaTheta, double rsqInt, double rsqSlp, double *X0, double *X1, double *like);

SEXP MCmargLogLikeAR_trendR(SEXP thetaR, SEXP NmissR, SEXP distMatR, SEXP g1R, SEXP g2R, SEXP yR, SEXP NR, SEXP alphaThetaR, SEXP betaThetaR, SEXP rsqIntR, SEXP rsqSlpR, SEXP X0R, SEXP X1R);

SEXP MCnullMargLogLikeAR_trend(SEXP thetaR, SEXP NmissR, SEXP distMatR, SEXP yR, SEXP NR, SEXP alphaThetaR, SEXP betaThetaR, SEXP X0R);
void quadformMatrix(double *X, double *S, int N, int p, double *Ans, double coeff, int invert);

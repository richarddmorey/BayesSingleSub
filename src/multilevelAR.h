


double max(double x, double y);
void rev(double *v, double *vrev, int n);
void levinsonDurban_step1(double *T, int N, double *backwardVecs);
void levinsonDurban_step2(double *T, int N, double *backwardVecs, double *solVecsN, int col);
void invertPsi(double theta, int N, double *distMat, double *invPsi, double thres, int *n);
double logPhiAminusPhiB(double a, double b);
double expXphiAminusphiB(double x, double a, double b, int returnlog);
int nrmatch(int x, int *y, int len); 
double log1mp(double logp);
double findlogPriorp(int iter, int N, int *test, double a, double b, double r1, double r2);
double findlogPriorpIT(int iter, int N, int *test, double aB1, double bB1, double aB3, double bB3, double r1, double r2, double r3, double r4); 



double funclogsig2e(double logsig2_e, double *Y, double B0, double B1,
                   double *t, double *invPsi, double mu_sig2_e, double sig2_sig2_e, int M);
                   
double funcMuSig2e(double mu_sig2_e, double mu_B1, double *logsig2_e,
                   double sig2_sig2_e, double mu_mu_sig2_e, double sig2_mu_sig2_e,
                 double g1, double g2, double *B1, int N);
                 
double funcTheta(double theta, int *M, int *Msqr, int NM, int NMsqr, int *cumM, int *cumMsqr, double *distMat, int N, double *Y, double *B0, double *B1, double *t,
                  double beta_theta, double *sig2_e, double thres, int *nr, double trunc);

double funcB0(double *Y, double B1, double *invPsi, double sig2_e, int M, double *t);

double funcB1(double *Y, double B0, double *invPsi, double *t,  double mu_B1,
              double sig2_e, double mu_sig2_e, double g2, int M, double *mu_B1i, double *sig_B1i);
              
double funcmuB1(double *B1, double mu_sig2_e, double g1, double g2, int N);

double funcsig2sig2e(double *logsig2_e, double mu_sig2_e, int N);

double funcg1(double r1, double mu_B1, double mu_sig2_e);

double funcg2(double *B1, double mu_B1, double mu_sig2_e, double r2, int N);

double metroplogsig2e(double logsig2_e, double *Y, double B0, double B1, double *t, double *invPsi,
                      double mu_sig2_e, double sig2_sig2_e, int M, double sdmet);
                      
double metropmusig2e(double mu_sig2_e, double mu_B1, double *logsig2_e, double sig2_sig2_e, double mu_mu_sig2_e, double sig2_mu_sig2_e, double g1, double g2, double *B1, int N, double sdmet);

double metropTheta(double theta, int *M, int *Msqr, int NM, int NMsqr, int *cumM, int *cumMsqr, double *distMat, int N, double *Y, double *B0, double *B1, double *t,
                  double beta_theta, double *sig2_e, double thres, int *nr, double trunc, double sdmet);
                  
void gibbsSampler(double *Y, int *M1, int *M2, int N, int iter, int thin, double r1, double r2, double beta_theta,
                    double a, double b, double thres,
                    double sdmetlogsig2e, double sdmetmusig2e, double sdmettheta,
                    double mu_mu_sig2_e, double sig2_mu_sig2_e, double *distMat, double trunc, int *test,
                    double *samples, double *BF, double *BFpatt, int progress, SEXP pBar, SEXP rho);
                    
SEXP gibbsSamplerCall(SEXP YR, SEXP M1R, SEXP M2R, SEXP NR, SEXP testR, SEXP iterR, SEXP thinR, SEXP r1R, SEXP r2R,
                      SEXP beta_thetaR, SEXP truncR, SEXP aR, SEXP bR, SEXP thresR,
                      SEXP sdmetlogsig2eR, SEXP sdmetmusig2eR, SEXP sdmetthetaR,
                      SEXP mu_mu_sig2_eR, SEXP sig2_mu_sig2_eR, SEXP distMatR,
                      SEXP progressR, SEXP pBar, SEXP rho);
                      
void BayesFactorsTrend(double aB1, double bB1, double aB3, double bB3, double *muB1, double *sig2B1, double *muB3, double *sig2B3, double *mu_sig2_e,
                        int *test, int iter, int N, double r1, double r2, double r3, double r4, double log_prioroddsIT,
                        double *BFint, double *BFtrend, double *BFit, double *BFpatt);
                        
SEXP BayesFactorsTrendCall(SEXP aB1R, SEXP bB1R, SEXP aB3R, SEXP bB3R, SEXP muB1R, SEXP sig2B1R, SEXP muB3R, SEXP sig2B3R,
                                SEXP mu_sig2_eR, SEXP testR, SEXP iterR, SEXP NR, SEXP r1R, SEXP r2R, SEXP r3R, SEXP r4R,
                                SEXP log_prioroddsITR);
                                
SEXP log1mpCall(SEXP logpR);

SEXP findlogPriorpCall(SEXP iterR, SEXP NR, SEXP testR, SEXP aR, SEXP bR, SEXP r1R, SEXP r2R);

SEXP findlogPriorpITCall(SEXP iterR, SEXP NR, SEXP testR, SEXP aB1R, SEXP bB1R, SEXP aB3R, SEXP bB3R, 
                            SEXP r1R, SEXP r2R, SEXP r3R, SEXP r4R);
                            
SEXP expXphiAminusphiBCall(SEXP xR, SEXP aR, SEXP bR, SEXP returnLogR);

SEXP matrixDetCall(SEXP AR, SEXP NR, SEXP returnLogR);


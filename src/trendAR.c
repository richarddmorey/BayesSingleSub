#include "utility.h"
#include "trendAR.h"


void sampleMisstrend(double *y, double B0, double B1, double B2, double B3, double *x1, double *x2, double *x3, double sig2_e, double *Psi,
                int Nobs, int Nmiss, int *miss, double *ySample){

  int i=0, j=0, l=0, m=0, a=0, b=0, c=0, d=0, one=1;
  int N = Nobs+Nmiss, Nyy = Nobs*Nobs, Nyym = Nobs*Nmiss, Nymym = Nmiss*Nmiss;
  double muym[Nmiss];
  double Syym[Nyym], Symym[Nymym], invSyy[Nyy], Symy[Nyym];
  double x[Nobs], C[Nyym], v[Nmiss], D[Nymym];
  double alpha = 1, beta=0;
  double SigmaCond[Nymym];


  // Construct partition of mu_y into mu_y_obs and mu_y_miss

  a=0;
  l=0;
  for(j=0;j<N;j++){

     if(j != miss[l]){
           x[a] =  y[j] - B0 - B1*x1[j] - B2*x2[j] - B3*x3[j];
           a+=1;
     }

     if(j == miss[l]){
         muym[l] =  B0 + B1*x1[j] + B2*x2[j] + B3*x3[j];
         l= minc((l+1),(Nmiss-1));
      }

  }


  //Construct partition of Sigma

  a=0;b=0; c=0;d=0; m=0;
  for(j=0;j<N;j++){
    l=0;
     for(i=0;i<N;i++){

       if(j == miss[m] && i != miss[l]){
      Syym[a] = Psi[N*j+i]*sig2_e;
       a+=1;
     }

      if(j == miss[m] && i == miss[l]){
      Symym[b] = Psi[N*j+i]*sig2_e;
        l= minc((l+1),(Nmiss-1));
        b+=1;
      }

      if(j != miss[m] && i != miss[l]){
      invSyy[c] = Psi[N*j+i]*sig2_e;
       c+=1;
      }

      if(j != miss[m] && i == miss[l]){
      Symy[d] = Psi[N*j+i]*sig2_e;
        l= minc((l+1),(Nmiss-1));
        d+=1;
      }

    }
       if(j == miss[m]){
       m = minc((m+1),(Nmiss-1));
      }
  }

 inverse(invSyy, Nobs);


 F77_NAME(dgemm)("n","n",&Nmiss,&Nobs,&Nobs,&alpha,Symy,&Nmiss,invSyy,&Nobs,&beta,C,&Nmiss);
 F77_NAME(dgemm)("n","n",&Nmiss,&Nmiss,&Nobs,&alpha,C,&Nmiss,Syym,&Nobs,&beta,D,&Nmiss);
 F77_NAME(dgemv)("n",&Nmiss,&Nobs,&alpha,C,&Nmiss,x,&one,&beta,v,&one);

 // Conditional mean of missings, written to ySample
 for(j=0;j<Nmiss;j++){
 ySample[j] = muym[j] + v[j];
 }

 // Conditional Sigma of missings
  for(j=0;j<Nymym;j++){
  SigmaCond[j] = Symym[j] - D[j];
  }

 rmvGaussianC(ySample, SigmaCond, Nmiss);

 }


void gibbsTwoSampleAR_trend(double *y, int N, double *X, int p, double rscaleInt, double rscaleSlp,
                            double alphaTheta, double betaTheta, double loInt, double upInt, double loSlp, double upSlp,
                            int iterations, int leftSidedInt, int leftSidedSlp, double sdmet, double *chains,
                            double *postp, double *postdens, int progress, int *miss, int Nmiss,
                            SEXP pBar, SEXP rho)
{
  int i=0, j=0, m=0, l=0, Nsqr=N*N, Nobs=N-Nmiss, iOne=1,iTwo=2,iThree=3;
	double aSig2, bSig2, ag, bg, rscIntsq=rscaleInt*rscaleInt, rscSlpsq=rscaleSlp*rscaleSlp, alpha1, alpha2;
	double sig2=0, theta = 0, g1 = 1, g2=1;
	double dZero=0,dOne=1,dNegOne=-1;
	double dOneOverSig2;
  double ySample[Nmiss];
  double postpInt[iterations], postpSlp[iterations];

	double ldensFullRestrict,ldensSlpRestrict,ldensIntRestrict,lAreaSlp,lAreaInt;

	double tempV[N];
	double tempV2[N];
	double invPsi[Nsqr], Psi[Nsqr];
	AZERO(invPsi,Nsqr);
	double beta[p];
	AZERO(beta,p);
	double tempNbyP[N*p];
	double Sigma[p*p];
	double tempNby2[2*N];
	double temp2by2[2*2];
	double temp2by2_2[2*2];
	double temp2by1[2];
	double temp2by1_2[2];


	double Xg[2*N], X2[2*N], X1[N], X3[3*N], beta2[2], beta3[3], betag[2];
	double betaVar;
	double betaMean;


	Memcpy(Xg, X, N);
	Memcpy(Xg + N, X + 2*N, N);
	Memcpy(X2, X + N, N);
	Memcpy(X2 + N, X + 3*N, N);

	int npars = p + 4;

	for(i=0;i<N;i++)
	{
		beta[0] += y[i]/(N*1.0);
		sig2 += y[i]*y[i];
	}


	sig2 = (sig2 - N*beta[0]*beta[0])/(N*1.0-1);


	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

	GetRNGstate();

	for(m=0;m<iterations;m++)
	{

		R_CheckUserInterrupt();

		//Check progress
		if(progress && !((m+1)%progress)){
			pSampCounter[0]=m+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}


		// Build invPsi matrix
		invPsi[0] = 1;
		invPsi[N*N-1] = 1;
		invPsi[1] = -theta;
		invPsi[N] = -theta;

		for(i=1;i<(N-1);i++)
		{
			invPsi[i + N*i] = (1 + theta*theta);
			invPsi[i + N*(i+1)] = -theta;
			invPsi[(i+1) + N*i] = -theta;
		}

     // Sample missing values

    if(Nmiss > 0){

    for(i=0;i<N;i++){
    for(j=0;j<N;j++){
     Psi[i*N + j] = 1/(1-theta*theta) * pow(theta,abs(i-j));
      }
    }

  sampleMisstrend(y, beta[0], beta[1], beta[2], beta[3], &X[N], &X[2*N], &X[3*N], sig2, Psi,
                Nobs, Nmiss, miss, ySample);

      l=0;
     for(i=0;i<N;i++){

      if(i == miss[l]){
        y[i] = ySample[l];
        l= minc((l+1),(Nmiss-1));
  }
  }

  }


		dOneOverSig2 = 1/sig2;
		F77_NAME(dgemm)("N", "N", &N, &p, &N, &dOneOverSig2, invPsi, &N, X, &N, &dZero, tempNbyP, &N);

		F77_NAME(dgemm)("T", "N", &p, &p, &N, &dOne, X, &N, tempNbyP, &N, &dZero, Sigma, &p);

		Sigma[1 + p] += 1/(sig2 * g1);
		Sigma[3 + p*3] += 1/(sig2 * g2);

		InvMatrixUpper(Sigma, p);
		internal_symmetrize(Sigma, p);

		F77_NAME(dgemv)("T", &N, &p, &dOne, tempNbyP, &N, y, &iOne, &dZero, tempV, &iOne);
		F77_NAME(dgemv)("N", &p, &p, &dOne, Sigma, &p, tempV, &iOne, &dZero, tempV2, &iOne);

		rmvGaussianC(tempV2, Sigma, p);
		Memcpy(beta, tempV2, p);

		//densities

		//slope restricted
		Memcpy(X1, X + 3*N, N);
		Memcpy(X3, Xg, 2*N);
		Memcpy(X3 + 2*N, X + N, N);
		beta3[0] = beta[0];
		beta3[1] = beta[2];
		beta3[2] = beta[1];
		Memcpy(tempV,y,N);
		F77_NAME(dgemv)("N", &N, &iThree, &dNegOne, X3, &N, beta3, &iOne, &dOne, tempV, &iOne);
		betaVar = 0;
		betaMean = 0;
		betaVar = quadform(X1,invPsi,N,1,N);
		F77_NAME(dgemv)("N", &N, &N, &dOne, invPsi, &N, X1, &iOne, &dZero, tempV2, &iOne);

		for(i=0;i<N;i++)
			betaMean += tempV2[i] * tempV[i];

		betaVar = 1 / (betaVar + 1/g2);
		betaMean = betaVar * betaMean;

		ldensSlpRestrict = dnorm(0,betaMean/sqrt(sig2),sqrt(betaVar),1);
		lAreaSlp = pnorm(upSlp,betaMean/sqrt(sig2),sqrt(betaVar),1,0) - pnorm(loSlp,betaMean/sqrt(sig2),sqrt(betaVar),1,0);
    postpSlp[m] = pnorm(0,betaMean/sqrt(sig2),sqrt(betaVar), leftSidedSlp, 1);
    
		if(m==0)
		{
			postdens[1] = ldensSlpRestrict;
			postdens[4] = lAreaSlp;
		}else{
			postdens[1] = LogOnePlusExpX(ldensSlpRestrict-postdens[1])+postdens[1];
			postdens[4] += lAreaSlp;
		}



		//intercept restricted
		Memcpy(X1, X + N, N);
		Memcpy(X3 + 2*N, X + 3*N, N);
		beta3[0] = beta[0];
		beta3[1] = beta[2];
		beta3[2] = beta[3];
		Memcpy(tempV,y,N);
		F77_NAME(dgemv)("N", &N, &iThree, &dNegOne, X3, &N, beta3, &iOne, &dOne, tempV, &iOne);
		betaVar = 0;
		betaMean = 0;
		betaVar = quadform(X1,invPsi,N,1,N);
		F77_NAME(dgemv)("N", &N, &N, &dOne, invPsi, &N, X1, &iOne, &dZero, tempV2, &iOne);

		for(i=0;i<N;i++)
			betaMean += tempV2[i] * tempV[i];

		betaVar = 1 / (betaVar + 1/g1);
		betaMean = betaVar * betaMean;
		ldensIntRestrict = dnorm(0,betaMean/sqrt(sig2),sqrt(betaVar),1);
		lAreaInt = pnorm(upInt,betaMean/sqrt(sig2),sqrt(betaVar),1,0) - pnorm(loInt,betaMean/sqrt(sig2),sqrt(betaVar),1,0);
    postpInt[m] = pnorm(0,betaMean/sqrt(sig2),sqrt(betaVar), leftSidedInt, 1);

		if(m==0)
		{
			postdens[2] = ldensIntRestrict;
			postdens[3] = lAreaInt;
		}else{
			postdens[2] = LogOnePlusExpX(ldensIntRestrict-postdens[2])+postdens[2];
			postdens[3] += lAreaInt;
		}


		//Both restricted
		betag[0] = beta[0];
		betag[1] = beta[2];
		Memcpy(tempV,y,N);

		//(y - Xg%*%betag)
		F77_NAME(dgemv)("N", &N, &iTwo, &dNegOne, Xg, &N, betag, &iOne, &dOne, tempV, &iOne);

		// invPsi%*%X2
		F77_NAME(dgemm)("N", "N", &N, &iTwo, &N, &dOne, invPsi, &N, X2, &N, &dZero, tempNby2, &N);
		//t(X2)%*%invPsi%*%(y - Xg%*%betag)
		F77_NAME(dgemv)("T", &N, &iTwo, &dOne, tempNby2, &N, tempV, &iOne, &dZero, temp2by1, &iOne);

		// t(X2)%*%invPsi%*%X2
		F77_NAME(dgemm)("T", "N", &iTwo, &iTwo, &N, &dOne, X2, &N, tempNby2, &N, &dZero, temp2by2, &iTwo);

		temp2by2[0] += 1/g1;
		temp2by2[3] += 1/g2;

		Memcpy(temp2by2_2,temp2by2,4);

		InvMatrixUpper(temp2by2, 2);
		internal_symmetrize(temp2by2, 2);

		dOneOverSig2 = 1/sqrt(sig2);

		//Sigma%*%t(X2)%*%invPsi%*%X2
		F77_NAME(dgemv)("N", &iTwo, &iTwo, &dOneOverSig2, temp2by2, &iTwo, temp2by1, &iOne, &dZero, temp2by1_2, &iOne);

		ldensFullRestrict = -log(2 * M_PI) - 0.5*matrixDet(temp2by2,2,2,1) - 0.5*quadform(temp2by1_2,temp2by2_2,2,1,2);

		if(m==0)
		{
			postdens[0] = ldensFullRestrict;
		}else{
			postdens[0] = LogOnePlusExpX(ldensFullRestrict-postdens[0])+postdens[0];
		}

		//sig2
		Memcpy(tempV,y,N);
		F77_NAME(dgemv)("N", &N, &p, &dNegOne, X, &N, beta, &iOne, &dOne, tempV, &iOne);
		aSig2 = 0.5*(N+2);
		bSig2 = 0.5*(quadform(tempV,invPsi,N,1,N) + beta[1]*beta[1]/g1 + beta[3]*beta[3]/g2);

		sig2 = 1/rgamma(aSig2,1/bSig2);

		//g1
		ag = 1;
		bg = 0.5*(beta[1]*beta[1]/sig2 + rscIntsq);
		g1 = 1/rgamma(ag,1/bg);

		//g2
		ag = 1;
		bg = 0.5*(beta[3]*beta[3]/sig2 + rscSlpsq);
		g2 = 1/rgamma(ag,1/bg);

		//theta
		theta = sampThetaAR_trend(theta, beta, X, sig2, y, N, p, alphaTheta, betaTheta, sdmet);

		// write chain
		Memcpy(chains + m*npars, beta, p);
		chains[p + m*npars] = sig2;
		chains[1 + p + m*npars] = g1;
		chains[2 + p + m*npars] =  g2;
		chains[3 + p + m*npars] = theta;
		//chains[4 + p + m*npars] = ldensFullRestrict;
		//chains[5 + p + m*npars] = ldensSlpRestrict;
		//chains[6 + p + m*npars] = ldensIntRestrict;
		//chains[7 + p + m*npars] = lAreaInt;
		//chains[8 + p + m*npars] = lAreaSlp;

	}
  
    postp[0] = logMeanExpLogs(postpInt, iterations);
    postp[1] = logMeanExpLogs(postpSlp, iterations);

	UNPROTECT(2);
	PutRNGstate();
}


SEXP RgibbsTwoSampleAR_trend(SEXP yR, SEXP NR, SEXP XR, SEXP pR, SEXP rscaleIntR, SEXP rscaleSlpR, SEXP alphaThetaR, SEXP betaThetaR, 
                            SEXP loIntR, SEXP upIntR, SEXP loSlpR, SEXP upSlpR, SEXP iterationsR, 
                            SEXP leftSidedIntR, SEXP leftSidedSlpR, SEXP sdmetR, SEXP missR, SEXP NmissR, SEXP progressR, SEXP pBar, SEXP rho)
{
  int iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscaleInt = NUMERIC_VALUE(rscaleIntR);
	double rscaleSlp = NUMERIC_VALUE(rscaleSlpR);
	double alphaTheta = NUMERIC_VALUE(alphaThetaR);
	double betaTheta = NUMERIC_VALUE(betaThetaR);
	double sdmet = NUMERIC_VALUE(sdmetR);
	double *y = REAL(yR);
	double *X = REAL(XR);
	int p = INTEGER_VALUE(pR);
	int npars = p + 4;
  int leftSidedInt = INTEGER_VALUE(leftSidedIntR);
  int leftSidedSlp = INTEGER_VALUE(leftSidedSlpR);

  int *miss = INTEGER_POINTER(missR);
  int Nmiss = INTEGER_VALUE(NmissR);

	double loInt = REAL(loIntR)[0];
	double upInt = REAL(upIntR)[0];
	double loSlp = REAL(loSlpR)[0];
	double upSlp = REAL(upSlpR)[0];


	SEXP chainsR;
  SEXP postpR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));
  PROTECT(postpR = allocVector(REALSXP, 2));

	// means
	SEXP postdensR;
	PROTECT(postdensR = allocMatrix(REALSXP, 1, 5));
	SEXP returnList;
	PROTECT(returnList = allocVector(VECSXP, 3));


	gibbsTwoSampleAR_trend(y, N, X, p, rscaleInt, rscaleSlp, alphaTheta, betaTheta, loInt, upInt, loSlp, upSlp, iterations, 
                        leftSidedInt, leftSidedSlp, sdmet, REAL(chainsR), REAL(postpR), REAL(postdensR), progress, miss, Nmiss, pBar, rho);

	SET_VECTOR_ELT(returnList, 0, chainsR);
  SET_VECTOR_ELT(returnList, 1, postdensR);
  SET_VECTOR_ELT(returnList, 2, postpR);

	UNPROTECT(4);

	return(returnList);
}

double sampThetaAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta, double sdmet)
{
	// sample theta with Metropolis-Hastings
	double cand, likeRat, b;

	cand = theta + rnorm(0,sdmet);

	if(cand<0 || cand>1)
	{
		return(theta);
	}
	b = log(runif(0,1));
	likeRat = thetaLogLikeAR_trend(cand, beta, X, sig2, y, N, p, alphaTheta, betaTheta) - thetaLogLikeAR_trend(theta, beta, X, sig2, y, N, p, alphaTheta, betaTheta);

	if(b>likeRat){
		return(theta);
	}else{
		return(cand);
	}
}

double thetaLogLikeAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta)
{
	int i,iOne=1;
	double loglike=0,tempV[N],invPsi[N*N],det,dNegOne=-1,dOne=1;

	Memcpy(tempV,y,N);
	F77_NAME(dgemv)("N", &N, &p, &dNegOne, X, &N, beta, &iOne, &dOne, tempV, &iOne);

	AZERO(invPsi,N*N);

	// Build invPsi matrix
	invPsi[0] = 1;
	invPsi[N*N-1] = 1;
	invPsi[1] = -theta;
	invPsi[N] = -theta;
	for(i=0;i<N;i++)
	{

		if(i>0 && i<(N-1)){
			invPsi[i + N*i] = (1 + theta*theta);
			invPsi[i + N*(i+1)] = -theta;
			invPsi[(i+1) + N*i] = -theta;
		}
	}
	det = log(1-theta*theta);

	loglike = 0.5 * det - 0.5/sig2 * quadform(tempV,invPsi,N,1,N) + (alphaTheta-1)*log(theta) + (betaTheta-1)*log(1-theta);

	return(loglike);
}

SEXP MCAR_trend(SEXP yR, SEXP NR, SEXP NmissR, SEXP distMatR, SEXP alphaThetaR, SEXP betaThetaR, SEXP rsqIntR, SEXP rsqSlpR, SEXP X0R, SEXP X1R, SEXP iterationsR, SEXP progressR, SEXP pBar, SEXP rho){

	int m,i, iterations = INTEGER_VALUE(iterationsR), N = INTEGER_VALUE(NR);
	int progress = INTEGER_VALUE(progressR);
  int Nmiss = INTEGER_VALUE(NmissR);
  int *distMat = INTEGER_POINTER(distMatR);
	double loglike[3], logsum[3], g1, g2, theta;
	double alphaTheta = REAL(alphaThetaR)[0];
	double betaTheta = REAL(betaThetaR)[0];
	double rsqInt = REAL(rsqIntR)[0];
	double rsqSlp = REAL(rsqSlpR)[0];

	SEXP returnR;
	PROTECT(returnR = allocVector(VECSXP, 2));

	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, 3, iterations));

	SEXP logmeanR;
	PROTECT(logmeanR = allocVector(REALSXP, 3));

	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

	GetRNGstate();

	for(m=0;m<iterations;m++){

		R_CheckUserInterrupt();

		//Check progress
		if(progress && !((m+1)%progress)){
			pSampCounter[0]=m+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}

		// sample
		g1 = 1/rgamma(0.5,2/rsqInt);
		g2 = 1/rgamma(0.5,2/rsqSlp);
		theta = rbeta(alphaTheta,betaTheta);

		MCmargLogLikeAR_trend(theta, Nmiss, distMat, g1, g2, REAL(yR), N, alphaTheta, betaTheta, rsqInt, rsqSlp, REAL(X0R), REAL(X1R), loglike);

		if(m==0)
		{
			Memcpy(logsum,loglike,3);
		}else{
			for(i=0;i<3;i++)
				logsum[i] = LogOnePlusExpX(loglike[i]-logsum[i])+logsum[i];
		}
		Memcpy(REAL(chainsR) + m*3, loglike,3);

	}

	for(i=0;i<3;i++)
		REAL(logmeanR)[i] = logsum[i] - log(iterations);

	SET_VECTOR_ELT(returnR, 0, logmeanR);
    SET_VECTOR_ELT(returnR, 1, chainsR);


	PutRNGstate();
	UNPROTECT(5);

	return(returnR);

}

SEXP MCmargLogLikeAR_trendR(SEXP thetaR, SEXP NmissR, SEXP distMatR, SEXP g1R, SEXP g2R, SEXP yR, SEXP NR, SEXP alphaThetaR, SEXP betaThetaR, SEXP rsqIntR, SEXP rsqSlpR, SEXP X0R, SEXP X1R){
	double theta = REAL(thetaR)[0], g1 = REAL(g1R)[0], g2 = REAL(g2R)[0];
	double alphaTheta = REAL(alphaThetaR)[0], betaTheta = REAL(betaThetaR)[0];
	int N = INTEGER_VALUE(NR);
  int Nmiss = INTEGER_VALUE(NmissR);
  int *distMat = INTEGER_POINTER(distMatR);
	double rsqInt = REAL(rsqIntR)[0],rsqSlp = REAL(rsqSlpR)[0];

	SEXP loglikeR;
	PROTECT(loglikeR = allocVector(REALSXP, 3));

	MCmargLogLikeAR_trend(theta, Nmiss, distMat, g1, g2, REAL(yR), N, alphaTheta, betaTheta, rsqInt, rsqSlp, REAL(X0R), REAL(X1R), REAL(loglikeR));

	UNPROTECT(1);
	return(loglikeR);
}

void MCmargLogLikeAR_trend(double theta,  int Nmiss, int *distMat, double g1, double g2, double *y, int N, double alphaTheta, double betaTheta, double rsqInt, double rsqSlp, double *X0, double *X1, double *like){

	int i,j,iOne=1,iTwo=2;
	double invPsi[N*N],detInvPsi,dOne=1, det1, det2, dZero=0, quad;

	double tempM1[2*2];
	double tempM2[2*N];
	double tempS1, tempS2;
	double tempV1[N];
	double tempV2[N];
	double Z1[N*N];
	double Z2[N*N];

	AZERO(invPsi,N*N);

  if(Nmiss == 0){
	// Build invPsi matrix
	invPsi[0] = 1;
	invPsi[N*N-1] = 1;
	invPsi[1] = -theta;
	invPsi[N] = -theta;
	for(i=0;i<N;i++)
	{

		if(i>0 && i<(N-1)){
			invPsi[i + N*i] = (1 + theta*theta);
			invPsi[i + N*(i+1)] = -theta;
			invPsi[(i+1) + N*i] = -theta;
		}
	}

	detInvPsi = log(1-theta*theta);
}

if(Nmiss > 0){
 
 for(i=0;i<N;i++)
{
invPsi[i + N*i] = 1/(1-theta*theta);
for(j=0;j<i;j++)
{
invPsi[i + N*j] = invPsi[i + N*i] * pow(theta,distMat[i + N*j]);  
invPsi[j + N*i] = invPsi[i + N*j];
}
}

InvMatrixUpper(invPsi, N);
internal_symmetrize(invPsi, N);
  
 detInvPsi = matrixDet(invPsi, N, N, 1);
}

	// Integral of beta0
	quadformMatrix(X0,invPsi,N,2,tempM1,1,1);
	//det1 = matrixDet(tempM1,2,2,1);
	det1 = log(tempM1[0]*tempM1[3] - tempM1[1]*tempM1[2]);


	F77_NAME(dgemm)("T", "N", &iTwo, &N, &N, &dOne, X0, &N, invPsi, &N, &dZero, tempM2, &iTwo);
	quadformMatrix(tempM2,tempM1,2,N,Z1,-1,0);
	for(i=0;i<(N*N);i++){
		Z1[i] += invPsi[i];
	}


	//Integral of beta1, slope restricted model
	tempS1 = quadform2(X1,Z1,N,1,N);
	tempS1 += 1/g1;
	tempS1 = 1/tempS1;

	det2 = log(tempS1);

	F77_NAME(dgemv)("T", &N, &N, &dOne, Z1, &N, X1, &iOne, &dZero, tempV1, &iOne);

	//outer product of tempV
	for(j=0;j<N;j++){
		Z2[j + j*N] = Z1[j + j*N] - tempV1[j]*tempV1[j]*tempS1;
		for(i=0;i<j;i++){
			Z2[i + j*N] = Z1[i + j*N] - tempV1[i]*tempV1[j]*tempS1;
			Z2[j + i*N] = Z2[i + j*N];
		}
	}
	quad = quadform2(y,Z2,N,1,N);


	like[1] = lgamma( (N-2) * 0.5) + 0.5 * detInvPsi +
	 -( (N-2) * 0.5) * log( quad ) +
	 0.5 * det1 + 0.5 * det2 + -0.5*log(g1);

	//Integral of beta1, intercept restricted model
	tempS1 = quadform2(X1 + N,Z1,N,1,N);
	tempS1 += 1/g2;
	tempS1 = 1/tempS1;

	det2 = log(tempS1);

	F77_NAME(dgemv)("T", &N, &N, &dOne, Z1, &N, X1+N, &iOne, &dZero, tempV1, &iOne);

	//outer product of tempV
	for(j=0;j<N;j++){
		Z2[j + j*N] = Z1[j + j*N] - tempV1[j]*tempV1[j]*tempS1;
		for(i=0;i<j;i++){
			Z2[i + j*N] = Z1[i + j*N] - tempV1[i]*tempV1[j]*tempS1;
			Z2[j + i*N] = Z2[i + j*N];
		}
	}

	quad = quadform2(y,Z2,N,1,N);

	like[2] = lgamma( (N-2) * 0.5) + 0.5 * detInvPsi +
	 -( (N-2) * 0.5) * log( quad ) +
	 0.5 * det1 + 0.5 * det2 + -0.5*log(g2);


	//Integral of beta1, general model
	quadformMatrix(X1,Z1,N,2,tempM1,1,0);
	tempM1[0] += 1/g1;
	tempM1[3] += 1/g2;

	InvMatrixUpper(tempM1, 2);
	internal_symmetrize(tempM1, 2);

	//det2 = matrixDet(tempM1,2,2,1);
	det2 = log(tempM1[0]*tempM1[3] - tempM1[1]*tempM1[2]);

	F77_NAME(dgemm)("T", "N", &iTwo, &N, &N, &dOne, X1, &N, Z1, &N, &dZero, tempM2, &iTwo);
	quadformMatrix(tempM2,tempM1,2,N,Z2,-1,0);
	for(i=0;i<(N*N);i++){
		Z2[i] += Z1[i];
	}

	/*int j;
	for(j=0;j<N;j++){
		for(i=0;i<N;i++)
			Rprintf("%f ",Z2[j*N+i]);
		Rprintf("\n");
	}*/

	quad = quadform2(y,Z2,N,1,N);

	like[0] = lgamma( (N-2) * 0.5) + 0.5 * detInvPsi +
	 -( (N-2) * 0.5) * log( quad ) +
	 0.5 * det1 + 0.5 * det2 + -0.5*log(g1) + -0.5*log(g2);


}

SEXP MCnullMargLogLikeAR_trend(SEXP thetaR, SEXP NmissR, SEXP distMatR, SEXP yR, SEXP NR, SEXP alphaThetaR, SEXP betaThetaR, SEXP X0R){

	double theta = REAL(thetaR)[0], *y = REAL(yR);
	int N = INTEGER_VALUE(NR);
  int Nmiss = INTEGER_VALUE(NmissR);
  int *distMat = INTEGER_POINTER(distMatR);
	double alphaTheta = REAL(alphaThetaR)[0];
	double betaTheta = REAL(betaThetaR)[0];
	double *X0 = REAL(X0R);

	int i, j=0,iOne=1,iTwo=2;
	double invPsi[N*N],detInvPsi,dOne=1, det1, dZero=0;

	double tempM1[2*2];
	double tempM2[2*N];
	double Z1[N*N];

	SEXP returnR;
	PROTECT(returnR = allocVector(REALSXP, 1));

  AZERO(invPsi,N*N);
  
	// Build invPsi matrix
  
  if(Nmiss == 0){
    
	invPsi[0] = 1;
	invPsi[N*N-1] = 1;
	invPsi[1] = -theta;
	invPsi[N] = -theta;
	for(i=0;i<N;i++)
	{

		if(i>0 && i<(N-1)){
			invPsi[i + N*i] = (1 + theta*theta);
			invPsi[i + N*(i+1)] = -theta;
			invPsi[(i+1) + N*i] = -theta;
		}
	}

	detInvPsi = log(1-theta*theta);
  
  }
  
  if(Nmiss > 0){
 
  for(i=0;i<N;i++)
  {
  invPsi[i + N*i] = 1/(1-theta*theta);
  for(j=0;j<i;j++)
  {
  invPsi[i + N*j] = invPsi[i + N*i] * pow(theta,distMat[i + N*j]);  
  invPsi[j + N*i] = invPsi[i + N*j];
  }
  }

  InvMatrixUpper(invPsi, N);
  internal_symmetrize(invPsi, N);
  
  detInvPsi = matrixDet(invPsi, N, N, 1);
  }


	// Integral of beta0
	quadformMatrix(X0,invPsi,N,2,tempM1,1,1);
	//det1 = matrixDet(tempM1,2,2,1);
	det1 = log(tempM1[0]*tempM1[3] - tempM1[1]*tempM1[2]);

	F77_NAME(dgemm)("T", "N", &iTwo, &N, &N, &dOne, X0, &N, invPsi, &N, &dZero, tempM2, &iTwo);

	quadformMatrix(tempM2,tempM1,2,N,Z1,-1,0);

	for(i=0;i<(N*N);i++){
		Z1[i] += invPsi[i];
	}

	REAL(returnR)[0] = lgamma((N-2) * 0.5) + 0.5 * detInvPsi +
	 -((N-2) * 0.5) * log( quadform2(y,Z1,N,1,N) ) +
	 0.5 * det1 + dbeta(theta,alphaTheta,betaTheta,1);

	UNPROTECT(1);
	return(returnR);

}


void quadformMatrix(double *X, double *S, int N, int p, double *Ans, double coeff, int invert){

	double tempM[N*p],dOne=1,dZero=0;

	F77_NAME(dgemm)("T", "N", &p, &N, &N, &dOne, X, &N, S, &N, &dZero, tempM, &p);
	F77_NAME(dgemm)("N", "N", &p, &p, &N, &coeff, tempM, &p, X, &N, &dZero, Ans, &p);

	if(invert)
	{
		InvMatrixUpper(Ans, p);
		internal_symmetrize(Ans, p);
	}

}

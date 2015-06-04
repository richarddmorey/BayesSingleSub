#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include <Rversion.h>
#include <Rconfig.h>
#include <R_ext/Constants.h>
#include <R_ext/Random.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include "utility.h"
#include "multilevelAR.h"



// max

double max(double x, double y){

  if((x) < (y)){
    return(y);
  }
  else{
    return(x);
  }
}

// reverse elements of vector v of length n; output written to vrev
void rev(double *v, double *vrev, int n){

   int i=0, j=(n);

   for(i=0;i<n;i++){
    j += -1;
    vrev[i] = v[j];
   }
}


// Compute inverse of Psi

// Step 1: Retrieve the backward vectors. Does NOT depend on y.
void levinsonDurban_step1(double *T, int N, double *backwardVecs){
  
  int i=1,j=0, k=0, m=0, q=0;
  double errori=0; 
  double tempV1[N], tempV2[N];
  
  backwardVecs[0] = 1;
  
  for(i=1;i<N;i++){
    
    errori = 0;
    for(j=0;j<i;j++){
      errori += T[j+1]*backwardVecs[k];
       k+=1;
      }
    
   tempV1[0] = 0;
    for(j=0;j<i;j++){
     tempV1[j+1] = backwardVecs[m];
     m+=1;
    }
  
    rev(tempV1, tempV2, (i+1));
  
   q=0;
    for(j=k;j<(k+i+1);j++){     
     backwardVecs[j] =  (1/(1-errori*errori))*tempV1[q] - (errori/(1-errori*errori))*tempV2[q];
     q+=1;
    }
    
  }

}

// Step 2: Find column 'col' of the inverse:
void levinsonDurban_step2(double *T, int N, double *backwardVecs, double *solVecsN, int col){
  
  int *y = Calloc(N,int);
  int i=1,j=0, k=0, q=0;
  double errori = 0;
  double tempV1[N];
  double solVecs[(N*(N+1)/2)];
  
  y[col] = 1;
  
  //Retrieve solution x for equation y = T*x:
  
  solVecs[0] = y[0];
  
  for(i=1;i<N;i++){
    
    errori = 0;
    tempV1[i] = 0;
    for(j=0;j<i;j++){
      errori += T[i*N+j]*solVecs[k];
      tempV1[j] = solVecs[k];
    
      //Rprintf("%d\n", solVecs[k]); 
      
       k+=1;
      }
    
   q=0;
    for(j=k;j<(k+i+1);j++){     
     solVecs[j] =  tempV1[q] + (y[i]-errori)*backwardVecs[j];
     q+=1;
    }    
  }
  
  for(j=0;j<N;j++){  
  solVecsN[j] = solVecs[(N*(N-1)/2)+j];
  }
  
  Free(y);
  
}


void invertPsi(double theta, int N, double *distMat, double *invPsi, double thres, int *n){
  
  int j=0, i=0, Nsqr=N*N;
  double backwardVecs[N*(N+1)/2], solVecsN[N], T[Nsqr];
  

  for(j=0;j<Nsqr;j++){
      T[j] = pow(theta,distMat[j]);       
    }

  
  levinsonDurban_step1(T, N, backwardVecs);
  
  for(i=0;i<N;i++){  
  
  levinsonDurban_step2(T, N, backwardVecs, solVecsN, i);
  
    for(j=0;j<N;j++){  
      
      invPsi[i*N+j] = solVecsN[j] * (1-theta*theta);
  
    }
  
  }
  
  
  // If thres > 0, check whether Psi %*% invPsi = I 
  if(thres>0){
    
    int i=0, j=0, *I = Calloc(Nsqr, int);
    double l, one=1, zero=0, PsiPsim1[Nsqr], Psi[Nsqr];
    
    l = (1-theta*theta);
    for(j=0;j<Nsqr;j++){
      Psi[j] = pow(theta,distMat[j])/l;
    }
    
    F77_CALL(dgemm)("N", "N", &N, &N, &N, &one, Psi, &N, invPsi, &N, &zero, PsiPsim1, &N);
    
    
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        I[i*N+j] = 1;
      }
    }
    
    for(j=0;j<Nsqr;j++){
      
       if((PsiPsim1[j] - I[j]) > thres){
        n[0]+=1;
        break;
      }
    }
    
    Free(I);
    
   }
  

  
}



// Compute log area under normal curve

double logPhiAminusPhiB(double a, double b){
  int i;
  double c0=.2316419;
  double c[5]={.319381530,-.356563782,1.781477937,-1.821255978,1.330274429};
  double da=0,db=0,ta,tb,fa,fb,g,lza,lzb;
  if(a<0||b<0){
    Rprintf("\n****Tried to use approximation inappropriately. a=%f, b=%f\n",a,b);
   }
  ta=1.0/(1.0+c0*a);
  tb=1.0/(1.0+c0*b);
  lza=-pow(a,2)/2;
  lzb=-pow(b,2)/2;
  for(i=0;i<5;i++){
    da+=c[i]*pow(ta,i+1);
    db+=c[i]*pow(tb,i+1);
  }
  fa=exp(lza+log(da));
  fb=exp(lzb+log(db));
  g=log(fb-fa);
  return(-.5*log(2*M_PI)+g);
}


double expXphiAminusphiB(double x, double a, double b, int returnlog){
  double ret;
  switch(abs(a)>5&&abs(b)>5){
  case 1:
    if(a>5&&b>5){
      ret=x+logPhiAminusPhiB(a,b);
    }else if(a<-5&&b<-5){
      ret=x+logPhiAminusPhiB(-b,-a);
    }else{
      ret=x+log(pnorm(a,0,1,1,0)-pnorm(b,0,1,1,0));
    }
    break;
  default:
    ret=x+log(pnorm(a,0,1,1,0)-pnorm(b,0,1,1,0));
    break;
  }
  if(!returnlog){
    return(exp(ret));
  }else{
    return(ret);
  }
}



// Match function; n indicates how often x occures in y

int nrmatch(int x, int *y, int len){

  int i=0, n=0;

  for(i=0;i<len;i++){

    if((x) == (y[i])){
      n += 1;
    }
  }

  return(n);

}



// Compute log(1-p)

double log1mp(double logp){

   if ((logp) > -M_LN2){
     return(log(-expm1(logp)));
   }

   else{
     return(log1p(-exp(logp)));

   }
}






double findlogPriorp(int iter, int N, int *test, double a, double b, double r1, double r2){
  
 int i=0, j=0, n1=0, n0=0;
 double clogprp0, clogprp1, *logprpV = Calloc(iter,double), logprp, g1, g2;
 
 for(i=0;i<N;i++){
  n1 += test[i];
 }
 
 n0 = N-n1;
 
 for(j=0;j<iter;j++){
  
  g1 = 1.0/rgamma(.5,2/(r1*r1));
  g2 = 1.0/rgamma(.5,2/(r2*r2));
  
  clogprp0 = max(log(pnorm(b,0,sqrt(g1+g2), 1, 0) - pnorm(a,0,sqrt(g1+g2), 1, 0)), -1000000);
  clogprp1 = max(log1mp(clogprp0), -1000000);
  
  logprpV[j] = n1*clogprp1 + n0*clogprp0;
}

logprp = logMeanExpLogs(logprpV, iter);

Free(logprpV);

return(logprp);

}

double findlogPriorpIT(int iter, int N, int *test, double aB1, double bB1, double aB3, double bB3, double r1, double r2, double r3, double r4){
  
 int i=0, j=0, n1=0, n0=0;
 double clogprp0, clogprp1, *logprpV = Calloc(iter,double), logprp, g1, g2, g3, g4;
 
 for(i=0;i<N;i++){
  n1 += test[i];
 }
 
 n0 = N-n1;
 
 for(j=0;j<iter;j++){
  
  g1 = 1.0/rgamma(.5,2/(r1*r1));
  g2 = 1.0/rgamma(.5,2/(r2*r2));
  g3 = 1.0/rgamma(.5,2/(r3*r3));
  g4 = 1.0/rgamma(.5,2/(r4*r4));
  
  clogprp0 = max(
            (
            log(pnorm(bB1,0,sqrt(g1+g2), 1, 0) - pnorm(aB1,0,sqrt(g1+g2), 1, 0)) + 
            log(pnorm(bB3,0,sqrt(g3+g4), 1, 0) - pnorm(aB3,0,sqrt(g3+g4), 1, 0))
            ),
            -1000000);
  clogprp1 = max(
            (
            log(pnorm(aB1,0,sqrt(g1+g2), 1, 0) + pnorm(bB1,0,sqrt(g1+g2), 0, 0)) + 
            log(pnorm(aB3,0,sqrt(g3+g4), 1, 0) + pnorm(bB3,0,sqrt(g3+g4), 0, 0))
            ),
            -1000000);
  
  logprpV[j] = n1*clogprp1 + n0*clogprp0;
}

logprp = logMeanExpLogs(logprpV, iter);

Free(logprpV);

return(logprp);

}


// /////////////////////////////////////////////////////////////////////////////////////



// Full conditional for logsig2e

double funclogsig2e(double logsig2_e, double *Y, double B0, double B1,
                   double *t, double *invPsi, double mu_sig2_e, double sig2_sig2_e, int M)
{
  int j=0, incx = 1;
  double deviations[M], logd;

         for(j=0;j<M;j++){

    deviations[j] = Y[j] - B0 - B1*t[j];

    }

 logd = -M*logsig2_e/2 - quadform2(deviations, invPsi, M, incx, M)/(2*exp(logsig2_e)) -
        (logsig2_e - mu_sig2_e)*(logsig2_e - mu_sig2_e)/(2*sig2_sig2_e);

return(logd);

}




// Full conditional for mu_sig2_e
double funcMuSig2e(double mu_sig2_e, double mu_B1, double *logsig2_e,
                   double sig2_sig2_e, double mu_mu_sig2_e, double sig2_mu_sig2_e,
                 double g1, double g2, double *B1, int N)
{
  double logd, sumlogsig2_e=0, sumsqdev = 0;
  int i=0;

  for(i=0;i<N;i++){
 sumlogsig2_e += logsig2_e[i];
 sumsqdev += (B1[i]-mu_B1)*(B1[i]-mu_B1);
  }

  logd = (mu_mu_sig2_e/sig2_mu_sig2_e + sumlogsig2_e/sig2_sig2_e  - (N+1)*.5 )*mu_sig2_e -
          (g1*sumsqdev + g2*mu_B1*mu_B1)/(2*g1*g2*exp(mu_sig2_e)) -
          .5*(N/sig2_sig2_e + 1/sig2_mu_sig2_e)*mu_sig2_e*mu_sig2_e;

  return(logd);
}



// Full conditional for theta
double funcTheta(double theta, int *M, int *Msqr, int NM, int NMsqr, int *cumM, int *cumMsqr, double *distMat, int N, double *Y, double *B0, double *B1, double *t,
                  double beta_theta, double *sig2_e, double thres, int *nr, double trunc)
{
  // declare types of all variables
  int i=0, j=0, inc=1, givelog=1;
  double invPsi[NMsqr], log_detInvPsi=0, power=0, logd, deviations[NM];


  if(theta < 0 || theta >= 1 ){

  	logd = log(0);

 	}


 	else{


  for(i=0;i<N;i++){

    invertPsi(theta, M[i], &distMat[cumMsqr[i]], &invPsi[cumMsqr[i]], thres, nr);

		for(j=0;j<M[i];j++){

	deviations[cumM[i]+j] = Y[cumM[i]+j] - B0[i] - B1[i]*t[cumM[i]+j];

		}
	}


 	for(i=0;i<N;i++){

  power += quadform2(&deviations[cumM[i]], &invPsi[cumMsqr[i]], M[i], inc, M[i]) / (2*sig2_e[i]);
	log_detInvPsi += matrixDet(&invPsi[cumMsqr[i]], M[i], M[i], givelog);
  }

  logd = (beta_theta-1)*log(trunc-theta) + .5*log_detInvPsi - power;

	}

  return(logd);

}


// Sample from full conditional for B0i

double funcB0(double *Y, double B1, double *invPsi, double sig2_e, int M, double *t)
{
  int i=0, inc=1;
  double zero=0, one = 1, mu_B0, sig_B0,
  ones[M], v[M], tdevinvPsiones = 0, onesinvPsiones, B0samples;

  for(i=0;i<M;i++){
  ones[i] = 1;
  }

  F77_NAME(dsymv)("U",&M,&one,invPsi,&M,ones,&inc,&zero,v,&inc);

   for(i=0;i<M;i++){
     tdevinvPsiones +=  (Y[i] - B1*t[i])*v[i];
   }

  onesinvPsiones = quadform2(ones, invPsi, M, inc, M);


    mu_B0 = tdevinvPsiones / onesinvPsiones;
    sig_B0 = sqrt(sig2_e / onesinvPsiones);


    B0samples = rnorm(mu_B0, sig_B0);

 return(B0samples);

}



// Sample from full conditional for B1i

double funcB1(double *Y, double B0, double *invPsi, double *t,  double mu_B1,
              double sig2_e, double mu_sig2_e, double g2, int M, double *mu_B1i, double *sig_B1i)
{
  int i=0, inc = 1;
  double v[M], zero=0, one=1, tdevinvPsit=0, tinvPsit, B1samples;

  F77_NAME(dsymv)("U",&M,&one,invPsi,&M,t,&inc,&zero,v,&inc);

   for(i=0;i<M;i++){
     tdevinvPsit +=  (Y[i] - B0)*v[i];
   }

   tinvPsit = quadform2(t, invPsi, M, inc, M);

   mu_B1i[0] =(tdevinvPsit + (sig2_e/(g2*exp(mu_sig2_e)))*mu_B1) / (tinvPsit + sig2_e/(g2*exp(mu_sig2_e)));
   sig_B1i[0] = sqrt(sig2_e / (tinvPsit + sig2_e/(g2*exp(mu_sig2_e))));

   B1samples = rnorm(mu_B1i[0], sig_B1i[0]);

   return(B1samples);

}



// Sample from full conditional for mu_B1

double funcmuB1(double *B1, double mu_sig2_e, double g1, double g2, int N)

{
  int i=0;
  double sumB1=0, mu_muB1, sig_muB1, muB1samples;

  for(i=0;i<N;i++){
    sumB1 += B1[i];
  }

  mu_muB1 = sumB1 / (N + g2/g1);
  sig_muB1 = sqrt((g2 * exp(mu_sig2_e)) / (N + g2/g1));

  muB1samples = rnorm(mu_muB1, sig_muB1);

  return(muB1samples);

}


// Sample from full conditional for sig2_sig2_e

double funcsig2sig2e(double *logsig2_e, double mu_sig2_e, int N)
{
  double a_sig2_sig2_e_pr =.5, b_sig2_sig2_e_pr = .5, a_sig2_sig2_e, b_sig2_sig2_e,
          sig2sig2esamples, sumdevlogsig2e=0;
  int i=0;

   for(i=0;i<N;i++){

     sumdevlogsig2e += (logsig2_e[i] - mu_sig2_e)*(logsig2_e[i] - mu_sig2_e);
   }

     a_sig2_sig2_e = .5*(N + 2*a_sig2_sig2_e_pr);
     b_sig2_sig2_e = .5*(sumdevlogsig2e + 2*b_sig2_sig2_e_pr);


     sig2sig2esamples = 1.0/rgamma(a_sig2_sig2_e,1/b_sig2_sig2_e);

     return(sig2sig2esamples);
}



// Sample from full conditional for g1

double funcg1(double r1, double mu_B1, double mu_sig2_e)
{
  double bg1, g1samples;

  bg1 = 0.5*r1*r1 + mu_B1*mu_B1 / (2*exp(mu_sig2_e));

  g1samples = 1.0/rgamma(1,1/bg1);

  return(g1samples);

}


// Sample from full conditional for g2

double funcg2(double *B1, double mu_B1, double mu_sig2_e, double r2, int N)
{
  double sumdevB1=0, g2samples, ag2, bg2;
  int i=0;

  for(i=0;i<N;i++){
  sumdevB1 += (B1[i]-mu_B1)*(B1[i]-mu_B1);
  }

  ag2 = (N+1)*.5;
  bg2 = 0.5*r2*r2 + sumdevB1/(2*exp(mu_sig2_e));

  g2samples = 1.0/rgamma(ag2,1/bg2);

  return(g2samples);
}






// /////////////////////////////////////////////////////////////////////////////////////

// Metropolis-Hastings algorithms

// sig2_e

double metroplogsig2e(double logsig2_e, double *Y, double B0, double B1, double *t, double *invPsi,
                      double mu_sig2_e, double sig2_sig2_e, int M, double sdmet)
{
	// sample log sig2_e with Metropolis-Hastings
	double cand, likeRat, b;

	cand = logsig2_e + rnorm(0,sdmet);

	likeRat =  funclogsig2e(cand, Y, B0, B1, t, invPsi, mu_sig2_e, sig2_sig2_e, M) -  funclogsig2e(logsig2_e, Y, B0, B1, t, invPsi, mu_sig2_e, sig2_sig2_e, M);

  if(likeRat >= 0){
    return(cand);
  }

  else{
  b = log(runif(0,1));

	if(b>likeRat){
		return(logsig2_e);
	}else{
		return(cand);
	}

  }
}



// mu_sig2_e

double metropmusig2e(double mu_sig2_e, double mu_B1, double *logsig2_e, double sig2_sig2_e, double mu_mu_sig2_e, double sig2_mu_sig2_e, double g1, double g2, double *B1, int N, double sdmet)
{
	// sample mu_sig2_e with Metropolis-Hastings
	double cand, likeRat, b;

	cand = mu_sig2_e + rnorm(0,sdmet);

	likeRat = funcMuSig2e(cand, mu_B1, logsig2_e, sig2_sig2_e, mu_mu_sig2_e, sig2_mu_sig2_e, g1, g2, B1, N) -  funcMuSig2e(mu_sig2_e, mu_B1, logsig2_e, sig2_sig2_e, mu_mu_sig2_e, sig2_mu_sig2_e, g1, g2, B1, N);

  if(likeRat >= 0){
    return(cand);
  }
  else{
  b = log(runif(0,1));

	if(b>likeRat){
		return(mu_sig2_e);
	}else{
		return(cand);
	}
}

}

                 
// theta

double metropTheta(double theta, int *M, int *Msqr, int NM, int NMsqr, int *cumM, int *cumMsqr, double *distMat, int N, double *Y, double *B0, double *B1, double *t,
                  double beta_theta, double *sig2_e, double thres, int *nr, double trunc, double sdmet)
{
	// sample theta with Metropolis-Hastings
	double cand, likeRat, b;

	cand = theta + rnorm(0,sdmet);

	if(cand<0 || cand>trunc)
	{
		return(theta);
  }

	likeRat = funcTheta(cand, M, Msqr, NM, NMsqr, cumM, cumMsqr, distMat, N, Y, B0, B1, t, beta_theta, sig2_e, thres, nr, trunc) -
            funcTheta(theta, M, Msqr, NM, NMsqr, cumM, cumMsqr, distMat, N, Y, B0, B1, t, beta_theta, sig2_e, thres, nr, trunc);

  if(likeRat >= 0){
   return(cand);
  }
  else{
  b = log(runif(0,1));

  if(b>likeRat){
		return(theta);
	}else{
		return(cand);
	}
}

}





// /////////////////////////////////////////////////////////////////////////////////////

void gibbsSampler(double *Y, int *M1, int *M2, int N, int iter, int thin, double r1, double r2, double beta_theta,
                    double a, double b, double thres,
                    double sdmetlogsig2e, double sdmetmusig2e, double sdmettheta,
                    double mu_mu_sig2_e, double sig2_mu_sig2_e, double *distMat, double trunc, int *test,
                    double *samples, double *BF, double *BFpatt, int progress, SEXP pBar, SEXP rho){


   int M[N], Msqr[N], NM=0, NMsqr=0, cumM[N+1],cumMsqr[N+1], nsave = iter/thin;
   int c=0, i=0, j=0, k=0;
   int alt[N];
  double B0[N], B1[N], sig2_e[N], logsig2_e[N],  muB1[N], sigB1[N];
  double mu_B1=0, mu_sig2_e=0, sumlsig2e = 0, sumlsig2esqr = 0, sig2_sig2_e, g1, g2, theta, mu_B1i[1], sig_B1i[1];
  double *sumY = Calloc(N,double), *sumY1sqr = Calloc(N,double), *sumY2sqr = Calloc(N,double),
         *sumY1 = Calloc(N,double), *sumY2 = Calloc(N,double), *logptestSamples = Calloc(nsave, double), *logpaltSamples = Calloc(nsave, double);
  double scale, log_priorpnull, log_priorpalt, log_priorodds, logpriorptest=0, logpriorpalt = 0, logptest, logpalt,
        log_postnull_samples[N*nsave], log_postalt_samples[N*nsave],
        left, right, log_postnull[N], log_postalt[N];

  int nr[1], nr2[1], excThres = 0;
  nr[0] = 0;
  nr2[0] = 0;



  cumM[0] = 0;
  cumMsqr[0] = 0;

  for(i=0;i<N;i++){

    M[i] = M1[i] + M2[i];
    Msqr[i] = M[i]*M[i];
    NM += M[i];
    NMsqr += Msqr[i];
    cumM[(i+1)] = cumM[i] + M[i];
    cumMsqr[(i+1)] = cumMsqr[i] + Msqr[i];
  }

   double invPsi[NMsqr], t[NM];


       for(i=0;i<N;i++){

     for(j=0;j<M[i];j++){
   sumY[i] += Y[cumM[i]+j];
     }

     for(j=0;j<M1[i];j++){
   sumY1[i] += Y[cumM[i]+j];
   sumY1sqr[i] += Y[cumM[i]+j]*Y[cumM[i]+j];
     }

     for(j=M1[i];j<M[i];j++){
   sumY2[i] += Y[cumM[i]+j];
   sumY2sqr[i] += Y[cumM[i]+j]*Y[cumM[i]+j];
     }
   }


  for(i=0;i<N;i++){
  B0[i] = sumY[i]/M[i];
  B1[i] = sumY2[i]/M2[i] - sumY1[i]/M1[i];
  sig2_e[i] = ((sumY1sqr[i] - sumY1[i]*sumY1[i]/M1[i])/M1[i] +
              (sumY2sqr[i] - sumY2[i]*sumY2[i]/M2[i])/M2[i])/2.0;
  logsig2_e[i] = log(sig2_e[i]);
  mu_B1 += B1[i]/N;
  mu_sig2_e += sig2_e[i]/N;
  }
  mu_sig2_e = log(mu_sig2_e);
      for(i=0;i<N;i++){
        sumlsig2e += logsig2_e[i];
        sumlsig2esqr += logsig2_e[i]*logsig2_e[i];
      }
  sig2_sig2_e = (sumlsig2esqr - sumlsig2e*sumlsig2e/N)/N;
  g1 = 1;
  g2 = 1;
  theta = .1;

  // Make t variable

  for(i=0;i<N;i++){

    for(j=0;j<M1[i];j++){
    t[cumM[i]+j] = -.5;
    }
    for(j=M1[i];j<M[i];j++){
    t[cumM[i]+j] = .5;
    }
  }


  // progress stuff
  SEXP sampCounter, R_fcall;
  int *pSampCounter;
  PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);


// Gibbs sampler
  for(k=0;k<iter;k++){

    R_CheckUserInterrupt();
    
        //Check progress
		if(progress && !((k+1)%progress)){
			pSampCounter[0]=k+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}

    //Rprintf("%d\n",k);


  for(i=0;i<N;i++){
    invertPsi(theta, M[i], &distMat[cumMsqr[i]], &invPsi[cumMsqr[i]], thres, &nr[0]);
  }


  // Sample B0
  for(i=0;i<N;i++){
  B0[i] = funcB0(&Y[cumM[i]], B1[i], &invPsi[cumMsqr[i]], sig2_e[i], M[i], &t[cumM[i]]);
   }

  // Sample B1
  for(i=0;i<N;i++){
  B1[i] = funcB1(&Y[cumM[i]], B0[i], &invPsi[cumMsqr[i]], &t[cumM[i]], mu_B1, sig2_e[i], mu_sig2_e, g2, M[i], mu_B1i, sig_B1i);
  muB1[i] = mu_B1i[0];
  sigB1[i] = sig_B1i[0];
  }

  // Sample mu_B1
  mu_B1 = funcmuB1(B1, mu_sig2_e, g1, g2, N);

  // Sample logsig2_e
  for(i=0;i<N;i++){
  logsig2_e[i] = metroplogsig2e(logsig2_e[i], &Y[cumM[i]], B0[i], B1[i], &t[cumM[i]], &invPsi[cumMsqr[i]], mu_sig2_e, sig2_sig2_e, M[i], sdmetlogsig2e);
  sig2_e[i] = exp(logsig2_e[i]);
    }

  // Sample mu_sig2_e
  mu_sig2_e = metropmusig2e(mu_sig2_e, mu_B1, logsig2_e, sig2_sig2_e, mu_mu_sig2_e, sig2_mu_sig2_e, g1, g2, B1, N, sdmetmusig2e);

  // Sample sig2_sig2_e
  sig2_sig2_e = funcsig2sig2e(logsig2_e, mu_sig2_e, N);

  // Sample g1
  g1 =  funcg1(r1, mu_B1, mu_sig2_e);

  // Sample g2
  g2 = funcg2(B1, mu_B1, mu_sig2_e, r2, N);

  // Sample theta
  theta = metropTheta(theta, M, Msqr, NM, NMsqr, cumM, cumMsqr, distMat, N, Y, B0, B1, t, beta_theta, sig2_e, thres, &nr2[0], trunc, sdmettheta);



  excThres += nr[0];

  // save parameters

  if(k == c*thin){
  c+=1;

  for(i=0;i<N;i++){
  samples[i*nsave+k/thin] = B0[i];
  samples[(i+N)*nsave+k/thin] = B1[i];
  samples[(i+2*N)*nsave+k/thin] = sig2_e[i];

  left = (a - muB1[i] / sqrt(exp(mu_sig2_e))) /
            sqrt(sigB1[i]*sigB1[i] /exp(mu_sig2_e));
  right = (b - muB1[i] / sqrt(exp(mu_sig2_e))) /
            sqrt(sigB1[i]*sigB1[i] /exp(mu_sig2_e));

  // Conditional posterior probabilities that subject i changed
  log_postnull_samples[i*nsave+k/thin] = max(expXphiAminusphiB(0, right, left, 1), -1000000);
  log_postalt_samples[i*nsave+k/thin] = max(log1mp(log_postnull_samples[i*nsave+k/thin]), -1000000);

  // Conditional posterior probabilities of the test-pattern and the alternative (all ones) pattern
  if(test[i] == 1){
    logptestSamples[k/thin] += log_postalt_samples[i*nsave+k/thin];
    }
  if(test[i] == 0){
    logptestSamples[k/thin] += log_postnull_samples[i*nsave+k/thin];
    }
  logpaltSamples[k/thin] += log_postalt_samples[i*nsave+k/thin];

  samples[(i+3*N+6)*nsave+k/thin] = log_postnull_samples[i*nsave+k/thin];
  samples[(i+4*N+6)*nsave+k/thin] = log_postalt_samples[i*nsave+k/thin];

  }
  samples[(3*N)*nsave+k/thin] =  mu_B1;
  samples[(3*N+1)*nsave+k/thin] = mu_sig2_e;
  samples[(3*N+2)*nsave+k/thin] = sig2_sig2_e;
  samples[(3*N+3)*nsave+k/thin] = g1;
  samples[(3*N+4)*nsave+k/thin] = g2;
  samples[(3*N+5)*nsave+k/thin] = theta;
  }

  }
  
  UNPROTECT(2);

  if(thres > 0){

     if(excThres > 0){

  Rprintf("Warning: at least one element of Psi x invPsi - I exceeds the treshold, in %d iterations \n", excThres);
     }

  }

    

  // Bayes factor per subject

 // Prior probability and odds that a subject changed
  scale = r1+r2;
   log_priorpalt = log(pcauchy(b,0,scale,0,0) + pcauchy(a,0,scale,1,0));
   log_priorpnull = log1mp(log_priorpalt);
  log_priorodds = log_priorpnull -log_priorpalt;

  // Prior probability of test pattern and alternative
  for(i=0;i<N;i++){
    alt[i] = 1;
  }
  logpriorptest = findlogPriorp(10000, N, test, a, b, r1, r2);
  logpriorpalt = findlogPriorp(10000, N, alt, a, b, r1, r2);
  

  // Marginal posterior probability and Bayes factor that subject i changed
  for(i=0;i<N;i++){
  log_postnull[i] = logMeanExpLogs(&log_postnull_samples[i*nsave], nsave);
  log_postalt[i] = log1mp(log_postnull[i]);
  BF[i] =  log_postnull[i] - log_postalt[i] - log_priorodds;
   }

  // Marginal posterior probability and Bayes factor for test-pattern and alternative (over subjects)

   logptest = logMeanExpLogs(logptestSamples, nsave);
   logpalt = logMeanExpLogs(logpaltSamples, nsave);

   BFpatt[0] = logptest - logpalt - logpriorptest + logpriorpalt;
   
   // REMOVE REMOVE REMOVE 
   //for(i=0;i<NMsqr;i++){
   //  samples[i] = invPsi[i];  
  // }

  Free(sumY);
  Free(sumY1);
  Free(sumY2);
  Free(sumY1sqr);
  Free(sumY2sqr);
  Free(logptestSamples);
  Free(logpaltSamples);

}



SEXP gibbsSamplerCall(SEXP YR, SEXP M1R, SEXP M2R, SEXP NR, SEXP testR, SEXP iterR, SEXP thinR, SEXP r1R, SEXP r2R,
                      SEXP beta_thetaR, SEXP truncR, SEXP aR, SEXP bR, SEXP thresR,
                      SEXP sdmetlogsig2eR, SEXP sdmetmusig2eR, SEXP sdmetthetaR,
                      SEXP mu_mu_sig2_eR, SEXP sig2_mu_sig2_eR, SEXP distMatR, SEXP progressR, SEXP pBar, SEXP rho){


  double *Y, r1, r2, beta_theta, a, b, thres, sdmetlogsig2e, sdmetmusig2e, sdmettheta,
          mu_mu_sig2_e, sig2_mu_sig2_e,*distMat, trunc;
  int K, *M1, *M2, N, iter, thin, nsave, *test, progress;

  Y = REAL(YR);
  r1 = NUMERIC_VALUE(r1R);
  r2 = NUMERIC_VALUE(r2R);
  M1 = INTEGER_POINTER(M1R);
  M2 = INTEGER_POINTER(M2R);
  N = INTEGER_VALUE(NR);
  test = INTEGER_POINTER(testR);
  iter = INTEGER_VALUE(iterR);
  thin = INTEGER_VALUE(thinR);
  beta_theta = NUMERIC_VALUE(beta_thetaR);
  trunc = NUMERIC_VALUE(truncR);
  a = NUMERIC_VALUE(aR);
  b = NUMERIC_VALUE(bR);
  thres = NUMERIC_VALUE(thresR);
  sdmetlogsig2e = NUMERIC_VALUE(sdmetlogsig2eR);
  sdmetmusig2e = NUMERIC_VALUE(sdmetmusig2eR);
  sdmettheta = NUMERIC_VALUE(sdmetthetaR);
  mu_mu_sig2_e = NUMERIC_VALUE(mu_mu_sig2_eR);
  sig2_mu_sig2_e = NUMERIC_VALUE(sig2_mu_sig2_eR);
  distMat = REAL(distMatR);
  progress = INTEGER_VALUE(progressR);



  nsave = iter/thin;
  K = nsave*(5*N+6);

  SEXP samplesR;
  SEXP BFR;
  SEXP BFpattR;
  SEXP returnR;

  PROTECT(samplesR = allocVector(REALSXP,K));
  PROTECT(BFR = allocVector(REALSXP,N));
  PROTECT(BFpattR = allocVector(REALSXP,1));
  PROTECT(returnR = allocVector(VECSXP, 3));

  GetRNGstate();
  gibbsSampler(Y, M1, M2, N, iter, thin, r1, r2, beta_theta, a,b, thres,
                        sdmetlogsig2e, sdmetmusig2e, sdmettheta,
                        mu_mu_sig2_e, sig2_mu_sig2_e, distMat, trunc, test, REAL(samplesR), REAL(BFR), REAL(BFpattR),
                        progress, pBar, rho);

  SET_VECTOR_ELT(returnR, 0, BFR);
  SET_VECTOR_ELT(returnR, 1, BFpattR);
  SET_VECTOR_ELT(returnR, 2, samplesR);

   UNPROTECT(4);
   PutRNGstate();
  return(returnR);

}
  
  


 
void BayesFactorsTrend(double aB1, double bB1, double aB3, double bB3, double *muB1, double *sig2B1, double *muB3, double *sig2B3, double *mu_sig2_e,
                        int *test, int iter, int N, double r1, double r2, double r3, double r4, double log_prioroddsIT,
                        double *BFint, double *BFtrend, double *BFit, double *BFpatt){
                          
  int k=0, i=0;
  
  int alt[N];
  double left1, right1, left2, right2;

  double *logptestSamplesInt = Calloc(iter, double), *logpaltSamplesInt = Calloc(iter, double),
         *logptestSamplesTrend = Calloc(iter, double), *logpaltSamplesTrend = Calloc(iter, double),
         *logptestSamplesIT = Calloc(iter, double), *logpaltSamplesIT = Calloc(iter, double);
  
  double scaleInt, scaleTrend,
        log_priorpnullInt, log_priorpaltInt, log_priorpnullTrend, log_priorpaltTrend, 
        log_prioroddsInt, log_prioroddsTrend, 
        logpriorptestInt=0, logpriorpaltInt = 0, logpriorptestTrend = 0, logpriorpaltTrend = 0,logpriorptestIT = 0, logpriorpaltIT = 0,
        logptestInt, logpaltInt, logptestTrend, logpaltTrend, logptestIT, logpaltIT,
        log_postnull_samplesInt[N*iter], log_postalt_samplesInt[N*iter], 
        log_postnull_samplesTrend[N*iter], log_postalt_samplesTrend[N*iter], 
        log_postnull_samplesIT[N*iter], log_postalt_samplesIT[N*iter],
        log_postnullInt[N], log_postaltInt[N], log_postnullTrend[N], log_postaltTrend[N], log_postnullIT[N], log_postaltIT[N];

 for(k=0;k<iter;k++){
   
  for(i=0;i<N;i++){
   

  left1 = (aB1 - muB1[i*iter + k] / sqrt(exp(mu_sig2_e[k]))) /
            sqrt(sig2B1[i*iter + k] /exp(mu_sig2_e[k]));
            
  right1 = (bB1 - muB1[i*iter + k] / sqrt(exp(mu_sig2_e[k]))) /
            sqrt(sig2B1[i*iter + k] /exp(mu_sig2_e[k]));
            
  left2 = (aB3 - muB3[i*iter + k] / sqrt(exp(mu_sig2_e[k]))) /
            sqrt(sig2B3[i*iter + k] /exp(mu_sig2_e[k]));
            
  right2 = (bB3 - muB3[i*iter + k] / sqrt(exp(mu_sig2_e[k]))) /
            sqrt(sig2B3[i*iter + k] /exp(mu_sig2_e[k]));

  // Conditional posterior probabilities that subject i changed
  // Intercept
  log_postnull_samplesInt[i*iter+k] = max(expXphiAminusphiB(0, right1, left1, 1), -1000000);
  log_postalt_samplesInt[i*iter+k] = max(log1mp(log_postnull_samplesInt[i*iter+k]), -1000000);
  
  //Trend
  log_postnull_samplesTrend[i*iter+k] = max(expXphiAminusphiB(0, right2, left2, 1), -1000000);
  log_postalt_samplesTrend[i*iter+k] = max(log1mp(log_postnull_samplesTrend[i*iter+k]), -1000000);
  
  //Intercept + trend
  log_postnull_samplesIT[i*iter+k] = log_postnull_samplesInt[i*iter+k] + log_postnull_samplesTrend[i*iter+k];
  log_postalt_samplesIT[i*iter+k] = log_postalt_samplesInt[i*iter+k] + log_postalt_samplesTrend[i*iter+k];
  

  // Conditional posterior probabilities of the test-pattern and the alternative (all ones) pattern
  // Intercept
  if(test[i] == 1){
    logptestSamplesInt[k] += log_postalt_samplesInt[i*iter+k];
    }
  if(test[i] == 0){
    logptestSamplesInt[k] += log_postnull_samplesInt[i*iter+k];
    }
  logpaltSamplesInt[k] += log_postalt_samplesInt[i*iter+k];
  
   // Trend
  if(test[i] == 1){
    logptestSamplesTrend[k] += log_postalt_samplesTrend[i*iter+k];
    }
  if(test[i] == 0){
    logptestSamplesTrend[k] += log_postnull_samplesTrend[i*iter+k];
    }
  logpaltSamplesTrend[k] += log_postalt_samplesTrend[i*iter+k];
  
  
  //intercept + trend
  if(test[i] == 1){
    logptestSamplesIT[k] += log_postalt_samplesIT[i*iter+k];
    }
  if(test[i] == 0){
    logptestSamplesIT[k] += log_postnull_samplesIT[i*iter+k];
    }
  logpaltSamplesIT[k] += log_postalt_samplesIT[i*iter+k];
  
    }

  }


  // Bayes factor per subject

 // Prior probability and odds that a subject changed
 // Intercept
  scaleInt = r1+r2;
   log_priorpaltInt = log(pcauchy(bB1,0,scaleInt,0,0) + pcauchy(aB1,0,scaleInt,1,0));
   log_priorpnullInt = log1mp(log_priorpaltInt);
  log_prioroddsInt = log_priorpnullInt -log_priorpaltInt;
  
  // Trend
  scaleTrend = r3+r4;
   log_priorpaltTrend = log(pcauchy(bB3,0,scaleTrend,0,0) + pcauchy(aB3,0,scaleTrend,1,0));
   log_priorpnullTrend = log1mp(log_priorpaltTrend);
  log_prioroddsTrend = log_priorpnullTrend -log_priorpaltTrend;
  
  // Intercept + trend computed in R
  

  // Prior probability of test pattern and alternative
    for(i=0;i<N;i++){
    alt[i] = 1;
  }
  // Intercept
  logpriorptestInt = findlogPriorp(10000, N, test, aB1, bB1, r1, r2);
  logpriorpaltInt = findlogPriorp(10000, N, alt, aB1, bB1, r1, r2);
  
  // Trend
  logpriorptestTrend = findlogPriorp(10000, N, test, aB3, bB3, r3, r4);
  logpriorpaltTrend = findlogPriorp(10000, N, alt, aB3, bB3, r3, r4);
  
  // Intercept + trend
  logpriorptestIT = findlogPriorpIT(10000, N, test, aB1, bB1, aB3, bB3, r1, r2, r3, r4);
  logpriorpaltIT = findlogPriorpIT(10000, N, alt, aB1, bB1, aB3, bB3, r1, r2, r3, r4);


  // Marginal posterior probability and Bayes factor that subject i changed
  // Intercept
  for(i=0;i<N;i++){
  log_postnullInt[i] = logMeanExpLogs(&log_postnull_samplesInt[i*iter], iter);
  log_postaltInt[i] = log1mp(log_postnullInt[i]);
  BFint[i] =  log_postnullInt[i] - log_postaltInt[i] - log_prioroddsInt;
   }
   
  // Trend
  for(i=0;i<N;i++){
  log_postnullTrend[i] = logMeanExpLogs(&log_postnull_samplesTrend[i*iter], iter);
  log_postaltTrend[i] = log1mp(log_postnullTrend[i]);
  BFtrend[i] =  log_postnullTrend[i] - log_postaltTrend[i] - log_prioroddsTrend;
   }
   
  // Intercept + trend
  for(i=0;i<N;i++){
  log_postnullIT[i] = logMeanExpLogs(&log_postnull_samplesIT[i*iter], iter);
  log_postaltIT[i] = log1mp(log_postnullIT[i]);
  BFit[i] =  log_postnullIT[i] - log_postaltIT[i] - log_prioroddsIT;
   }

  // Marginal posterior probability and Bayes factor for test-pattern and alternative (over subjects)
  // Intercept  
  logptestInt = logMeanExpLogs(logptestSamplesInt, iter);
  logpaltInt = logMeanExpLogs(logpaltSamplesInt, iter);
  BFpatt[0] = logptestInt - logpaltInt - logpriorptestInt + logpriorpaltInt;

  // Trend  
  logptestTrend = logMeanExpLogs(logptestSamplesTrend, iter);
  logpaltTrend = logMeanExpLogs(logpaltSamplesTrend, iter);
  BFpatt[1] = logptestTrend - logpaltTrend - logpriorptestTrend + logpriorpaltTrend;
  
  // Intercept + trend  
  logptestIT = logMeanExpLogs(logptestSamplesIT, iter);
  logpaltIT = logMeanExpLogs(logpaltSamplesIT, iter);
  BFpatt[2] = logptestIT - logpaltIT - logpriorptestIT + logpriorpaltIT;

  Free(logptestSamplesInt);
  Free(logpaltSamplesInt);
  Free(logptestSamplesTrend);
  Free(logpaltSamplesTrend);
  Free(logptestSamplesIT);
  Free(logpaltSamplesIT);

}


SEXP BayesFactorsTrendCall(SEXP aB1R, SEXP bB1R, SEXP aB3R, SEXP bB3R, SEXP muB1R, SEXP sig2B1R, SEXP muB3R, SEXP sig2B3R,
                                SEXP mu_sig2_eR, SEXP testR, SEXP iterR, SEXP NR, SEXP r1R, SEXP r2R, SEXP r3R, SEXP r4R,
                                SEXP log_prioroddsITR){


  double r1, r2, r3, r4, aB1, bB1, aB3, bB3, log_prioroddsIT, *muB1, *sig2B1, *muB3, *sig2B3, *mu_sig2_e;
  int N, iter, *test;

  r1 = NUMERIC_VALUE(r1R);
  r2 = NUMERIC_VALUE(r2R);
  r3 = NUMERIC_VALUE(r3R);
  r4 = NUMERIC_VALUE(r4R);
  N = INTEGER_VALUE(NR);
  test = INTEGER_POINTER(testR);
  iter = INTEGER_VALUE(iterR);
  aB1 = NUMERIC_VALUE(aB1R);
  bB1 = NUMERIC_VALUE(bB1R);
  aB3 = NUMERIC_VALUE(aB3R);
  bB3 = NUMERIC_VALUE(bB3R);
  muB1 = REAL(muB1R);
  sig2B1 = REAL(sig2B1R);
  muB3 = REAL(muB3R);
  sig2B3 = REAL(sig2B3R);
  mu_sig2_e = REAL(mu_sig2_eR);
  log_prioroddsIT  = NUMERIC_VALUE(log_prioroddsITR);


  SEXP BFintR;
  SEXP BFtrendR;
  SEXP BFitR;
  SEXP BFpattR;
  SEXP returnR;


  PROTECT(BFintR = allocVector(REALSXP,N));

  PROTECT(BFtrendR = allocVector(REALSXP,N));
 
  PROTECT(BFitR = allocVector(REALSXP,N));
  PROTECT(BFpattR = allocVector(REALSXP,3));


  
  PROTECT(returnR = allocVector(VECSXP, 4));

  GetRNGstate();
 
 
 
  BayesFactorsTrend(aB1, bB1, aB3, bB3, muB1, sig2B1, muB3, sig2B3, mu_sig2_e,
                        test, iter, N, r1, r2, r3, r4, log_prioroddsIT, 
                        REAL(BFintR), REAL(BFtrendR), REAL(BFitR), REAL(BFpattR));
  
  
                        
  SET_VECTOR_ELT(returnR, 0, BFintR);
 SET_VECTOR_ELT(returnR, 1, BFtrendR);
 SET_VECTOR_ELT(returnR, 2, BFitR);
  SET_VECTOR_ELT(returnR, 3, BFpattR);

   PutRNGstate();

 UNPROTECT(5);
   
   
  return(returnR);

}

                    

SEXP log1mpCall(SEXP logpR){

  double logp, *log1mpC;
  SEXP log1mpR;

  logp = NUMERIC_VALUE(logpR);

  PROTECT(log1mpR = allocVector(REALSXP,1));

  log1mpC = REAL(log1mpR);

  log1mpC[0] = log1mp(logp);

  UNPROTECT(1);
  return(log1mpR);

}

SEXP findlogPriorpCall(SEXP iterR, SEXP NR, SEXP testR, SEXP aR, SEXP bR, SEXP r1R, SEXP r2R){
  
  int iter, N, *test;
  double a, b, r1, r2, *logprpC;
  SEXP logprpR;
  
  iter = INTEGER_VALUE(iterR);
  N = INTEGER_VALUE(NR);
  test = INTEGER_POINTER(testR);
  a = NUMERIC_VALUE(aR);
  b = NUMERIC_VALUE(bR);
  r1= NUMERIC_VALUE(r1R);
  r2 = NUMERIC_VALUE(r2R);
  
  PROTECT(logprpR = allocVector(REALSXP,1));
  
  logprpC = REAL(logprpR);
  
  GetRNGstate();
  logprpC[0] = findlogPriorp(iter, N, test, a, b, r1, r2);
  
  UNPROTECT(1);
  PutRNGstate();
  return(logprpR);
  
}

SEXP findlogPriorpITCall(SEXP iterR, SEXP NR, SEXP testR, SEXP aB1R, SEXP bB1R, SEXP aB3R, SEXP bB3R, 
                            SEXP r1R, SEXP r2R, SEXP r3R, SEXP r4R){
  
  int iter, N, *test;
  double aB1, bB1, aB3, bB3, r1, r2, r3, r4, *logprpC;
  SEXP logprpR;
  
  iter = INTEGER_VALUE(iterR);
  N = INTEGER_VALUE(NR);
  test = INTEGER_POINTER(testR);
  aB1 = NUMERIC_VALUE(aB1R);
  bB1 = NUMERIC_VALUE(bB1R);
  aB3 = NUMERIC_VALUE(aB3R);
  bB3 = NUMERIC_VALUE(bB3R);
  r1 = NUMERIC_VALUE(r1R);
  r2 = NUMERIC_VALUE(r2R);
  r3 = NUMERIC_VALUE(r3R);
  r4 = NUMERIC_VALUE(r4R);
  
  PROTECT(logprpR = allocVector(REALSXP,1));
  
  logprpC = REAL(logprpR);
  
  GetRNGstate();
  logprpC[0] = findlogPriorpIT(iter, N, test, aB1, bB1, aB3, bB3, r1, r2, r3, r4);
  
  UNPROTECT(1);
  PutRNGstate();
  return(logprpR);
  
}


SEXP expXphiAminusphiBCall(SEXP xR, SEXP aR, SEXP bR, SEXP returnLogR){
  
  double x, a, b;
  int returnLog;
  
  SEXP logpR;
  
  x = NUMERIC_VALUE(xR);
  a = NUMERIC_VALUE(aR);
  b = NUMERIC_VALUE(bR);
  returnLog = INTEGER_VALUE(returnLogR);
  
  PROTECT(logpR = allocVector(REALSXP,1));
  
  REAL(logpR)[0] = expXphiAminusphiB(x, a, b, returnLog);
 
  UNPROTECT(1);
  
  return(logpR);
}

SEXP matrixDetCall(SEXP AR, SEXP NR, SEXP returnLogR){
  
  double *A;
  int N, returnLog;
  
  SEXP detR;
  
  A = REAL(AR);
  N = INTEGER_VALUE(NR);
  returnLog = INTEGER_VALUE(returnLogR);
  
  PROTECT(detR = allocVector(REALSXP,1));
  
  REAL(detR)[0] = matrixDet(A, N, N, returnLog);
  
  UNPROTECT(1);
  
  return(detR);

}




















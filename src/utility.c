#include "utility.h"

int minc(int x, int y){

  if( (x)<(y)){
  return(x);
  }

	else{
	return(y);
	}
}


void inverse(double *A, int N)
{
  int IPIV[N+1], INFO, LWORK = N*N;
  double WORK[LWORK];
  
  F77_NAME(dgetrf)(&N, &N, A, &N, IPIV, &INFO);
  F77_NAME(dgetri)(&N, A, &N, IPIV, WORK, &LWORK, &INFO);
  
}

void debugPrintMatrix(double *X, int rows, int cols)
{
	int i=0,j=0;
	
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			Rprintf("%f ",X[j*rows+i]);
		}
		Rprintf("\n");
	}
}

void debugPrintVector(double *x, int len)
{
	int i=0;
	
	for(i=0;i<len;i++){
		Rprintf("%f ",x[i]);
	}
	Rprintf("\n");
}

SEXP RLogMeanExpLogs(SEXP Rv, SEXP Rlen)
{
	double *v;
	int len;
	SEXP ret;
	double *retp;
	
	len = INTEGER_VALUE(Rlen);
	v = REAL(Rv);
	PROTECT(ret = allocVector(REALSXP,1));
	retp = REAL(ret);
	
	retp[0] = logMeanExpLogs(v,len);
	
	UNPROTECT(1);
	return(ret);
}


double logSumExpLogs(double *v, int len)
{
	int i=0;
	double sum=v[0];
	
	for(i=1;i<len;i++)
	{
		sum = LogOnePlusExpX(v[i]-sum)+sum;
	}
	
	return(sum);
}

double logMeanExpLogs(double *v, int len)
{	
	double sum=0;
	sum = logSumExpLogs(v,len);
	return(sum - log(len));
}



// Calculate log(1 + x), preventing loss of precision for small values of x.
// The input x must be larger than -1 so that log(1 + x) is real.
// Taken from http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
        
        error("Attempt to compute log(1+x) on x<=-1!");
    }

	if (fabs(x) > 0.375)
    {
        // x is sufficiently large that the obvious evaluation is OK
        return log(1.0 + x);
    }

	// For smaller arguments we use a rational approximation
	// to the function log(1+x) to avoid the loss of precision
	// that would occur if we simply added 1 to x then took the log.

    const double p1 =  -0.129418923021993e+01;
    const double p2 =   0.405303492862024e+00;
    const double p3 =  -0.178874546012214e-01;
    const double q1 =  -0.162752256355323e+01;
    const double q2 =   0.747811014037616e+00;
    const double q3 =  -0.845104217945565e-01;
    double t, t2, w;

    t = x/(x + 2.0);
    t2 = t*t;
    w = (((p3*t2 + p2)*t2 + p1)*t2 + 1.0)/(((q3*t2 + q2)*t2 + q1)*t2 + 1.0);
    return 2.0*t*w;
}

// return log(1 + exp(x)), preventing cancellation and overflow */
// From http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
double LogOnePlusExpX(double x)
{
    const double LOG_DBL_EPSILON = log(DBL_EPSILON);
    const double LOG_ONE_QUARTER = log(0.25);

    if (x > -LOG_DBL_EPSILON)
    {
        // log(exp(x) + 1) == x to machine precision
        return x;
    }
    else if (x > LOG_ONE_QUARTER)
    {
        return log( 1.0 + exp(x) );
    }
    else
    {
        // Prevent loss of precision that would result from adding small argument to 1.
        return LogOnePlusX( exp(x) );
    }
}


// Compute determinant of an N by N matrix A
double matrixDet(double *A, int N, int LDA, int doLog)
{
//SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
	int i=0, info=0;
	double B[N*N], logDet=0;
	
	//Memcpy(B,A,N*N);
	for(i=0;i<N;i++){
		Memcpy(&B[i*N],&A[i*LDA],N);
	}
	
	F77_CALL(dpotrf)("U", &N, B, &N, &info);
	if(info){
		Rprintf("Cholesky decomposition in matrixDet() returned nonzero info %d.\n",info);
	}
	
	for(i=0;i<N;i++)
	{
		logDet += 2 * log(B[i*N+i]);
	}
	
	if(doLog){
		return(logDet);
	}else{
		return(exp(logDet));
	}
}

double quadform(double *x, double *A, int N, int incx, int LDA)
{
  
  int Nsqr = N*N,info,i=0,j=0;
  double *B = Calloc(Nsqr,double);
  //double one=1;
  //double zero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;

  for(i=0;i<N;i++){
    y[i] = x[i*incx];
  }
  for(i=0;i<N;i++){
	Memcpy(&B[i*N],&A[i*LDA],N);
  }

  F77_NAME(dpotrf)("U", &N, B, &N, &info);
  F77_NAME(dtrmv)("U","N","N", &N, B, &N, y, &iOne);
  
  for(i=0;i<N;i++){
    sumSq += y[i]*y[i];
  }
  
  Free(B);
  
  return(sumSq);
}

//version without cholesky decomp for non positive definite matrix A
double quadform2(double *x, double *A, int N, int incx, int LDA)
{
  
  int i=0;
  double dOne=1;
  double dZero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;
  
  F77_NAME(dgemv)("N", &N, &N, &dOne, A, &LDA, x, &incx, &dZero, y, &iOne);
 
  for(i=0;i<N;i++){
    sumSq += y[i]*x[i];
  }
    
  return(sumSq);
}



int InvMatrixUpper(double *A, int p)
{
      int info1, info2;
      F77_NAME(dpotrf)("U", &p, A, &p, &info1);
      F77_NAME(dpotri)("U", &p, A, &p, &info2);      
      //make sure you make it symmetric later...
      return(info1);
}



void rmvGaussianC(double *mu, double *Sigma, int p)
{
  double ans[p];
  int info, psqr,j=0, intOne=1;
  double *scCp, one = 1; //zero = 0;
  
  psqr = p * p;
  scCp = Memcpy(Calloc(psqr,double), Sigma, psqr);

  F77_NAME(dpotrf)("L", &p, scCp, &p, &info);
  if (info){
	error("Nonzero info from dpotrf: Sigma matrix is not positive-definite");
  }
  //GetRNGstate();
  for(j=0;j<p;j++)
    {
      ans[j] = rnorm(0,1);
    }
  F77_NAME(dtrmv)("L","N","N", &p, scCp, &p, ans, &intOne);
  F77_NAME(daxpy)(&p, &one, ans, &intOne, mu, &intOne);
  //PutRNGstate();
  Free(scCp);
}

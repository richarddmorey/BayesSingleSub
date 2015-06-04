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

int minc(int x, int y);
void inverse(double *A, int N);

void debugPrintMatrix(double *X, int rows, int cols);
void debugPrintVector(double *x, int len);

SEXP RLogMeanExpLogs(SEXP Rv, SEXP Rlen);
double logSumExpLogs(double *v, int len);
double logMeanExpLogs(double *v, int len);

double LogOnePlusX(double x);
double LogOnePlusExpX(double x);

double matrixDet(double *A, int N, int LDA, int doLog);
double quadform(double *x, double *A, int N, int incx, int LDA);
double quadform2(double *x, double *A, int N, int incx, int LDA);
int InvMatrixUpper(double *A, int p);
void rmvGaussianC(double *mu, double *Sigma, int p);




#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static R_INLINE double*
internal_symmetrize(double *a, int nc)
{
    int i,j;
    for (i = 1; i < nc; i++)
  for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

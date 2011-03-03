#ifndef __glm_h__
#define __glm_h__

// C++ headers
#include <vector>

// External library headers
#include <gsl/gsl_cdf.h> 

// Local headers
#include "fortran_matrix.h"
#include "svd.h"
#include "rank.h" 

using namespace std;

//#include "inner.h" // Inner product of two vectors

// This file contains data structures/routines for handling the
// internals of the GLM function originally written in Matlab.


// Return values for the GLM, as per the original Matlab routine
struct GLMData
{
  // Solution to the pseudoinverse problem X * beta = y
  vector<double> beta;

  // 1 - fcdf(F, V1, V2), as computed in the glm function
  double p;

  double F; // input to fcdf

  // y'*y - beta'*Xy, as computed in the glm function
  double ErrorSS;

  // length(y) - rank(X), as computed in the glm function
  int V2;
};


// This version of the glm function assumes Kt is a vector.  It makes
// the code a good deal cleaner (but not really any faster when Kt is
// a 1-by-N matrix) than the more general case where Kt is a matrix.
/*
void glm(const FortranMatrix &Xf, 
const FortranMatrix &XftXf, 
const FortranMatrix &XftXfti, 
const vector<double> &snp, 
const double snptsnp, 
const double snpty, 
const double yty, 
const vector<double> &Kt, 
const vector<double> &Xfty, 
const double rK, 
GLMData& glm_data)
*/
void glm(const FortranMatrix& X,
	 const vector<double>& y,
	 const vector<double>& Kt,
	 GLMData& glm_data /*return value, to be filled in*/)
{  
  int n  = y.size();
  int V1 = 1; // Kt is assumed to be a row vector in this version



  // Compute "V2" now that we know rX == rank(X) == rank(X^T X)
  glm_data.V2 = n - rX;

  // Compute ErrorSS for return in the GLMData data structure
  //! @todo cblas_ddot()
  //glm_data.ErrorSS = inner(y,y) - inner(glm_data.beta, XTy);
  double yty = cblas_ddot(y.size(), &y[0], 1, &y[0], 1);
  glm_data.ErrorSS = yty - cblas_ddot(glm_data.beta.size(), &glm_data.beta[0], 1, &XTy[0], 1);

  // Debugging: print ErrorSS
  // cout << "glm_data.ErrorSS=" << glm_data.ErrorSS << endl;

  ////////////////////////////////////////////////////////////////////////////////
  // Up to here, the matrix-Kt and vector-Kt versions of the glm
  // function are identical...
  
  // Compute the matrix-vector product, Kb = Kt * beta.
  // Note if Kt is actually a vector, Kb will be a scalar...
  //double Kb = inner(Kt, glm_data.beta);
  double Kb = cblas_ddot(Kt.size(), &Kt[0], 1, &glm_data.beta[0], 1);

  // We need the rank of the Kt "matrix" now.  Since Kt is
  // assumed to be a vector, we can just assume it has rank 1.
  int rK = 1;

  // F = Kb' * inv(Kt * G * Kt') * Kb * V2 / (rK * ErrorSS);
  //
  // Note Kb is a scalar in this version

  // 1.) Re-use the SVD of X^T * X to compute temp = G * Kt', for the
  // temporary result vector, temp.  
  vector<double> temp;
  svd_apply(U, S, VT, /*result=*/temp, Kt);
  
  // 2.) Compute the inner product c = (Kt,temp), using 'temp'
  // from the step above.  
  //double c = inner(Kt,temp);
  double c = cblas_ddot(Kt.size(), &Kt[0], 1, &temp[0], 1);

  // 3.) Note: for F we need to compute Kb' * inv(c) * Kb, but since c
  // is a scalar, we just divide!
  double F = Kb * Kb / c * static_cast<double>(glm_data.V2) / static_cast<double>(rK) / glm_data.ErrorSS; 

  // Debugging: print F
  // cout << "F=" << F << endl;

  // F must be positive otherwise the F distribution will return nan
  if (F < 0)
    {
      cerr << "Error! F<0 obtained, F distribution can only args >= 0." << endl;
      exit(1);
    }
  
  // Compute probability, p.  This part requires the GSL.
  glm_data.p = 1. - gsl_cdf_fdist_P(F, static_cast<double>(V1), static_cast<double>(glm_data.V2));
}


#endif

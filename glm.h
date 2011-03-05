#ifndef __glm_h__
#define __glm_h__

// C++ headers
#include <vector>

// External library headers
#include <gsl/gsl_cdf.h> 
#include <cblas.h>

// Local headers
#include "fortran_matrix.h"

#define swap(item1, item2, itemT) (itemT = item1; item1 = item2; item2 = item1)

using namespace std;

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
void glm(const FortranMatrix &X, 
	 const FortranMatrix &XtX, 
	 const FortranMatrix &XtXi, 
	 const vector<double> &XtSNP,
	 const double SNPtSNP, 
	 const double SNPty, 
	 const double yty, 
	 const vector<double> &Kt, 
	 const vector<double> &Xty, 
	 double rX,
	 GLMData& glm_data)
{  
  int m  = X.get_n_rows(), n = X.get_n_cols();
  int V1 = 1; // Kt is assumed to be a row vector in this version

  /*
  vector<double> SNPtX(m); // SNPtX = (XtSNP)t
  for(unsigned i = 0; i < m; i++)
    SNPtX[i] = XtSNP[m - 1 - i];
  */

  // G = XtXi
  // compute transpose of SNPtXG: nx1
  vector<double> GtXtSNP(n, 0.0);

  //! @todo use cblas_dsymv for this
  cblas_dgemv(CblasColMajor,
	      CblasTrans,
	      n,
	      n,
	      1.0,
	      &XtXi.values[0],
	      n,
	      &XtSNP[0],
	      1,
	      0.0,
	      &GtXtSNP[0],
	      1);

  // compute SNPtXGXtSNP (scalar)
  double SNPtXGXtSNP = cblas_ddot(n, &XtSNP[0], 1, &GtXtSNP[0], 1);

  // compute S = Schur complement of partitioned matrix to invert
  double S = SNPtSNP - SNPtXGXtSNP;
  if(!S){
    // bad news
    glm_data.F = 0.0;
    return;
  }

  S = 1.0 / S;

  // compute G' = (X'tX')i, X' = [X SNP]
  // Gn = [G + S*(snptXG'*snptXG), -S*snptXG'; -S*snptXG, S]; % n + 1 x n + 1
  FortranMatrix Gn = XtXi;

  // compute G + S *(SNPtXG'*SNPtXG) (n x n)
  // = G + S * (GtXtSNP*GtXtSNP')?= to within machine precision
  /*! @todo cblas_dspr(CblasColMajor, ); ? requires packed symmetric matrix,
    future operations using this result must also assume Gn is packed
   */
  //cblas_dsyr(CblasColMajor, CblasUpper, n, S, &GtXtSNP[0], 1, &Gn.values[0], n);
  cblas_dger(CblasColMajor, n,
	     n, S, &GtXtSNP[0], 1, &GtXtSNP[0], 1,
	     &Gn.values[0], n);

  //Gn.print("Gn before resize");

  Gn.resize_retain(n + 1, n + 1);

  //Gn.print("Gn after resize");


  // compute right and bottom edges of Gn
  // right edge
  cblas_daxpy(n, S, &GtXtSNP[0], 1, &Gn.values[n * (n + 1)], 1);
  cblas_dcopy(n, &Gn.values[n * (n + 1)], 1, &Gn.values[n], n + 1);

  // store bottom right element of Gn
  Gn(n, n) = S;
 
  // Xtyn = [Xty; snpty]; % n + 1 x 1
  vector<double> Xtyn(n+1);
  Xtyn.assign(Xty.begin(), Xty.end()); // sets size to n
  Xtyn.push_back(SNPty); // append 1
  
  // compute beta (n + 1 x 1)
  // beta = Gn * Xtyn
  glm_data.beta.resize(n + 1);

  //! @todo use cblas_dsymv()
  cblas_dgemv(CblasColMajor, CblasNoTrans, n + 1, n + 1, 1.0, 
	      &Gn.values[0], n + 1, &Xtyn[0], 1, 0.0, &glm_data.beta[0], 1);

  //! @todo beta is wrong
  
  // Compute ErrorSS for return in the GLMData data structure
  glm_data.ErrorSS = yty - 
    cblas_ddot(n + 1, &glm_data.beta[0], 1, &Xtyn[0], 1);

  ////////////////////////////////////////////////////////////////////////////////
  // Up to here, the matrix-Kt and vector-Kt versions of the glm
  // function are identical...
  
  // Compute the matrix-vector product, Kb = Kt * beta.
  // Note if Kt is actually a vector, Kb will be a scalar...
  //double Kb = inner(Kt, glm_data.beta);
  double Kb = glm_data.beta.back();

  // We need the rank of the Kt "matrix" now.  Since Kt is
  // assumed to be a vector, we can just assume it has rank 1.
  int rK = 1;

  //! @todo update rX from trace
  double trn = S * (SNPtSNP - SNPtXGXtSNP);
  rX++;

  // Compute "V2" now that we know rX == rank(X) == rank(X^T X)
  glm_data.V2 = m - rX;


  // 3.) Note: for F we need to compute Kb' * inv(c) * Kb, but since c
  // is a scalar, we just divide!
  // F = Kb' * inv(Kt * G * Kt') * Kb * V2 / (rK * ErrorSS);
  glm_data.F = (1.0 / S) * Kb * Kb * glm_data.V2 / (rK * glm_data.ErrorSS);

  //Kb * Kb / c * static_cast<double>(glm_data.V2) / static_cast<double>(rK) / glm_data.ErrorSS; 
}


#endif

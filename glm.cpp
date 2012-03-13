/*


Copyright (c) 2011, The Arizona Board of Regents on behalf of 
The University of Arizona

All rights reserved.

Developed by Peter Bailey, Tapasya Patki, and Greg Striemer with
support from the iPlant Collaborative as a collaboration between
participants at BIO5 at The University of Arizona (the primary hosting
institution), Cold Spring Harbor Laboratory, The University of Texas
at Austin, and individual contributors. Find out more at
http://www.iplantcollaborative.org/.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright 
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution.
 * Neither the name of the iPlant Collaborative, BIO5, The University 
   of Arizona, Cold Spring Harbor Laboratory, The University of Texas at 
   Austin, nor the names of other contributors may be used to endorse or 
   promote products derived from this software without specific prior 
   written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


*/
#include <gsl/gsl_cdf.h> 
extern "C"{
#include <cblas.h>
}
#include <sstream>

#include "glm.h"

using namespace std;

// This version of the glm function assumes Kt is a vector.  It makes
// the code a good deal cleaner (but not really any faster when Kt is
// a 1-by-N matrix) than the more general case where Kt is a matrix.
void glm(unsigned id, unsigned iteration,
	 int n,
	 FortranMatrix &XtXi, // updated
	 const double *XtSNP,
	 const double SNPtSNP, 
	 const double SNPty, 
	 const double yty, 
	 vector<double> &Xty, // updated
	 double rX,
	 // output
	 GLMData& glm_data)
{  
  //int 
    //m  = X.get_n_rows(),
    //n = X.get_n_cols();
  //int V1 = 1; // Kt is assumed to be a row vector in this version

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
#ifdef _DEBUG
  if(!id){
    stringstream ss;
    ss << "GtXtSNP_" << iteration << "p.dat";
    writeD(ss.str(), GtXtSNP);
  }
#endif


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
  // XtXi = [G + S*(snptXG'*snptXG), -S*snptXG'; -S*snptXG, S]; % n + 1 x n + 1

  // compute G + S *(SNPtXG'*SNPtXG) (n x n)
  // = G + S * (GtXtSNP*GtXtSNP')?= to within machine precision
  /*! @todo cblas_dspr(CblasColMajor, ); ? requires packed symmetric matrix,
    future operations using this result must also assume XtXi is packed
   */
  //cblas_dsyr(CblasColMajor, CblasUpper, n, S, &GtXtSNP[0], 1, &XtXi.values[0], n);
  cblas_dger(CblasColMajor, n,
	     n, S, &GtXtSNP[0], 1, &GtXtSNP[0], 1,
	     &XtXi.values[0], n);

  //! @todo this could be avoided by use of lda in computation of XtXi;
  // just allocate XtXi as n+1xn+1 and use lda=n+1, M = n, N = n
  XtXi.resize_retain(n + 1, n + 1);

  // compute right and bottom edges of XtXi
  // right edge
  cblas_daxpy(n, -S, &GtXtSNP[0], 1, &XtXi.values[n * (n + 1)], 1);
  cblas_dcopy(n, &XtXi.values[n * (n + 1)], 1, &XtXi.values[n], n + 1);

  // store bottom right element of XtXi
  XtXi(n, n) = S;
 
  // Xtyn = [Xty; snpty]; % n + 1 x 1
  Xty.push_back(SNPty); // append 1

  // compute beta (n + 1 x 1)
  // beta = XtXi * Xty
  glm_data.beta.resize(n + 1);

  //! @todo use cblas_dsymv()
  cblas_dgemv(CblasColMajor, CblasNoTrans, n + 1, n + 1, 1.0, 
	      &XtXi.values[0], n + 1, &Xty[0], 1, 0.0, &glm_data.beta[0], 1);

  // Compute ErrorSS for return in the GLMData data structure
  glm_data.ErrorSS = yty - 
    cblas_ddot(n + 1, &glm_data.beta[0], 1, &Xty[0], 1);

  // Compute the matrix-vector product, Kb = Kt * beta.
  // Note if Kt is actually a vector, Kb will be a scalar...
  //double Kb = inner(Kt, glm_data.beta);
  double Kb = glm_data.beta.back();

  // We need the rank of the Kt "matrix" now.  Since Kt is
  // assumed to be a vector, we can just assume it has rank 1.

  //! @todo update rX from trace
  //double trn = S * (SNPtSNP - SNPtXGXtSNP);
  rX++;

  // Compute "V2" now that we know rX == rank(X) == rank(X^T X)
  glm_data.V2--;


  // F = Kb' * inv(Kt * G * Kt') * Kb * V2 / (rK * ErrorSS);
  glm_data.F = (1.0 / S) * Kb * Kb * glm_data.V2 / glm_data.ErrorSS;
}

void plm(
	 // inputs
	 const FortranMatrix &XtXi, 
	 const double *XtSNP,
	 const double SNPtSNP, 
	 const double SNPty, 
	 const double yty, 
	 const vector<double> &Xty, 
	 const unsigned rX,
	 float *F,
	 double ErrorSS,
	 unsigned V2
	 //const GLMData &glm_data,
	 // output
	 //GLMData& glm_data_new
	 )
{
  // previous V2 in glm_data

  //int m  = X.get_n_rows(), n = X.get_n_cols();
  int n = XtXi.get_n_rows();

  // G = XtXi
  // compute transpose of SNPtXG: 1xn
  vector<double> GtXtSNP(n, 0.0);

  //! @todo use cblas_dsymv for this
  cblas_dgemv(CblasColMajor,
	      CblasTrans, //! @todo G is symmetric; don't need transpose
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

  writeD("GtXtSNP.dat", GtXtSNP); // ok


  // compute SNPtXGXtSNP (scalar)
  // <SNPtX GtXtSNP> == <XtSNP GtXtSNP>
  double SNPtXGXtSNP = cblas_ddot(n, &XtSNP[0], 1, &GtXtSNP[0], 1);

  // compute S = Schur complement of partitioned matrix to invert
  double S = SNPtSNP - SNPtXGXtSNP;
  if(!S){ //! @todo if zero within tolerance
    // bad news
    *F = 0.0;
    return;
  }

  S = 1.0 / S;

  // compute snpty - snptXGXty = snptMy == scalar
  // already know snpty, snptXG', Xty
  double SNPtMy = -cblas_ddot(n, &GtXtSNP[0], 1, &Xty[0], 1);
  SNPtMy += SNPty;

  double SSM = SNPtMy * SNPtMy * S;
  V2--;
  ErrorSS = ErrorSS - SSM;
  *F = V2 * SSM / ErrorSS;
}

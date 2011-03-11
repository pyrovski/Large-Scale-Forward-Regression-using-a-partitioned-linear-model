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
#include "inner.h" // Inner product of two vectors

using namespace std;

// This file contains data structures/routines for handling the
// internals of the GLM function originally written in Matlab.


// Return values for the GLM, as per the original Matlab routine
struct GLMData
{
  // Solution to the pseudoinverse problem X * beta = y
  std::vector<double> beta;

  // 1 - fcdf(F, V1, V2), as computed in the glm function
  double p;

  // y'*y - beta'*Xy, as computed in the glm function
  double ErrorSS;

  // length(y) - rank(X), as computed in the glm function
  int V2;
};






// The glm function declares the input matrix operators X,
// y, and Kt constant to indicate to the user that they will
// not be changed internally by the routine. 
//
// Note: There are 2 versions of this function, one which assumes Kt
// is a vector, and one which assumes it is a matrix
void glm(const FortranMatrix& X,
	 const std::vector<double>& y,
	 const FortranMatrix& Kt,
	 GLMData& glm_data /*return value, to be filled in*/)
{
  // Debugging: identify ourself.
  // std::cout << "Calling glm function with Kt matrix." << std::endl;
  
  // Debugging: print the Kt
  // Kt.print("Kt at start of glm");
  
  int n  = y.size();
  int V1 = Kt.get_n_rows();

  // Debugging
  // std::cout << "V1=" << V1 << std::endl;
  
  // Compute the matrix-vector product, XTy := X' * y.  
  std::vector<double> XTy = matvec(X, y, /*transX=*/true);
  
  // Create the matrix X^T * X  via multiplication.  X^T * X will
  // actually be used in the rest of the computations, X is no longer
  // needed...
  FortranMatrix XTX = matmat(X, X, /*transA=*/true, /*transB=*/false);
  
  // Solve (X^T * X)*beta = X^T*y for beta.  Note that X and X^T * X
  // have the same rank.

  // Initialize SVD components, A = U * S * V^T
  std::vector<double> S;
  FortranMatrix U, VT;

  // Create the SVD of X^T * X 
  svd_create(XTX, U, S, VT);

  // Apply the SVD to obtain beta
  int rX = svd_apply(U, S, VT, /*result=*/glm_data.beta, XTy);

  // Debugging
//  std::cout << "rows(X^T * X)="
//	    << XTX.get_n_rows()
//	    << ", cols(X^T * X)="
//	    << XTX.get_n_cols()
//	    << ", rX=" << rX << std::endl;

  // Debugging: print beta
  //std::cout << "beta=" << std::endl;
  //for (unsigned i=0; i<glm_data.beta.size(); ++i)
  //  std::cout << glm_data.beta[i] << std::endl;
  
  // Compute "V2" now that we know rX == rank(X) == rank(X^T X)
  glm_data.V2 = n - rX;

  // Debugging: print V2
  // std::cout << "glm_data.V2=" << glm_data.V2 << std::endl;

  // Compute ErrorSS for return in the GLMData data structure
  glm_data.ErrorSS = inner(y,y) - inner(glm_data.beta, XTy);

  // Debugging: print ErrorSS
  // std::cout << "glm_data.ErrorSS=" << glm_data.ErrorSS << std::endl;

  ////////////////////////////////////////////////////////////////////////////////
  // Up to here, the matrix-Kt and vector-Kt versions of the glm
  // function are identical...
  
  // Compute the matrix-vector product, Kb = Kt * beta.
  std::vector<double> Kb = matvec(Kt, glm_data.beta);

  
  // Debugging:
  //std::cout << "Kb=" << std::endl;
  //for (unsigned i=0; i<Kb.size(); ++i)
  //  std::cout << Kb[i] << std::endl;

  // Debugging:
  //Kt.print("before rank(Kt)");
  
  // We need the rank of the Kt "matrix" now.  In this version,
  // we assume Kt is *not* a vector, and therefore a full rank
  // (involving the SVD) calculation is needed...
  int rK = rank(Kt);
  
  // Debugging:
  // Kt.print("after rank(Kt)");
  
  // Debugging (this should be=1 if Kt is a vector...)
  // std::cout << "rK=" << rK << std::endl;
  if (rK==0)
    {
      std::cerr << "Error! Invalid Kt matrix passed with rank 0!" << std::endl;
      exit(1);
    }

  // F = Kb' * inv(Kt * G * Kt') * Kb * V2 / (rK * ErrorSS);
  //
  // Note Kb is always a vector, so F is always a scalar, regardless
  // of the size of Kt, G, etc.

  // 1.) Re-use the SVD of X^T * X to compute B = G * Kt', for the
  // temporary matrix B.  
  FortranMatrix B;
  svd_apply(U, S, VT, /*result=*/B, Kt, /*transKt=*/true);

  // Debugging: print information about B
  //  B.print("B");
  //  std::cout << "rows(B)="
  //	    << B.get_n_rows()
  //	    << ", cols(B)="
  //	    << B.get_n_cols()
  //	    << std::endl;
  

  // 2.) Compute the matrix-matrix product C = Kt * B, using the B
  // from the step above.  
  FortranMatrix C = matmat(Kt, B);

  // Debugging:
  // C.print("C");

  // 3.) Solve for gamma in the equation inv(C) * Kb = gamma, for 'gamma'.
  //
  // We should also verify whether C is always of full rank,
  // always rank-deficient, or sometimes one/sometimes the other.
  // If C is always of full rank, it will probably be faster to use LU
  // or even Cholesky decompositions here.
  std::vector<double> gamma;
  int rC = svd_solve(C, /*result=*/gamma, Kb);
      
  // Debugging/FYI:
//  std::cout << "rows(C)="
//	    << C.get_n_rows()
//	    << ", cols(C)="
//	    << C.get_n_cols()
//	    << ", rC="
//	    << rC
//	    << std::endl;
      
  // .) Compute the inner product of Kb and gamma
  double inner_Kb_gamma = inner(Kb, gamma);

  // Debugging: print scalar result of the inner product between Kb, gamma
  // std::cout << "inner_Kb_gamma=" << inner_Kb_gamma << std::endl;

  // .) Compute F
  double F = inner_Kb_gamma * static_cast<double>(glm_data.V2) / static_cast<double>(rK) / glm_data.ErrorSS; 

  // Debugging: print F, must be positive otherwise the F distribution will return nan
  // std::cout << "F=" << F << std::endl;

  if (F < 0)
    {
      std::cerr << "Error! F<0 obtained, F distribution can only args >= 0." << std::endl;
      exit(1);
    }
  
  // Compute probability, p.  This part requires the GSL.
  glm_data.p = 1. - gsl_cdf_fdist_P(F, static_cast<double>(V1), static_cast<double>(glm_data.V2));
}







// This version of the glm function assumes Kt is a vector.  It makes
// the code a good deal cleaner (but not really any faster when Kt is
// a 1-by-N matrix) than the more general case where Kt is a matrix.
void glm(const FortranMatrix& X,
	 const std::vector<double>& y,
	 const std::vector<double>& Kt,
	 GLMData& glm_data /*return value, to be filled in*/)
{
  // Debugging: identify ourself.
  // std::cout << "Calling glm function with Kt vector." << std::endl;
  
  // Debugging: print the Kt
  // Kt.print("Kt at start of glm");
  
  int n  = y.size();
  int V1 = 1; // Kt is assumed to be a row vector in this version

  // Debugging
  // std::cout << "V1=" << V1 << std::endl;
  
  // Compute the matrix-vector product, XTy := X' * y.  
  std::vector<double> XTy = matvec(X, y, /*transX=*/true);

  // Create the matrix X^T * X  via multiplication.  X^T * X will
  // actually be used in the rest of the computations, X is no longer
  // needed...
  FortranMatrix XTX = matmat(X, X, /*transA=*/true, /*transB=*/false);
  
  // Solve (X^T * X)*beta = X^T*y for beta.  Note that X and X^T * X
  // have the same rank.

  // Initialize SVD components, A = U * S * V^T
  std::vector<double> S;
  FortranMatrix U, VT;

  // Create the SVD of X^T * X 
  svd_create(XTX, U, S, VT);

  // Apply the SVD to obtain beta
  int rX = svd_apply(U, S, VT, /*result=*/glm_data.beta, XTy);

  // Debugging
//  std::cout << "rows(X^T * X)="
//	    << XTX.get_n_rows()
//	    << ", cols(X^T * X)="
//	    << XTX.get_n_cols()
//	    << ", rX=" << rX << std::endl;

  // Debugging: print beta
  //std::cout << "beta=" << std::endl;
  //for (unsigned i=0; i<glm_data.beta.size(); ++i)
  //  std::cout << glm_data.beta[i] << std::endl;
  
  // Compute "V2" now that we know rX == rank(X) == rank(X^T X)
  glm_data.V2 = n - rX;

  // Debugging: print V2
  // std::cout << "glm_data.V2=" << glm_data.V2 << std::endl;

  // Compute ErrorSS for return in the GLMData data structure
  glm_data.ErrorSS = inner(y,y) - inner(glm_data.beta, XTy);

  // Debugging: print ErrorSS
  // std::cout << "glm_data.ErrorSS=" << glm_data.ErrorSS << std::endl;

  ////////////////////////////////////////////////////////////////////////////////
  // Up to here, the matrix-Kt and vector-Kt versions of the glm
  // function are identical...
  
  // Compute the matrix-vector product, Kb = Kt * beta.
  // Note if Kt is actually a vector, Kb will be a scalar...
  double Kb = inner(Kt, glm_data.beta);

  // We need the rank of the Kt "matrix" now.  Since Kt is
  // assumed to be a vector, we can just assume it has rank 1.
  int rK = 1;

  // F = Kb' * inv(Kt * G * Kt') * Kb * V2 / (rK * ErrorSS);
  //
  // Note Kb is a scalar in this version

  // 1.) Re-use the SVD of X^T * X to compute temp = G * Kt', for the
  // temporary result vector, temp.  
  std::vector<double> temp;
  svd_apply(U, S, VT, /*result=*/temp, Kt);
  
  // 2.) Compute the inner product c = (Kt,temp), using 'temp'
  // from the step above.  
  double c = inner(Kt,temp);

  // 3.) Note: for F we need to compute Kb' * inv(c) * Kb, but since c
  // is a scalar, we just divide!
  double F = Kb * Kb / c * static_cast<double>(glm_data.V2) / static_cast<double>(rK) / glm_data.ErrorSS; 

  // Debugging: print F
  // std::cout << "F=" << F << std::endl;

  // F must be positive otherwise the F distribution will return nan
  if (F < 0)
    {
      std::cerr << "Error! F<0 obtained, F distribution can only args >= 0." << std::endl;
      exit(1);
    }
  
  // Compute probability, p.  This part requires the GSL.
  glm_data.p = 1. - gsl_cdf_fdist_P(F, static_cast<double>(V1), static_cast<double>(glm_data.V2));
}


#endif

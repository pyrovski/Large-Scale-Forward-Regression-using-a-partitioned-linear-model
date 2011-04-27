#ifndef __glm_h__
#define __glm_h__

// C++ headers
#include <vector>

// External library headers

// Local headers
#include "fortran_matrix.h"


// This file contains data structures/routines for handling the
// internals of the GLM function originally written in Matlab.


// Return values for the GLM, as per the original Matlab routine
struct GLMData
{
  // Solution to the pseudoinverse problem X * beta = y
  std::vector<double> beta;

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
void glm(unsigned id, unsigned iteration, 
	 int n,
	 FortranMatrix &XtXi, // updated
	 const double *XtSNP,
	 const double SNPtSNP, 
	 const double SNPty, 
	 const double yty, 
	 std::vector<double> &Xty, // updated
	 double rX,
	 // output
	 GLMData& glm_data);


void plm(
	 // inputs
	 const FortranMatrix &X, 
	 const FortranMatrix &XtXi, 
	 const std::vector<double> &XtSNP,
	 const double SNPtSNP, 
	 const double SNPty, 
	 const double yty, 
	 const std::vector<double> &Xty, 
	 const unsigned rX,
	 const GLMData &glm_data,
	 // output
	 GLMData& glm_data_new);

#endif

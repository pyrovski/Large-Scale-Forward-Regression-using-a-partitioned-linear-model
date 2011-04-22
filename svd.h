#ifndef __svd_solve_h__
#define __svd_solve_h__

// C++ headers
#include <vector>
#include <iostream>
#include <limits>

// Local headers
#include "fortran_matrix.h"
#include "my_lapack_blas.h"
#include "print_matrix.h"


void svd_create(/*in */FortranMatrix& A,
		/*out*/FortranMatrix& U,
		/*out*/std::vector<double>& S,
		/*out*/FortranMatrix& VT);

int svd_apply(/*in */const FortranMatrix& U,
	      /*in */const std::vector<double>& S,
	      /*in */const FortranMatrix& VT,
	      /*out*/std::vector<double>& x,
	      /*in */const std::vector<double>& b,
	      /*in,default*/double tol=-1.0);

int svd_apply(/*in */ FortranMatrix& U,
	      /*in */ std::vector<double>& S,
	      /*in */ FortranMatrix& VT,
	      /*out*/ FortranMatrix& X,
	      /*in */ const FortranMatrix& B,
	      /*in */ bool transB=false, // if true, apply the SVD to the transpose of B rather than B itself.
	      /*in,default*/double tol = -1.0);

int svd_solve(FortranMatrix& A,
	      std::vector<double>& x, // solution
	      const std::vector<double>& b, // rhs
	      double tol = -1.0);

#endif

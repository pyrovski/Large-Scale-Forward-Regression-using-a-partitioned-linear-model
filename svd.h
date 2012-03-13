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


void svd_create(/*in */FortranMatrix A,
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

int svd_solve(FortranMatrix A,
	      std::vector<double>& x, // solution
	      const std::vector<double>& b, // rhs
	      double tol = -1.0);

#endif

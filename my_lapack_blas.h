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
#ifndef __my_lapack_blas_h__
#define __my_lapack_blas_h__


// Prototypes for Lapack/BLAS functions we will need.  Assume this will
// be compiled with a C++ compiler

// Use this #define to set the precision and underscore convention of your system. 
// This is the PETSc technique for doing this. 
#define LAPACKgesdd_ dgesdd_
#define BLASgemv_ dgemv_
#define BLASgemm_ dgemm_

// Don't mangle the Lapack function name
extern "C"
{
  // Lapack SVD algorithm using "domain decomposition" algorithm.
  void LAPACKgesdd_( const char* jobz, int* m, int* n, double* a,
		     int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
		     double* work, int* lwork, int* iwork, int* info );

  void BLASgemv_(const char*, int* m, int* n,
		 double* alpha, double* a, int* lda,
		 double* x, int* incx, double* beta,
		 double* y, int* incy);

  // In general, GEMM computes
  // C <- alpha*op( A )*op( B ) + beta*C,
  void BLASgemm_(const char* transa, const char* transb,
		 int* m, int* n, int* k,
		 double* alpha, double* a, int* lda,
		 double* b, int* ldb, double* beta,
		 double* c, int* ldc);  

}



#endif

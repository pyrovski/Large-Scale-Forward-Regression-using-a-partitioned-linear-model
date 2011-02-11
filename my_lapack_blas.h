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

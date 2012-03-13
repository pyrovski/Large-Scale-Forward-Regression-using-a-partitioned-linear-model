#ifndef __rank_h__
#define __rank_h__

// C++ headers
#include <vector>
#include <iostream>

// Local headers
#include "fortran_matrix.h"
#include "my_lapack_blas.h"
#include "print_matrix.h"

// Function for computing the rank of a FortranMatrix.  Uses adaptive
// tolerance if tol==-1. is passed for the tol parameter.  Tolerance
// determines if a singular value is nonzero, this is similar to
// Matlab/Octave.  Note that for computing the rank we need not
// compute the singular vectors (U, V^T) and so we can avoid fully
// allocating those arrays (they are not referenced if the first
// argument to dgesdd is 'N'.)
//
// Note that because it is undefined behavior to take the address of
// the first element of an empty vector, we simply make the u and
// vt vectors here have a single entry...
//
// Because Lapack overwrites the input A vector (at least, my implementation
// of Lapack does) even if we pass 'N' as the first argument to gesdd,
// we make a copy of 'a' before performing the SVD.  This is in line with
// the behavior of the similarly named Matlab function, which does not
// modify 'a'.
int rank(const FortranMatrix& a, double tol=-1.)
{
  // Make the copy!
  FortranMatrix mutable_a = a;
  
  // Number of rows/cols in a
  int m = mutable_a.get_n_rows();
  int n = mutable_a.get_n_cols();

  // "Leading dimension" variables
  int lda=m, ldu=m, ldvt=n;

  // Work array sizes and info variables
  int info, lwork;
  double wkopt;

  // iwork dimension should be at least 8*min(m,n) per Lapack manual
  std::vector<int> iwork(8*n);

  // Note: U is nominally M-by-M, V^T is nominally N-by-N, here we are
  // not requesting eigenvectors and therefore we only allocate them to
  // have a size of 1 each to save memory.
  std::vector<double> s(n), u(1 /*ldu*m*/), vt(1 /*ldvt*n*/), work;

  // if -1 on input, gesdd returns the optimal lwork value in "wkopt" 
  lwork = -1; 
  
  LAPACKgesdd_( "No vectors", // 'A', 'S', 'O', 'N': all, first min(m,n), overwrite, or none 
		&m, // number of rows of A (pointer) 
		&n, // number of cols of A (pointer) 
		&(mutable_a.values[0]),  // array, dimension lda-by-n, may be overwritten if 'O' is passed as 1st arg 
		&lda,   // leading dim. of A aka number of rows (pointer) 
		&s[0],      // singular values array 
		&u[0],      // m-by-m matrix of orthonormal "output" vectors, eigenvectors of A A^T 
		&ldu,   // leading dim of U (pointer) 
		&vt[0],     // n-by-n matrix of orthonormal "input" vectors, eigenvectors of A^T A 
		&ldvt,  // leading dimension of V, aka number of rows (pointer) 
		&wkopt, // double (pointer) to the lwork array, or if lwork==-1, on output contains the optimal lwork 
		&lwork, // if -1 on input, the optimal lwork value is returned, in general, larger values==better performance  
		&iwork[0],  // INTEGER array, dimension (8*min(M,N)) 
		&info   // 0==success, -i==argument i was illegal, >0 == failure to converge 
		);

  //  if (info == 0)
  //    {
  //      std::cout << " Optimal lwork size returned in wkopt=" << wkopt << std::endl;
  //    }
  
  lwork = static_cast<int>(wkopt);
  work.resize(lwork);
  
  /* Compute SVD, this time we actually do the computation */
  LAPACKgesdd_( "No vectors",
		&m,
		&n,
		&(mutable_a.values[0]),
		&lda,
		&s[0],
		&u[0],
		&ldu,
		&vt[0],
		&ldvt,
		&work[0],   /* pointer to the now optimally-allocated work array */
		&lwork, /* integer (pointer) The (optimal) current size of the work array */
		&iwork[0],
		&info
		);
  
  /* Check for convergence */
  if( info > 0 ) {
    std::cout << "The algorithm computing SVD failed to converge." << std::endl;
    exit( 1 );
  }

  // Debugging
  // print_matrix( "Singular values", 1, n, &s[0], 1 );

  // Compute adaptive tolerance ala Matlab/Octave:
  // tol = max (M,N) * s[0] * eps;
  // where eps is machine precision
  // Note that the singular values are real, non-negative,
  // and in descending order, so we don't have to do things like
  // check absolute values.
  if (tol == -1.)
    tol = std::max(m,n) * s[0] * std::numeric_limits<double>::epsilon();

  // std::cout << "Using tol=" << tol << std::endl;

  // Loop, count singular values larger than tol
  int rank=0;
  for (unsigned i=0; i<s.size(); ++i)
    if (s[i] > tol)
      rank++;

  return rank;
}

#endif

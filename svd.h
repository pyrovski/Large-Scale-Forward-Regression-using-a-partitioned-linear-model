#ifndef __svd_solve_h__
#define __svd_solve_h__

// C++ headers
#include <vector>
#include <iostream>

// Local headers
#include "fortran_matrix.h"
#include "my_lapack_blas.h"
#include "print_matrix.h"

// Use this routine to create the U, S, and V^T matrices.  The results
// can then be used with the svd_apply routine, to solve a
// system of equations Ax=b, or to compute the action of A^+ M on another
// matrix M.  
void svd_create(/*in */FortranMatrix& A,
		/*out*/FortranMatrix& U,
		/*out*/std::vector<double>& S,
		/*out*/FortranMatrix& VT)
{
  // Number of rows/cols in A
  int m = A.get_n_rows();
  int n = A.get_n_cols();

  // "Leading dimension" variables
  int lda=m, ldu=m, ldvt=n;

  // Work array sizes and info variables
  int info, lwork;
  double wkopt;

  // iwork dimension should be at least 8*min(m,n) per Lapack manual
  std::vector<int> iwork(8*n);

  // A = U         S      V^T
  //     (mxm) * (mxn) * (nxn)
  // U is M-by-M, V^T is N-by-N, these are stored as FortranMatrix objects,
  // so that we can conveniently perform operations on them subsequently...
  // FortranMatrix U(ldu,m), VT(ldvt,n);
  U.resize(ldu,m);
  VT.resize(ldvt,n);
  
  // Array of singular values
  S.resize(n);

  // if -1 on input, gesdd returns the optimal lwork value in "wkopt" 
  lwork = -1; 
  
  LAPACKgesdd_( "All vectors", // 'A', 'S', 'O', 'N': all, first min(m,n), overwrite, or none 
		&m, // number of rows of A (pointer) 
		&n, // number of cols of A (pointer) 
		&(A.values[0]),  // array, dimension lda-by-n, may be overwritten if 'O' is passed as 1st arg 
		&lda,   // leading dim. of A aka number of rows (pointer) 
		&S[0],      // singular values array 
		&(U.values[0]),      // m-by-m matrix of orthonormal "output" vectors, eigenvectors of A A^T 
		&ldu,   // leading dim of U (pointer) 
		&(VT.values[0]),     // n-by-n matrix of orthonormal "input" vectors, eigenvectors of A^T A 
		&ldvt,  // leading dimension of V, aka number of rows (pointer) 
		&wkopt, // double (pointer) to the lwork array, or if lwork==-1, on output contains the optimal lwork 
		&lwork, // if -1 on input, the optimal lwork value is returned, in general, larger values==better performance  
		&iwork[0],  // INTEGER array, dimension (8*min(M,N)) 
		&info   // 0==success, -i==argument i was illegal, >0 == failure to converge 
		);

  if (info < 0)
    {
      std::cerr << "Argument " << -info << " was illegal." << std::endl;
      exit(1);
    }
    
    //  if (info == 0)
    //    {
    //      /* For this example problem, wkopt should be 332 */
    //      std::cout << " Optimal lwork size returned in wkopt=" << wkopt << std::endl;
    //    }
    
  lwork = static_cast<int>(wkopt);

  // Work array
  std::vector<double> work(lwork);
  
  /* Compute SVD, this time we actually do the computation */
  LAPACKgesdd_( "All vectors",
		&m,
		&n,
		&(A.values[0]),
		&lda,
		&S[0],
		&(U.values[0]),
		&ldu,
		&(VT.values[0]),
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
  // print_matrix( "Singular values", 1, S.size(), &S[0], 1 );
}



// Use this routine with data created by the svd_create() routine.
// The output x will be the solution to the linear system Ax=b,
// including singular values up to the tolerance, tol.  As a
// by-product of applyin the SVD, the approximate rank (to within tol)
// of the matrix is returned.
int svd_apply(/*in */FortranMatrix& U,
	      /*in */std::vector<double>& S,
	      /*in */FortranMatrix& VT,
	      /*out*/std::vector<double>& x,
	      /*in */const std::vector<double>& b,
	      /*in,default*/double tol=-1.0)
{
  int m = U.get_n_rows(); // U is mxm
  int n = VT.get_n_rows(); // V^T is nxn

  // Debugging:
  // std::cout << "m=" << m << std::endl;
  // std::cout << "n=" << n << std::endl;
  
  // Check the size of the rhs for compatibility
  //   A   *   x   =   b
  // (mxn)   (nx1)   (mx1)
  if (b.size() != m)
    {
      std::cerr << "Error! Number of rows in A must match the length of b!" << std::endl;
      exit (1);
    }
  
  if (tol == -1.)
    tol = std::max(m,n) * S[0] * std::numeric_limits<double>::epsilon();

  // Debugging
  //  std::cout << "Using tol=" << tol << std::endl;


  // Compute y     = U^T   * b
  //         (mx1) = (mxm) * (mx1)
  std::vector<double> y = matvec(U, b, /*transU=*/true);

  // Scale by *inverse* singular values, unless they are below the
  // tolerance, in which case, zero the entries.  Also, count up the
  // rank of the matrix.
  //
  // S always has size n.  In cases where n<m, we must only loop up to
  // n (although y has size m).  If n<m the additional entries of y will be
  // "chopped" off by the resize step below...
  int rank=0;
  for (unsigned i=0; i<std::min(m,n); ++i)
    {
      // Debugging
      // std::cout << "S[" << i << "]=" << S[i] << std::endl;
      
      if (S[i] <= tol)
	y[i] = 0.;
      else
	{
	  y[i] /= S[i];
	  rank++;
	}
    }

  // Chop/off add entries to y so that it now has size n.
  // This is a consequence of multiplying by S^-1, which is
  // a diagonal matrix of size (nxm)
  y.resize(n);

  // Debugging
  // print_matrix( "U^T * b", 1, y.size(), &y[0], 1 );

  // Now compute   x   =   V   *   y
  //             (nx1) = (nxn) * (nx1)
  // Note that we have V^T, so we have to use the transposed matvec again...
  x = matvec(VT, y, /*transVT=*/true);

  return rank;
}



// Generalized version of the svd_apply routine which operates on an
// m-by-r matrix B and returns an n-by-r matrix X as the result.  The
// svd_apply function which operates on vectors could be replaced by
// this generalized function.  Also computes the tol-dependent rank
// as a by-product of applying the SVD.
int svd_apply(/*in */ FortranMatrix& U,
	      /*in */ std::vector<double>& S,
	      /*in */ FortranMatrix& VT,
	      /*out*/ FortranMatrix& X,
	      /*in */ const FortranMatrix& B,
	      /*in */ bool transB=false, // if true, apply the SVD to the transpose of B rather than B itself.
	      /*in,default*/double tol=-1.0)
{
  int m = U.get_n_rows(); // U is mxm
  int n = VT.get_n_rows(); // V^T is nxn

  // The number of rows in op(B), must be m to continue
  int q = (transB ? B.get_n_cols() : B.get_n_rows());
  
  // The number of columns in op(B)
  int r = (transB ? B.get_n_rows() : B.get_n_cols());

  // Debugging:  In all cases I've seen, S.size() == n, NOT m
  // std::cout << "S.size()=" << S.size() << std::endl;
  // std::cout << "m=" << m << std::endl;
  // std::cout << "n=" << n << std::endl;
  
  // Check the size of B for compatibility
  //   X   =    V   *  S^+  *  U^T  *   B
  // (nxr)    (nxn)   (nxm)   (mxm)   (mxr)
  if ( q != m )
    {
      std::cerr << "Error! The number of rows of op(B), " << q
		<< ", does not match the number of cols of U^T, " << m << "!" << std::endl;
      exit (1);
    }

  
  // Compute U^T * B, Result will be m-by-r
  FortranMatrix UTB = (transB
		       ? matmat(U, B, /*transA=*/true, /*transB=*/true)
		       : matmat(U, B, /*transA=*/true, /*transB=*/false)
		       );

  // Debugging
  // print_matrix( "U^T * B", UTB.get_n_rows(), UTB.get_n_cols(), &(UTB.values[0]), UTB.get_n_rows() );

  // Compute tolerance if necessary
  if (tol == -1.)
    tol = std::max(m,n) * S[0] * std::numeric_limits<double>::epsilon();

  // Debugging:
  // std::cout << "Tolerance = " << tol << std::endl;
  
  // Next, scale each row of U^T B by 1/S[i] if S[i] is larger than
  // the tolerance, otherwise, zero it.  Also compute the approximate
  // rank at this point.  In all cases I've tested, S has n entries.
  // If m > n, we must only loop up to n to avoid going past the end
  // of the S array.  If m < n, then the last n-m entries of S are
  // identically zero.

  // Debugging:
  //if (m < n)
  //  for (unsigned i=0; i<S.size(); ++i)
  //    std::cout << "S[" << i << "]=" << S[i] << std::endl;
  
  int rank = 0;
  for (unsigned i=0; i<std::min(m,n); ++i)
    {
      // Debugging:
      //std::cout << "S[" << i << "]=" << S[i] << std::endl;
      
      if (S[i] <= tol)
	{
	  for (unsigned j=0; j<r; ++j)
	    {
	      UTB(i,j) = 0.;
	    }
	}
      else
	{
	  rank++;
	  for (unsigned j=0; j<r; ++j)
	    {
	      UTB(i,j) /= S[i];
	    }
	}
    }

  // Debugging
  // print_matrix( "S^+ * U^T * B", UTB.get_n_rows(), UTB.get_n_cols(), &(UTB.values[0]), UTB.get_n_rows() );

  // Note: multiplication by S^+ leads to the nxr matrix S^+ U^T B,
  // our matrix currently still has m logical rows....so we'd like to
  // call matmat(VT, UTB, true, false) directly, but the sizes won't be compatible.
  //
  // If n>m, we can just add extra zero "rows" to UTB to make the
  // sizes match.  If n <= m then only the first m rows of V will
  // take part in the multiplication.  We can handle this by looping the
  // middle index (in this case, j) only to min(m,n).
  //
  // FIXME: Do this in BLAS GEMM with LDA tricks?
  X.resize(n,r);

  // Compute "middle" index as the smaller of n (n. rows of S^+) and m (n. cols of S^+)
  unsigned jmax = std::min(m,n);
  
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<jmax; ++j)
      for (unsigned int k=0; k<r; ++k)
	{
	  X(i,k) += VT(j,i) * UTB(j,k); // Note: we have VT, want V, hence VT(j,i)
	}
  
  
//  if (n>m)
//    {
//      std::cout << "n>m"  << std::endl;
//      //exit(1);
//      for (unsigned int i=0; i<n; ++i)
//	for (unsigned int j=0; j<m/*min(m,n)*/; ++j)
//	  for (unsigned int k=0; k<r; ++k)
//	    {
//		X(i,k) += VT(j,i) * UTB(j,k); // Note: we have VT, want V, hence VT(j,i)
//	    }
//    }
//  else
//    {
//      std::cout << "n<=m"  << std::endl;
//      for (unsigned int i=0; i<n; ++i)
//	for (unsigned int j=0; j<n; ++j)
//	  for (unsigned int k=0; k<r; ++k)
//	    X(i,k) += VT(j,i) * UTB(j,k); // Note: we have VT, want V, hence VT(j,i)
//    }

  
  return rank;
}



// Function for solving the system Ax = b using the pseudoinverse, as
// determined by the SVD.  Uses adaptive tolerance if tol==-1.0 is
// passed for the tol parameter.  Tolerance determines if a singular
// value is nonzero, this is similar to Matlab/Octave.  As a by-product
// of the SVD solve, the approximate rank (to within tol) of the
// matrix is returned.  Use this routine if you don't need to
// reuse the SVD for any additional vectors and don't want to manage
// the storage of the U, S, and V^T matrices yourself...
int svd_solve(FortranMatrix& A,
	      std::vector<double>& x, // solution
	      const std::vector<double>& b, // rhs
	      double tol=-1.0)
{
  // SVD components
  std::vector<double> S;
  FortranMatrix U, VT;

  // Create the SVD...
  svd_create(A, U, S, VT);

  // And use it to solve Ax=b
  return svd_apply(U, S, VT, x, b, tol);
}


#endif

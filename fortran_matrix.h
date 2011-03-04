#ifndef __fortran_matrix_h__
#define __fortran_matrix_h__

#include <stdlib.h>

// C++ includes
#include <iostream>
#include <vector>

// Local project includes
#include "my_lapack_blas.h"
#include "print_matrix.h"

// This Matrix is stored in column-major order.  That way it can be
// created as normal in C++, and passed to Fortran LAPACK/BLAS routines.  Note:
// it will be more efficient to fill this matrix column-by-column...
class FortranMatrix
{
public:
  // Constructor, by default, an empty matrix is constructed
  FortranMatrix(unsigned nr=0, unsigned nc=0) : n_rows(nr), n_cols(nc), values(nr*nc) {}

  //! @todo operator = or FortranMatrix(const FortranMatrix &)

  // Returns a writable reference to the (i,j) entry of the matrix
  double& operator()(unsigned i, unsigned j)
  {
    return values[i + n_rows*j];
  }

  // The print function is logically const.
  void print(std::string title="") const
  {
    // Uses this for backwards-compatibility, since I started off using print_matrix
    // before creating the FortranMatrix object....
    print_matrix( title.c_str(),
		  this->get_n_rows(),
		  this->get_n_cols(),
		  &(this->values[0]),
		  this->get_n_rows() );    
  }



  // Replace this matrix with its transpose.  Here we simply
  // use n_rows*n_cols temporary storage.  In-place transposition
  // of a non-square matrix is a non-trivial algorithm.
  // http://en.wikipedia.org/wiki/In-place_matrix_transposition
  void transpose_self()
  {
    std::vector<double> transpose(n_cols*n_rows);

    for (unsigned int j=0; j<n_cols; ++j)
      for (unsigned int i=0; i<n_rows; ++i)
	{
	  // transpose(j,i) = values(i,j)
	  // The new matrix is still column major, but has
	  // "n_cols" rows.
	  transpose[j + n_cols*i] = 
	    values[i + n_rows*j];
	}

    values.swap(transpose);
    std::swap(n_rows, n_cols);
  }


  // Change this matix to have new_n_rows rows and new_n_cols columns
  // Do not rely on any previous values in the matrix for this routine.
  void resize(unsigned new_n_rows, unsigned new_n_cols)
  {
    n_rows = new_n_rows;
    n_cols = new_n_cols;
    values.resize(n_rows*n_cols);
  }
  
  // retain old values in correct coordinates
  void resize_retain(unsigned new_n_rows, unsigned new_n_cols){
    this->resize(new_n_rows, new_n_cols);
    //! @todo finish
    //for(unsigned 
  }

  // Return a writable reference to the number of rows/cols of the matrix.
  // You should only change this if you know what you are doing and have also
  // changed the "values" array in a consistent manner.
  unsigned& get_n_rows() { return n_rows; }
  unsigned& get_n_cols() { return n_cols; }

  // Const version of the above, just returns a copy of the number of
  // rows/cols, does not allow anything to be modified.
  unsigned get_n_rows() const { return n_rows; }
  unsigned get_n_cols() const { return n_cols; }


  //
  // Data
  //
  
  // Users have direct access to the array of values
  std::vector<double> values;
  
private: 
  unsigned n_rows, n_cols;
};






// Computes and returns (by value) the matrix containing the product
// A * B, possibly with transposing either or both matrices.  By default
// neither are transposed.  Although BLAS knows nothing about const arguments,
// A and B are not overwritten by dgemm and this function may therefore be
// called with const operators.
FortranMatrix matmat(const FortranMatrix& A, const FortranMatrix& B,
		     bool transA=false, bool transB=false)
{

  // number of rows of op(A)
  int M = transA ? A.get_n_cols() : A.get_n_rows();

  // number of cols of op(B)
  int N = transB ? B.get_n_rows() : B.get_n_cols();

  // number of cols of op(A)
  int K = transA ? A.get_n_rows() : A.get_n_cols();

  // number of rows of op(B) [not needed by BLAS]
  int L = transB ? B.get_n_cols() : B.get_n_rows();
    
  // For compatibility, the number of cols of op(A) must match the number of rows of op(B)
  if (K != L)
    {
      std::cerr << "Error! Matrix inner dimensions must match!" << std::endl;
      exit(1);
    }
      
  // Routine actually computes C <- alpha*op( A )*op( B ) + beta*C
  double alpha = 1., beta = 0.;
  
  FortranMatrix result(M,N);

  // Not sure if Fortran needs a NULL-terminated string here or if we can just
  // do the ternary operator "in place" or if it has to be a separate variable...
  char
    char_transA = transA ? 'T' : 'N',
    char_transB = transB ? 'T' : 'N';

  // Call BLAS routine
  BLASgemm_(&char_transA,
	    &char_transB,
	    &M, // number of rows of op(A)
	    &N, // number of cols of op(B)
	    &K, // number of cols of op(A)
	    &alpha,
	    const_cast<double*>(&(A.values[0])), // Cast away const-ness, BLAS will NOT overwrite
	    transA ? &K : &M, // leading dimension of A (see BLAS docs)
	    const_cast<double*>(&(B.values[0])), // Cast away const-ness, BLAS will NOT overwrite
	    transB ? &N : &K, // leading dimension of B (see BLAS docs)
	    &beta,
	    &(result.values[0]),
	    &M  // at least max(1,M)
	    );

  // Any extra copy will hopefully be removed by RVO.  
  return result;
}




// Computes y     =   A     * x
//          (mx1) = (mxn)   * (nx1)
// Or
//          y     =   A^T * x
//          (nx1) = (nxm) * (mx1)
//
// depending on the transA paramter (default false).  We use the BLAS
// dgemv routine for this.  This 'free' version of the matvec function
// returns its result by value, but hopefully this is taken care of by
// RVO.
std::vector<double> matvec(const FortranMatrix& A, const std::vector<double>& x,
			   bool transA=false)
{
  // BLAS needs ints, not unsigneds
  int m = A.get_n_rows(), n = A.get_n_cols();

  // Check to be sure that the size of x matches either the number of rows or
  // columns of A, depending on transA.
  if ( ((transA) && (x.size() != m)) || ((!transA) && (x.size() != n)) )
    {
      std::cerr << "Error! Sizes of A and x do not match!" << std::endl;
      exit(1);
    }

  // Create solution vector
  std::vector<double> y ( transA ? n : m );

  // Routine actually computes y <- (alpha) * A * x + (beta) * y
  double alpha = 1., beta = 0.;

  // The "increment" for the elements of x and y
  int incx = 1, incy = 1;

  // Temporary char we can take the address of to pass to BLAS.
  char
    char_transA = transA ? 'T' : 'N';

  BLASgemv_(&char_transA,
	    &m,
	    &n,
	    &alpha,
	    const_cast<double*>(&A.values[0]),
	    &m, // leading dimension of A
	    const_cast<double*>(&x[0]),
	    &incx,
	    &beta,
	    &y[0],
	    &incy);

  return y;
}


#endif

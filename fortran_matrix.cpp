
#include <stdlib.h>
#include <string.h>
#include <strings.h>

// C++ includes
#include <iostream>
#include <vector>
#include <iostream>

extern "C"{
#include <cblas.h>
}

// Local project includes
#include "fortran_matrix.h"
#include "my_lapack_blas.h"
#include "print_matrix.h"

using namespace std;

// This Matrix is stored in column-major order.  That way it can be
// created as normal in C++, and passed to Fortran LAPACK/BLAS routines.  Note:
// it will be more efficient to fill this matrix column-by-column...

  // Constructor, by default, an empty matrix is constructed
FortranMatrix::FortranMatrix(uint64_t nr, uint64_t nc) : 
  n_rows(nr), n_cols(nc), values(nr*nc) 
{}


FortranMatrix::FortranMatrix() : 
  n_rows(0), n_cols(0), values(0) 
{}

const FortranMatrix & FortranMatrix::operator = (const FortranMatrix &rhs){
    this->values = rhs.values;
    this->n_rows = rhs.n_rows;
    this->n_cols = rhs.n_cols;
    return *this;
  }

  void FortranMatrix::add(const FortranMatrix &rhs){
    if(rhs.n_rows != n_rows || rhs.n_cols != n_cols){
      cerr << "invalid matrix dimensions" << endl;
      exit(1);
    }
    for(uint64_t col = 0; col < n_cols; col++)
      for(uint64_t row = 0; row < n_rows; row++)
	values[row + n_rows * col] += rhs.values[row + n_rows * col];
  }


  // The print function is logically const.
  void FortranMatrix::print(std::string title="") const
  {
    // Uses this for backwards-compatibility, since I started off using print_matrix
    // before creating the FortranMatrix object....
    print_matrix( title.c_str(),
		  this->get_n_rows(),
		  this->get_n_cols(),
		  &(this->values[0]),
		  this->get_n_rows() );    
  }

  int FortranMatrix::write(string filename){
    write_matrix(filename.c_str(),
		 n_rows,
		 n_cols,
		 &values[0],
		 n_rows);
    return 0;
  }


#ifdef _DEBUG
  int FortranMatrix::writeD(string filename){return this->write(filename);}
#else
  int FortranMatrix::writeD(string filename){return 0;}
#endif

void FortranMatrix::transpose_dims(){
    std::swap(n_rows, n_cols);
}

  // Replace this matrix with its transpose.  Here we simply
  // use n_rows*n_cols temporary storage.  In-place transposition
  // of a non-square matrix is a non-trivial algorithm.
  // http://en.wikipedia.org/wiki/In-place_matrix_transposition
  void FortranMatrix::transpose_self()
  {
    std::vector<double> transpose(n_cols*n_rows);

    for (uint64_t j=0; j<n_cols; ++j)
      for (uint64_t i=0; i<n_rows; ++i)
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
  void FortranMatrix::resize(uint64_t new_n_rows, uint64_t new_n_cols)
  {
    n_rows = new_n_rows;
    n_cols = new_n_cols;
    values.resize(n_rows*n_cols);
  }
  
  // retain old values in correct coordinates
  void FortranMatrix::resize_retain(uint64_t new_n_rows, uint64_t new_n_cols){
    uint64_t old_n_rows = n_rows, old_n_cols = n_cols;
    if(new_n_rows < old_n_rows || new_n_cols < old_n_cols){
      cerr << "invalid resize; fixme" << endl;
      exit(1);
    } else if (new_n_rows == old_n_rows && new_n_cols == old_n_cols)
      return;
    
    this->resize(new_n_rows, new_n_cols);
    // first column is already in the right place
    for(uint64_t col = old_n_cols - 1; col > 0; col--){
      //! @todo this could be faster; memmove uses temporary storage
      memmove(&values[new_n_rows * col], &values[old_n_rows * col], old_n_rows * sizeof(double));
      /*
      double *base = &values[old_n_rows * col];
      double *dest = &values[new_n_rows * col];
      for(uint64_t row = old_n_rows - 1; row >= 0; row--)
	dest[row] = base[row];
      */
      // clear newly resized space (rest of column)
      bzero(&values[new_n_rows * col + old_n_rows], 
	    sizeof(double) * (new_n_rows - old_n_rows));
    }
    // zero first column
    bzero(&values[old_n_rows], 
	  sizeof(double) * (new_n_rows - old_n_rows));
  }

  // Return a writable reference to the number of rows/cols of the matrix.
  // You should only change this if you know what you are doing and have also
  // changed the "values" array in a consistent manner.

/*
  uint64_t& FortranMatrix::get_n_rows() { return n_rows; }
  uint64_t& FortranMatrix::get_n_cols() { return n_cols; }
*/

  //
  // Data
  //
  



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
  
  FortranMatrix resultasdf;
  resultasdf.resize(M, N);

  // Not sure if Fortran needs a NULL-terminated string here or if we can just
  // do the ternary operator "in place" or if it has to be a separate variable...
  CBLAS_TRANSPOSE
    cblas_transA = transA ? CblasTrans : CblasNoTrans,
    cblas_transB = transB ? CblasTrans : CblasNoTrans;

  // Call BLAS routine
  cblas_dgemm(CblasColMajor,
	      cblas_transA,
	      cblas_transB,
	      M, // number of rows of op(A)
	      N, // number of cols of op(B)
	      K, // number of cols of op(A)
	      alpha,
	      const_cast<double*>(&(A.values[0])), // Cast away const-ness, BLAS will NOT overwrite
	      transA ? K : M, // leading dimension of A (see BLAS docs)
	      const_cast<double*>(&(B.values[0])), // Cast away const-ness, BLAS will NOT overwrite
	      transB ? N : K, // leading dimension of B (see BLAS docs)
	      beta,
	      &resultasdf.values[0],
	      M  // at least max(1,M)
	      );

  // Any extra copy will hopefully be removed by RVO.  
  return resultasdf;
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

#ifdef _DEBUG
void writeD(string filename, const vector<double> &v){
  write_matrix(filename.c_str(), v.size(), 1, &v[0], 1);
}
void writeD(string filename, const double *v, unsigned length){
  write_matrix(filename.c_str(), length, 1, &v[0], 1);
}
#else
void writeD(string filename, const vector<double> &v){
}
void writeD(string filename, const double *v, unsigned length){
}
#endif
void write(string filename, const vector<double> &v){
  write_matrix(filename.c_str(), v.size(), 1, &v[0], 1);
}


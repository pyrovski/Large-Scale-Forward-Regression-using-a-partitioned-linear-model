#ifndef __fortran_matrix_h__
#define __fortran_matrix_h__

#include <string>
#include <vector>
#include <stdint.h>
extern "C"{
#include <cblas.h>
}

#include "print_matrix.h"

class FortranMatrix
{
public:
  // Constructor, by default, an empty matrix is constructed
  FortranMatrix(uint64_t nr, uint64_t nc);
  FortranMatrix();

  void add(const FortranMatrix &rhs);

  const FortranMatrix & operator = (const FortranMatrix &rhs);
  
  // Returns a writable reference to the (i,j) entry of the matrix
  inline double& operator()(uint64_t i, uint64_t j)
  {
    return values[i + n_rows*j];
  }

  void print(std::string title) const;
  int write(std::string filename);
  int writeD(std::string filename);

  // Replace this matrix with its transpose.  Here we simply
  // use n_rows*n_cols temporary storage.  In-place transposition
  // of a non-square matrix is a non-trivial algorithm.
  // http://en.wikipedia.org/wiki/In-place_matrix_transposition
  void transpose_self();
  void transpose_dims();

  // Change this matix to have new_n_rows rows and new_n_cols columns
  // Do not rely on any previous values in the matrix for this routine.
  void resize(uint64_t new_n_rows, uint64_t new_n_cols);

  // retain old values in correct coordinates
  void resize_retain(uint64_t new_n_rows, uint64_t new_n_cols);

  // Return a writable reference to the number of rows/cols of the matrix.
  // You should only change this if you know what you are doing and have also
  // changed the "values" array in a consistent manner.

  /*
  uint64_t& get_n_rows();
  uint64_t& get_n_cols();
  */

  // Const version of the above, just returns a copy of the number of
  // rows/cols, does not allow anything to be modified.
  inline uint64_t get_n_rows() const { return n_rows; }
  inline uint64_t get_n_cols() const { return n_cols; }


  //
  // Data
  //
  
  // Users have direct access to the array of values
  std::vector<double> values;
  
 private:
  uint64_t n_rows, n_cols;
};

#ifdef _DEBUG
template <class T> void writeD(std::string filename, const std::vector<T> &v){
  write_matrix(filename.c_str(), v.size(), 1, &v[0], 1);
}
template <class T> void writeD(std::string filename, const T *v, unsigned length){
  write_matrix(filename.c_str(), length, 1, &v[0], 1);
}
#else
template <class T> void writeD(std::string filename, const std::vector<T> &v){
}
template <class T> void writeD(std::string filename, const T *v, unsigned length){
}
#endif
template <class T> void write(std::string filename, const std::vector<T> &v){
  write_matrix(filename.c_str(), v.size(), 1, &v[0], 1);
}

FortranMatrix matmat(const FortranMatrix& A, const FortranMatrix& B,
		     bool transA, bool transB);

std::vector<double> matvec(const FortranMatrix& A, 
			   const std::vector<double>& x,
			   bool transA);
#endif

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
  inline double& operator()(uint64_t row, uint64_t col)
  {
    // storage is column-major
    return values[row + n_rows*col];
  }

  void print(std::string title) const;
  int write(std::string filename);
  int writeD(std::string filename);
  bool checkNeg();

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

template <class T> bool checkNeg(T* data, uint64_t length){
  for(uint64_t i = 0; i < length; i++)
    if(data[i] < 0)
      return true;
  return false;
}
#else
template <class T> void writeD(std::string filename, const std::vector<T> &v){}
template <class T> void writeD(std::string filename, const T *v, unsigned length){}
template <class T> bool checkNeg(T* data, uint64_t length){}
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

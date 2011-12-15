#ifndef __print_matrix_h__
#define __print_matrix_h__
#include <fstream>
#include <iostream>
#include <iomanip>

void print_matrix( const char* desc, int m, int n, const double* a, int lda);
template <class T> void write_matrix(const char *filename, int m, int n, const T *a, int lda){
  std::fstream file;
  file.open(filename, std::fstream::out);
  if(file.fail()){
    std::cout << "failed to open file: " << filename << std::endl;
    exit(1);
  }
  file << std::setprecision(15);
  //file << std::scientific << std::setprecision(5);
  //  file << std::scientific << std::setprecision(15);
  
  for( int i = 0; i < m; i++ ) {
    for( int j = 0; j < n; j++ )
      file /*<< std::setw(13)*/ << " " << a[i+j*lda];
    file << std::endl;
  }
  file.close();
}

#endif

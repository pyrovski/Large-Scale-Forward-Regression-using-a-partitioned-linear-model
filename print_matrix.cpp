
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

// Print function: actually uses C++ for printing even though it has a C interface.
void print_matrix( const char* desc, int m, int n, const double* a, int lda)
{
  std::ios::fmtflags flags = std::cout.flags();
  std::cout << std::setprecision(2);
  //std::cout << std::scientific << std::setprecision(5);
  //  std::cout << std::scientific << std::setprecision(15);
  
  std::cout << "\n" << desc << "\n";
  for( int i = 0; i < m; i++ ) {
    for( int j = 0; j < n; j++ )
      std::cout /*<< std::setw(13)*/ << " " << a[i+j*lda];
    std::cout << std::endl;
  }
  std::cout.flags(flags);
}

void write_matrix(const char *filename, int m, int n, const double *a, int lda){
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


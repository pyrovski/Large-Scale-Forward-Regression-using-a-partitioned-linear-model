#ifndef __print_matrix_h__
#define __print_matrix_h__

#include <iostream>
#include <iomanip>

// Print function: actually uses C++ for printing even though it has a C interface.
void print_matrix( const char* desc, int m, int n, const double* a, int lda )
{
  std::ios::fmtflags flags = std::cout.flags();
  //std::cout << std::fixed << std::setprecision(2);
  // std::cout << std::scientific << std::setprecision(5);
  std::cout << std::scientific << std::setprecision(15);
  
  std::cout << "\n" << desc << "\n";
  for( int i = 0; i < m; i++ ) {
    for( int j = 0; j < n; j++ )
      std::cout /*<< std::setw(13)*/ << " " << a[i+j*lda];
    std::cout << std::endl;
  }
  std::cout.flags(flags);
}

#endif

#ifndef __inner_h__
#define __inner_h__

#include <vector>

// Compute the (real-valued) inner product of 2 vectors
double inner(const std::vector<double>& a, const std::vector<double>& b)
{
  if (a.size() != b.size())
    {
      std::cerr << "Error! Input vectors must have the same length!" << std::endl;
      exit(1);
    }
  
  double result = 0.;
  for (unsigned i=0; i<a.size(); ++i)
    result += a[i]*b[i];
  
  return result;
}

#endif

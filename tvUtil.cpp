#include "tvUtil.h"

// subtract timeval b from a
struct timeval operator-(const struct timeval &lhs, const struct timeval &rhs){
  struct timeval dest;
  dest.tv_sec = lhs.tv_sec - rhs.tv_sec - (rhs.tv_usec > lhs.tv_usec);
  dest.tv_usec = rhs.tv_usec > lhs.tv_usec ? 
    1000000 + lhs.tv_usec - rhs.tv_usec: 
    lhs.tv_usec - rhs.tv_usec;
  return dest;
}

double tvDouble(const struct timeval &tv){
  double result = tv.tv_sec + tv.tv_usec / 1000000.0;
  return result;
}

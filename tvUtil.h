#ifndef TVUTIL_H
#define TVUTIL_H
#include <sys/time.h>

double tvDouble(const struct timeval &ts);
struct timeval operator-(const struct timeval &lhs, const struct timeval &rhs);
#endif

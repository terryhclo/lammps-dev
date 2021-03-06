#include <iso646.h>
#include "erf.h"
#include "direct.h"
#include "math.h"
// LAMMPS uses usleep with 100 ms arguments, no microsecond precision needed
#include "sleep.h"

// some symbols have different names in Windows

#undef ATOBIGINT
#define ATOBIGINT _atoi64

#define pclose _pclose
#define __restrict__ __restrict

// the following functions ared defined to get rid of 
// 'ambiguous call to overloaded function' error in VSS for mismathched type arguments

inline double pow(int i, int j){
  return pow((double)i,(double) j);
}

inline double sqrt(int i){
  return sqrt((double) i);
}

inline double fabs(int i){
  return fabs((double) i);
}

inline double trunc(double x) {
  return x > 0 ? floor(x) : ceil(x);
}

// Windows version of mkdir function does not have permission flags
# define S_IRWXU 0
# define S_IRGRP 0
# define S_IXGRP 0
inline int mkdir(const char *path, int flags){
  return _mkdir(path);
}

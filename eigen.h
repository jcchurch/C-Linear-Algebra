#ifndef __ISOMAP_EIGEN
#define __ISOMAP_EIGEN
#include "matrix.h"

double powerMethod(matrix* a);
matrix* francisQRstep(matrix* a);
matrix* eigenvector(matrix* a, double eigenvalue);

#endif

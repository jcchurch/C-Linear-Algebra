#ifndef __ISOMAP_MATRIXADV
#define __ISOMAP_MATRIXADV
#include "matrix.h"

void LUdecomposition(matrix* a, matrix** l, matrix** u);
double determinantMatrix(matrix* a);
matrix* matrixInverse(matrix* a);
matrix* solver(matrix* a, matrix* b);

#endif

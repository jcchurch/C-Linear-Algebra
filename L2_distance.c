#include "matrix.h"
#include <stdlib.h>

/*===========================================================================
 * l2_disatnce
 * Given a two matrices with identical widths (i.e. same dimensionality),
 * this method returns an a->height by b->height matrix of the l2-norm
 * distance measures (without the square root calculation)
 *
 * if a and b are vectors, the l2-norm is given by:
 *  ||a - b|| = sqrt( ||a||^2 + ||b||^2 - 2 * a * b )
 *=========================================================================*/
matrix* L2_distance(matrix* a, matrix* b) {
    matrix* out;
    matrix* aa;
    matrix* bb;
    matrix* bt;
    matrix* ab;
    int i, j;
    double* ptrAB;
    double* ptrOut;

    assert(a->width == b->width, "Matrices a and b must be the same dimensionality.");

    out = makeMatrix(b->height, a->height);

    aa = dotDiagonalMatrix(a, NULL);
    bb = dotDiagonalMatrix(b, NULL);
    bt = transposeMatrix(b);
    ab = multiplyMatrix(a, bt);

    ptrOut = out->data;
    ptrAB = ab->data;
    for (i = 0; i < a->height; i++) {   
        for (j = 0; j < b->height; j++) {
            *ptrOut = aa->data[i] + bb->data[j] - 2 * *ptrAB;
            ptrOut++;
            ptrAB++;
        }
    }

    freeMatrix(aa);
    freeMatrix(bb);
    freeMatrix(bt);
    freeMatrix(ab);
    return out;
}

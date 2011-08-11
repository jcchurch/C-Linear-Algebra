#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "eigen.h"
#include "svd.h"

void svd(matrix* a) {
    matrix* aat;
    matrix* ata;
    matrix* at;
    matrix* eigen_aat;
    matrix* eigen_ata;

    assert(a->height == a->width, "Matrix must be square.");

    at = transposeMatrix(a);
    ata = multiplyMatrix(at, a);
    aat = multiplyMatrix(a, at);

    eigen_aat = francisQRstep(aat);
    eigen_ata = francisQRstep(ata);

    printMatrix(eigen_aat);
    printMatrix(eigen_ata);

    freeMatrix(eigen_aat);
    freeMatrix(eigen_ata);
    freeMatrix(at);
    freeMatrix(ata);
    freeMatrix(aat);
}

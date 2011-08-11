#include <stdio.h>
#include "matrix.h"
#include "matrixadv.h"

int main(int argc, char** argv) {
    matrix* l = NULL;
    matrix* u = NULL;
    matrix* x;
    matrix* i;
    matrix* p;
    matrix* q;
    matrix* a = readMatrix(argv[1]);

    printf("Original -----------------\n");
    printMatrix(a);

    i = eyeMatrix(a->width);
    printf("Eye ----------------------\n");
    printMatrix(i);
    LUdecomposition(a, &l, &u);
    p = multiplyMatrix(l, u);
    x = solver(a, i);
    q = multiplyMatrix(a, x);

    printf("Lower --------------------\n");
    printMatrix(l);
    printf("Upper --------------------\n");
    printMatrix(u);
    printf("Product ------------------\n");
    printMatrix(p);
    printf("Solver with Identity -----\n");
    printMatrix(x);
    printf("Identity -----------------\n");
    printMatrix(q);

    freeMatrix(i);
    freeMatrix(q);
    freeMatrix(x);
    freeMatrix(a);
    freeMatrix(p);
    freeMatrix(l);
    freeMatrix(u);
    return 0;
}

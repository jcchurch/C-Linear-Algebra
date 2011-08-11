#include <stdio.h>
#include "matrix.h"
#include "qr.h"

int main(int argc, char** argv) {
    matrix* b;
    matrix* c;
    matrix* a = readMatrix(argv[1]);
    b = unitVectorRows(a);
    c = unitVectorColumns(a);

    printf("Original -----------------\n");
    printMatrix(a);

    printf("Row Norm -----------------\n");
    printMatrix(b);
    printf("Col Norm -----------------\n");
    printMatrix(c);

    freeMatrix(a);
    freeMatrix(b);
    freeMatrix(c);
    return 0;
}

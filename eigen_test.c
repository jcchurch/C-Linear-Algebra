#include <stdio.h>
#include "matrix.h"
#include "eigen.h"

int main(int argc, char** argv) {
    matrix* a = readMatrix(argv[1]);
    printf("Read.\n");
    matrix* e = francisQRstep(a);
    matrix* evec;
    double eigenvalue;

    printMatrix(e);

    eigenvalue = e->data[ e->height - 1 ];

    printf("Largest eigenvalue: %f (francis qr step)\n", eigenvalue);

    eigenvalue = powerMethod(a);

    printf("Largest eigenvalue: %f (power method) \n", eigenvalue);

    evec = eigenvector(a, eigenvalue);

    printMatrix(evec);

    freeMatrix(evec);
    freeMatrix(e);
    freeMatrix(a);
    return 0;
}

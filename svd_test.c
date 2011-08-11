#include <stdio.h>
#include "matrix.h"
#include "svd.h"

int main(int argc, char** argv) {
    matrix* a = readMatrix(argv[1]);

    svd(a);

    freeMatrix(a);
    return 0;
}

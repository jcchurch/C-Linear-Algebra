#include <stdio.h>
#include "qsort.h"

int main(int argc, char** argv) {
    int i;
    double a[] = {7,6,5,4,3,2,1};

    quicksort(a, 0, 6);

    for (i = 0; i < 7; i++)
        printf("%f\n", a[i]);
    return 0;
}

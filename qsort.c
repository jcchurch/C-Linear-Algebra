#include "qsort.h"

/*===========================================================================
 * quicksort
 * This algorithm performs the quicksort algorithm given an array of
 * doubles, a start index, and an end index.
 *=========================================================================*/
void quicksort(double *a, int start, int end) {
    int pivot;

    if (start < end) {
        pivot = partition (a, start, end);
        quicksort (a, start, pivot-1);
        quicksort (a, pivot+1, end);
    }
} 

/*===========================================================================
 * partition
 * This algorithm partitions an array so that all of the low values are
 * near the beginning and high values are near the end.
 * The partitioning value is the first element.
 *=========================================================================*/
int partition(double *a, int start, int end) {
    double val = a[start];
    int pivot = start;
    double temp;
    int i;

    for (i = start+1; i <= end; i++) {
        if (a[i] < val) {
            pivot++;
            temp = a[pivot];
            a[pivot] = a[i];
            a[i] = temp;
        }
    }

    temp = a[start];
    a[start] = a[pivot];
    a[pivot] = temp;
    return pivot;
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "matrixadv.h"
#include "qr.h"
#include "qsort.h"

/*===========================================================================
 * francisQRstep
 * This algorithm performs the Francis QR Step to find the eigen values
 * of a square matrix.
 *
 * I don't know what this algorithm should produce and I just read this
 * technique on wikipedia. Currently I have this method to test the
 * approach.
 *=========================================================================*/
matrix* francisQRstep(matrix* a) {
    int i, j, k = 0;
    int upper_triangle;
    const double EPSILON = 0.001;
    matrix* a_k;   // Matrix A_k
    matrix* a_kp1; // Matrix A_(k+1)
    matrix* q = NULL;
    matrix* r = NULL;
    matrix* out;
    double* ptr;
    double* outPtr;

    assert(a->width == a->height, "Matrix must be square.");

    out = makeMatrix(1, a->height);

    a_k = copyMatrix(a);

    while (1) {
        // Perform the QR decomposition of a_k
        printf("Loop.\n");
        gram_schmidt(a_k, &q, &r);

        // a_kp1 is the product of r * q
        a_kp1 = multiplyMatrix(r, q);

        // And the cleanup
        freeMatrix(a_k);
        freeMatrix(q);
        freeMatrix(r);

        a_k = copyMatrix(a_kp1);

        freeMatrix(a_kp1);

        q = NULL;
        r = NULL;
        printMatrix(a_k);

        // Test for completion.
        // We need to make sure this is an upper triangular matrix
        upper_triangle = 1; // TRUE
        for (i = 1; i < a_k->height; i++) {
            for (j = 0; j < i; j++) {
                if (fabs(a_k->data[i * a_k->width + j]) >= EPSILON) {
                    upper_triangle = 0; // FALSE
                    break;
                }
            }
            if (upper_triangle == 0) {
                break;
            }
        }

        if (upper_triangle == 1) {
            break;
        }

        k++;
    }

    // Gather up all of the elements along the diagonal.
    ptr = a_k->data;
    outPtr = out->data;
    for (i = 0; i < a_k->width; i++) {
        *outPtr = *ptr;
        outPtr++;
        ptr += a_k->width + 1;
    }

    // Sort the eigen values
    quicksort(out->data, 0, out->width * out->height - 1);

    freeMatrix(a_k);
    return out;
}

/*===========================================================================
 * powerMethod
 * This algorithm determines the largest eigenvalue of a matrix using the
 * power method.
 *
 * This was described to me in a Randomized Algoirthms course.
 *=========================================================================*/
double powerMethod(matrix* a) {
    matrix* x;
    matrix* xp1; // x plus 1
    const double EPSILON = 0.001;
    double sum;
    int i;
    int k = 0;
    int converge;

    assert(a->width == a->height, "Matrix must be square.");

    srand(time(0)); // Initalize our RNG

    // Let's initalize x to a random vector
    x = makeMatrix(1, a->height);
    for (i = 0; i < a->height; i++) {
        x->data[i] = (double) rand() / RAND_MAX;
    }

    // Iterate until the x vector converges.
    while (1) {
        k++;

        // Multiply A * x
        xp1 = multiplyMatrix(a, x);

        // Add up all of the values in xp1
        sum = 0;
        for (i = 0; i < a->height; i++) {
            sum += xp1->data[i];
        }

        // Divide each value in xp1 by sum. (Normalize)
        for (i = 0; i < a->height; i++) {
            xp1->data[i] /= sum;
        }

        // Test to see if we need to quit.
        converge = 1; // Converged
        for (i = 0; i < a->height; i++) {
            if (fabs(x->data[i] - xp1->data[i]) >= EPSILON) {
                converge = 0; // Not converged.
                break;
            }
        }

        // Set up for the next loop.
        freeMatrix(x);
        x = copyMatrix(xp1);
        freeMatrix(xp1);

        // Really test for quit.
        if (converge == 1) {
            break;
        }
    }

    freeMatrix(x);
    return sum;
}

/*===========================================================================
 * eigenvector
 * This algorithm determines the eigenvector of a matrix given an eigenvalue.
 *=========================================================================*/
matrix* eigenvector(matrix* a, double eigenvalue) {
    matrix* b; // This matrix will store A-eI
    matrix* zero; // This matrix will store a column vector of zeros
    matrix* out;
    double* ptr;
    int i;

    assert(a->width == a->height, "Matrix must be square.");

    // Create our column vector of zeros
    zero = makeMatrix(1, a->height);

    // Copy A
    b = copyMatrix(a);

    // Subtract eigenvalue from the diagonal elements
    ptr = b->data;
    for (i = 0; i < b->height; i++) {
        *ptr -= eigenvalue;
        ptr += b->width + 1;
    }

    // Find the eigenvector
    out = solver(b, zero);

    freeMatrix(b);
    freeMatrix(zero);
    return out;
}

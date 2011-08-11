#include "matrix.h"
#include <stdio.h>

/*===========================================================================
 * LUdecomposition
 * Given a matrix A and two more null values, the function uses Crout's
 * algorithm (page 44 of Numerical Recipes in C) to perform the
 * LU decomposition of matrix A.
 *=========================================================================*/
void LUdecomposition(matrix* a, matrix** l, matrix** u) {
    int i, j, k;
    double* ptrA;
    double* ptrL;
    double* ptrU;
    double sum;

    assert(a->width == a->height, "Matrix A must be square");

    assert(*l == NULL && *u == NULL, "Matricies L and U must be null");

    *l = makeMatrix(a->width, a->height);
    *u = makeMatrix(a->width, a->height);

    // Step 1: Assign 1 to the diagonal of the lower matrix.
    ptrL = (*l)->data;
    for (i = 0; i < a->width; i++) {
        *ptrL = 1.0;
        ptrL += a->width + 1;
    }

    // Step 2
    for (j = 0; j < a->width; j++) {

        // Part A: Solve for the upper matrix.
        for (i = 0; i <= j; i++) {

            sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += (*l)->data[i * a->width + k] * (*u)->data[k * a->width + j];
            }

            (*u)->data[i * a->width + j] = a->data[i * a->width + j] - sum;
        }

        // Part B: Solve fpr the lower matrix
        for (i = j+1; i < a->width; i++) {

            sum = 0.0;
            for (k = 0; k < j; k++) {
                sum += (*l)->data[i * a->width + k] * (*u)->data[k * a->width + j];
            }

            (*l)->data[i * a->width + j] = 1.0/(*u)->data[j * a->width + j] * ( a->data[i * a->width + j] - sum);
        }
    }

    return;
}

/*===========================================================================
 * determinantMatrix
 * Given a matrix A, returns the determinant.
 * Based on the formula on page 49 of Numerical Recipes in C.
 *=========================================================================*/
double determinantMatrix(matrix* a) {
    double product = 0.0;
    matrix* l = NULL;
    matrix* u = NULL;
    int i;

    assert(a->width == a->height, "Matrix A must be square.");
    LUdecomposition(a, &l, &u);

    // Get the product of upper matrix diagonal
    // We don't need the lower matrix for this calculation.
    for (i = 0; i < a->width; i++) {
        product *= u->data[i * a->width + i];
    }

    freeMatrix(l);
    freeMatrix(u);
    return product;
}

/*===========================================================================
 * solver
 * Solves a system of liner equations using LU decomposition. This requires
 * an A: n by n matrix and B: n by p matrix, where
 * A * X = b
 * The algorithm returns X, which will be a n by p matrix.
 *
 * This algorithm is described on page 121 of the book "Matrix Computations"
 * by Golub and Loan.
 *=========================================================================*/
matrix* solver(matrix* a, matrix* b) {
    int i, j, k;
    double sum;
    matrix *l = NULL;
    matrix *u = NULL;
    matrix *y;
    matrix *x;
    double* row;

    assert(a->height == b->height, "In the solver, both matrices should have the same height.");

    LUdecomposition(a, &l, &u);
    y = makeMatrix(1, a->height);
    x = makeMatrix(b->width, b->height);

    for (k = 0; k < b->width; k++) {
        // Perform backward subsitituion with L
        // L * y = B_k
        for (i = 0; i < a->height; i++) {
            row = l->data + i * a->width;
            sum = 0;
            for (j = 0; j < i; j++) {
                sum += y->data[j] * (*row++);
            }

            y->data[i] = (b->data[i * b->width + k] - sum) / *row;
        }

        // Perform backward subsitituion again with U
        // U * x = y
        for (i = a->height - 1; i >= 0; i--) {
            row = u->data + i * a->width + (a->width - 1);
            sum = 0;
            for (j = a->width - 1; j > i; j--) {
                sum += x->data[j * b->width + k] * (*row--);
            }

            x->data[i * b->width + k] = (y->data[i] - sum) / *row;
        }
    }

    freeMatrix(l);
    freeMatrix(u);
    freeMatrix(y);
    return x;
}

/*===========================================================================
 * inverseMatrix
 * Given a matrix A, returns the determinant.
 *
 * This algorithm is described on page 121 of the book "Matrix Computations"
 * by Golub and Loan.
 *=========================================================================*/
matrix* inverseMatrix(matrix* a) {
    matrix* eye;
    matrix* inv;

    assert(a->width == a->height, "Matrix A must be square.");
    eye = eyeMatrix(a->width);

    inv = solver(a, eye);
    freeMatrix(eye);
    return inv;
}

/*===========================================================================
 * row_echelon_form
 * Given a matrix A, returns the same matrix in row-echelon form
 *=========================================================================*/
matrix* row_echelon_form(matrix* a) {
    matrix* out = copyMatrix(a);
    double* outPtr = out->data;
    double* prevPtr;
    int i, j;

    // First, divide each element in the first row by the first element.
    double scale = *outPtr;
    for (i = 0; i < out->width; i++) {
        *outPtr++ /= scale;
    }

    // Next, for each row beyond the first, we scale the previous
    // row by the first non-zero value in the current row, then
    // subtract that result from the current row. That should make
    // The first non-zero value now zero. Then we divide each element
    // on the current row by the next non-zero value.
    for (j = 1; j < out->height; j++) {
        outPtr = out->data + (j * out->width) + j;    
        prevPtr = out->data + (j-1 * out->width) + j;    

        scale = *outPtr;
        for (i = j; i < out->width; i++) {
            *outPtr = *outPtr - (*prevPtr * scale);
        }

       
    }
}

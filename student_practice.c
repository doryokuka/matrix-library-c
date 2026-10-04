/*******************************************************************************
* File Name   : student_practice.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : 第七节课课内练习程序。
*******************************************************************************/

#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_solve.h"
#include "matrix_lu.h"

#include <math.h>
#include <stdio.h>

static int AlmostEqual(REAL a, REAL b)
{
    return fabs(a - b) <= 1e-12;
}

int main(void)
{
    Matrix A, b, x_gauss, L, U, x_lu;
    MatrixInit(&A); MatrixInit(&b); MatrixInit(&x_gauss);
    MatrixInit(&L); MatrixInit(&U); MatrixInit(&x_lu);

    MatrixCreate(&A, 2, 2);
    MatrixCreate(&b, 2, 1);
    MatrixCreate(&x_gauss, 2, 1);
    MatrixCreate(&L, 2, 2);
    MatrixCreate(&U, 2, 2);
    MatrixCreate(&x_lu, 2, 1);

    /* 同一个线性系统分别用高斯消元和 LU 分解求解。精确解为 x=[1,2]^T。 */
    MatrixSet(&A, 0, 0, 4.0); MatrixSet(&A, 0, 1, 3.0);
    MatrixSet(&A, 1, 0, 6.0); MatrixSet(&A, 1, 1, 3.0);
    MatrixSet(&b, 0, 0, 10.0); MatrixSet(&b, 1, 0, 12.0);

    printf("Practice: Gaussian elimination and LU solve\n\n");
    MatrixPrint(&A, "A");
    MatrixPrint(&b, "b");

    MatrixError eg = GaussianSolvePartialPivot(&A, &b, &x_gauss, 1e-12);
    if (eg == MATRIX_SUCCESS) {
        MatrixPrint(&x_gauss, "x from Gaussian elimination");
    } else {
        printf("Gaussian solve failed: %s\n", MatrixErrorMessage(eg));
    }

    MatrixError elu = LUDecomposeNoPivot(&A, &L, &U, 1e-12);
    if (elu == MATRIX_SUCCESS) {
        elu = LUSolve(&L, &U, &b, &x_lu, 1e-12);
    }
    if (elu == MATRIX_SUCCESS) {
        MatrixPrint(&L, "L");
        MatrixPrint(&U, "U");
        MatrixPrint(&x_lu, "x from LU solve");
    } else {
        printf("LU solve failed: %s\n", MatrixErrorMessage(elu));
    }

    if (eg == MATRIX_SUCCESS && elu == MATRIX_SUCCESS) {
        int ok = AlmostEqual(x_gauss.data[0], x_lu.data[0]) &&
                 AlmostEqual(x_gauss.data[1], x_lu.data[1]);
        printf("Compare two solutions: %s\n", ok ? "PASS" : "FAIL");
    }

    MatrixFree(&A); MatrixFree(&b); MatrixFree(&x_gauss);
    MatrixFree(&L); MatrixFree(&U); MatrixFree(&x_lu);
    return 0;
}

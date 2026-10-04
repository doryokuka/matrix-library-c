/*******************************************************************************
* File Name   : test_runner.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : 单元测试、错误路径测试和线性求解验证程序。
*******************************************************************************/

#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_solve.h"
#include "matrix_lu.h"
#include "matrix_diag.h"

#include <math.h>
#include <stdio.h>

static int total_tests = 0;
static int passed_tests = 0;

#define CHECK_TRUE(condition, name) \
    do { \
        total_tests++; \
        if (condition) { \
            passed_tests++; \
            printf("[PASS] %s\n", name); \
        } else { \
            printf("[FAIL] %s\n", name); \
        } \
    } while (0)
// equal to: #define CHECK_TRUE(condition, name) do { total_tests++; if (condition) { passed_tests++; printf("[PASS] %s\n", name); } else { printf("[FAIL] %s\n", name); } } while (0)
static int AlmostEqual(REAL a, REAL b, REAL tol)
{
    return fabs(a - b) <= tol;
}

static int CreateMatrixForTest(Matrix *A, int row, int column, const char *name)
{
    MatrixError error = MatrixCreate(A, row, column);
    if (error != MATRIX_SUCCESS) {
        printf("Create %s failed: %s\n", name, MatrixErrorMessage(error));
        return 0;
    }
    return 1;
}

static void TestMatrixAdd(void)
{
    Matrix A, B, C;
    MatrixInit(&A); MatrixInit(&B); MatrixInit(&C);

    if (!CreateMatrixForTest(&A, 2, 2, "A") ||
        !CreateMatrixForTest(&B, 2, 2, "B") ||
        !CreateMatrixForTest(&C, 2, 2, "C")) {
        MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
        return;
    }

    /* 手算样例：A=[1,2;3,4]，B=[10,20;30,40]，A+B=[11,22;33,44]。 */
    MatrixFillSequence(&A, 1.0, 1.0);
    MatrixFillSequence(&B, 10.0, 10.0);

    CHECK_TRUE(MatrixAdd(&A, &B, &C) == MATRIX_SUCCESS,
               "MatrixAdd returns MATRIX_SUCCESS");
    CHECK_TRUE(AlmostEqual(C.data[0], 11.0, 1e-12),
               "MatrixAdd C(0,0) = 11");
    CHECK_TRUE(AlmostEqual(C.data[3], 44.0, 1e-12),
               "MatrixAdd C(1,1) = 44");

    MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
}

static void TestMatrixMultiply(void)
{
    Matrix A, B, C;
    MatrixInit(&A); MatrixInit(&B); MatrixInit(&C);

    if (!CreateMatrixForTest(&A, 2, 3, "A") ||
        !CreateMatrixForTest(&B, 3, 2, "B") ||
        !CreateMatrixForTest(&C, 2, 2, "C")) {
        MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
        return;
    }

    /* A=[1,2,3;4,5,6]，B=[1,2;3,4;5,6]，AB=[22,28;49,64]。 */
    MatrixFillSequence(&A, 1.0, 1.0);
    MatrixFillSequence(&B, 1.0, 1.0);

    CHECK_TRUE(MatrixMultiply(&A, &B, &C) == MATRIX_SUCCESS,
               "MatrixMultiply returns MATRIX_SUCCESS");
    CHECK_TRUE(AlmostEqual(C.data[0], 22.0, 1e-12),
               "MatrixMultiply C(0,0) = 22");
    CHECK_TRUE(AlmostEqual(C.data[3], 64.0, 1e-12),
               "MatrixMultiply C(1,1) = 64");

    MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
}

static void TestErrorPaths(void)
{
    Matrix A, B, C;
    MatrixInit(&A); MatrixInit(&B); MatrixInit(&C);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&B, 3, 3, "B");
    CreateMatrixForTest(&C, 2, 2, "C");

    CHECK_TRUE(MatrixAdd(&A, &B, &C) == MATRIX_ERROR_SIZE_MISMATCH,
               "MatrixAdd detects size mismatch");
    CHECK_TRUE(MatrixSet(&A, 5, 0, 1.0) == MATRIX_ERROR_INDEX_OUT_OF_RANGE,
               "MatrixSet detects index out of range");

    MatrixFree(&A); MatrixFree(&B); MatrixFree(&C);
}

static void TestGaussianNoPivot(void)
{
    Matrix A, b, x;
    MatrixInit(&A); MatrixInit(&b); MatrixInit(&x);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&b, 2, 1, "b");
    CreateMatrixForTest(&x, 2, 1, "x");

    /* 线性系统：2x1+x2=5，x1+3x2=7，精确解为 x=[1.6,1.8]^T。 */
    MatrixSet(&A, 0, 0, 2.0); MatrixSet(&A, 0, 1, 1.0);
    MatrixSet(&A, 1, 0, 1.0); MatrixSet(&A, 1, 1, 3.0);
    MatrixSet(&b, 0, 0, 5.0); MatrixSet(&b, 1, 0, 7.0);

    CHECK_TRUE(GaussianSolveNoPivot(&A, &b, &x, 1e-12) == MATRIX_SUCCESS,
               "GaussianSolveNoPivot returns MATRIX_SUCCESS");
    CHECK_TRUE(AlmostEqual(x.data[0], 1.6, 1e-12),
               "GaussianSolveNoPivot x(0) = 1.6");
    CHECK_TRUE(AlmostEqual(x.data[1], 1.8, 1e-12),
               "GaussianSolveNoPivot x(1) = 1.8");

    MatrixFree(&A); MatrixFree(&b); MatrixFree(&x);
}

static void TestGaussianPartialPivot(void)
{
    Matrix A, b, x;
    MatrixInit(&A); MatrixInit(&b); MatrixInit(&x);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&b, 2, 1, "b");
    CreateMatrixForTest(&x, 2, 1, "x");

    /* 第一个对角元为 0。无主元算法失败，列主元算法换行后得到 x=[1,2]^T。 */
    MatrixSet(&A, 0, 0, 0.0); MatrixSet(&A, 0, 1, 2.0);
    MatrixSet(&A, 1, 0, 1.0); MatrixSet(&A, 1, 1, 1.0);
    MatrixSet(&b, 0, 0, 4.0); MatrixSet(&b, 1, 0, 3.0);

    CHECK_TRUE(GaussianSolveNoPivot(&A, &b, &x, 1e-12) == MATRIX_ERROR_SINGULAR,
               "GaussianSolveNoPivot detects zero pivot");
    MatrixFillZero(&x);
    CHECK_TRUE(GaussianSolvePartialPivot(&A, &b, &x, 1e-12) == MATRIX_SUCCESS,
               "GaussianSolvePartialPivot returns MATRIX_SUCCESS");
    CHECK_TRUE(AlmostEqual(x.data[0], 1.0, 1e-12),
               "GaussianSolvePartialPivot x(0) = 1");
    CHECK_TRUE(AlmostEqual(x.data[1], 2.0, 1e-12),
               "GaussianSolvePartialPivot x(1) = 2");

    MatrixFree(&A); MatrixFree(&b); MatrixFree(&x);
}

static void TestGaussianSingular(void)
{
    Matrix A, b, x;
    MatrixInit(&A); MatrixInit(&b); MatrixInit(&x);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&b, 2, 1, "b");
    CreateMatrixForTest(&x, 2, 1, "x");

    MatrixSet(&A, 0, 0, 1.0); MatrixSet(&A, 0, 1, 2.0);
    MatrixSet(&A, 1, 0, 2.0); MatrixSet(&A, 1, 1, 4.0);
    MatrixSet(&b, 0, 0, 3.0); MatrixSet(&b, 1, 0, 6.0);

    CHECK_TRUE(GaussianSolvePartialPivot(&A, &b, &x, 1e-12) == MATRIX_ERROR_SINGULAR,
               "GaussianSolvePartialPivot detects singular matrix");

    MatrixFree(&A); MatrixFree(&b); MatrixFree(&x);
}

static void TestMatrixSwapRows(void)
{
    Matrix A;
    MatrixInit(&A);
    CreateMatrixForTest(&A, 2, 2, "A");
    MatrixFillSequence(&A, 1.0, 1.0);

    CHECK_TRUE(MatrixSwapRows(&A, 0, 1) == MATRIX_SUCCESS,
               "MatrixSwapRows returns MATRIX_SUCCESS");
    CHECK_TRUE(AlmostEqual(A.data[0], 3.0, 1e-12) && AlmostEqual(A.data[2], 1.0, 1e-12),
               "MatrixSwapRows exchanges two rows");
    CHECK_TRUE(MatrixSwapRows(&A, 0, 5) == MATRIX_ERROR_INDEX_OUT_OF_RANGE,
               "MatrixSwapRows detects index out of range");

    MatrixFree(&A);
}

static void TestLUDecompose(void)
{
    Matrix A, L, U, R;
    MatrixInit(&A); MatrixInit(&L); MatrixInit(&U); MatrixInit(&R);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&L, 2, 2, "L");
    CreateMatrixForTest(&U, 2, 2, "U");
    CreateMatrixForTest(&R, 2, 2, "R");

    MatrixSet(&A, 0, 0, 4.0); MatrixSet(&A, 0, 1, 3.0);
    MatrixSet(&A, 1, 0, 6.0); MatrixSet(&A, 1, 1, 3.0);

    CHECK_TRUE(LUDecomposeNoPivot(&A, &L, &U, 1e-12) == MATRIX_SUCCESS,
               "LUDecomposeNoPivot returns MATRIX_SUCCESS");
    CHECK_TRUE(MatrixMultiply(&L, &U, &R) == MATRIX_SUCCESS,
               "MatrixMultiply forms L * U");

    for (int i = 0; i < A.row; ++i) {
        for (int j = 0; j < A.column; ++j) {
            int id = MatrixIndex(&A, i, j);
            char name[80];
            snprintf(name, sizeof(name), "LU reconstruction A(%d,%d)", i, j);
            CHECK_TRUE(AlmostEqual(R.data[id], A.data[id], 1e-12), name);
        }
    }

    MatrixFree(&A); MatrixFree(&L); MatrixFree(&U); MatrixFree(&R);
}

static void TestLUSolve(void)
{
    Matrix A, L, U, b, x;
    MatrixInit(&A); MatrixInit(&L); MatrixInit(&U); MatrixInit(&b); MatrixInit(&x);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&L, 2, 2, "L");
    CreateMatrixForTest(&U, 2, 2, "U");
    CreateMatrixForTest(&b, 2, 1, "b");
    CreateMatrixForTest(&x, 2, 1, "x");

    MatrixSet(&A, 0, 0, 4.0); MatrixSet(&A, 0, 1, 3.0);
    MatrixSet(&A, 1, 0, 6.0); MatrixSet(&A, 1, 1, 3.0);
    MatrixSet(&b, 0, 0, 10.0); MatrixSet(&b, 1, 0, 12.0);

    CHECK_TRUE(LUDecomposeNoPivot(&A, &L, &U, 1e-12) == MATRIX_SUCCESS,
               "LU decomposition for solve example");
    CHECK_TRUE(LUSolve(&L, &U, &b, &x, 1e-12) == MATRIX_SUCCESS,
               "LUSolve returns MATRIX_SUCCESS");
    CHECK_TRUE(AlmostEqual(x.data[0], 1.0, 1e-12), "LUSolve x(0) = 1");
    CHECK_TRUE(AlmostEqual(x.data[1], 2.0, 1e-12), "LUSolve x(1) = 2");

    MatrixFree(&A); MatrixFree(&L); MatrixFree(&U); MatrixFree(&b); MatrixFree(&x);
}

static void TestSingularLU(void)
{
    Matrix A, L, U;
    MatrixInit(&A); MatrixInit(&L); MatrixInit(&U);

    CreateMatrixForTest(&A, 2, 2, "A");
    CreateMatrixForTest(&L, 2, 2, "L");
    CreateMatrixForTest(&U, 2, 2, "U");

    MatrixSet(&A, 0, 0, 1.0); MatrixSet(&A, 0, 1, 2.0);
    MatrixSet(&A, 1, 0, 2.0); MatrixSet(&A, 1, 1, 4.0);

    CHECK_TRUE(LUDecomposeNoPivot(&A, &L, &U, 1e-12) == MATRIX_ERROR_SINGULAR,
               "LU detects singular matrix");

    MatrixFree(&A); MatrixFree(&L); MatrixFree(&U);
}

static void TestMatrixDiagonalizeSymmetricJacobi(void)
{
    Matrix A, Q, D, QT, AQ, QTAQ;
    int i, j;
    MatrixInit(&A); MatrixInit(&Q); MatrixInit(&D);
    MatrixInit(&QT); MatrixInit(&AQ); MatrixInit(&QTAQ);

    if (!CreateMatrixForTest(&A, 3, 3, "A") ||
        !CreateMatrixForTest(&Q, 3, 3, "Q") ||
        !CreateMatrixForTest(&D, 3, 3, "D") ||
        !CreateMatrixForTest(&QT, 3, 3, "QT") ||
        !CreateMatrixForTest(&AQ, 3, 3, "AQ") ||
        !CreateMatrixForTest(&QTAQ, 3, 3, "QTAQ")) {
        MatrixFree(&A); MatrixFree(&Q); MatrixFree(&D);
        MatrixFree(&QT); MatrixFree(&AQ); MatrixFree(&QTAQ);
        return;
    }

    /*
     * A is a 3 x 3 real symmetric matrix.
     * The Jacobi method should produce Q and D such that Q^T * A * Q is nearly diagonal.
     */
    MatrixSet(&A, 0, 0, 4.0); MatrixSet(&A, 0, 1, 1.0); MatrixSet(&A, 0, 2, 1.0);
    MatrixSet(&A, 1, 0, 1.0); MatrixSet(&A, 1, 1, 3.0); MatrixSet(&A, 1, 2, 1.0);
    MatrixSet(&A, 2, 0, 1.0); MatrixSet(&A, 2, 1, 1.0); MatrixSet(&A, 2, 2, 2.0);

    CHECK_TRUE(MatrixIsSymmetric(&A, 1e-12) == MATRIX_SUCCESS,
               "MatrixIsSymmetric detects symmetric matrix");
    CHECK_TRUE(MatrixDiagonalizeSymmetricJacobi(&A, &Q, &D, 200, 1e-12) == MATRIX_SUCCESS,
               "MatrixDiagonalizeSymmetricJacobi returns MATRIX_SUCCESS");

    for (i = 0; i < D.row; i++) {
        for (j = 0; j < D.column; j++) {
            if (i != j) {
                char name[100];
                snprintf(name, sizeof(name), "Jacobi D(%d,%d) is nearly zero", i, j);
                CHECK_TRUE(AlmostEqual(D.data[MatrixIndex(&D, i, j)], 0.0, 1e-8), name);
            }
        }
    }

    CHECK_TRUE(MatrixTranspose(&Q, &QT) == MATRIX_SUCCESS,
               "Jacobi test computes Q^T");
    CHECK_TRUE(MatrixMultiply(&A, &Q, &AQ) == MATRIX_SUCCESS,
               "Jacobi test computes A * Q");
    CHECK_TRUE(MatrixMultiply(&QT, &AQ, &QTAQ) == MATRIX_SUCCESS,
               "Jacobi test computes Q^T * A * Q");

    for (i = 0; i < D.row; i++) {
        for (j = 0; j < D.column; j++) {
            char name[100];
            snprintf(name, sizeof(name), "Q^T A Q matches D at (%d,%d)", i, j);
            CHECK_TRUE(AlmostEqual(QTAQ.data[MatrixIndex(&QTAQ, i, j)],
                                   D.data[MatrixIndex(&D, i, j)],
                                   1e-8), name);
        }
    }

    MatrixFree(&A); MatrixFree(&Q); MatrixFree(&D);
    MatrixFree(&QT); MatrixFree(&AQ); MatrixFree(&QTAQ);
}

int main(void)
{
    printf("Lab 7: tests for the matrix solver library\n\n");

    TestMatrixAdd();
    TestMatrixMultiply();
    TestErrorPaths();
    TestMatrixSwapRows();
    TestGaussianNoPivot();
    TestGaussianPartialPivot();
    TestGaussianSingular();
    TestLUDecompose();
    TestLUSolve();
    TestSingularLU();
    TestMatrixDiagonalizeSymmetricJacobi();

    printf("\nPassed %d / %d tests.\n", passed_tests, total_tests);
    return (passed_tests == total_tests) ? 0 : 1;
}

/*******************************************************************************
* File Name   : benchmark.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : 第七节课性能测试程序。比较矩阵乘法循环顺序、高斯消元
*               和 LU 分解在不同矩阵规模下的运行时间。
*******************************************************************************/

#include "matrix_core.h"
#include "matrix_ops.h"
#include "matrix_solve.h"
#include "matrix_lu.h"

#include <stdio.h>
#include <time.h>

/* 将 clock() 的计数差转换为秒。 */
static double ElapsedSeconds(clock_t start, clock_t end)
{
    return (double)(end - start) / CLOCKS_PER_SEC;
}

/*
 * 对同一个矩阵乘法函数重复计时，返回平均时间。
 * multiply 可以指向 MatrixMultiply(ijk) 或 MatrixMultiplyIKJ(i-k-j)。
 */
static double AverageMultiplyTime(MatrixError (*multiply)(const Matrix *, const Matrix *, Matrix *),
                                  const Matrix *A,
                                  const Matrix *B,
                                  Matrix *C,
                                  int repeats)
{
    double total = 0.0;
    for (int r = 0; r < repeats; ++r) {
        clock_t start = clock();
        MatrixError error = multiply(A, B, C);
        clock_t end = clock();
        if (error != MATRIX_SUCCESS) {
            printf("Matrix multiplication failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
        total += ElapsedSeconds(start, end);
    }
    return total / (double)repeats;
}

/*
 * 生成对角占优矩阵。该矩阵在课堂测试中具有较稳定的主元，
 * 适合用于高斯消元和无主元 LU 分解的性能计时。
 */
static void FillWellConditionedMatrix(Matrix *A)
{
    MatrixFillZero(A);
    for (int i = 0; i < A->row; ++i) {
        for (int j = 0; j < A->column; ++j) {
            REAL value = (i == j) ? (REAL)(A->row + 2) : 1.0 / (REAL)(1 + i + j);
            A->data[MatrixIndex(A, i, j)] = value;
        }
    }
}

/* 生成右端项向量 b。 */
static void FillRightHandSide(Matrix *b)
{
    for (int i = 0; i < b->row; ++i) {
        b->data[i] = 1.0 + 0.01 * (REAL)i;
    }
}

/* 对列主元高斯消元重复计时，返回平均时间。 */
static double AverageGaussianTime(const Matrix *A, const Matrix *b, Matrix *x, int repeats)
{
    double total = 0.0;
    for (int r = 0; r < repeats; ++r) {
        clock_t start = clock();
        MatrixError error = GaussianSolvePartialPivot(A, b, x, 1e-12);
        clock_t end = clock();
        if (error != MATRIX_SUCCESS) {
            printf("Gaussian elimination failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
        total += ElapsedSeconds(start, end);
    }
    return total / (double)repeats;
}

/* 对无主元 LU 分解重复计时，返回平均时间。 */
static double AverageLUTime(const Matrix *A, Matrix *L, Matrix *U, int repeats)
{
    double total = 0.0;
    for (int r = 0; r < repeats; ++r) {
        clock_t start = clock();
        MatrixError error = LUDecomposeNoPivot(A, L, U, 1e-12);
        clock_t end = clock();
        if (error != MATRIX_SUCCESS) {
            printf("LU decomposition failed: %s\n", MatrixErrorMessage(error));
            return -1.0;
        }
        total += ElapsedSeconds(start, end);
    }
    return total / (double)repeats;
}

/*
 * 性能测试一：矩阵乘法。
 * 比较 ijk 三重循环与 i-k-j 三重循环的平均运行时间。
 */
static void BenchmarkMatrixMultiplication(void)
{
    int sizes[] = {100, 200, 400};
    int count = (int)(sizeof(sizes) / sizeof(sizes[0]));
    int repeats = 3;

    printf("[Benchmark 1] Matrix multiplication\n");
    printf("repeats = %d\n", repeats);
    for (int s = 0; s < count; ++s) {
        int n = sizes[s];
        Matrix A, B, C;
        MatrixInit(&A);
        MatrixInit(&B);
        MatrixInit(&C);

        MatrixCreate(&A, n, n);
        MatrixCreate(&B, n, n);
        MatrixCreate(&C, n, n);
        MatrixFillSequence(&A, 1.0, 0.001);
        MatrixFillSequence(&B, 2.0, 0.001);

        double tij = AverageMultiplyTime(MatrixMultiply, &A, &B, &C, repeats);
        double tikj = AverageMultiplyTime(MatrixMultiplyIKJ, &A, &B, &C, repeats);
        printf("n = %4d | ijk = %.6f s | ikj = %.6f s\n", n, tij, tikj);

        MatrixFree(&A);
        MatrixFree(&B);
        MatrixFree(&C);
    }
}

/*
 * 性能测试二：线性方程组相关算法。
 * 比较列主元高斯消元求解与无主元 LU 分解的平均运行时间。
 */
static void BenchmarkLinearSolvers(void)
{
    int sizes[] = {50, 100, 200};
    int count = (int)(sizeof(sizes) / sizeof(sizes[0]));
    int repeats = 3;

    printf("\n[Benchmark 2] Linear solver kernels\n");
    printf("repeats = %d\n", repeats);
    for (int s = 0; s < count; ++s) {
        int n = sizes[s];
        Matrix A, b, x, L, U;
        MatrixInit(&A);
        MatrixInit(&b);
        MatrixInit(&x);
        MatrixInit(&L);
        MatrixInit(&U);

        MatrixCreate(&A, n, n);
        MatrixCreate(&b, n, 1);
        MatrixCreate(&x, n, 1);
        MatrixCreate(&L, n, n);
        MatrixCreate(&U, n, n);
        FillWellConditionedMatrix(&A);
        FillRightHandSide(&b);

        double tg = AverageGaussianTime(&A, &b, &x, repeats);
        double tlu = AverageLUTime(&A, &L, &U, repeats);
        printf("n = %4d | Gaussian(partial pivot) = %.6f s | LU = %.6f s\n", n, tg, tlu);

        MatrixFree(&A);
        MatrixFree(&b);
        MatrixFree(&x);
        MatrixFree(&L);
        MatrixFree(&U);
    }
}

int main(void)
{
    BenchmarkMatrixMultiplication();
    BenchmarkLinearSolvers();
    return 0;
}

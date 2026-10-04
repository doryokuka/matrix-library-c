/*******************************************************************************
* File Name   : matrix_solve.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : Gaussian elimination 线性方程组求解函数实现。
*******************************************************************************/

#include "matrix_solve.h"
#include <math.h>

MatrixError MatrixSwapRows(Matrix *A, int r1, int r2)
{
    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (r1 < 0 || r1 >= A->row || r2 < 0 || r2 >= A->row) {
        return MATRIX_ERROR_INDEX_OUT_OF_RANGE;
    }
    if (r1 == r2) {
        return MATRIX_SUCCESS;
    }

    /* 对矩阵的两行逐列交换。右端项 b 也可以看成 n x 1 的矩阵使用该函数。 */
    for (int j = 0; j < A->column; ++j) {
        int idx1 = MatrixIndex(A, r1, j);
        int idx2 = MatrixIndex(A, r2, j);
        REAL tmp = A->data[idx1];
        A->data[idx1] = A->data[idx2];
        A->data[idx2] = tmp;
    }
    return MATRIX_SUCCESS;
}

static MatrixError CheckLinearSystem(const Matrix *A, const Matrix *b, Matrix *x)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(b) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }
    if (b->row != A->row || b->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    if (x->row != A->column || x->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

MatrixError GaussianSolveNoPivot(const Matrix *A, const Matrix *b, Matrix *x, REAL tol)
{
    MatrixError error = CheckLinearSystem(A, b, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;
    Matrix U, rhs;
    MatrixInit(&U);
    MatrixInit(&rhs);

    error = MatrixCreate(&U, n, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    error = MatrixCreate(&rhs, n, 1);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        return error;
    }
    MatrixCopy(A, &U);
    MatrixCopy(b, &rhs);

    /* 消元阶段：第 k 步使用 U(k,k) 消去其下方元素。 */
    for (int k = 0; k < n - 1; ++k) {
        REAL pivot = U.data[MatrixIndex(&U, k, k)];
        if (fabs(pivot) < tol) {
            MatrixFree(&U);
            MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }
        for (int i = k + 1; i < n; ++i) {
            REAL factor = U.data[MatrixIndex(&U, i, k)] / pivot;
            U.data[MatrixIndex(&U, i, k)] = 0.0;
            for (int j = k + 1; j < n; ++j) {
                U.data[MatrixIndex(&U, i, j)] -= factor * U.data[MatrixIndex(&U, k, j)];
            }
            rhs.data[i] -= factor * rhs.data[k];
        }
    }

    /* 回代阶段：上三角系统 U x = rhs 从最后一行向第一行求解。 */
    for (int i = n - 1; i >= 0; --i) {
        REAL sum = rhs.data[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= U.data[MatrixIndex(&U, i, j)] * x->data[j];
        }
        REAL diag = U.data[MatrixIndex(&U, i, i)];
        if (fabs(diag) < tol) {
            MatrixFree(&U);
            MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }
        x->data[i] = sum / diag;
    }

    MatrixFree(&U);
    MatrixFree(&rhs);
    return MATRIX_SUCCESS;
}

MatrixError GaussianSolvePartialPivot(const Matrix *A, const Matrix *b, Matrix *x, REAL tol)
{
    MatrixError error = CheckLinearSystem(A, b, x);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    int n = A->row;
    Matrix U, rhs;
    MatrixInit(&U);
    MatrixInit(&rhs);

    error = MatrixCreate(&U, n, n);
    if (error != MATRIX_SUCCESS) {
        return error;
    }
    error = MatrixCreate(&rhs, n, 1);
    if (error != MATRIX_SUCCESS) {
        MatrixFree(&U);
        return error;
    }
    MatrixCopy(A, &U);
    MatrixCopy(b, &rhs);

    for (int k = 0; k < n - 1; ++k) {
        int pivot_row = k;
        REAL pivot_abs = fabs(U.data[MatrixIndex(&U, k, k)]);

        /* 在第 k 列的当前列段中寻找绝对值最大的主元。 */
        for (int i = k + 1; i < n; ++i) {
            REAL value_abs = fabs(U.data[MatrixIndex(&U, i, k)]);
            if (value_abs > pivot_abs) {
                pivot_abs = value_abs;
                pivot_row = i;
            }
        }
        if (pivot_abs < tol) {
            MatrixFree(&U);
            MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }

        /* 系数矩阵换行时，右端项必须同步换行。 */
        MatrixSwapRows(&U, k, pivot_row);
        MatrixSwapRows(&rhs, k, pivot_row);

        REAL pivot = U.data[MatrixIndex(&U, k, k)];
        for (int i = k + 1; i < n; ++i) {
            REAL factor = U.data[MatrixIndex(&U, i, k)] / pivot;
            U.data[MatrixIndex(&U, i, k)] = 0.0;
            for (int j = k + 1; j < n; ++j) {
                U.data[MatrixIndex(&U, i, j)] -= factor * U.data[MatrixIndex(&U, k, j)];
            }
            rhs.data[i] -= factor * rhs.data[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        REAL sum = rhs.data[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= U.data[MatrixIndex(&U, i, j)] * x->data[j];
        }
        REAL diag = U.data[MatrixIndex(&U, i, i)];
        if (fabs(diag) < tol) {
            MatrixFree(&U);
            MatrixFree(&rhs);
            return MATRIX_ERROR_SINGULAR;
        }
        x->data[i] = sum / diag;
    }

    MatrixFree(&U);
    MatrixFree(&rhs);
    return MATRIX_SUCCESS;
}

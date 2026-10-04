/*******************************************************************************
* File Name   : matrix_lu.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : LU 分解、前代和回代求解函数实现。
*******************************************************************************/

#include "matrix_lu.h"
#include <math.h>

MatrixError LUDecomposeNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol)
{
    if (!MatrixIsValid(A) || !MatrixIsValid(L) || !MatrixIsValid(U)) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    int n = A->row;
    if (!MatrixHasShape(L, n, n) || !MatrixHasShape(U, n, n)) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    /* Doolittle 型 LU 分解：L 的对角线取 1，U 保存主元行。 */
    MatrixFillZero(L);
    MatrixFillZero(U);
    for (int i = 0; i < n; ++i) {
        L->data[MatrixIndex(L, i, i)] = 1.0;
    }

    for (int k = 0; k < n; ++k) {
        /* 第 k 步：由 A=LU 的第 k 行关系计算 U(k,j), j>=k。 */
        for (int j = k; j < n; ++j) {
            REAL sum = 0.0;
            for (int s = 0; s < k; ++s) {
                sum += L->data[MatrixIndex(L, k, s)] * U->data[MatrixIndex(U, s, j)];
            }
            U->data[MatrixIndex(U, k, j)] = A->data[MatrixIndex(A, k, j)] - sum;
        }

        REAL pivot = U->data[MatrixIndex(U, k, k)];
        if (fabs(pivot) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }

        /* 第 k 步：由 A=LU 的第 k 列关系计算 L(i,k), i>k。 */
        for (int i = k + 1; i < n; ++i) {
            REAL sum = 0.0;
            for (int s = 0; s < k; ++s) {
                sum += L->data[MatrixIndex(L, i, s)] * U->data[MatrixIndex(U, s, k)];
            }
            L->data[MatrixIndex(L, i, k)] = (A->data[MatrixIndex(A, i, k)] - sum) / pivot;
        }
    }
    return MATRIX_SUCCESS;
}

MatrixError ForwardSubstitution(const Matrix *L, const Matrix *b, Matrix *y, REAL tol)
{
    if (!MatrixIsValid(L) || !MatrixIsValid(b) || !MatrixIsValid(y)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = L->row;
    if (L->row != L->column || b->row != n || b->column != 1 ||
        y->row != n || y->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    /* 求解 L y = b。第 i 行只依赖 y(0),...,y(i-1)。 */
    for (int i = 0; i < n; ++i) {
        REAL sum = b->data[i];
        for (int j = 0; j < i; ++j) {
            sum -= L->data[MatrixIndex(L, i, j)] * y->data[j];
        }
        REAL diag = L->data[MatrixIndex(L, i, i)];
        if (fabs(diag) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }
        y->data[i] = sum / diag;
    }
    return MATRIX_SUCCESS;
}

MatrixError BackSubstitution(const Matrix *U, const Matrix *y, Matrix *x, REAL tol)
{
    if (!MatrixIsValid(U) || !MatrixIsValid(y) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = U->row;
    if (U->row != U->column || y->row != n || y->column != 1 ||
        x->row != n || x->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    /* 求解 U x = y。计算顺序从最后一行向第一行推进。 */
    for (int i = n - 1; i >= 0; --i) {
        REAL sum = y->data[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= U->data[MatrixIndex(U, i, j)] * x->data[j];
        }
        REAL diag = U->data[MatrixIndex(U, i, i)];
        if (fabs(diag) < tol) {
            return MATRIX_ERROR_SINGULAR;
        }
        x->data[i] = sum / diag;
    }
    return MATRIX_SUCCESS;
}

MatrixError LUSolve(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol)
{
    if (!MatrixIsValid(L) || !MatrixIsValid(U) || !MatrixIsValid(b) || !MatrixIsValid(x)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    int n = L->row;
    if (L->row != L->column || U->row != U->column || U->row != n ||
        b->row != n || b->column != 1 || x->row != n || x->column != 1) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    Matrix y;
    MatrixInit(&y);
    MatrixError error = MatrixCreate(&y, n, 1);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    /* 先前代求 y，再回代求 x。临时矩阵 y 在函数结束前释放。 */
    error = ForwardSubstitution(L, b, &y, tol);
    if (error == MATRIX_SUCCESS) {
        error = BackSubstitution(U, &y, x, tol);
    }

    MatrixFree(&y);
    return error;
}

MatrixError LUDeterminant(const Matrix *U, REAL *det)
{
    if (!MatrixIsValid(U) || det == NULL) {
        return MATRIX_ERROR_NULL_POINTER;
    }
    if (U->row != U->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    REAL value = 1.0;
    for (int i = 0; i < U->row; ++i) {
        value *= U->data[MatrixIndex(U, i, i)];
    }
    *det = value;
    return MATRIX_SUCCESS;
}

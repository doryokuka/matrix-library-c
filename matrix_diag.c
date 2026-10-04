

/*******************************************************************************
* File Name   : matrix_diag.c
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Description : Jacobi diagonalization for real symmetric matrices.
*******************************************************************************/

#include <math.h>
#include "matrix_diag.h"
#include "matrix_core.h"

static REAL AbsReal(REAL x)
{
    return x < 0.0 ? -x : x;
}

static MatrixError MatrixSetIdentity(Matrix *A)
{
    int i;

    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    MatrixFillZero(A);

    for (i = 0; i < A->row; i++) {
        A->data[MatrixIndex(A, i, i)] = 1.0;
    }

    return MATRIX_SUCCESS;
}

MatrixError MatrixIsSymmetric(const Matrix *A, REAL tol)
{
    int i, j;

    if (!MatrixIsValid(A)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    if (tol < 0.0) {
        tol = -tol;
    }

    for (i = 0; i < A->row; i++) {
        for (j = i + 1; j < A->column; j++) {
            REAL aij = A->data[MatrixIndex(A, i, j)];
            REAL aji = A->data[MatrixIndex(A, j, i)];

            if (AbsReal(aij - aji) > tol) {
                return MATRIX_ERROR_NOT_SYMMETRIC;
            }
        }
    }

    return MATRIX_SUCCESS;
}

static void FindMaxOffDiagonal(const Matrix *D, int *p, int *q, REAL *max_value)
{
    int i, j;
    int n = D->row;

    *p = 0;
    *q = 1;
    *max_value = 0.0;

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            REAL value = AbsReal(D->data[MatrixIndex(D, i, j)]);

            if (value > *max_value) {
                *max_value = value;
                *p = i;
                *q = j;
            }
        }
    }
}

MatrixError MatrixDiagonalizeSymmetricJacobi(
    const Matrix *A,
    Matrix *Q,
    Matrix *D,
    int max_iter,
    REAL tol
)
{
    int n;
    int iter;
    MatrixError error;

    if (!MatrixIsValid(A) || !MatrixIsValid(Q) || !MatrixIsValid(D)) {
        return MATRIX_ERROR_NULL_POINTER;
    }

    if (A->row != A->column) {
        return MATRIX_ERROR_NOT_SQUARE;
    }

    n = A->row;

    if (!MatrixHasShape(Q, n, n) || !MatrixHasShape(D, n, n)) {
        return MATRIX_ERROR_SIZE_MISMATCH;
    }

    if (max_iter <= 0) {
        return MATRIX_ERROR_INVALID_SIZE;
    }

    if (tol < 0.0) {
        tol = -tol;
    }

    error = MatrixIsSymmetric(A, tol);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixCopy(A, D);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    error = MatrixSetIdentity(Q);
    if (error != MATRIX_SUCCESS) {
        return error;
    }

    if (n == 1) {
        return MATRIX_SUCCESS;
    }

    for (iter = 0; iter < max_iter; iter++) {
        int p, q, k;
        REAL max_offdiag;
        REAL app, aqq, apq;
        REAL theta, c, s;
        REAL dpp, dqq;

        FindMaxOffDiagonal(D, &p, &q, &max_offdiag);

        if (max_offdiag < tol) {
            return MATRIX_SUCCESS;
        }

        app = D->data[MatrixIndex(D, p, p)];
        aqq = D->data[MatrixIndex(D, q, q)];
        apq = D->data[MatrixIndex(D, p, q)];

        if (AbsReal(apq) < tol) {
            D->data[MatrixIndex(D, p, q)] = 0.0;
            D->data[MatrixIndex(D, q, p)] = 0.0;
            continue;
        }

        theta = 0.5 * atan2(2.0 * apq, aqq - app);
        c = cos(theta);
        s = sin(theta);

        for (k = 0; k < n; k++) {
            if (k != p && k != q) {
                REAL dkp = D->data[MatrixIndex(D, k, p)];
                REAL dkq = D->data[MatrixIndex(D, k, q)];

                REAL new_dkp = c * dkp - s * dkq;
                REAL new_dkq = s * dkp + c * dkq;

                D->data[MatrixIndex(D, k, p)] = new_dkp;
                D->data[MatrixIndex(D, p, k)] = new_dkp;
                D->data[MatrixIndex(D, k, q)] = new_dkq;
                D->data[MatrixIndex(D, q, k)] = new_dkq;
            }
        }

        dpp = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        dqq = s * s * app + 2.0 * s * c * apq + c * c * aqq;

        D->data[MatrixIndex(D, p, p)] = dpp;
        D->data[MatrixIndex(D, q, q)] = dqq;
        D->data[MatrixIndex(D, p, q)] = 0.0;
        D->data[MatrixIndex(D, q, p)] = 0.0;

        for (k = 0; k < n; k++) {
            REAL qkp = Q->data[MatrixIndex(Q, k, p)];
            REAL qkq = Q->data[MatrixIndex(Q, k, q)];

            Q->data[MatrixIndex(Q, k, p)] = c * qkp - s * qkq;
            Q->data[MatrixIndex(Q, k, q)] = s * qkp + c * qkq;
        }
    }

    return MATRIX_ERROR_NOT_CONVERGED;
}

#ifndef MATRIX_DIAG_H
#define MATRIX_DIAG_H

#include "matrix_core.h"

/*
 * Check whether a matrix is symmetric within a given tolerance.
 *
 * Return:
 * MATRIX_SUCCESS if A is symmetric.
 * MATRIX_ERROR_NOT_SQUARE if A is not square.
 * MATRIX_ERROR_NULL_POINTER if A is NULL.
 * MATRIX_ERROR_NOT_SYMMETRIC if A is not symmetric.
 */
MatrixError MatrixIsSymmetric(const Matrix *A, REAL tol);

/*
 * Diagonalize a real symmetric matrix using the Jacobi method.
 *
 * For a symmetric matrix A, this function computes matrices Q and D such that:
 *
 *     Q^T * A * Q ≈ D
 *
 * where D is approximately diagonal and Q is approximately orthogonal.
 *
 * Parameters:
 * A        input n x n symmetric matrix
 * Q        output n x n orthogonal matrix
 * D        output n x n diagonal matrix
 * max_iter maximum number of Jacobi iterations
 * tol      tolerance for convergence
 *
 * Return:
 * MATRIX_SUCCESS if the algorithm converges.
 * MATRIX_ERROR_NULL_POINTER if any pointer is NULL.
 * MATRIX_ERROR_NOT_SQUARE if A is not square.
 * MATRIX_ERROR_SIZE_MISMATCH if matrix sizes do not match.
 * MATRIX_ERROR_NOT_SYMMETRIC if A is not symmetric.
 * MATRIX_ERROR_NOT_CONVERGED if the algorithm does not converge.
 */
MatrixError MatrixDiagonalizeSymmetricJacobi(
    const Matrix *A,
    Matrix *Q,
    Matrix *D,
    int max_iter,
    REAL tol
);

#endif

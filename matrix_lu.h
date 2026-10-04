/*******************************************************************************
* File Name   : matrix_lu.h
* Course      : C Program and Algorithm Analysis Lab -- Matrix Solver Library
* Author      : Gong Helin
* Year        : 2026
* Description : LU 分解与三角求解接口声明。
*******************************************************************************/

#ifndef MATRIX_LU_H
#define MATRIX_LU_H

#include "matrix_core.h"

MatrixError LUDecomposeNoPivot(const Matrix *A, Matrix *L, Matrix *U, REAL tol);
MatrixError ForwardSubstitution(const Matrix *L, const Matrix *b, Matrix *y, REAL tol);
MatrixError BackSubstitution(const Matrix *U, const Matrix *y, Matrix *x, REAL tol);
MatrixError LUSolve(const Matrix *L, const Matrix *U, const Matrix *b, Matrix *x, REAL tol);
MatrixError LUDeterminant(const Matrix *U, REAL *det);

#endif

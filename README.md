# Matrix Computation Library in C

A course-based C99 matrix library extended with **Jacobi diagonalization for real symmetric matrices**.

This repository is based on starter code provided for the **C Program and Algorithm Analysis Lab — Matrix Solver Library** (2026). The original course files credit **Gong Helin**; their author and course notices are preserved.

中文简介：本项目基于课程老师提供的矩阵计算库框架，新增实对称矩阵的 Jacobi 正交对角化。本人确认的贡献仅为 Jacobi 及其相关接口、错误处理、测试和构建接入；LU 分解、转置等已有功能不列为本人贡献。

## My Contributions

The confirmed contribution is the Jacobi extension, supported by the source changes relative to the supplied `lec7code.zip`:

- **Jacobi implementation and API:** `matrix_diag.c` and `matrix_diag.h`, including `MatrixIsSymmetric` and `MatrixDiagonalizeSymmetricJacobi`.
- **Error handling:** symmetry and convergence checks, with `MATRIX_ERROR_NOT_SYMMETRIC` and `MATRIX_ERROR_NOT_CONVERGED` added to the core error codes and messages.
- **Testing:** the Jacobi test in `test_runner.c` checks a symmetric 3 × 3 example, successful convergence, near-zero off-diagonal entries of `D`, and `Q^T A Q ≈ D`.
- **Build integration:** `matrix_diag.o` and the relevant header dependencies are included in the Makefile.

Matrix transpose, LU decomposition, Gaussian elimination, and other existing functionality are included as part of the supplied course code. They are not claimed as my contributions.

## Features

| Module | Functionality | Origin in this version |
| --- | --- | --- |
| `matrix_core` | Matrix storage, allocation, access, copying, printing, error messages | Course starter; extended with Jacobi error codes/messages |
| `matrix_ops` | Addition, subtraction, scaling, transpose, Frobenius norm, multiplication | Supplied course code |
| `matrix_solve` | Gaussian elimination with and without partial pivoting | Supplied course code |
| `matrix_lu` | LU decomposition without pivoting, triangular solves, determinant | Supplied course code |
| `matrix_diag` | Symmetry check and Jacobi orthogonal diagonalization | Confirmed extension |

For a real symmetric matrix `A`, the Jacobi method computes `Q` and `D` such that:

```text
Q^T A Q ≈ D
```

`Q` contains the corresponding eigenvectors as columns, and the diagonal entries of `D` approximate the eigenvalues. Eigenvalues are not sorted. Each iteration selects the largest off-diagonal entry in magnitude and applies a plane rotation; convergence is checked against the absolute tolerance.

## Project Structure

The original flat source layout is retained.

```text
matrix-library-c/
├── README.md
├── README_Lab7.md          Original course documentation
├── .gitignore
├── Makefile
├── matrix_core.c / .h
├── matrix_ops.c / .h
├── matrix_solve.c / .h
├── matrix_lu.c / .h
├── matrix_diag.c / .h      Jacobi extension
├── test_runner.c
├── benchmark.c
└── student_practice.c
```

## Build and Run

Requirements: a C99 compiler (GCC or Clang), Make, and the standard C math library. Run the following commands from this directory:

```sh
make
make run-tests
make run-practice
```

`make` builds `lab7_tests`, `lab7_benchmark`, and `lab7_practice`. To select Clang explicitly, use `make CC=clang`.

Optional performance benchmarks for the existing multiplication and linear solver routines:

```sh
make run-bench
```

Remove generated object files and executables:

```sh
make clean
```

## Jacobi API

```c
#include "matrix_diag.h"

MatrixError MatrixDiagonalizeSymmetricJacobi(
    const Matrix *A,
    Matrix *Q,
    Matrix *D,
    int max_iter,
    REAL tol
);
```

Initialize and allocate `A`, `Q`, and `D` with the core API before calling the function. All three matrices must be `n × n`; use separate matrices for the input and outputs. Supply a positive iteration budget and a positive finite tolerance (for example, `200` and `1e-12` for the included 3 × 3 test). Inspect the returned `MatrixError` before using the result, and release allocations with `MatrixFree` afterwards.

See `TestMatrixDiagonalizeSymmetricJacobi` in `test_runner.c` for a complete allocation, input, call, verification, and cleanup example.

## Validation

The included test runner reports **50 / 50 checks passed** on macOS with Apple Clang during repository preparation. This includes 30 existing course checks and 20 checks in the Jacobi test. All three programs build with the Makefile's C99 warning flags, and the practice program reports matching Gaussian and LU solutions.

The Jacobi test verifies convergence and the diagonalization relation for one 3 × 3 matrix. It does not directly test `Q^T Q ≈ I`, non-symmetric rejection, or the iteration-limit error path; those checks should not be inferred from the passing count.

## Acknowledgement and Source Provenance

- **Starter code:** provided by the course instructor, credited to Gong Helin (2026) in the supplied files. `README_Lab7.md` is retained unchanged from the archive as the original course reference.
- **Confirmed extension:** Jacobi diagonalization and its related API, error handling, tests, and Makefile integration.
- **Preparation sources:** the supplied `lec7code.zip` contains the starter version. The sibling local `lec7code` folder contains the Jacobi extension; its two new module files and four modified integration files were copied into this repository. The other source files match the archive.
- **Preparation edits:** this README and `.gitignore` were added. Final newlines were added to the two Jacobi module files to remove compiler warnings. No algorithm or source layout was changed.

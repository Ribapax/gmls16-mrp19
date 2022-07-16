#ifndef LU_FACTORIZATION_HEADER
#define LU_FACTORIZATION_HEADER

#include "../common/common.h"

// invertMatrix returns the inverted matrix of the matrix A of size n
// 1) Get L and U by solving -> L = generateMatrixL(A, B, U, n);
// 2) Get the y arrays by solving -> Ly = b;
// 3) Get the inverted matrix x by solving -> Ux = y for each y.
// 4) Returns x.
RealNumber **InvertMatrix(RealNumber **A, int n);

// generateMatrixL returns the multiplier's matrix L of a matrix A of size n
// using Gauss Elimination with partial pivoting. The partial pivoting should
// also replace lines in the B matrix (identity matrix).
// Also stores the echelon form of A inside U.
RealNumber **generateMatrixL(RealNumber **A, RealNumber **B, RealNumber **U, int n);

// CalculateResidue returns the L2 Norm of the residue for the current result.
// R = I - A * A^-1
// L2Norm = sqrt(sum(R[i][j]^2)), ∀ 1 <= i, j <= n
RealNumber CalculateResidue(RealNumber **A, RealNumber **invertedA, int n);

// finPivotIndex finds the line index of the pivot element of the given columnIndex
unsigned int findPivotIndex(double** Matrix, unsigned int columnIndex, unsigned int systemSize);

// replaceLines replaces the line of index 'index' with the line of index 'pivotIndex', and vice versa.
void replaceLines(
    double **Matrix,
    double **identityMatrix,
    unsigned int index,
    unsigned int pivotIndex
);

#endif
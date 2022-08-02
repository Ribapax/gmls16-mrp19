/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#ifndef LU_FACTORIZATION_HEADER
#define LU_FACTORIZATION_HEADER

#include "../common/common.h"

#define RESIDUE_THRESHOLD 0.00000000000001

// Returns the inverted matrix of the matrix `A` of size `n`.
//  - It also stores the average time of solving the Linear Systems
//    Ly = b and Ux = y inside `averageLinearSystemTime`;
//  - And stores the L and U matrix calculation time inside `LUTime`.
// Algorithm:
// 1) Get L and U by solving -> L = generateMatrixL(A, B, U, n);
// 2) Get the y arrays by solving -> Ly = b;
// 3) Get the inverted matrix x by solving -> Ux = y for each y.
// 4) Return x.
RealNumber **InvertMatrix(
    RealNumber **A,
    int n,
    Time *averageLinearSystemTime,
    Time *LUTime
);

// Returns the multiplier's matrix L of a matrix `A` of size `n`
// using Gaussian Elimination with partial pivoting. The partial pivoting should
// also replace lines in the `B` matrix (identity matrix).
// Also stores the echelon form of A inside `U`.
RealNumber **generateMatrixL(RealNumber **A, RealNumber **B, RealNumber **U, int n);

// Returns the L2 Norm of the residue for the current result.
// R = I - A * A^-1
// L2Norm = sqrt(sum(R[i][j]^2)), ∀ 1 <= i, j <= n
// Params:
// - A: matrix to be inverted;
// - invertedA: matrix A inverted;
// - n: matrix dimension.
RealNumber CalculateResidue(RealNumber **A, RealNumber **invertedA, int n);

// Finds the line index of the pivot element of the given column index.
// Params:
// - Matrix: matrix that will be pivoted;
// - columnIndex: index of the column that will be pivoted;
// - systemSize: matrix dimension.
unsigned int findPivotIndex(double** Matrix, unsigned int columnIndex, unsigned int systemSize);

// Replaces the line of index 'index' with the line of index 'pivotIndex', and vice versa.
// This operation is made in both `Matrix` and `identityMatrix`.
void replaceLinesWithIdentityMatrix(
    double **Matrix,
    double **identityMatrix,
    unsigned int index,
    unsigned int pivotIndex
);

#endif
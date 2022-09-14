/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#ifndef LU_FACTORIZATION_HEADER
#define LU_FACTORIZATION_HEADER

#include "../common/common.h"

#define RESIDUE_THRESHOLD 0.0000000000000001

#define ENABLE_PARTIAL_PIVOTING 0

/*
 * Returns the solution of `n` Linear Systems by solving AX = B using the
 *    matrices L and U, previously calculated by the LU decomposition.
 * Params:
 * - B: the independent term's matrix;
 * - n: the matrices dimension;
 * - averageLinearSystemTime: pointer to the variable which accumulates the time
 *   to solve all linear systems in this program;
 *  - L: Lower matrix from LU decomposition;
 *  - U: Upper matrix from LU decomposition;
 */
RealNumber *SolveLinearSystems(
    const RealNumber *B,
    int n,
    Time *averageLinearSystemTime,
    const RealNumber *L,
    const RealNumber *U
);

/*
 *  Decomposes the Matrix A of size n into L (Lower) and U (Upper).
 *  Params:
 *  - A: Matrix to be decomposed;
 *  - B: identity matrix pointer, in case some line is changed by partial pivoting;
 *  - U: Upper Matrix pointer, should be allocated;
 *  - L: Lower Matrix pointer, should be allocated;
 *  - n: matrix dimension;
 *  Returns:
 *  - 0 (success);
 *  - 1 (error);
 */
int LUDecomposition(RealNumber *A, RealNumber *U, RealNumber *L, int n);

/*
 * Returns the L2 Norm of the residue for the current result.
 *  R = I - A * A^-1
 *  L2Norm = sqrt(sum(R[i][j]^2)), ∀ 1 <= i, j <= n
 * Params:
 * - A: matrix to be inverted;
 * - B: identity matrix;
 * - invertedA: matrix A inverted;
 * - n: matrix dimension.
 */
RealNumber CalculateResidueL2Norm(RealNumber *A, RealNumber *B, RealNumber *invertedA, int n);

/*
 * Returns the line index of the pivot element of the given column index.
 * Params:
 * - Matrix: matrix that will be pivoted;
 * - columnIndex: index of the column that will be pivoted;
 * - systemSize: matrix dimension.
 */
unsigned int findPivotIndex(double* Matrix, unsigned int columnIndex, unsigned int systemSize);

// Replaces the line of index 'index' with the line of index 'pivotIndex', and vice versa.
void replaceLinesWithIdentityMatrix(
    double *Matrix,
    unsigned int index,
    unsigned int pivotIndex,
    unsigned int n
);

#endif
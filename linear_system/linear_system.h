/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#ifndef LINEAR_SYSTEM_HEADER
#define LINEAR_SYSTEM_HEADER
#include "../common/common.h"
#include "../lu_factorization/lu_factorization.h"

// Helper to access unidimensional arrays as they were matrices
#define Index(i, j, n) (((i) * (n)) + (j))

// Max value of the generated Linear System coefficients
#define COEFFICIENT_LIMIT 32.0

// TODO: description
int FillMatrix(RealNumber *A, RealNumber coefficientLimit, unsigned int n);

// Copies matrix `A` into matrix `B`, both of dimension `n`.
void copyMatrix(const RealNumber *A, RealNumber *B, int n);

// Multiplies the matrices `A` and `B`, both of dimension `n`.
// Returns a new matrix of size `n`.
RealNumber *multiplyMatricesOfEqualSize(const RealNumber *A, const RealNumber *B, int n);

// Returns the result of `A` - `B`, both of dimension `n`.
RealNumber *subtractMatrices(const RealNumber *A, const RealNumber *B, int n);

// Returns the identity matrix of dimension `n`.
RealNumber *GenerateIdentityMatrix(int n);

// Returns 1 if the matrix `A` of dimension `n` is invertible, and 0 if not.
int MatrixIsInvertible(const RealNumber *A, int n);

// TODO: description
RealNumber* AllocateMatrix(unsigned int n);


// TODO: update description
/*
 * Given a first solution, returns a more refined solution.
 * Params:
 *  - A: matrix to be inverted;
 *  - B: identity matrix;
 *  - invertedMatrix: first solution, inverted A;
 *  - L: Lower matrix from LU decomposition;
 *  - U: Upper matrix from LU decomposition;
 *  - averageLinearSystemTime: pointer to the variable which accumulates the time
 *    to solve all linear systems in this program;
 *  - n: matrix dimension;
 * Returns a matrix of dimension `n` with the refined solution.
 */

 RealNumber * RefineSolution(
    RealNumber * restrict A,
    RealNumber * restrict B,
    RealNumber * restrict invertedMatrix,
    RealNumber * restrict L,
    RealNumber * restrict U,
    PivotArray *P,
    Time *averageLinearSystemTime,
    int n
);

/*
 * Asserts if the solution has reached the stopping criteria.
 * Also, test if the last residue L2 norm is smaller than the current one.
 * If the residue is getting bigger, this function returns 0, and stops the execution.
 * Params:
 *  - iteration: current iteration;
 *  - iterationsLimit: maximum amount of iterations;
 *  - currentResidueL2Norm: the current L2 norm of the residue of the current iteration solution;
 *  - lastResidueL2Norm: the last L2 norm of the residue of the solution of the last iteration;
 * Returns:
 *  - 1 (yes);
 *  - 0 (no).
 */
int HasNotReachedStoppingCriteria(
    int iteration,
    int iterationsLimit,
    RealNumber currentResidueL2Norm,
    RealNumber lastResidueL2Norm
);

// Verifies if the L2 Norm of the Residue has increased based on the last iteration
// Params:
// - currentResidueL2Norm: the L2 Norm of the Residue of the current iteration;
// - lastResidueL2Norm: the L2 Norm of the Residue of the last iteration;
// Returns:
// - 1 (true);
// - 0 (false).
int ResidueIsIncreasing(RealNumber currentResidueL2Norm, RealNumber lastResidueL2Norm);

#endif


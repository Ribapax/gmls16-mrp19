/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#ifndef LINEAR_SYSTEM_HEADER
#define LINEAR_SYSTEM_HEADER
#include "../common/common.h"

// Max value of the generated Linear System coefficients
#define COEFFICIENT_LIMIT 32.0

// Matrix allocation type
typedef enum {
    PointerToPointer = 0, // Array of N pointers to arrays of size N
    PointerToArray       // Array of N pointers for a single array of size N * N
} MatrixAllocationType;

// Linear System struct
typedef struct {
    RealNumber **A; // coefficients
    RealNumber *b; // independent terms
    unsigned int n; // size
    MatrixAllocationType matrixAllocationType; // allocation type of the coefficients matrix
} LinearSystem;

// Matrix types
typedef enum {
    GenericMatrix = 0,
    HilbertMatrix,
    DominantDiagonal,
    NullEquation,
    ProportionalEquation,
    LinearCombinationEquation
} MatrixType;

/*
 * Allocate a whole linear system.
 * Params:
 *  - n: linear system dimension;
 *  - type: matrix allocation type (enumerator).
 * Returns:
 * - NULL (if an error occurs);
 * - LinearSystem* (a pointer to a linear system correctly allocated).
 */
LinearSystem* AllocateLinearSystem(unsigned int n, MatrixAllocationType type);

/*
 * Free all memory allocated to the given Linear System.
 * Params:
 *  - SL: a linear system pointer;
 * Returns nothing.
 */
void freeLinearSystem(LinearSystem *SL);

/*
 * Fills a Linear System with random data.
 * Params:
 *  - SL: a linear system pointer;
 *  - type: type of matrix that will be allocated (enumerator);
 *  - coefficientLimit: maximum value of the matrix's coefficients.
 * Returns:
 * - 0 (success);
 * - -1 (error).
 */
int FillLinearSystem(LinearSystem *SL, MatrixType type, RealNumber coefficientLimit);

/*
 * Prints an array to stdout. Useful for debugging.
 * Params:
 *  - arr: array to be printed;
 *  - n: array size;
 * Returns nothing.
 */
void printArray(RealNumber *arr, unsigned int n);

// Copies matrix `A` into matrix `B`, both of dimension `n`.
void copyMatrix(RealNumber **A, RealNumber **B, int n);

// Multiplies the matrices `A` and `B`, both of dimension `n`.
// Returns a new matrix of size `n`.
RealNumber **multiplyMatrixOfEqualSize(RealNumber **A, RealNumber **B, int n);

// Multiplies the matrix `A` with the array `B`, both of dimension `n`.
// Returns a new array of size `n`.
RealNumber *multiplyMatrixWithArray(RealNumber **A, const RealNumber *B, int n);

// Returns the result of `A` - `B`, both of dimension `n`.
RealNumber **subtractMatrix(RealNumber **A, RealNumber **B, int n);

// Returns the result of `A` - `B`, both of dimension `n`.
RealNumber *subtractArrays(const RealNumber *A, const RealNumber *B, int n);

// Returns the identity matrix of dimension `n`.
RealNumber **GetIdentityMatrix(int n);

// Returns 1 if the matrix `A` of dimension `n` is invertible, and 0 if not.
int MatrixIsInvertible(RealNumber **A, int n);

/*
 * Executes a Gaussian Elimination the given linear system.
 * Params:
 *  - A: coefficients matrix;
 *  - B: independent terms;
 *  - n: linear system dimension;
 * Returns an array with the solution of the system.
 * TODO: error handling.
 */
RealNumber *GaussElimination(RealNumber **A, RealNumber *B, int n);

/*
 * Given a first solution, returns a more refined solution.
 * Params:
 *  - A: matrix to be inverted;
 *  - B: identity matrix;
 *  - invertedMatrix: first solution, inverted A;
 *  - n: matrix dimension;
 * Returns a matrix of dimension `n` with the refined solution.
 */
RealNumber **refineSolution(RealNumber **A, RealNumber **B, RealNumber **invertedMatrix, int n);

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
int hasNotReachedStoppingCriteria(
    int iteration,
    int iterationsLimit,
    RealNumber currentResidueL2Norm,
    RealNumber lastResidueL2Norm
);

// Executes a forward substitution on the given triangular Linear System.
// Params:
// - A: coefficient matrix;
// - b: independent terms;
// - x: result array;
// - n: linear system dimension.
// Returns nothing.
void forwardSubstitution(RealNumber **A, const RealNumber *b, RealNumber *x, unsigned int n);

// Replace lines of the matrix and independent terms.
// Params:
// - Matrix: coefficient matrix;
// - independentTerms: array of independent terms;
// - index: index of the line that will be replaced;
// - pivotIndex: index of the line replacing.
// Returns nothing.
void replaceLines(
    double **Matrix,
    double *independentTerms,
    unsigned int index,
    unsigned int pivotIndex
);

#endif


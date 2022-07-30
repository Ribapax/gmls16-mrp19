#ifndef LINEAR_SYSTEM_HEADER
#define LINEAR_SYSTEM_HEADER
#include "../common/common.h"

#define hasNotReachedStoppingCriteria(it, limit, res) (it <= limit && res > RESIDUE_THRESHOLD)

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


LinearSystem* AllocateLinearSystem(unsigned int n, MatrixAllocationType type);
void freeLinearSystem(LinearSystem *SL);
int FillLinearSystem(LinearSystem *SL, MatrixType type, RealNumber coefficientLimit);

LinearSystem *readLinearSystem();
void printLinearSystem(LinearSystem *SL);
void printArray(RealNumber *arr, unsigned int n);

// copyMatrix copies A into B
void copyMatrix(RealNumber **A, RealNumber **B, int n);

// multiplyMatrixOfEqualSize returns the result of A x B
RealNumber **multiplyMatrixOfEqualSize(RealNumber **A, RealNumber **B, int n);

RealNumber *multiplyMatrixWithArray(RealNumber **A, const RealNumber *B, int n);

// subtractMatrix returns the result of A - B
RealNumber **subtractMatrix(RealNumber **A, RealNumber **B, int n);

RealNumber *subtractArrays(const RealNumber *A, const RealNumber *B, int n);

// GetIdentityMatrix returns an identity matrix of size n
RealNumber **GetIdentityMatrix(int n);

// MatrixIsInvertible returns 1 if the matrix is invertible and 0 if not
int MatrixIsInvertible(RealNumber **A, int n);

RealNumber *GaussElimination(RealNumber **A, RealNumber *B, int n);

#endif


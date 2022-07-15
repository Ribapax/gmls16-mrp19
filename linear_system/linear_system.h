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


LinearSystem* allocateLinearSystem(unsigned int n, MatrixAllocationType type);
void freeLinearSystem(LinearSystem *SL);
void fillLinearSystem(LinearSystem *SL, MatrixType type, RealNumber coefficientLimit);

LinearSystem *readLinearSystem();
void printLinearSystem(LinearSystem *SL);
void printArray(RealNumber *arr, unsigned int n);

#endif


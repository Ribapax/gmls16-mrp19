#include <stdio.h>
#include <math.h>
#include "lu_factorization.h"
#include "linear_system.h"

unsigned int findPivotIndex(double** Matrix, unsigned int columnIndex, unsigned int systemSize) {
    RealNumber greatestValue = Matrix[columnIndex][columnIndex];
    unsigned int pivotIndex = columnIndex;
    for (unsigned int i = columnIndex; i < systemSize; i++) {
        // FIXME: we probably need a decent float comparison here
        if (fabs(Matrix[i][columnIndex]) > fabs(greatestValue)) {
            greatestValue = Matrix[i][columnIndex];
            pivotIndex = i;
        }
    }
    return pivotIndex;
}

void replaceLines(
    double **Matrix,
    double **identityMatrix,
    unsigned int index,
    unsigned int pivotIndex
) {
    RealNumber *lineToBeReplaced = Matrix[index];
    Matrix[index] = Matrix[pivotIndex];
    Matrix[pivotIndex] = lineToBeReplaced;

    RealNumber *identityLineToBeReplaced = identityMatrix[index];
    identityMatrix[index] = identityMatrix[pivotIndex];
    identityMatrix[pivotIndex] = identityLineToBeReplaced;
}

RealNumber **generateMatrixL(RealNumber **A, RealNumber **B, RealNumber **U, int n) {
    RealNumber **L = allocateLinearSystem(n, PointerToPointer)->A;
    copyMatrix(A, U, n);
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;

        // Partial pivoting
        unsigned int pivotIndex = findPivotIndex(U, i, n);
        if (i != pivotIndex) {
            replaceLines(U, B, i, pivotIndex);
        }

        // Gauss Elimination
        for (int k = i + 1; k < n; k++) {
            if (U[k][k] == 0) {
                fprintf(stderr, "%s\n", "gaussian elimination error: division by zero");
                return NULL; // TODO: return error codes
            }
            double m = U[k][i] / U[i][i];

            // Filling multiplier's matrix
            L[k][i] = m;

            U[k][i] = 0.0;
            for (int j = i + 1; j < n; j++) {
                U[k][j] -= U[i][j] * m;
            }
        }
    }

    return L;
}

RealNumber **InvertMatrix(RealNumber **A, int n) {
    // 1) Get L and U by solving -> L = generateMatrixL(A, B, U, n);
    RealNumber **invertedMatrix = allocateLinearSystem(n, PointerToPointer)->A;
    RealNumber **B = getIdentityMatrix(n);
    RealNumber **U = allocateLinearSystem(n, PointerToPointer)->A;
    RealNumber **L = generateMatrixL(A, B, U, n);

    // 2) Get the y arrays by solving -> Ly = b;
        // TODO

    // 3) Get the inverted matrix x by solving -> Ux = y for each y.
        // TODO

    // 4) Return x.
    return invertedMatrix;
}

RealNumber CalculateResidue(RealNumber **A, RealNumber **invertedA, int n) {
    RealNumber **multiplication = multiplyMatrix(A, invertedA, n);
    RealNumber **identityMatrix = getIdentityMatrix(n);
    RealNumber **R = subtractMatrix(identityMatrix, multiplication, n);
    RealNumber sum = 0.;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += pow(R[i][j], 2);
        }
    }
    return sqrt(sum);
}
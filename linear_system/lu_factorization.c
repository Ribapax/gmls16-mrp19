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
    RealNumber **L = AllocateLinearSystem(n, PointerToPointer)->A;
    if (L == NULL) {
        return NULL;
    }
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
                return NULL;
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

// TODO: handle errors
RealNumber **InvertMatrix(
    RealNumber **A,
    int n,
    Time *averageLinearSystemTime,
    Time *LUTime
) {
    Time time; // Stores the time of solving the linear systems
    RealNumber **invertedMatrix = AllocateLinearSystem(n, PointerToPointer)->A;
    if (invertedMatrix == NULL) {
        fprintf(stderr, "could not allocate inverted matrix\n");
        return NULL;
    }

    // 1) Get L and U by solving -> L = generateMatrixL(A, B, U, n);
    RealNumber **B = getIdentityMatrix(n);
    if (B == NULL) {
        fprintf(stderr, "could not allocate identity matrix\n");
        return NULL;
    }
    RealNumber **U = AllocateLinearSystem(n, PointerToPointer)->A;
    if (U == NULL) {
        fprintf(stderr, "could not allocate \"U\" matrix\n");
        return NULL;
    }
    *LUTime = timestamp();
    RealNumber **L = generateMatrixL(A, B, U, n);
    if (L == NULL) {
        fprintf(stderr, "could not generate \"L\" matrix\n");
        return NULL;
    }
    *LUTime = timestamp() - *LUTime;

    // 2) Get the y arrays by solving -> Ly = b;
    *averageLinearSystemTime = 0;
    for (int i = 0; i < n; ++i) {
        time = timestamp();
        // TODO: solve Ly = b here
        time = timestamp() - time;
        *averageLinearSystemTime += time;
    }

    // 3) Get the inverted matrix x by solving -> Ux = y for each y.
    for (int i = 0; i < n; ++i) {
        time = timestamp();
        // TODO: solve Ux = y here
        time = timestamp() - time;
        *averageLinearSystemTime += time;
    }

    *averageLinearSystemTime /= (n + n);

    // 4) Return x.
    //return invertedMatrix;

    // TODO: remove it later
    RealNumber **stubMatrix = AllocateLinearSystem(n, PointerToPointer)->A;
    stubMatrix[0][0] = 25;
    stubMatrix[0][1] = 5;
    stubMatrix[0][2] = 1;

    stubMatrix[1][0] = 64;
    stubMatrix[1][1] = 8;
    stubMatrix[1][2] = 1;

    stubMatrix[2][0] = 144;
    stubMatrix[2][1] = 12;
    stubMatrix[2][2] = 1;
    return stubMatrix;
}

RealNumber CalculateResidue(RealNumber **A, RealNumber **invertedA, int n) {
    RealNumber **multiplication = multiplyMatrix(A, invertedA, n);
    RealNumber **identityMatrix = getIdentityMatrix(n);
    if (identityMatrix == NULL) {
        fprintf(stderr, "could not calculate residue: could not allocate identity matrix\n");
        exit(-1);
    }
    RealNumber **R = subtractMatrix(identityMatrix, multiplication, n);
    RealNumber sum = 0.;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += pow(R[i][j], 2);
        }
    }
    return sqrt(sum);
}
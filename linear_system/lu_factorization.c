/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <math.h>


#include "lu_factorization.h"
#include "linear_system.h"

unsigned int findPivotIndex(double** Matrix, unsigned int columnIndex, unsigned int systemSize) {
    RealNumber greatestValue = fabs(Matrix[columnIndex][columnIndex]);
    unsigned int pivotIndex = columnIndex;
    for (unsigned int i = columnIndex + 1; i < systemSize; i++) {
        RealNumber v = fabs(Matrix[i][columnIndex]);
        if (v > greatestValue) {
            greatestValue = v;
            pivotIndex = i;
        }
    }
    return pivotIndex;
}

void replaceLinesWithIdentityMatrix(
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

RealNumber **LUDecomposition(
    RealNumber **A,
    RealNumber **B,
    RealNumber **U,
    RealNumber **L,
    int n
) {
    copyMatrix(A, U, n);
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;

        // Partial pivoting
        if (ENABLE_PARTIAL_PIVOTING) {
            unsigned int pivotIndex = findPivotIndex(U, i, n);
            if (i != pivotIndex) {
                replaceLinesWithIdentityMatrix(U, B, i, pivotIndex);
            }
        }

        // Triangularization
        for (int k = i + 1; k < n; k++) {
            if (U[i][i] == 0) {
                fprintf(stderr, "%s\n", "error: division by zero");
                return NULL;
            }
            double m = U[k][i] / U[i][i];

            L[k][i] = m; // Filling multiplier's matrix

            U[k][i] = 0.0;
            for (int j = i + 1; j < n; j++) {
                U[k][j] -= U[i][j] * m;
            }
        }
    }

    return L;
}

// TODO: calculate times
RealNumber **InvertMatrix(
    RealNumber **A,
    RealNumber **B,
    int n,
    Time *averageLinearSystemTime,
    Time *LUTime,
    RealNumber **L,
    RealNumber **U
) {
    Time time; // Stores the time of solving the linear systems
    RealNumber **invertedMatrix = AllocateLinearSystem(n, PointerToPointer)->A;
    if (invertedMatrix == NULL) {
        fprintf(stderr, "could not allocate inverted matrix\n");
        return NULL;
    }

    // 1) Get L and U by solving -> L = LUDecomposition(A, B, U, n);
    RealNumber **Y = AllocateLinearSystem(n, PointerToPointer)->A;
    if (Y == NULL) {
        fprintf(stderr, "could not allocate \"Y\" matrix\n");
        return NULL;
    }
    RealNumber **X = AllocateLinearSystem(n, PointerToPointer)->A;
    if (X == NULL) {
        fprintf(stderr, "could not allocate \"X\" matrix\n");
        return NULL;
    }

    *LUTime = timestamp();
    LUDecomposition(A, B, U, L, n);
    if (L == NULL) {
        fprintf(stderr, "could not generate \"L\" matrix\n");
        return NULL;
    }
    *LUTime = timestamp() - *LUTime;

    // 2) Get the y arrays by solving -> Ly = b;
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; i++) {
            Y[i][k] = B[i][k];
            for (unsigned int j = 0; j < i; j++) {
                Y[i][k] -= L[i][j] * Y[j][k];
            }
            Y[i][k] /= L[i][i];
        }
    }

    // 3) Get the inverted matrix x by solving -> Ux = y for each y.
    for (int k = 0; k < n; ++k) {
        for (int i = n - 1; i >= 0; i--) {
            X[i][k] = Y[i][k];
            for (unsigned int j = i + 1; j < n; j++) {
                X[i][k] -= U[i][j] * X[j][k];
            }
            X[i][k] /= U[i][i];
        }
    }

    return X;
}

RealNumber CalculateResidue(RealNumber **A, RealNumber **B, RealNumber **invertedA, int n) {
    RealNumber **multiplication = multiplyMatrixOfEqualSize(A, invertedA, n);
    RealNumber **R = subtractMatrix(B, multiplication, n);
    RealNumber sum = 0.;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += pow(R[i][j], 2);
        }
    }
    return sqrt(sum);
}
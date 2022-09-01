/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <math.h>
#include "lu_factorization.h"
#include "../linear_system/linear_system.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

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

int LUDecomposition(
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
                return -1;
            }
            double m = U[k][i] / U[i][i];

            L[k][i] = m; // Filling multiplier's matrix

            U[k][i] = 0.0;
            for (int j = i + 1; j < n; j++) {
                U[k][j] -= U[i][j] * m;
            }
        }
    }
    return 0;
}

RealNumber **SolveLinearSystems(
    RealNumber **B,
    int n,
    Time *averageLinearSystemTime,
    RealNumber **L,
    RealNumber **U
) {
    RealNumber **Y = AllocateLinearSystem(n, PointerToPointer)->A;
    if (Y == NULL) {
        fprintf(stderr, "could not allocate \"Y\" matrix\n");
        return NULL;
    }
    // 2) Get the Y matrix by solving -> LY = B for each y and b;
    Time linearSystemTime;
    for (int k = 0; k < n; ++k) {
        linearSystemTime = GetTimestamp();
        for (int i = 0; i < n; i++) {
            Y[i][k] = B[i][k];
            for (unsigned int j = 0; j < i; j++) {
                Y[i][k] -= L[i][j] * Y[j][k];
            }
            Y[i][k] /= L[i][i];
        }
        linearSystemTime = GetTimestamp() - linearSystemTime;
        *averageLinearSystemTime += linearSystemTime;
    }

    // 3) Get the inverted matrix X by solving -> UX = Y for each y and x
    RealNumber **X = AllocateLinearSystem(n, PointerToPointer)->A;
    if (X == NULL) {
        fprintf(stderr, "could not allocate \"X\" matrix\n");
        return NULL;
    }
    for (int k = 0; k < n; ++k) {
        linearSystemTime = GetTimestamp();
        for (int i = n - 1; i >= 0; i--) {
            X[i][k] = Y[i][k];
            for (unsigned int j = i + 1; j < n; j++) {
                X[i][k] -= U[i][j] * X[j][k];
            }
            X[i][k] /= U[i][i];
        }
        linearSystemTime = GetTimestamp() - linearSystemTime;
        *averageLinearSystemTime += linearSystemTime;
    }

    return X;
}

RealNumber CalculateResidueL2Norm(RealNumber **A, RealNumber **B, RealNumber **invertedA, int n) {
    LIKWID_MARKER_START("RESIDUE_CALCULATION");
    RealNumber **multiplication = multiplyMatricesOfEqualSize(A, invertedA, n);
    RealNumber **R = subtractMatrices(B, multiplication, n);
    LIKWID_MARKER_STOP("RESIDUE_CALCULATION");
    RealNumber sum = 0.;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += pow(R[i][j], 2);
        }
    }
    return sqrt(sum);
}
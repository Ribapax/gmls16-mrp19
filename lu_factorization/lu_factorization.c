/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include "lu_factorization.h"
#include "../io/io.h"
#include "../linear_system/linear_system.h"
#include <math.h>
#include <stdio.h>

// This block enables to compile the code with and without the LIKWID header in
// place
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

PivotArray *AllocatePivotamento(unsigned int n) {
    PivotArray *aux = malloc(sizeof(PivotArray));
    aux->plinha = (int *) malloc(n * sizeof(int));
    aux->olinha = (int *) malloc(n * sizeof(int));
    aux->tam = 0;
    return aux;
}

void addLineToPivotArray(PivotArray* restrict P, int i, int pivotIndex) {
    P->olinha[P->tam] = i;
    P->plinha[P->tam] = pivotIndex;
    P->tam++;
}

unsigned int findPivotIndex(
    RealNumber* restrict Matrix,
    unsigned int columnIndex,
    unsigned int systemSize
) {
    RealNumber greatestValue = fabs(Matrix[Index(columnIndex, columnIndex, systemSize)]);
    unsigned int pivotIndex = columnIndex;
    for (register unsigned int i = columnIndex + 1; i < systemSize; i++) {
        RealNumber v = fabs(Matrix[Index(i, columnIndex, systemSize)]);
        if (v > greatestValue) {
            greatestValue = v;
            pivotIndex = i;
        }
    }
    return pivotIndex;
}

void replaceLinesWithIdentityMatrix(
    RealNumber* restrict Matrix,
    unsigned int index,
    unsigned int pivotIndex,
    unsigned int n
) {
    double *aux = malloc(n * sizeof(double));
    for (register unsigned int i = 0; i < n; ++i) {
        aux[i] = Matrix[Index(index, i, n)];
        Matrix[Index(index, i, n)] = Matrix[Index(pivotIndex, i, n)];
    }
    for (register unsigned int i = 0; i < n; ++i) {
        Matrix[Index(pivotIndex, i, n)] = aux[i];
    }
}

int LUDecomposition(
    RealNumber* A,
    RealNumber* U,
    RealNumber* L,
    PivotArray* P,
    int n
) {
    copyMatrix(A, U, n);
    for (register unsigned int i = 0; i < n; i++) {
        if (ENABLE_PARTIAL_PIVOTING) {
            unsigned int pivotIndex = findPivotIndex(U, i, n);
            if (i != pivotIndex) {
                replaceLinesWithIdentityMatrix(U, i, pivotIndex, n);
                replaceLinesWithIdentityMatrix(L, i, pivotIndex, n);
                addLineToPivotArray(P, i, pivotIndex);
            }
        }

        L[Index(i, i, n)] = 1;

        // Triangularization
        for (register unsigned int k = i + 1; k < n; k++) {
            if (fabs(U[Index(i, i, n)] - 0.0) < RESIDUE_THRESHOLD) {
                fprintf(stderr, "%s\n", "error: division by zero");
                return -1;
            }
            double m = U[Index(k, i, n)] / U[Index(i, i, n)];

            L[Index(k, i, n)] = m; // Filling multiplier's matrix

            U[Index(k, i, n)] = 0.0;
            for (int j = i + 1; j < n; j++) {
                U[Index(k, j, n)] -= U[Index(i, j, n)] * m;
            }
        }
    }
    return 0;
}

RealNumber *SolveLinearSystems(
    const RealNumber* B,
    int n,
    Time *averageLinearSystemTime,
    const RealNumber* L,
    const RealNumber* U
) {
    RealNumber *Y = AllocateMatrix(n);
    if (Y == NULL) {
        fprintf(stderr, "could not allocate \"Y\" matrix\n");
        return NULL;
    }

    // Get the Y matrix by solving -> LY = B for each y and b;
    Time linearSystemTime;
    for (register unsigned int k = 0; k < n; ++k) {
        linearSystemTime = GetTimestamp();
        for (register unsigned int i = 0; i < n; i++) {
            Y[Index(i, k, n)] = B[Index(i, k, n)];
            for (register unsigned int j = 0; j < i; j++) {
                Y[Index(i, k, n)] -= L[Index(i, j, n)] * Y[Index(j, k, n)];
            }
            Y[Index(i, k, n)] /= L[Index(i, i, n)];
        }
        linearSystemTime = GetTimestamp() - linearSystemTime;
        *averageLinearSystemTime += linearSystemTime;
    }

    //  Get the inverted matrix X by solving -> UX = Y for each y and x
    RealNumber *X = AllocateMatrix(n);
    if (X == NULL) {
        fprintf(stderr, "could not allocate \"X\" matrix\n");
        return NULL;
    }
    for (register unsigned int k = 0; k < n; ++k) {
        linearSystemTime = GetTimestamp();
        for (register unsigned int i = n - 1; i >= 0; i--) {
            X[Index(i, k, n)] = Y[Index(i, k, n)];
            for (register unsigned int j = i + 1; j < n; j++) {
                X[Index(i, k, n)] -= U[Index(i, j, n)] * X[Index(j, k, n)];
            }
            X[Index(i, k, n)] /= U[Index(i, i, n)];
        }
        linearSystemTime = GetTimestamp() - linearSystemTime;
        *averageLinearSystemTime += linearSystemTime;
    }

    return X;
}

RealNumber CalculateResidueL2Norm(
    RealNumber* A,
    RealNumber* B,
    RealNumber* invertedA,
    int n
) {
    LIKWID_MARKER_START("RESIDUE_CALCULATION");
    RealNumber *multiplication = multiplyMatricesOfEqualSize(A, invertedA, n);
    RealNumber *R = subtractMatrices(B, multiplication, n);
    LIKWID_MARKER_STOP("RESIDUE_CALCULATION");
    RealNumber sum = 0.;
    for (register unsigned int i = 0; i < n; ++i) {
        for (register unsigned int j = 0; j < n; ++j) {
            sum += R[Index(i, j, n)] * R[Index(i, j, n)];
        }
    }

    return sqrt(sum);
}
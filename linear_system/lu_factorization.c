/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <math.h>


#include "lu_factorization.h"
#include "linear_system.h"
#include "../io/io.h"

unsigned int findPivotIndex(double** Matrix, unsigned int columnIndex, unsigned int systemSize) {
    RealNumber greatestValue = fabs(Matrix[columnIndex][columnIndex]);
    unsigned int pivotIndex = columnIndex;
    for (unsigned int i = columnIndex + 1; i < systemSize; i++) {
        // FIXME: we probably need a decent float comparison here
        RealNumber v = fabs(Matrix[i][columnIndex]);
        if (v > greatestValue) {
            greatestValue = v;
            pivotIndex = i;
        }
//        if (fabs(Matrix[i][columnIndex]) > fabs(greatestValue)) {
//            greatestValue = Matrix[i][columnIndex];
//            pivotIndex = i;
//        }
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

RealNumber **generateMatrixL(RealNumber **A, RealNumber **B, RealNumber **U, int n) {
    RealNumber **L = AllocateLinearSystem(n, PointerToPointer)->A;
    if (L == NULL) {
        return NULL;
    }
    copyMatrix(A, U, n);
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;

        // Partial pivoting
//        unsigned int pivotIndex = findPivotIndex(U, i, n);
//        if (i != pivotIndex) {
//            replaceLinesWithIdentityMatrix(U, B, i, pivotIndex);
//        }

        // Triangularization
        for (int k = i + 1; k < n; k++) {
            if (U[i][i] == 0) {
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

RealNumber **InvertMatrix(
    RealNumber **A,
    RealNumber **B,
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

    printf("\nL MATRIX\n");
    PrintMatrix(stdout, L, n);
    printf("\nL MATRIX\n");

    printf("\nU MATRIX\n");
    PrintMatrix(stdout, U, n);
    printf("\nU MATRIX\n");

    printf("\nB MATRIX\n");
    PrintMatrix(stdout, B, n);
    printf("\nB MATRIX\n");

    // TODO: calculate times


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
    printf("\nY MATRIX\n");
    PrintMatrix(stdout, Y, n);
    printf("\nY MATRIX\n");

    // Old version of (2) - Ribamar
//    *averageLinearSystemTime = 0;
//    for (int i = 0; i < n; ++i) {
//        time = timestamp();
//        for (int k = 0; k < n; k++) { // Iterate over Y columns
//            for (int line = 0; line < n; line++ ) {
//                if (line == 0) {
//                    Y[line][k] = B[line][k];
//                } else {
//                    RealNumber sum = 0.0;
//                    for (int j = 0; j < line; j++) {
//                        sum += Y[j][k] * L[line][j];
//                    }
//                    Y[line][k] = B[line][k] - sum;
//                }
//            }
//        }
//        time = timestamp() - time;
//        *averageLinearSystemTime += time;
//    }

    // 3) Get the inverted matrix x by solving -> U=x = y for each y.
    for (int k = 0; k < n; ++k) {
        for (int i = n - 1; i >= 0; i--) {
            X[i][k] = Y[i][k];
            for (unsigned int j = i + 1; j < n; j++) {
                X[i][k] -= U[i][j] * X[j][k];
            }
            X[i][k] /= U[i][i];
        }
    }

    printf("\nX MATRIX\n");
    PrintMatrix(stdout, X, n);
    printf("\nX MATRIX\n");

    // Old version of (3) - Ribamar
//    for (int i = 0; i < n; ++i) {
//        time = timestamp();
//        for (int k = n - 1; k >= 0; k--) { // iterates over X columns
//            RealNumber retro;
//            for (int i = n - 1; i >= 0 ; --i) {
//                retro = 0;
//                for (int j = i + 1; j < n; j++) {
//                    retro += U[i][j] * X[j][k];
//                }
//                X[i][k] = (Y[i][k] - retro) / U[i][i];
//            }
//        }
//        time = timestamp() - time;
//        *averageLinearSystemTime += time;
//    }
//    *averageLinearSystemTime /= (n + n);
    // 4) Return x.
    //return invertedMatrix the X matrix;
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
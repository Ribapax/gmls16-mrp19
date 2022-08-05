/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <stdlib.h>

#include "../common/common.h"
#include "linear_system.h"
#include "lu_factorization.h"

LinearSystem* AllocateLinearSystem(unsigned int n, MatrixAllocationType type) {
    LinearSystem *SL = (LinearSystem *) malloc(sizeof(LinearSystem));
    if (!SL) {
        return SL;
    }

    SL->n = n;
    SL->matrixAllocationType = type;
    SL->A = (RealNumber **) calloc(n, sizeof(RealNumber *));
    SL->b = (RealNumber *) malloc(n * sizeof(RealNumber));
    if (!(SL->A) || !(SL->b)) {
        freeLinearSystem(SL);
        return NULL;
    }

    if (type == PointerToArray) {
        SL->A[0] = (RealNumber *) malloc(n * n * sizeof(RealNumber));
        if (!(SL->A[0])) {
            freeLinearSystem(SL);
            return NULL;
        }
        for (int i = 1; i < n; ++i) {
            SL->A[i] = SL->A[i - 1] + n;
        }
    } else if (type == PointerToPointer) {
        for (int i=0; i < n; ++i) {
            SL->A[i] = (RealNumber *) calloc(n, sizeof(RealNumber));
        }
    }

    return SL;
}

void freeLinearSystem(LinearSystem *SL) {
    if (!SL) {
        return;
    }
    if (SL->b) {
        free(SL->b);
    }
    if (!SL->A) {
        free(SL);
        return;
    }
    if (SL->matrixAllocationType == PointerToArray && SL->A[0]) {
        free(SL->A[0]);
    } else if (SL->matrixAllocationType == PointerToPointer) {
        for (int i=0; i < SL->n; ++i) {
            free(SL->A[i]);
        }
    }
    free(SL->A);
    free(SL);
}

int FillLinearSystem(LinearSystem *SL, MatrixType type, RealNumber coefficientLimit) {
    unsigned int n = SL->n;
    RealNumber invRandMax = ((RealNumber) coefficientLimit / (RealNumber) RAND_MAX);
    for (unsigned int i = 0; i < n; ++i) {
        SL->b[i] = (RealNumber) rand() * invRandMax;
    }
    if (type == HilbertMatrix) {
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                SL->A[i][j] = 1.0 / (RealNumber) (i + j + 1);
            }
        }
        return 0;
    }
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[i][j] = (RealNumber) rand() * invRandMax;
        }
    }
    if (type == NullEquation) {
        unsigned int nula = rand() % n;
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[nula][j] = 0.0;
        }
        SL->b[nula] = 0.0;
        return 0;
    } else if (type == ProportionalEquation) {
        unsigned int propDst = rand() % n;
        unsigned int propSrc = (propDst + 1) % n;
        RealNumber mult = (RealNumber) rand() * invRandMax;
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[propDst][j] = SL->A[propSrc][j] * mult;
        }
        SL->b[propDst] = SL->b[propSrc] * mult;
        return 0;
    } else if (type == LinearCombinationEquation) {
        unsigned int combDst = rand() % n;
        unsigned int combSrc1 = (combDst + 1) % n;
        unsigned int combSrc2 = (combDst + 2) % n;
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[combDst][j] = SL->A[combSrc1][j] + SL->A[combSrc2][j];
        }
        SL->b[combDst] = SL->b[combSrc1] + SL->b[combSrc2];
        return 0;
    } else if (type == DominantDiagonal) {
        for (unsigned int i = 0; i < n; ++i) {
            RealNumber soma = 0.0;
            for (unsigned int j = 0; j < i; ++j) soma += SL->A[i][j];
            for (unsigned int j = i + 1; j < n; ++j) soma += SL->A[i][j];
            SL->A[i][i] += soma;
        }
        return 0;
    } else if (type == GenericMatrix) {
        return 0;
    } else {
        return -1;
    }
}

void copyMatrix(RealNumber **A, RealNumber **B, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = A[i][j];
        }
    }
}

RealNumber **multiplyMatricesOfEqualSize(RealNumber **A, RealNumber **B, int n) {
    RealNumber **Result = AllocateLinearSystem(n, PointerToPointer)->A;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                Result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return Result;
}

RealNumber **subtractMatrices(RealNumber **A, RealNumber **B, int n) {
    RealNumber **Result = AllocateLinearSystem(n, PointerToPointer)->A;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Result[i][j] = A[i][j] - B[i][j];
        }
    }
    return Result;
}

RealNumber **GenerateIdentityMatrix(int n) {
    RealNumber **I = AllocateLinearSystem(n, PointerToPointer)->A;
    if (I == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                I[i][j] = 1;
            } else {
                I[i][j] = 0;
            }
        }
    }
    return I;
}

int MatrixIsInvertible(RealNumber **A, int n) {
    RealNumber acc = 1;
    for (int i = 0; i < n; ++i) {
        acc *= A[i][i];
    }
    if (acc == 0) {
        return 0;
    }
    return 1;
}

RealNumber **RefineSolution(
    RealNumber **A,
    RealNumber **B,
    RealNumber **invertedMatrix,
    RealNumber **L,
    RealNumber **U,
    int n
) {
    RealNumber **residue;
    // 1) A x A^-1
    residue = multiplyMatricesOfEqualSize(A, invertedMatrix, n);

    // 2) B - (A x A^-1)
    residue = subtractMatrices(B, residue, n);

    // 3) AW = B - (A x A^-1) -> for each W column
    // 3.1) Get the Y matrix by solving -> Ly = B - (A x A^-1); for each y.
    RealNumber **Y = AllocateLinearSystem(n, PointerToPointer)->A;
    if (Y == NULL) {
        fprintf(stderr, "could not allocate \"Y\" matrix\n");
        return NULL;
    }
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; i++) {
            Y[i][k] = residue[i][k];
            for (unsigned int j = 0; j < i; j++) {
                Y[i][k] -= L[i][j] * Y[j][k];
            }
            Y[i][k] /= L[i][i];
        }
    }

    // 3.2) Get the inverted matrix X by solving -> Ux = y for each y and x.
    RealNumber **X = AllocateLinearSystem(n, PointerToPointer)->A;
    if (X == NULL) {
        fprintf(stderr, "could not allocate \"X\" matrix\n");
        return NULL;
    }
    for (int k = 0; k < n; ++k) {
        for (int i = n - 1; i >= 0; i--) {
            X[i][k] = Y[i][k];
            for (unsigned int j = i + 1; j < n; j++) {
                X[i][k] -= U[i][j] * X[j][k];
            }
            X[i][k] /= U[i][i];
        }
    }

    // 4) X(1) = X(0) + W
    RealNumber **refinedSolution = AllocateLinearSystem(n, PointerToPointer)->A;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            refinedSolution[i][j] = invertedMatrix[i][j] + X[i][j];
        }
    }

    return refinedSolution;
}

int HasNotReachedStoppingCriteria(
    int iteration,
    int iterationsLimit,
    RealNumber currentResidueL2Norm,
    RealNumber lastResidueL2Norm
) {
    if (ResidueIsIncreasing(currentResidueL2Norm, lastResidueL2Norm)) {
        return 0;
    }
    if (iteration <= iterationsLimit && currentResidueL2Norm > RESIDUE_THRESHOLD) {
        return 1;
    }
    return 0;
}

int ResidueIsIncreasing(RealNumber currentResidueL2Norm, RealNumber lastResidueL2Norm) {
    if (
        lastResidueL2Norm != 1 + RESIDUE_THRESHOLD &&
        (currentResidueL2Norm - lastResidueL2Norm) > RESIDUE_THRESHOLD
    ) {
        return 1;
    }
    return 0;
}

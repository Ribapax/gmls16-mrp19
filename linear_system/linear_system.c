#include <stdio.h>
#include <stdlib.h>

#include "../common/common.h"
#include "linear_system.h"
#include "lu_factorization.h"

// TODO: error handling
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


/*!
  \brief Generate coefficients and independent terms
  \param SL
  \param type
  \param coefficientLimit
*/
int FillLinearSystem(LinearSystem *SL, MatrixType type, RealNumber coefficientLimit) {
    unsigned int n = SL->n;
    RealNumber invRandMax = ((RealNumber) coefficientLimit / (RealNumber) RAND_MAX);
    for (unsigned int i=0; i<n; ++i) {
        SL->b[i] = (RealNumber) rand() * invRandMax;
    }
    if (type == HilbertMatrix) {
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j)  {
                SL->A[i][j] = 1.0 / (RealNumber)(i + j + 1);
            }
        }
        return 0;
    }
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < n; ++j)  {
            SL->A[i][j] = (RealNumber) rand() * invRandMax;
        }
    }
    if (type == NullEquation) {
        unsigned int nula = rand() % n;
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[nula][j] = 0.0;
        }
        SL->b[nula] = 0.0;
    } else if (type == ProportionalEquation) {
        unsigned int propDst = rand() % n;
        unsigned int propSrc = (propDst + 1) % n;
        RealNumber mult = (RealNumber) rand() * invRandMax;
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[propDst][j] = SL->A[propSrc][j] * mult;
        }
        SL->b[propDst] = SL->b[propSrc] * mult;
    } else if (type == LinearCombinationEquation) {
        unsigned int combDst = rand() % n;
        unsigned int combSrc1 = (combDst + 1) % n;
        unsigned int combSrc2 = (combDst + 2) % n;
        for (unsigned int j = 0; j < n; ++j) {
            SL->A[combDst][j] = SL->A[combSrc1][j] + SL->A[combSrc2][j];
        }
        SL->b[combDst] = SL->b[combSrc1] + SL->b[combSrc2];
    } else if (type == DominantDiagonal) {
        for (unsigned int i = 0; i < n; ++i) {
            RealNumber soma = 0.0;
            for (unsigned int j = 0; j < i; ++j) soma += SL->A[i][j];
            for (unsigned int j = i + 1; j < n; ++j) soma += SL->A[i][j];
            SL->A[i][i] += soma;
        }
    } else {
        return -1;
    }
    return 0;
}


LinearSystem *readLinearSystem(MatrixAllocationType type) {
    unsigned int n;
    LinearSystem *SL;
    scanf("%d", &n);
    SL = AllocateLinearSystem(n, type);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scanf("%lg", &SL->A[i][j]);
        }
    }
    for(int i=0; i < n; ++i) {
        scanf ("%lg", &SL->b[i]);
    }

    return SL;
}


void printLinearSystem(LinearSystem *SL) {
    int n = SL->n;
    for(int i = 0; i < n; ++i) {
        printf("\n  ");
        for(int j = 0; j < n; ++j) {
            printf ("%10g", SL->A[i][j]);
        }
        printf ("   |   %g", SL->b[i]);
    }
    printf("\n\n");
}

void printArray(RealNumber *arr, unsigned int n) {
    int i;
    printf ("\n");
    for(i = 0; i < n; ++i) {
        printf ("%10g ", arr[i]);
    }
    printf ("\n\n");
}

void copyMatrix(RealNumber **A, RealNumber **B, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = A[i][j];
        }
    }
}

RealNumber **multiplyMatrixOfEqualSize(RealNumber **A, RealNumber **B, int n) {
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

RealNumber **subtractMatrix(RealNumber **A, RealNumber **B, int n) {
    RealNumber **Result = AllocateLinearSystem(n, PointerToPointer)->A;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Result[i][j] = A[i][j] - B[i][j];
        }
    }
    return Result;
}

RealNumber *subtractArrays(const RealNumber *A, const RealNumber *B, int n) {
    RealNumber *result = malloc(sizeof(RealNumber)*n);
    for (int i = 0; i < n; ++i){
        result[i]= A[i] - B[i];
    }
    return result;
}

RealNumber **GetIdentityMatrix(int n) {
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

void forwardSubstitution(RealNumber **A, const RealNumber *b, RealNumber *x, unsigned int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (unsigned int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

void replaceLines(
    double **Matrix,
    double *independentTerms,
    unsigned int index,
    unsigned int pivotIndex
) {
    RealNumber *lineToBeReplaced = Matrix[index];
    Matrix[index] = Matrix[pivotIndex];
    Matrix[pivotIndex] = lineToBeReplaced;

    RealNumber termToBeReplaced = independentTerms[index];
    independentTerms[index] = independentTerms[pivotIndex];
    independentTerms[pivotIndex] = termToBeReplaced;
}

RealNumber *GaussElimination(RealNumber **A, RealNumber *B, int n) {
    for (int i = 0; i < n; i++) {
        // Partial Pivoting
//        unsigned int pivotIndex = findPivotIndex(A, i, n);
//        if (i != pivotIndex) {
//            replaceLines(A, B, i, pivotIndex);
//        }
        for (int k = i + 1; k < n; k++) {
            if (A[k][k] == 0) {
                fprintf(stderr, "%s\n", "gaussian elimination error: division by zero");
                exit(-1);
            }
            double m = A[k][i] / A[i][i];
            A[k][i] = 0.0;
            for (int j = i + 1; j < n; j++) {
                A[k][j] -= A[i][j] * m;
            }
            B[k] -= B[i] * m;
        }
    }
    RealNumber *x = malloc(sizeof(RealNumber) * n);
    forwardSubstitution(A, B, x, n);
    return x;
}


RealNumber *multiplyMatrixWithArray(RealNumber **A, const RealNumber *B, int n) {
    RealNumber *solution = malloc(sizeof(RealNumber)*n);
    RealNumber sum = 0.;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * B[j];
        }
        solution[i] = sum;
        sum = 0.;
    }
    return solution;
}

RealNumber **refineSolution(RealNumber **A, RealNumber **invertedMatrix, int n) {
    RealNumber **refinedSolutionCOLUMNS = AllocateLinearSystem(n, PointerToPointer)->A;
    RealNumber **partialRefinedSolution = AllocateLinearSystem(n, PointerToPointer)->A;
    RealNumber **refinedSolution = AllocateLinearSystem(n, PointerToPointer)->A;
    RealNumber **identityMatrix = GetIdentityMatrix(n);
    RealNumber **residue;
    residue = multiplyMatrixOfEqualSize(A, invertedMatrix, n);
    residue = subtractMatrix(identityMatrix, residue, n);
    for (int i = 0; i < n; ++i) {
        refinedSolutionCOLUMNS[i] = GaussElimination(A, residue[i], n);
        for (int j = 0; j < n; ++j) {
            partialRefinedSolution[j][i] = refinedSolutionCOLUMNS[i][j];
        }
    }
    // TODO: sum matrix function
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            refinedSolution[i][j] = partialRefinedSolution[i][j] + invertedMatrix[i][j];
        }
    }
    return refinedSolution;
}

int hasNotReachedStoppingCriteria(
    int iteration,
    int iterationsLimit,
    RealNumber currentResidueL2Norm,
    RealNumber lastResidueL2Norm
) {
    if (
        lastResidueL2Norm != 1 + RESIDUE_THRESHOLD &&
        (currentResidueL2Norm - lastResidueL2Norm) > RESIDUE_THRESHOLD
    ) {
        fprintf(
            stderr,
            "\nerror: residue increasing, solution does not converge, stopping. Try to decrease the number of iterations.\n"
        );
        return 0;
    }
    if (iteration <= iterationsLimit && currentResidueL2Norm > RESIDUE_THRESHOLD) {
        return 1;
    }
    return 0;
}

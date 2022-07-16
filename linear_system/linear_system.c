#include <stdio.h>
#include <stdlib.h>

#include "../common/common.h"
#include "linear_system.h"

// TODO: error handling
LinearSystem* allocateLinearSystem(unsigned int n, MatrixAllocationType type) {
    LinearSystem *SL = (LinearSystem *) malloc(sizeof(LinearSystem));
    if (!SL) {
        return SL;
    }

    SL->n = n;
    SL->matrixAllocationType = type;
    SL->A = (RealNumber **) malloc(n * sizeof(RealNumber *));
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
            SL->A[i] = (RealNumber *) malloc(n * sizeof(RealNumber));
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
void fillLinearSystem(LinearSystem *SL, MatrixType type, RealNumber coefficientLimit) {
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
        return;
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
        // TODO: handle error
        return;
    }
}


LinearSystem *readLinearSystem(MatrixAllocationType type) {
    unsigned int n;
    LinearSystem *SL;
    scanf("%d", &n);
    SL = allocateLinearSystem(n, type);
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


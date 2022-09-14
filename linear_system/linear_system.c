/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <stdlib.h>

#include "../common/common.h"
#include "../io/io.h"
#include "../lu_factorization/lu_factorization.h"
#include "linear_system.h"

RealNumber *AllocateMatrix(unsigned int n) {
  return (RealNumber *)calloc(n * n, sizeof(RealNumber));
}

int FillMatrix(RealNumber *A, RealNumber coefficientLimit, unsigned int n) {
  RealNumber invRandMax = ((RealNumber)coefficientLimit / (RealNumber)RAND_MAX);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      A[Index(i, j, n)] = (RealNumber)rand() * invRandMax;
    }
  }
  return 0;
}

void copyMatrix(const RealNumber *A, RealNumber *B, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      B[Index(i, j, n)] = A[Index(i, j, n)];
    }
  }
}

RealNumber *multiplyMatricesOfEqualSize(const RealNumber *A,
                                        const RealNumber *B, int n) {
  RealNumber *Result = AllocateMatrix(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        Result[Index(i, j, n)] += A[Index(i, k, n)] * B[Index(k, j, n)];
      }
    }
  }
  return Result;
}

RealNumber *subtractMatrices(const RealNumber *A, const RealNumber *B, int n) {
  RealNumber *Result = AllocateMatrix(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      Result[Index(i, j, n)] = A[Index(i, j, n)] - B[Index(i, j, n)];
    }
  }
  return Result;
}

RealNumber *GenerateIdentityMatrix(int n) {
  RealNumber *I = AllocateMatrix(n);
  if (I == NULL) {
    return NULL;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        I[Index(i, j, n)] = 1;
      } else {
        I[Index(i, j, n)] = 0;
      }
    }
  }
  return I;
}

int MatrixIsInvertible(const RealNumber *A, int n) {
  RealNumber acc = 1;
  for (int i = 0; i < n; ++i) {
    acc *= A[Index(i, i, n)];
  }
  if (acc == 0) {
    return 0;
  }
  return 1;
}

RealNumber *RefineSolution(RealNumber *A, RealNumber *B,
                           RealNumber *invertedMatrix, RealNumber *L,
                           RealNumber *U, PivotArray *P,
                           Time *averageLinearSystemTime, int n) {
  RealNumber *residue;
  // 1) A x A^-1
  residue = multiplyMatricesOfEqualSize(A, invertedMatrix, n);

  if (ENABLE_PARTIAL_PIVOTING) {

    fprintf(stdout, "\n%d\n", P->tam);
    for (int i = 0; i < P->tam; i++) {
      fprintf(stdout, "%d %d\n", P->olinha[i], P->plinha[i]);
      replaceLinesWithIdentityMatrix(residue, P->olinha[i], P->plinha[i], n);
      // fprintf(stdout,"%d %d",P->olinha[i],P->plinha[i]);
    }
    fprintf(stdout, "\nResidue - SL\n");
    PrintMatrix(stdout, residue, n);
  }

  // 2) B - (A x A^-1)

  residue = subtractMatrices(B, residue, n);

  

  // 3) AW = B - (A x A^-1)
  RealNumber *X =
      SolveLinearSystems(residue, n, averageLinearSystemTime, L, P, U);

  // 4) X(1) = X(0) + W
  RealNumber *refinedSolution = AllocateMatrix(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      refinedSolution[Index(i, j, n)] =
          invertedMatrix[Index(i, j, n)] + X[Index(i, j, n)];
    }
  }

  return refinedSolution;
}

int HasNotReachedStoppingCriteria(int iteration, int iterationsLimit,
                                  RealNumber currentResidueL2Norm,
                                  RealNumber lastResidueL2Norm) {
  if (ResidueIsIncreasing(currentResidueL2Norm, lastResidueL2Norm)) {
    return 0;
  }
  if (iteration <= iterationsLimit &&
      currentResidueL2Norm > RESIDUE_THRESHOLD) {
    return 1;
  }
  return 0;
}

int ResidueIsIncreasing(RealNumber currentResidueL2Norm,
                        RealNumber lastResidueL2Norm) {
  if (lastResidueL2Norm != 1 + RESIDUE_THRESHOLD &&
      (currentResidueL2Norm - lastResidueL2Norm) > RESIDUE_THRESHOLD) {
    return 1;
  }
  return 0;
}
